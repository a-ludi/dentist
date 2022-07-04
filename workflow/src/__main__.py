import json
import logging
import os
import shlex
from itertools import chain
from os import environ
from pathlib import Path

import util.dentist as dentist
from dentist.workflow.engine import (
    FileList,
    ShellScript,
    Workflow,
    cli_parser,
    python_code,
    safe,
)
from util.dazzler import *

log = logging.getLogger(__name__)


class DentistGapClosing(Workflow):
    dust_mask = "dust"
    tandem_mask = "tan"

    def __init__(
        self,
        *args,
        workflow_config,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        with workflow_config.open() as config:
            self.config = json.load(config)
            self.config["__path__"] = workflow_config.resolve()

        log.info(f"changing directory to {workflow_config.parent}")
        os.chdir(workflow_config.parent)

        self.reference_fasta = Path(self.config["inputs"]["reference"])
        self.reads_fasta = Path(self.config["inputs"]["reads"])
        self.workdir = Path(self.config["workdir"])
        self.logdir = Path(self.config["logdir"])
        self.reads_type = self.config["inputs"]["reads_type"]
        self.dbsplit_options = {
            "reference": self.config["reference_dbsplit"],
            "reads": self.config["reads_dbsplit"],
        }
        self.lamerge_opts = ["-v"]
        if "TMPDIR" in environ:
            self.lamerge_opts.append(f"-P{environ['TMPDIR']}")
        self.tanmask_opts = self.config.get("tanmask", [])
        self.tanmask_opts = ensure_flags(
            self.tanmask_opts,
            {
                "-v": True,
                "-n": self.tandem_mask,
            },
        )
        self.self_mask = self.config.get("self_mask", "dentist-self")
        self.reads_mask = self.config.get("reads_mask", "dentist-reads")
        self.dentist_flags = shlex.split(environ.get("DENTIST_FLAGS", ""))

    def run(self):
        self.create_dirs("create_workdirs", [self.workdir, self.logdir])
        self.execute_jobs()

        self.dentist_config = self.create_dentist_config()
        self.reference = self.fasta2dazzler(self.reference_fasta, role="reference")
        self.reads = self.fasta2dazzler(self.reads_fasta, role="reads")
        self.execute_jobs()

        dentist.validate_config(self.dentist_config)

        self.mask_dust(self.reference)
        self.execute_jobs()
        self.tandem_alignment(self.reference)
        self.execute_jobs()
        self.mask_tandem(self.reference)
        self.execute_jobs()
        self.self_alignment(self.reference)
        self.execute_jobs()
        self.mask_self(self.reference)
        self.execute_jobs()
        self.ref_vs_reads_alignment(self.reference, self.reads)
        self.execute_jobs()
        self.mask_reads()
        self.execute_jobs()

    def create_dentist_config(self):
        @self.collect_job(
            inputs=[self.config["__path__"]],
            outputs=[self.workdir / "dentist.json"],
            exec_local=True,
        )
        @python_code
        def create_dentist_config(inputs, outputs):
            with inputs[0].open() as cf:
                dentist_config = json.load(cf).get("dentist", {})
            with outputs[0].open("w") as cf:
                json.dump(dentist_config, cf)

        return self.jobs["create_dentist_config"].outputs[0]

    def create_dirs(self, name, dirs):
        @python_code
        def _create_dirs(outputs):
            for directory in outputs:
                directory.mkdir(parents=True, exist_ok=True)

        self.collect_job(
            name=name,
            inputs=[],
            outputs=dirs,
            exec_local=True,
            action=_create_dirs,
        )

    def fasta2dazzler(self, fasta, role):
        if role == "reads" and self.reads_type == "PACBIO_SMRT":
            db_ext = ".db"
            fasta2dazz = "fasta2DB"
        else:
            db_ext = ".dam"
            fasta2dazz = "fasta2DAM"
        workdb = fasta_to_workdb(fasta, db_ext, self.workdir)
        self.collect_job(
            name=f"{role}2dazzler",
            inputs=[fasta],
            outputs=db_files(workdb),
            action=lambda inputs, outputs: ShellScript(
                # make sure no DB is present because sequences would be
                # appended instead of creating a fresh one
                ("rm", "-f", *outputs),
                # make sure FASTA lines are not too long for fasta2dazz command
                ("fold", "-w9998", inputs[0], safe("|"), fasta2dazz, "-i", outputs[0]),
                # splitting DBs is required by the pipeline
                ("DBsplit", *self.dbsplit_options[role], outputs[0]),
            ),
        )

        return workdb

    def mask_dust(self, db):
        dustcmd = dentist.generate_options_for("DBdust reference", self.dentist_config)
        dustcmd = deduplicate_flags(dustcmd)

        return self.collect_job(
            name=f"mask_dust_{db.stem}",
            inputs=db_files(db),
            outputs=mask_files(db, self.dust_mask),
            action=lambda inputs: ShellScript(
                (*dustcmd, inputs[0].with_suffix("")),
            ),
        )

    def tandem_alignment(self, db):
        with self.grouped_jobs(f"{__name__}.{db.stem}"):
            for i in range(get_num_blocks(db)):
                self.tandem_alignment_block(db, block=i + 1)
            self.execute_jobs()
            self.lamerge(
                self.workdir / alignment_file("TAN", db),
                [
                    self.workdir / las
                    for las in block_alignments("TAN", db, block_a=FULL_DB)
                ],
                job=f"tandem_alignment_{db.stem}",
                log=self.log_file(f"tandem-alignment.{db.stem}"),
            )

    def tandem_alignment_block(self, db, block):
        aligncmd = dentist.generate_options_for("tandem", self.dentist_config)
        aligncmd = ensure_threads_flag(aligncmd, threads=1)

        return self.collect_job(
            name=f"tandem_alignment_block_{db.stem}_{block}",
            inputs=FileList(
                db=db_files(db),
                config=self.dentist_config,
            ),
            outputs=[self.workdir / alignment_file("TAN", db, block_b=block)],
            log=self.log_file(f"tandem-alignment.{db.stem}.{block}"),
            action=lambda inputs: ShellScript(
                ("cd", self.workdir), (*aligncmd, f"{inputs.db[0].stem}.{block}")
            ),
        )

    def mask_tandem(self, db):
        with self.grouped_jobs(f"{__name__}.{db.stem}"):
            for i in range(get_num_blocks(db)):
                self.mask_tandem_block(db, block=i + 1)
            self.execute_jobs()
            self.catrack(
                db,
                self.tandem_mask,
                job=f"mask_tandem_{db.stem}",
                log=self.log_file(f"mask-tandem.{db.stem}"),
            )

    def mask_tandem_block(self, db, block):
        return self.collect_job(
            name=f"mask_tandem_block_{db.stem}_{block}",
            inputs=FileList(
                db=db_files(db),
                las=self.workdir / alignment_file("TAN", db, block_b=block),
                config=self.dentist_config,
            ),
            outputs=mask_files(db, block_mask(self.tandem_mask, block)),
            log=self.log_file(f"mask-tandem.{db.stem}.{block}"),
            action=lambda inputs: ShellScript(
                ("TANmask", *self.tanmask_opts, inputs.db[0], inputs.las),
            ),
        )

    def self_alignment(self, db):
        with self.grouped_jobs(f"{__name__}.{db.stem}"):
            num_blocks = get_num_blocks(db)
            for i in range(num_blocks):
                for j in range(i, num_blocks):
                    self.self_alignment_block(db, block_a=i + 1, block_b=j + 1)
            self.execute_jobs()
            self.lamerge(
                self.workdir / alignment_file(db),
                [self.workdir / las for las in block_alignments(db)],
                job=f"self_alignment_{db.stem}",
                log=self.log_file(f"self-alignment.{db.stem}"),
            )

    def self_alignment_block(self, db, block_a, block_b):
        aligncmd = dentist.generate_options_for(
            "self",
            self.dentist_config,
            ensure_masks([], self.dust_mask, self.tandem_mask),
        )
        aligncmd = ensure_threads_flag(aligncmd, threads=1)
        aligncmd = deduplicate_flags(aligncmd)

        return self.collect_job(
            name=f"self_alignment_block_{db.stem}_{block_a}_{block_b}",
            inputs=FileList(
                db=db_files(db),
                dust_mask=mask_files(db, self.dust_mask),
                tandem_mask=mask_files(db, self.tandem_mask),
                config=self.dentist_config,
            ),
            outputs=[
                self.workdir / alignment_file(db, block_a=block_a, block_b=block_b),
                self.workdir / alignment_file(db, block_a=block_b, block_b=block_a),
            ],
            log=self.log_file(f"self-alignment.{db.stem}.{block_a}.{block_b}"),
            action=lambda inputs: ShellScript(
                ("cd", self.workdir),
                (
                    *aligncmd,
                    f"{inputs.db[0].stem}.{block_a}",
                    f"{inputs.db[0].stem}.{block_b}",
                ),
            ),
        )

    def mask_self(self, db):
        self.collect_job(
            name=f"mask_self_{db.stem}",
            inputs=FileList(
                db=db_files(db),
                las=self.workdir / alignment_file(db),
                config=self.dentist_config,
            ),
            outputs=mask_files(db, self.self_mask),
            log=self.log_file(f"mask-self.{db.stem}"),
            action=lambda inputs: ShellScript(
                (
                    "dentist",
                    "mask",
                    f"--config={self.dentist_config}",
                    *self.dentist_flags,
                    inputs.db[0],
                    inputs.las,
                    self.self_mask,
                )
            ),
        )

    def ref_vs_reads_alignment(self, refdb, readsdb):
        with self.grouped_jobs(f"{__name__}.{refdb.stem}.{readsdb.stem}"):
            reads_blocks = get_num_blocks(readsdb)
            for j in range(reads_blocks):
                self.ref_vs_reads_alignment_block(refdb, readsdb, j + 1)
            self.execute_jobs()
            self.lamerge(
                self.workdir / alignment_file(refdb, readsdb),
                [
                    self.workdir / las
                    for las in block_alignments(refdb, readsdb, block_a=FULL_DB)
                ],
                job=f"ref_vs_reads_alignment_{refdb.stem}_{readsdb.stem}",
                log=self.log_file(
                    f"ref-vs-reads-alignment.{refdb.stem}.{readsdb.stem}"
                ),
            )

    def ref_vs_reads_alignment_block(self, refdb, readsdb, block_reads):
        aligncmd = dentist.generate_options_for(
            "reads",
            self.dentist_config,
            ensure_masks([], self.dust_mask, self.tandem_mask, self.self_mask),
        )
        aligncmd = ensure_threads_flag(aligncmd, threads=1)
        aligncmd = deduplicate_flags(aligncmd)

        self.collect_job(
            name=f"ref_vs_reads_alignment_block_{refdb.stem}_{readsdb.stem}_{block_reads}",
            inputs=FileList(
                refdb=db_files(refdb),
                readsdb=db_files(readsdb),
                dust_mask=mask_files(refdb, self.dust_mask),
                tandem_mask=mask_files(refdb, self.tandem_mask),
                self_mask=mask_files(refdb, self.self_mask),
                config=self.dentist_config,
            ),
            outputs=[
                self.workdir / alignment_file(refdb, readsdb, block_b=block_reads),
                self.workdir / alignment_file(readsdb, refdb, block_a=block_reads),
            ],
            log=self.log_file(
                f"ref-vs-reads-alignment.{refdb.stem}.{readsdb.stem}.{block_reads}"
            ),
            action=lambda inputs: ShellScript(
                ("cd", self.workdir),
                (
                    *aligncmd,
                    inputs.refdb[0].stem,
                    f"{inputs.readsdb[0].stem}.{block_reads}",
                ),
            ),
        )

    def mask_reads(self):
        self.collect_job(
            name=f"mask_reads",
            inputs=FileList(
                refdb=db_files(self.reference),
                readsdb=db_files(self.reads),
                las=self.workdir / alignment_file(self.reference, self.reads),
                config=self.dentist_config,
            ),
            outputs=mask_files(self.reference, self.reads_mask),
            log=self.log_file("mask-reads"),
            action=lambda inputs: ShellScript(
                (
                    "dentist",
                    "mask",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    inputs.refdb[0],
                    inputs.readsdb[0],
                    inputs.las,
                    self.reads_mask,
                )
            ),
        )

    def lamerge(self, merged, parts, job, log):
        self.collect_job(
            name=job,
            inputs=parts,
            outputs=[merged],
            log=log,
            action=lambda inputs, outputs: ShellScript(
                ("LAmerge", *self.lamerge_opts, outputs[0], *inputs)
            ),
        )

    def catrack(self, db, mask, job, log):
        self.collect_job(
            name=job,
            inputs=FileList(
                db=db_files(db),
                mask_files=chain.from_iterable(
                    mask_files(db, bm) for bm in block_masks(mask, db)
                ),
            ),
            outputs=mask_files(db, mask),
            log=log,
            action=lambda inputs, outputs: ShellScript(
                ("rm", "-f", *outputs),
                ("Catrack", "-v", inputs.db[0], mask),
            ),
        )

    def log_file(self, id):
        return self.logdir / f"{id}.log"


def main():
    import logging
    import sys

    logging.basicConfig(level=logging.DEBUG)
    parser = cli_parser()
    parser.add_argument(
        "workflow_config",
        metavar="<config:json>",
        type=Path,
        help="workflow configuraiton file in JSON format",
        default="workflow.json",
    )
    params = vars(parser.parse_args())

    DentistGapClosing(**params)()


if __name__ == "__main__":
    main()
