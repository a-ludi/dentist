import json
import logging
import os
from os import environ
from pathlib import Path

import util.dentist as dentist
from dentist.workflow.engine import ShellScript, Workflow, cli_parser, python_code, safe
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

    def run(self):
        self.create_dirs("create_workdirs", [self.workdir, self.logdir])
        self.execute_jobs()

        self.dentist_config = self.create_dentist_config()
        self.reference = self.fasta2dazzler(self.reference_fasta, role="reference")
        self.reads = self.fasta2dazzler(self.reads_fasta, role="reads")
        self.execute_jobs()

        dentist.validate_config(self.dentist_config)

        self.mask_dust(self.reference)
        self.tandem_alignment(self.reference)

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
                    self.workdir / a
                    for a in block_alignments("TAN", db, block_a=FULL_DB)
                ],
                job=f"tandem_alignment_{db.stem}",
                log=self.log_file(f"tandem-alignment.{db.stem}"),
            )

    def tandem_alignment_block(self, db, block):
        aligncmd = dentist.generate_options_for("tandem", self.dentist_config)
        aligncmd = ensure_threads_flag(aligncmd, threads=1)

        return self.collect_job(
            name=f"tandem_alignment_block_{db.stem}_{block}",
            inputs=[
                *db_files(db),
                self.dentist_config,
            ],
            outputs=[self.workdir / alignment_file("TAN", db.stem, block_b=block)],
            action=lambda inputs: ShellScript(
                (safe("{")),
                ("cd", self.workdir),
                (*aligncmd, f"{db.stem}.{block}"),
                (
                    safe("}"),
                    safe("&>"),
                    self.log_file(f"tandem-alignment.{db.stem}.{block}"),
                ),
            ),
        )

    def self_alignment_block(self, db, block_a, block_b):
        aligncmd = dentist.generate_options_for(
            "self", self.dentist_config, ensure_masks([], dust_mask, tandem_mask)
        )

        return self.collect_job(
            name=f"self_alignment_block_{db.stem}_{block}",
            inputs=[
                *db_files(db),
                self.dentist_config,
                *mask_files(db, self.dust_mask),
                *mask_files(db, self.tandem_mask),
            ],
            outputs=[
                alignment_file(db.stem, block_a=block_a, block_b=block_b),
                alignment_file(db.stem, block_b=block_a, block_a=block_b),
            ],
            action=lambda inputs: ShellScript(
                (safe("{")),
                ("cd", self.workdir),
                (*aligncmd, f"{db.stem}.{block_a}", f"{db.stem}.{block_b}"),
                (
                    safe("}"),
                    safe("&>"),
                    self.log_file(f"self-alignment.{db.stem}.{block_a}.{block_b}"),
                ),
            ),
        )

    def lamerge(self, merged, parts, job, log):
        self.collect_job(
            name=job,
            inputs=parts,
            outputs=[merged],
            action=lambda inputs: ShellScript(
                ("LAmerge", *self.lamerge_opts, merged, *parts, safe("&>"), log)
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
