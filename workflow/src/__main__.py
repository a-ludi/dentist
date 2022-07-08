import json
import logging
import os
import shlex
from argparse import Namespace
from itertools import chain
from os import environ
from pathlib import Path

import util.dentist as dentist
from dentist.workflow.engine import (
    FileList,
    MultiIndex,
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
    batched_jobs = (
        "tandem_alignment",
        "mask_tandem",
        "self_alignment",
        "ref_vs_reads_alignment",
        "homogenize_mask",
    )

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
        self.masks = (
            self.self_mask,
            self.tandem_mask,
            self.reads_mask,
        )
        batch_size = dict((job, 1) for job in self.batched_jobs)
        batch_size |= self.config.get("batch_size", dict())
        self.batch_size = Namespace(**batch_size)
        self.dentist_flags = shlex.split(environ.get("DENTIST_FLAGS", ""))
        self.pile_ups = self.workdir / "pile-ups.db"

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
        for mask in self.masks:
            self.homogenize_mask(self.reference, self.reads, mask)
        self.execute_jobs()
        self.collect()
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
        workdb = self.workdir / f"{role}{db_ext}"
        self.collect_job(
            name=f"{role}2dazzler",
            inputs=[fasta],
            outputs=db_files(workdb),
            exec_local=True,
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
            resources="mask_dust",
            action=lambda inputs: ShellScript(
                (*dustcmd, inputs[0].with_suffix("")),
            ),
        )

    def tandem_alignment(self, db):
        with self.grouped_jobs(f"tandem_alignment_{db.stem}"):
            for blocks in get_blocks(db, self.batch_size.tandem_alignment):
                self.tandem_alignment_block(db, blocks=blocks)
            self.execute_jobs()
            self.lamerge(
                self.workdir / alignment_file("TAN", db),
                [
                    self.workdir / las
                    for las in block_alignments("TAN", db, block_a=FULL_DB)
                ],
                job=f"tandem_alignment_{db.stem}",
                resources="tandem_alignment",
                log=self.log_file(f"tandem-alignment.{db.stem}"),
            )

    def tandem_alignment_block(self, db, blocks):
        aligncmd = dentist.generate_options_for("tandem", self.dentist_config)
        index = MultiIndex((blocks[0], blocks[-1]))

        return self.collect_job(
            name=f"tandem_alignment_block_{db.stem}",
            index=index,
            inputs=FileList(
                db=db_files(db),
                config=self.dentist_config,
            ),
            outputs=[
                self.workdir / alignment_file("TAN", db, block_b=block)
                for block in index.values()
            ],
            log=self.log_file(f"tandem-alignment.{db.stem}.{index}"),
            resources="tandem_alignment_block",
            action=lambda inputs, index, resources: ShellScript(
                ("cd", self.workdir),
                (
                    *ensure_threads_flag(aligncmd, resources["threads"]),
                    *(f"{inputs.db[0].stem}.{block}" for block in index.values()),
                ),
            ),
        )

    def mask_tandem(self, db):
        with self.grouped_jobs(f"mask_tandem_{db.stem}"):
            for blocks in get_blocks(db, self.batch_size.mask_tandem):
                self.mask_tandem_block(db, blocks=blocks)
            self.execute_jobs()
            self.catrack(
                db,
                self.tandem_mask,
                job=f"mask_tandem_{db.stem}",
                log=self.log_file(f"mask-tandem.{db.stem}"),
            )

    def mask_tandem_block(self, db, blocks):
        index = MultiIndex((blocks[0], blocks[-1]))

        return self.collect_job(
            name=f"mask_tandem_block_{db.stem}",
            index=index,
            inputs=FileList(
                db=db_files(db),
                las=[
                    self.workdir / alignment_file("TAN", db, block_b=block)
                    for block in index.values()
                ],
                config=self.dentist_config,
            ),
            outputs=[
                mask_files(db, block_mask(self.tandem_mask, block)) for block in blocks
            ],
            log=self.log_file(f"mask-tandem.{db.stem}.{index}"),
            resources="mask_tandem_block",
            action=lambda inputs, index: ShellScript(
                ("TANmask", *self.tanmask_opts, inputs.db[0], *inputs.las),
            ),
        )

    def self_alignment(self, db):
        with self.grouped_jobs(f"self_alignment_{db.stem}"):
            num_blocks = get_num_blocks(db)
            for i in range(1, num_blocks + 1):
                for blocks in get_blocks(
                    num_blocks, first=i, batch_size=self.batch_size.self_alignment
                ):
                    self.self_alignment_block(db, block_a=i, blocks_b=blocks)
            self.execute_jobs()
            self.lamerge(
                self.workdir / alignment_file(db),
                [self.workdir / las for las in block_alignments(db)],
                job=f"self_alignment_{db.stem}",
                resources="self_alignment",
                log=self.log_file(f"self-alignment.{db.stem}"),
            )

    def self_alignment_block(self, db, block_a, blocks_b):
        aligncmd = dentist.generate_options_for(
            "self",
            self.dentist_config,
            ensure_masks([], self.dust_mask, self.tandem_mask),
        )
        aligncmd = deduplicate_flags(aligncmd)
        index = MultiIndex(block_a, (blocks_b[0], blocks_b[-1]))

        return self.collect_job(
            name=f"self_alignment_block_{db.stem}",
            index=index,
            inputs=FileList(
                db=db_files(db),
                dust_mask=mask_files(db, self.dust_mask),
                tandem_mask=mask_files(db, self.tandem_mask),
                config=self.dentist_config,
            ),
            outputs=[
                (
                    self.workdir / alignment_file(db, block_a=block_a, block_b=block_b),
                    self.workdir / alignment_file(db, block_a=block_b, block_b=block_a),
                )
                for block_a, block_b in index.values()
            ],
            log=self.log_file(f"self-alignment.{db.stem}.{index}"),
            resources="self_alignment_block",
            action=lambda inputs, index, resources: ShellScript(
                ("cd", self.workdir),
                (
                    *ensure_threads_flag(aligncmd, resources["threads"]),
                    f"{inputs.db[0].stem}.{index[0]}",
                    *(
                        f"{inputs.db[0].stem}.{block_b}"
                        for block_a, block_b in index.values()
                    ),
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
            resources="mask_self",
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
        with self.grouped_jobs(f"ref_vs_reads_alignment_{refdb.stem}_{readsdb.stem}"):
            for blocks in get_blocks(readsdb, self.batch_size.ref_vs_reads_alignment):
                self.ref_vs_reads_alignment_block(refdb, readsdb, blocks)
            self.execute_jobs()
            self.lamerge(
                self.workdir / alignment_file(refdb, readsdb),
                [
                    self.workdir / las
                    for las in block_alignments(refdb, readsdb, block_a=FULL_DB)
                ],
                job=f"ref_vs_reads_alignment_{refdb.stem}_{readsdb.stem}",
                resources="ref_vs_reads_alignment",
                log=self.log_file(
                    f"ref-vs-reads-alignment.{refdb.stem}.{readsdb.stem}"
                ),
            )

    def ref_vs_reads_alignment_block(self, refdb, readsdb, blocks_reads):
        aligncmd = dentist.generate_options_for(
            "reads",
            self.dentist_config,
            ensure_masks([], self.dust_mask, self.tandem_mask, self.self_mask),
        )
        aligncmd = deduplicate_flags(aligncmd)
        index = MultiIndex((blocks_reads[0], blocks_reads[-1]))

        self.collect_job(
            name=f"ref_vs_reads_alignment_block_{refdb.stem}_{readsdb.stem}",
            index=index,
            inputs=FileList(
                refdb=db_files(refdb),
                readsdb=db_files(readsdb),
                dust_mask=mask_files(refdb, self.dust_mask),
                tandem_mask=mask_files(refdb, self.tandem_mask),
                self_mask=mask_files(refdb, self.self_mask),
                config=self.dentist_config,
            ),
            outputs=[
                (
                    self.workdir / alignment_file(refdb, readsdb, block_b=block_reads),
                    self.workdir / alignment_file(readsdb, refdb, block_a=block_reads),
                )
                for block_reads in index.values()
            ],
            log=self.log_file(
                f"ref-vs-reads-alignment.{refdb.stem}.{readsdb.stem}.{index}"
            ),
            resources="ref_vs_reads_alignment_block",
            action=lambda inputs, index, resources: ShellScript(
                ("cd", self.workdir),
                (
                    *ensure_threads_flag(aligncmd, resources["threads"]),
                    inputs.refdb[0].stem,
                    *(
                        f"{inputs.readsdb[0].stem}.{block_reads}"
                        for block_reads in index.values()
                    ),
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

    def homogenize_mask(self, refdb, readsdb, mask):
        with self.grouped_jobs(f"homogenize_mask_{refdb.stem}_{readsdb.stem}_{mask}"):
            for blocks in get_blocks(readsdb, self.batch_size.homogenize_mask):
                self.homogenize_mask_block(refdb, readsdb, mask, blocks)
            self.execute_jobs()
            self.merge_masks(
                refdb,
                homogenized_mask(mask),
                [
                    pseudo_block_mask(
                        homogenized_mask(mask), MultiIndex((blocks[0], blocks[-1]))
                    )
                    for blocks in get_blocks(readsdb, self.batch_size.homogenize_mask)
                ],
                job=f"homogenize_mask_{refdb.stem}_{readsdb.stem}_{mask}",
                log=self.log_file(f"homogenize-mask.{refdb.stem}.{readsdb.stem}"),
            )

    def homogenize_mask_block(self, refdb, readsdb, mask, blocks_reads):
        index = MultiIndex((blocks_reads[0], blocks_reads[-1]))

        self.collect_job(
            name=f"homogenize_mask_block_{refdb.stem}_{readsdb.stem}_{mask}",
            index=index,
            inputs=FileList(
                refdb=db_files(refdb),
                readsdb=db_files(readsdb),
                mask=mask_files(refdb, mask),
                las=FileList(
                    ref2reads=(
                        self.workdir
                        / alignment_file(refdb, readsdb, block_b=block_reads)
                        for block_reads in index.values()
                    ),
                    reads2ref=(
                        self.workdir
                        / alignment_file(readsdb, refdb, block_a=block_reads)
                        for block_reads in index.values()
                    ),
                ),
                config=self.dentist_config,
            ),
            outputs=[mask_files(refdb, homogenized_mask(mask), pseudo_block=index)],
            log=self.log_file(f"homogenize-mask.{refdb.stem}.{readsdb.stem}.{index}"),
            resources="homogenize_mask_block",
            action=lambda inputs, index, resources: ShellScript(
                (
                    "dentist",
                    "propagate-mask",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    f"--mask={mask}",
                    inputs.refdb[0],
                    inputs.readsdb[0],
                    *inputs.las.ref2reads,
                    pseudo_block_mask(mask, index),
                ),
                (
                    "dentist",
                    "propagate-mask",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    f"--mask={pseudo_block_mask(mask, index)}",
                    inputs.readsdb[0],
                    inputs.refdb[0],
                    *inputs.las.reads2ref,
                    pseudo_block_mask(homogenized_mask(mask), index),
                ),
            ),
        )

    def collect(self):
        self.collect_job(
            name="collect",
            inputs=FileList(
                refdb=db_files(self.reference),
                readsdb=db_files(self.reads),
                las=self.workdir / alignment_file(self.reference, self.reads),
                masks=dict(
                    (mask, mask_files(self.reference, homogenized_mask(mask)))
                    for mask in self.masks
                ),
                config=self.dentist_config,
            ),
            outputs=[self.pile_ups],
            log=self.log_file("collect"),
            action=lambda inputs, outputs, resources: ShellScript(
                (
                    "dentist",
                    "collect",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    f"--threads={dentist.main_threads(resources['threads'])}",
                    f"--auxiliary-threads={dentist.auxiliary_threads(resources['threads'])}",
                    f"--mask={','.join(inputs.masks.keys())}",
                    inputs.refdb[0],
                    inputs.readsdb[0],
                    inputs.las,
                    *outputs,
                )
            ),
        )

    def lamerge(self, merged, parts, job, log, resources=None):
        self.collect_job(
            name=job,
            inputs=parts,
            outputs=[merged],
            log=log,
            resources=resources,
            action=lambda inputs, outputs: ShellScript(
                ("LAmerge", *self.lamerge_opts, outputs[0], *inputs)
            ),
        )

    def catrack(self, db, mask, job, log):
        self.collect_job(
            name=job,
            inputs=FileList(
                db=db_files(db),
                mask=[mask_files(db, bm) for bm in block_masks(mask, db)],
            ),
            outputs={mask: mask_files(db, mask)},
            log=log,
            exec_local=True,
            action=lambda inputs, outputs: ShellScript(
                ("rm", "-f", *outputs),
                ("Catrack", "-v", inputs.db[0], *outputs.keys()),
            ),
        )

    def merge_masks(self, db, merged_mask, masks, job, log):
        self.collect_job(
            name=job,
            inputs=FileList(
                db=db_files(db),
                masks=dict((mask, mask_files(db, mask)) for mask in masks),
                config=self.dentist_config,
            ),
            outputs={merged_mask: mask_files(db, merged_mask)},
            log=log,
            exec_local=True,
            action=lambda inputs, outputs: ShellScript(
                (
                    "dentist",
                    "merge-masks",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    inputs.db[0],
                    *outputs.keys(),
                    *inputs.masks.keys(),
                )
            ),
        )

    def log_file(self, id):
        return self.logdir / f"{id}.log"


def main():
    import logging
    import sys

    config_exts = ("json", "yaml", "yml")
    override_defaults = dict(
        workflow_root="directory of <config:json>",
        resources=f"resources.({'|'.join(config_exts)}) in --workflow-root",
    )
    parser = cli_parser(log_level=True, **override_defaults)
    parser.add_argument(
        "workflow_config",
        metavar="<config:json>",
        type=Path,
        help="workflow configuration file in JSON format",
        default="workflow.json",
    )
    params = vars(parser.parse_args())
    logging.basicConfig(level=params.pop("log_level"))

    if str(params["workflow_root"]) == override_defaults["workflow_root"]:
        params["workflow_root"] = params["workflow_config"].parent

    if str(params["resources"]) == override_defaults["resources"]:
        params["resources"] = None
        for ext in config_exts:
            p = params["workflow_root"] / f"resources.{ext}"
            if p.exists():
                params["resources"] = p.name
                break

    DentistGapClosing(**params)()


if __name__ == "__main__":
    main()
