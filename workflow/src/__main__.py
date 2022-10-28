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
from util.json import multi_json_load

log = logging.getLogger("dentist.workflow.DentistGapClosing")


class DentistGapClosing(Workflow):
    dust_mask = "dust"
    tandem_mask = "tan"
    batched_jobs = (
        "tandem_alignment",
        "mask_tandem",
        "self_alignment",
        "ref_vs_reads_alignment",
        "homogenize_mask",
        "process_pile_ups",
        "lasplit_by_reference",
        "lamerge_by_reference",
        "validate_regions",
    )
    preliminary_output_revert_options = [
        "scaffolding",
        "skip-gaps",
        "skip-gaps-file",
        "agp-dazzler",
    ] + (["cache-contig-alignments"] if dentist.is_testing() else [])
    closed_gaps_mask = "closed-gaps"
    weak_coverage_mask = "dentist-weak-coverage"

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
        self.config["__chdir__"] = workflow_config.parent

        log.info(f"changing directory to {workflow_config.parent}")
        os.chdir(workflow_config.parent)

        self.flags = self.config.get("flags", {})
        self.flags = Namespace(
            purge_output=self.flags.get("purge_output", True),
            force_validation=self.flags.get("force_validation", False),
        )

        self.reference_fasta = Path(self.config["inputs"]["reference"])
        self.reads_fasta = Path(self.config["inputs"]["reads"])
        self.workdir = Path(self.config["workdir"])
        self.logdir = Path(self.config["logdir"])
        self.reads_type = self.config["inputs"]["reads_type"]
        self.dbsplit_options = {
            "reference": self.config["reference_dbsplit"],
            "gap-closed-preliminary": self.config["reference_dbsplit"],
            "reads": self.config["reads_dbsplit"],
        }
        self.lamerge_opts = ["-v"]
        self.lasplit_opts = ["-v"]
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
        default_batch_size = self.config.get("batch_size", dict()).get("__default__", 1)
        batch_size = dict((job, default_batch_size) for job in self.batched_jobs)
        batch_size |= self.config.get("batch_size", dict())
        self.batch_size = Namespace(**batch_size)
        self.dentist_flags = shlex.split(environ.get("DENTIST_FLAGS", ""))
        self.pile_ups = self.workdir / "pile-ups.db"
        self.insertion_batches_dir = self.workdir / "insertions"
        self.insertions = self.workdir / "insertions.db"
        self.gap_closed_fasta = Path(self.config["outputs"]["output_assembly"])
        self.validation_report_block = str(
            self.workdir / "validation-report.{block}.json"
        )

    def run(self):
        self.create_dirs("create_workdirs", [self.workdir, self.logdir])
        self.execute_jobs()

        self.dentist_config = self.create_dentist_config()
        self.reference = self.fasta2dazzler(self.reference_fasta, role="reference")
        self.reads = self.fasta2dazzler(self.reads_fasta, role="reads")
        self.execute_jobs()

        dentist.validate_config(self.dentist_config)

        self.until_reads_alignment(self.reference, self.reads)
        self.execute_jobs()
        self.mask_reads()
        self.execute_jobs()
        for mask in self.masks:
            self.homogenize_mask(self.reference, self.reads, mask)
        self.execute_jobs()
        self.collect_pile_ups()
        self.execute_jobs()
        self.process_pile_ups()
        self.execute_jobs()

        if self.flags.purge_output or self.flags.force_validation:
            self.assembly_output(self.reference, self.reads, preliminary=True)
            self.execute_jobs()
            self.preliminary_gap_closed = self.fasta2dazzler(
                self.jobs["preliminary_output"].outputs.fasta, "gap-closed-preliminary"
            )
            self.execute_jobs()
            self.bed2mask(
                self.preliminary_gap_closed,
                self.jobs["preliminary_output"].outputs.bed,
                self.closed_gaps_mask,
                "closed_gaps_bed2mask",
                self.log_file("closed-gaps-bed2mask"),
            )
            self.until_reads_alignment(self.preliminary_gap_closed, self.reads)
            self.execute_jobs()
            self.split_and_merge_by_reference(self.preliminary_gap_closed, self.reads)
            self.execute_jobs()
            self.validate_regions(self.preliminary_gap_closed, self.reads)
            self.execute_jobs()
            self.skip_gaps()
            self.execute_jobs()

        if self.flags.purge_output:
            self.assembly_output(
                self.reference,
                self.reads,
                skip_gaps=self.jobs["skip_gaps"].outputs.skip_gaps,
            )
        else:
            self.assembly_output(self.reference, self.reads)

    def on_finished(self):
        job_names = []
        if self.flags.purge_output:
            job_names.append("purged_output")
        else:
            job_names.append("unpurged_output")

        if self.flags.purge_output or self.flags.force_validation:
            job_names.append(
                f"validate_regions_report_{self.preliminary_gap_closed.stem}_{self.reads.stem}"
            )
            job_names.append("skip_gaps")

        outputs = [
            str(self.config["__chdir__"] / output)
            for job_name in job_names
            for output in self.jobs[job_name].outputs
        ]
        outputs = "\n  " + "\n  ".join(outputs)

        log.info(f"Check out these result files:{outputs}")

    def until_reads_alignment(self, refdb, readsdb):
        self.mask_dust(refdb)
        self.execute_jobs()
        self.tandem_alignment(refdb)
        self.execute_jobs()
        self.mask_tandem(refdb)
        self.execute_jobs()
        self.self_alignment(refdb)
        self.execute_jobs()
        self.mask_self(refdb)
        self.execute_jobs()
        self.ref_vs_reads_alignment(refdb, readsdb)
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
            action=lambda inputs, outputs: ShellScript(
                ("rm", "-f", *outputs),
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
                [mask]
                + [
                    pseudo_block_mask(
                        homogenized_mask(mask), MultiIndex((blocks[0], blocks[-1]))
                    )
                    for blocks in get_blocks(readsdb, self.batch_size.homogenize_mask)
                ],
                job=f"homogenize_mask_{refdb.stem}_{readsdb.stem}_{mask}",
                log=self.log_file(
                    f"homogenize-mask.{refdb.stem}.{readsdb.stem}.{mask}"
                ),
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
            log=self.log_file(
                f"homogenize-mask.{refdb.stem}.{readsdb.stem}.{mask}.{index}"
            ),
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

    def collect_pile_ups(self):
        self.collect_job(
            name="collect_pile_ups",
            inputs=FileList(
                refdb=db_files(self.reference),
                readsdb=db_files(self.reads),
                las=self.workdir / alignment_file(self.reference, self.reads),
                masks=dict(
                    (
                        homogenized_mask(mask),
                        mask_files(self.reference, homogenized_mask(mask)),
                    )
                    for mask in self.masks
                ),
                config=self.dentist_config,
            ),
            outputs=[self.pile_ups],
            log=self.log_file("collect_pile_ups"),
            action=lambda inputs, outputs, resources: ShellScript(
                (
                    "dentist",
                    "collect-pile-ups",
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

    def process_pile_ups(self):
        self.create_dirs("create_insertion_batches_dir", self.insertion_batches_dir)
        self.execute_jobs()

        with self.grouped_jobs("process_pile_ups"):
            num_pile_ups = dentist.get_num_pile_ups(self.pile_ups)
            for batch in dentist.batch_ranges(
                num_pile_ups, self.batch_size.process_pile_ups
            ):
                self.process_pile_ups_batch(batch)
            self.execute_jobs()
            self.merge_insertions()

    def process_pile_ups_batch(self, batch):
        index = MultiIndex((batch[0], batch[-1] + 1))
        self.collect_job(
            name="process_pile_ups",
            index=index,
            inputs=FileList(
                refdb=db_files(self.reference),
                readsdb=db_files(self.reads),
                pile_ups=self.pile_ups,
                masks=dict(
                    (
                        homogenized_mask(mask),
                        mask_files(self.reference, homogenized_mask(mask)),
                    )
                    for mask in self.masks
                ),
                config=self.dentist_config,
            ),
            outputs=[self.insertion_batches_dir / f"batch.{index}.db"],
            log=self.log_file(f"process_pile_ups.{index}"),
            action=lambda inputs, outputs, index, resources: ShellScript(
                (
                    "dentist",
                    "process-pile-ups",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    f"--threads={dentist.main_threads(resources['threads'])}",
                    f"--auxiliary-threads={dentist.auxiliary_threads(resources['threads'])}",
                    f"--mask={','.join(inputs.masks.keys())}",
                    f"--batch={index.to_str(range_sep='..')}",
                    inputs.refdb[0],
                    inputs.readsdb[0],
                    inputs.pile_ups,
                    *outputs,
                )
            ),
        )

    def merge_insertions(self):
        merge_config = self.workdir / "dentist.merge.json"
        insertion_batches = list(
            chain.from_iterable(
                job.outputs for job in self.jobs["process_pile_ups"].values()
            )
        )

        @self.collect_job(
            inputs=insertion_batches,
            outputs=[merge_config],
            log=self.log_file("merge-insertions"),
            exec_local=True,
        )
        @python_code
        def make_merge_config(inputs, outputs):
            dentist_config = self.config.get("dentist", {})
            merge_config = {
                "__default__": dentist_config.get("__default__", {}),
                "merge-insertions": dentist_config.get("merge-insertions", {}),
            }
            for cmd in ("__default__", "merge-insertions"):
                for arg in ("insertions", "partitioned-insertions"):
                    if arg in merge_config[cmd]:
                        del merge_config[cmd][arg]
            merge_config["merge-insertions"]["partitioned-insertions"] = list(
                str(path) for path in inputs
            )

            with outputs[0].open("w") as merge_config_file:
                json.dump(merge_config, merge_config_file)

        self.execute_jobs()
        self.collect_job(
            name="merge_insertions",
            inputs=FileList(
                config=merge_config,
                insertion_batches=insertion_batches,
            ),
            outputs=[self.insertions],
            log=self.log_file("merge-insertions"),
            action=lambda inputs, outputs: ShellScript(
                (
                    "dentist",
                    "merge-insertions",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    *outputs,
                    "-",
                )
            ),
        )

    def split_and_merge_by_reference(self, refdb, readsdb):
        with self.grouped_jobs(
            f"split_and_merge_by_reference_{refdb.stem}_{readsdb.stem}"
        ):
            for blocks_b in get_blocks(readsdb, self.batch_size.lasplit_by_reference):
                self.lasplit(
                    refdb,
                    readsdb,
                    blocks_b,
                    job_stem="lasplit_by_reference",
                    resources="lasplit_by_reference",
                )
            self.execute_jobs()
            for blocks_a in get_blocks(refdb, self.batch_size.lamerge_by_reference):
                self.lamerge_by_reference(refdb, readsdb, blocks_a)

    def validate_regions(self, refdb, readsdb):
        with self.grouped_jobs(f"validate_regions_{refdb.stem}_{readsdb.stem}"):
            block_jobs = [
                self.validate_regions_block(refdb, readsdb, blocks)
                for blocks in get_blocks(refdb, self.batch_size.validate_regions)
            ]
            weak_coverage_masks = list(
                chain.from_iterable(job.outputs.masks.keys() for job in block_jobs)
            )
            reports = list(
                chain.from_iterable(job.outputs.reports for job in block_jobs)
            )
            self.execute_jobs()
            self.merge_masks(
                refdb,
                self.weak_coverage_mask,
                weak_coverage_masks,
                job=f"weak_coverage_mask_{refdb.stem}_{readsdb.stem}",
                log=self.log_file(f"weak-coverage-mask.{refdb.stem}.{readsdb.stem}"),
            )
            self.concatenate(
                self.workdir / f"validation-report.json",
                reports,
                job=f"validate_regions_report_{refdb.stem}_{readsdb.stem}",
                resources="validate_regions_report",
                log=self.log_file(
                    f"validate-regions-report.{refdb.stem}.{readsdb.stem}"
                ),
            )

    def validate_regions_block(self, refdb, readsdb, blocks):
        index = MultiIndex((blocks[0], blocks[-1]))
        weak_coverage_masks = [
            pseudo_block_mask(self.weak_coverage_mask, block)
            for block in index.values()
        ]

        return self.collect_job(
            name=f"validate_regions_{refdb.stem}_{readsdb.stem}",
            index=index,
            inputs=FileList(
                refdb=db_files(refdb),
                readsdb=db_files(readsdb),
                las=[
                    self.workdir / alignment_file(refdb, readsdb, block_a=block)
                    for block in index.values()
                ],
                mask={self.closed_gaps_mask: mask_files(refdb, self.closed_gaps_mask)},
                config=self.dentist_config,
            ),
            outputs=FileList(
                reports=[
                    self.validation_report_block.format(block=block)
                    for block in index.values()
                ],
                masks=dict(
                    (mask, mask_files(refdb, mask)) for mask in weak_coverage_masks
                ),
            ),
            log=self.log_file(f"validate-regions.{index}"),
            action=lambda inputs, outputs, resources: ShellScript(
                *(
                    (
                        "dentist",
                        "validate-regions",
                        f"--config={inputs.config}",
                        f"--threads={resources['threads']}",
                        f"--weak-coverage-mask={weak_coverage_mask}",
                        inputs.refdb[0],
                        inputs.readsdb[0],
                        las,
                        *inputs.mask.keys(),
                        safe(">"),
                        report,
                    )
                    for las, report, weak_coverage_mask in zip(
                        inputs.las, outputs.reports, outputs.masks.keys()
                    )
                )
            ),
        )

    def skip_gaps(self):
        @self.collect_job(
            inputs=FileList(report=self.workdir / f"validation-report.json"),
            outputs=FileList(skip_gaps=self.workdir / f"skip-gaps.txt"),
            log=self.log_file("skip-gaps"),
            exec_local=True,
        )
        @python_code
        def skip_gaps(inputs, outputs):
            validations = None
            with inputs.report.open("r") as validation_report_file:
                validations = multi_json_load(validation_report_file)

            with outputs.skip_gaps.open("w") as skip_gaps_file:
                for validation in validations:
                    if not validation.get("isValid", False):
                        gap_spec = "-".join(
                            (str(cid) for cid in validation["contigIds"])
                        )
                        print(gap_spec, file=skip_gaps_file)

        return skip_gaps

    def lasplit(self, dba, dbb, blocks_b, job_stem="lasplit", resources=None):
        index = MultiIndex((blocks_b[0], blocks_b[-1]))
        las_batch = [
            self.workdir / alignment_file(dba, dbb, block_b=block_b)
            for block_b in index.values()
        ]
        split_files = dict(
            (
                str(
                    self.workdir
                    / alignment_file(dba, dbb, block_a="@", block_b=block_b)
                ),
                [
                    self.workdir
                    / alignment_file(dba, dbb, block_a=block_a, block_b=block_b)
                    for block_a in range(1, get_num_blocks(dba) + 1)
                ],
            )
            for block_b in index.values()
        )

        self.collect_job(
            name=f"{job_stem}_{dba.stem}_{dbb.stem}",
            index=index,
            inputs=FileList(
                las_batch=las_batch,
                db=db_files(dba),
            ),
            outputs=split_files,
            log=self.log_file(
                f"{job_stem.replace('_', '-')}.{dba.stem}.{dbb.stem}.{index}"
            ),
            resources=resources,
            action=lambda inputs, outputs: ShellScript(
                *(
                    (
                        "LAsplit",
                        *self.lasplit_opts,
                        target,
                        inputs.db[0],
                        safe("<"),
                        las,
                    )
                    for las, target in zip(inputs.las_batch, outputs.keys())
                ),
                ("file", *outputs),
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

    def lamerge_by_reference(self, refdb, readsdb, batch):
        index = MultiIndex((batch[0], batch[-1]))

        self.collect_job(
            name=f"lamerge_by_reference_{refdb.stem}_{readsdb.stem}",
            index=index,
            inputs=[
                [
                    self.workdir / las
                    for las in block_alignments(refdb, readsdb, block_a=block)
                ]
                for block in index.values()
            ],
            outputs=[
                self.workdir / alignment_file(refdb, readsdb, block_a=block)
                for block in index.values()
            ],
            resources="lamerge_by_reference",
            log=self.log_file(f"lamerge-by-reference.{refdb.stem}.{readsdb.stem}"),
            action=lambda inputs, outputs: ShellScript(
                *(
                    ("LAmerge", *self.lamerge_opts, outputs[0], *inputs)
                    for las_inputs, las_output in zip(inputs, outputs)
                )
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

    def bed2mask(self, db, bed, mask, job, log):
        self.collect_job(
            name=job,
            inputs=FileList(
                db=db_files(db),
                bed=bed,
                config=self.dentist_config,
            ),
            outputs={mask: mask_files(db, mask)},
            log=log,
            exec_local=True,
            action=lambda inputs, outputs: ShellScript(
                (
                    "dentist",
                    "bed2mask",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    "--data-comments",
                    f"--bed={inputs.bed}",
                    inputs.db[0],
                    *outputs.keys(),
                )
            ),
        )

    def concatenate(self, output, inputs, job, resources, log):
        @python_code
        def concatenate(inputs, outputs):
            with outputs[0].open("w") as outfile:
                for input_ in inputs:
                    with input_.open() as infile:
                        outfile.write(infile.read())

        return self.collect_job(
            name=job,
            inputs=inputs,
            outputs=[output],
            resources=resources,
            log=log,
            exec_local=True,
            action=concatenate,
        )

    def assembly_output(self, refdb, readsdb, skip_gaps=None, preliminary=False):
        if preliminary:
            job_name = "preliminary_output"
            options = (
                f"--revert={','.join(self.preliminary_output_revert_options)}",
                "--agp-dazzler",
            )
            gap_closed_fasta = self.workdir / "gap-closed-preliminary.fasta"
        else:
            job_name = "unpurged_output" if skip_gaps is None else "purged_output"
            options = tuple()
            gap_closed_fasta = self.gap_closed_fasta

        self.collect_job(
            name=job_name,
            inputs=FileList(
                refdb=db_files(refdb),
                readsdb=db_files(readsdb),
                insertions=self.insertions,
                config=self.dentist_config,
                skip_gaps=[] if skip_gaps is None else skip_gaps,
            ),
            outputs=FileList(
                fasta=gap_closed_fasta,
                bed=gap_closed_fasta.with_suffix(f".{self.closed_gaps_mask}.bed"),
                agp=gap_closed_fasta.with_suffix(f".agp"),
            ),
            log=self.log_file(job_name.replace("_", "-")),
            action=lambda inputs, outputs: ShellScript(
                (
                    "dentist",
                    "output",
                    f"--config={inputs.config}",
                    *self.dentist_flags,
                    *options,
                    f"--closed-gaps-bed={outputs.bed}",
                    f"--agp={outputs.agp}",
                    inputs.refdb[0],
                    inputs.readsdb[0],
                    inputs.insertions,
                    outputs.fasta,
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
