DENTIST
=======

![DENTIST](./docs/logo.png?sanitize=true&raw=true)

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat)](https://github.com/RichardLitt/standard-readme)
[![DUB](https://img.shields.io/dub/v/dentist)](https://code.dlang.org/packages/dentist)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/aludi/dentist?logo=docker&sort=semver)](https://hub.docker.com/repository/docker/aludi/dentist)

> Close assembly gaps using long-reads with focus on correctness.

Today, many genome sequencing project have been conducted using
second-generation sequencers which produce short reads. Such assemblies have
many gaps. `dentist` closes these gaps using a (small) set of long reads.
Furthermore, it can be used to scaffold contigs freely using a set of long
reads. This can be used to fix known scaffolding errors or to further scaffold
output of a long-read assembly pipeline.


Table of Contents
-----------------

- [Install](#install)
- [Usage](#usage)
- [Configuration](#configuration)
- [Citation](#citation)
- [Maintainer](#maintainer)
- [Contributing](#contributing)
- [License](#license)


Install
--------

### Use a Singularity Container (recommended)

Make sure [Singularity][singularity] is installed on your system. You can then use the container like so:

```sh
# launch an interactive shell
singularity shell docker://aludi/dentist:latest

# execute a single command inside the container
singularity exec docker://aludi/dentist:latest dentist --version

# run the whole workflow on a cluster using Singularity
snakemake --configfile=snakemake.yml --use-singularity --profile=slurm
```

The last command is explained in more detail below in
[the usage section](#usage).


[singularity]: https://sylabs.io/guides/3.7/user-guide/index.html


### Use Pre-Built Binaries

Download the latest pre-built binaries from the [releases section][release]
and extract the contents. The tarball contains a `dentist` binary as well as
the Snakemake workflow, example config files and this README. In short, everything you to run DENTIST.


[release]: https://github.com/a-ludi/dentist/releases


### Build from Source

Be sure to install the D package manager [DUB][DUB]. Install using either

```sh
dub install dentist
```

or

```sh
git clone https://github.com/a-ludi/dentist.git
cd dentist
dub build
```

### Runtime Dependencies

The following software packages are required to run `dentist`:

- [The Dazzler Data Base][DAZZ_DB] (>=2020-07-27)
  > Manage sequences (reads and assemblies) in 4bit encoding alongside 
  > auxiliary information such as masks or QV tracks
- [DALIGNER][daligner] (=2020-01-15)
  > Find significant local alignments.
- [DAMAPPER][damapper] (>=2020-03-10)
  > Find alignment chains, i.e. sequences of significant local alignments
  > possibly with unaligned gaps.
- [DAMASKER][damasker] (>=2020-01-15)
  > Discover tandem repeats.
- [DASCRUBBER][dascrubber] (>=2020-07-26)
  > Estimate coverage and compute QVs.
- [daccord][daccord] (>=v0.0.17)
  > Compute reference-based consensus sequence for gap filling.

Please see their own documentation for installation instructions. Note, the
available packages on Bioconda are outdated and should not be used at the
moment.

Please use the following versions in your dependencies in case you experience
troubles. These should be the same versions used in the [Dockerfile](./blob/develop/Dockerfile):

- [DENTIST@1.0.0](https://github.com/a-ludi/dentist/tree/v1.0.0)
- [snakemake@5.32.1](https://snakemake.readthedocs.io/en/v5.32.1/getting_started/installation.html)
- [DAZZ_DB@d22ae58](https://github.com/thegenemyers/DAZZ_DB/tree/d22ae58d32a663d09325699f17373ccf8c6f93a0)
- [DALIGNER@c2b47da](https://github.com/thegenemyers/DALIGNER/tree/c2b47da6b3c94ed248a6be395c5b96a4e63b3f63)
- [DAMAPPER@b2c9d7f](https://github.com/thegenemyers/DAMAPPER/tree/b2c9d7fd64bb4dd2dde7c69ff3cc8a04cbeeebbc)
- [DAMASKER@22139ff](https://github.com/thegenemyers/DAMASKER/tree/22139ff1c2b2c0ff2589fbc9cc948370be799827)
- [DASCRUBBER@a53dbe8](https://github.com/thegenemyers/DASCRUBBER/tree/a53dbe879a716e7b08338f397de5a0403637641e)
- [daccord@0.0.17](https://gitlab.com/german.tischler/daccord/tree/d54b10bb863b14103cf8e03c07efe4f93c5772d8)


[DUB]: https://code.dlang.org/download "Download DUB"
[DAZZ_DB]: https://github.com/thegenemyers/DAZZ_DB
[daligner]: https://github.com/thegenemyers/DALIGNER
[damapper]: https://github.com/thegenemyers/DAMAPPER
[damasker]: https://github.com/thegenemyers/DAMASKER
[dascrubber]: https://github.com/thegenemyers/DASCRUBBER
[daccord]: https://gitlab.com/german.tischler/daccord

Usage
-----

Suppose we have the genome assembly `reference.fasta` that is to be updated
and a set of reads `reads.fasta` with 25× coverage.


### Quick execution with Snakemake (and Singularity)

> TL;DR
>
>     # edit dentist.json and snakemake.yml
>     snakemake --configfile=snakemake.yml --use-singularity --profile=slurm

Install [Snakemake][snakemake] version >=5.32.1 and copy these files into your
working directory:

- `./snakemake/Snakefile`
- `./snakemake/snakemake.yml`
- `./snakemake/dentist.json`

Next edit `snakemake.yml` and `dentist.json` to fit your needs and optionally
test your configuration with

    snakemake --configfile=snakemake.yml --use-singularity --cores=1 -f -- validate_dentist_config

If no errors occurred the whole workflow can be executed using

    snakemake --configfile=snakemake.yml --use-singularity --cores=all

For small genomes of a few 100 Mbp this should run on a regular workstation.
One may use Snakemake's `--jobs` to run independent jobs in parallel. Larger
data sets may require a cluster in which case you can use Snakemake's
[cloud][snakemake-cloud] or [cluster][snakemake-cluster] facilities.


[snakemake]: https://snakemake.readthedocs.io/en/stable/index.html
[snakemake-cloud]: https://snakemake.readthedocs.io/en/stable/executable.html#cloud-support
[snakemake-cluster]: https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution


#### Executing on a Cluster

To make execution on a cluster easy DENTIST comes with examples files to make
Snakemake use SLURM via DRMAA. Please read the [documentation of
Snakemake][snakemake-cluster] if this does not suit your needs. Another good
starting point is [the Snakemake-Profiles project][smp-project].

Start by copying these files to your working/home directory:
    
- `./snakemake/Snakefile`
- `./snakemake/snakemake.yml`
- `./snakemake/cluster.yml`
- `./snakemake/profile-slurm.yml` → `~/.config/snakemake/slurm/config.yaml`

Next [adjust the profile][snakemake-profiles] according to your cluster. This
should enable Snakemake to submit and track jobs on your cluster. You may use
the configuration values specified in `cluster.yml` to configure job names and
resource allocation for each step of the pipeline. Now, submit the workflow
to your cluster by

    snakemake --configfile=snakemake.yml --profile=slurm --use-singularity

Note, parameters specified in the profile provide default values and can be
overridden by specififying different value on the CLI.


[smp-project]: https://github.com/snakemake-profiles/doc
[snakemake-profiles]: https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles


### Manual execution

Please inspect the Snakemake workflow to get all the details. It might be
useful to execute Snakemake with the `-p` switch which causes Snakemake to
print the shell commands. If you plan to write your own workflow management
for DENTIST please feel free to contact the maintainer!


Configuration
-------------

DENTIST comprises a complex pipeline of with many options for tweaking. This
section points out some important parameters and their effect on the result or
performance.


#### How to Choose DENTIST Parameters

The following list comprises the important/influential parameters for DENTIST
itself. Please keep in mind that the alignments generated by daligner/damapper
have immense influence on the performance of DENTIST.

- `--max-insertion-error`: Strong influence on quality and sensitivity. Lower
  values lead to lower sensitivity but higher quality. The maximum recommended value is `0.05`.

- `--min-anchor-length`: Higher values results in higher accuracy but lower
  sensitivity. Especially, large gaps cannot be closed if the value is too 
  high. Usually the value should be at least `500` and up to `10_000`.

- `--min-reads-per-pile-up`: Choosing higher values for the minimum number of
  reads drastically reduces sensitivity but has little effect on the quality.
  Small values may be chosen to get the maximum sensitivity in _de novo_
  assemblies. Make sure to throughly validate the results though. 

- `--min-spanning-reads`:  Higher values give more confidence on the
  correctness of closed gaps but reduce sensitivity. The value must be well
  below the expected coverage.

- `--allow-single-reads`: May be used under careful consideration. This is
  intended for one of the following scenarios:

  1. DENTIST is meant to close as many gaps as possible in a _de novo_
     assembly. Then the closed gaps must be validated by other means
     afterwards.
  2. DENTIST is used not with real reads but with an independent assembly.

- `--existing-gap-bonus`: If DENTIST finds evidence to join two contigs that
  are already consecutive in the input assembly (i.e. joined by `N`s) then it
  will preferred over conflicting joins (if present) with this bonus. The
  default value is rather conservative, i.e. the preferred join almost always
  wins over other joins in case of a conflict.

- `--join-policy`: Choose according to your needs:
  
  - `scaffoldGaps`: Closes only gaps that are marked by `N`s in the assembly.
    This is the default mode of operation. Use this if you do not want to alter
    the scaffolding of the assembly. See also `--existing-gap-bonus`.
  - `scaffolds`: Allows whole scaffolds to be joined in addition to the effects
    of `scaffoldGaps`. Use this if you have (many) scaffolds that are not yet
    full chromosome-scale.
  - `contigs`: Allows contigs to be rearranged freely. This is especially
    useful in _de novo_ assemblies **before** applying any other scaffolding
    methods as it increases the contiguity thus increasing the chance that
    large-scale scaffolding (e.g. Bionano or Hi-C) finds proper joins.


#### Choosing the Read Type

In the examples PacBio long reads are assumed but DENTIST can be run using any
kind of long reads. Currently, this is either PacBio or Oxford Nanopore reads.
For using none-PacBio reads, the `reads_type` in `snakemake.yml` must be set
to anything other than `PACBIO_SMRT`. The recommendation is to use
`OXFORD_NANOPORE` for Oxford Nanopore. These names are borrowed from the NCBI.
Further details on the rationale can found in [this issue][issue-nanopore].


[issue-nanopore]: https://github.com/a-ludi/dentist/issues/1#issuecomment-610764625


#### Cluster/Cloud Execution

Cluster job schedulers can become unresponsive or even crash if too many jobs
with short running time are submitted to the cluster. It is therefore
advisable to adjust the workflow accordingly. We tried to provide a default
configuration that works in most cases as is but the application scenarios can
be very diverse and manual adjustments may become necessary. Here is a small
guide which config parameters influence the number of jobs and how much
resources they consume.

- `max_threads`: Sets the maximum number of threads/cores a single job may
  use. A single-threaded job will always allocate a single core but
  thread-parallel steps, e.g. the sequence alignments, will use up to
  `max_threads` if snakemake has been provided enough cores via `--cores`.
- `-s<block_size:uint>`: The assembly and reads FAST/A files are converted into
  Dazzler DBs. These DBs store the sequence in a 2-bit encoding and have
  additional features like tracks (similar to BED files). Also they are split
  into blocks of `<block_size>`Mb. Alignments are calculated on the basis of
  these blocks which enables easy distribution onto the cluster. The larger the
  block size the longer are the alignment jobs and the more memory they require
  but also the number of jobs is reduced. Experience shows that the block size
  should be between 200Mb and 500Mb.
- `propagate_batch_size`: The repeat masks are homogenized by propagating them
  from the assembly to the reads and back again. Usually these jobs are very
  short because the propagation is parallelized over the blocks of the reads
  DB. To reduce the number of jobs both propagation directions are grouped
  together and submitted in batches of `propagate_batch_size` read blocks.
  Increasing `propagate_batch_size` reduces the number of submitted jobs and
  increases the run time per job. It has no effect on the memory requirements.
- `batch_size`: In the `collect` step DENTIST identifies candidates for gap
  closing each consisting of a pile up of reads. From these pile ups
  consensus sequences are computed and validated in the `process` step. Each
  job process `batch_size` pile ups. Increasing `batch_size` reduces the
  number of submitted jobs and increases the run time per job. It has no
  effect on the memory requirements.
- `validation_blocks`: The preliminarily closed gaps are validated by analyzing
  how the reads align to each closed gap. The validation is conducted in
  independent jobs for `validation_blocks` many blocks of the gap-closed
  assembly. Decreasing `validation_blocks` reduces the number of submitted
  jobs and increases the run time and memory requirements per job. The memory requirement is proportional to the size of the read alignment blocks. 


Citation
--------

> Arne Ludwig, Martin Pippel, Gene Myers, Michael Hiller. DENTIST – close
> assembly gaps with high confidence. _In preparation._


Maintainer
----------

DENTIST is being developed by Arne Ludwig &lt;<ludwig@mpi-cbg.de>&gt; at
the Planck Institute of Molecular Cell Biology and Genetics, Dresden, Germany.


Contributing
------------

Contributions are warmly welcome. Just create an [issue][gh-issues] or [pull request][gh-pr] on GitHub. If you submit a pull request please make sure that:

- the code compiles on Linux using the current release of [dmd][dmd-download],
- your code is covered with unit tests (if feasible) and
- `dub test` runs successfully.

It is recommended to install the Git hooks included in the repository to avoid premature pull requests. You can enable all shipped hooks with this command:

```sh
git config --local core.hooksPath .githooks/
```

If you do not want to enable just a subset use `ln -s .githooks/{hook} .git/hooks`. If you want to audit code changes before they get executed on your machine you can you `cp .githooks/{hook} .git/hooks` instead.


[gh-issues]: https://github.com/a-ludi/dentist/issues
[gh-pr]: https://github.com/a-ludi/dentist/pulls
[dmd-download]: https://dlang.org/download.html#dmd


License
-------

This project is licensed under MIT License (see [LICENSE](./LICENSE)).
