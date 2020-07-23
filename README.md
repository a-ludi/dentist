dentist
=========

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

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
- [Maintainer](#maintainer)
- [Contributing](#contributing)
- [License](#license)


Install
--------

### Use Pre-Built Binaries

Download the latest pre-built binaries from the [releases section][release]
and extract the contents. The tarball contains a `dentist` binary as well as
the snakemake workflow, example config files and this README. In short, everything you to run DENTIST.


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

- [The Dazzler Data Base][DAZZ_DB] (>=2020-01-15)
- [`daligner`][daligner] (>=2019-07-21 && <=2020-01-15)
- [`damapper`][damapper] (>=2020-03-10)
- [`TANmask`][damasker] (>=2020-01-15)
- [`daccord`][daccord] (>=v0.0.17)

Please see their own documentation for installtion instructions. Note, the
available packages on Bioconda are outdated and should not be used at the
moment.


[DUB]: https://code.dlang.org/download "Download DUB"
[DAZZ_DB]: https://github.com/thegenemyers/DAZZ_DB
[daligner]: https://github.com/thegenemyers/DALIGNER
[damapper]: https://github.com/thegenemyers/DAMAPPER
[damasker]: https://github.com/thegenemyers/DAMASKER
[daccord]: https://gitlab.com/german.tischler/daccord

Usage
-----

Suppose we have the genome assembly `reference.fasta` that is to be updated
and a set of reads `reads.fasta` with 25× coverage.


### Quick execution with `snakemake`

Install [snakemake][snakemake] version >=5.10.0 and copy these files into your
working directory:

- `./snakemake/Snakefile`
- `./snakemake/workflow_helper.py`
- `./snakemake/snakemake.example.yml` → `./snakemake/snakemake.yml`

Next edit `snakemake.yml` to fit your needs and test your configuration with

    snakemake --configfile=snakemake.yml extend_dentist_config

If no errors occurred the whole workflow can be executed using

    snakemake --configfile=snakemake.yml

For small genomes of a few 100 Mbp this should run on a regular workstation.
Larger data sets may require a cluster in which case you can use Snakemake's
[cloud][snakemake-cloud] or [cluster][snakemake-cluster] facilities.


[snakemake]: https://snakemake.readthedocs.io/en/stable/index.html
[snakemake-cloud]: https://snakemake.readthedocs.io/en/stable/executable.html#cloud-support
[snakemake-cluster]: https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution


#### Executing on a Cluster

To make execution on a cluster easy DENTIST comes with examples files to make
Snakemake use SLURM via DRMAA. Please read the [documentation of
Snakemake][snakemake-cluster] if this does not suit your needs. Another good
starting point is [the Snakemake-Profiles project][smp-project].

Start by copying these files to your working directory:
    
- `./snakemake/profile-slurm.yml` → `~/.config/snakemake/<profile>/config.yaml`
- `./snakemake/cluster.example.yml` → `./snakemake/cluster.yml`

Next [adjust the profile][snakemake-profiles] according to your cluster. This should enable
Snakemake to submit and track jobs on your cluster. You may use the
configuration values specified in `cluster.yml` to configure job names and
resource allocation for each step of the pipeline. Now, submit the workflow
to your cluster by

    snakemake --configfile=snakemake.yml --profile=<profile>

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
section points out some important parameters and their effect on the result.


#### How to Choose DENTIST Parameters

The following list comprises the important/influential parameters for DENTIST
itself. Please keep in mind that the alignments generated by daligner/damapper
have immense influence on the performance of DENTIST.

- `--max-insertion-error`: Strong influence on quality and sensitivity. Lower
  values lead to lower sensitivity but higher quality. The maximum recommended value is `0.05`.

- `--min-anchor-length`: Higher values results in higher accuracy but lower
  sensitivity. Especially, large gaps cannot be closed if the value is too 
  high. Usually the value should be at least `500` and up to `10_000`.

- `--reference-error`, `--reads-error`: Determines the `-e` parameter for
  daligner/damapper. Use `dentist generate-dazzler-options` to see the effect
  of these parameters or consult [`commandline.d`][src-commandline] in the
  source code.

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
     assembly. Then the closed gaps must validated by other means afterwards.
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


[src-commandline]: ./blob/master/source/dentist/commandline.d
[issue-nanopore]: https://github.com/a-ludi/dentist/issues/1#issuecomment-610764625


Maintainer
----------

Arne Ludwig &lt;<arne.ludwig@posteo.de>&gt;


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

This project is licensed under MIT License (see license in [LICENSE](./LICENSE).
