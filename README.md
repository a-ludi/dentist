dentist
=========

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

> Close assembly gaps using long-reads with focus on correctness.

Today, many genome sequencing project have been conducted using second-generation sequencers which produce short reads. Such assemblies have many gaps. `dentist` closes these gaps using a (small) set of long reads. Furthermore, it can be used to scaffold contigs freely using a set of long reads. This can be used to fix known scaffolding errors or to further scaffold output of a long-read assembly pipeline.


Table of Contents
-----------------

- [Install](#install)
- [Usage](#usage)
- [Maintainer](#maintainer)
- [Contributing](#contributing)
- [License](#license)


Install
--------

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

- [The Dazzler Data Base][DAZZ_DB]
- [`daligner`][daligner]
- [`damapper`][damapper]
- [`daccord`][daccord]

Please see their own documentation for installtion instructions. Note, the
available packages on Bioconda are outdated and should not be used at the
moment.


[DUB]: https://code.dlang.org/download "Download DUB"
[DAZZ_DB]: https://github.com/thegenemyers/DAZZ_DB
[daligner]: https://github.com/thegenemyers/DALIGNER
[damapper]: https://github.com/thegenemyers/DAMAPPER
[daccord]: https://gitlab.com/german.tischler/daccord

Usage
-----

Suppose we have the genome assembly `reference.fasta` that is to be updated and a set of reads `reads.fasta` with 25× coverage.


### Quick setup with `snakemake`

Install [snakemake][snakemake] version >=5.4.0 and copy these files next to your data:

    - `./snakemake/Snakefile`
    - `./snakemake/snakemake.yml`

Next execute `snakemake`. For small genomes of a few 100 Mbp this should run
on a regular workstation. Larger data sets may require a cluster in which case
you can use Snakemake's [cloud][snakemake-cloud] or
[cluster][snakemake-cluster] facilities. The cluster config `cluster.yml` and
Snakemake profile `profile-slurm.yml` under `./snakemake` provide a starting
point for a cluster setup.


[snakemake]: https://snakemake.readthedocs.io/en/stable/index.html
[snakemake-cloud]: https://snakemake.readthedocs.io/en/stable/executable.html#cloud-support
[snakemake-cluster]: https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution


### Manual execution

Create a directory structure like this:

```
project/
|-- workdir/
|-- reference.fasta
`-- reads.fasta
```

A typical sequence of commands to run `dentist` is:

```sh
# Create DB files for fast access to sequences
fasta2DAM workdir/reference.dam input/reference.fasta
fasta2DB workdir/reads.db input/reads.fasta

# Clip sequences shorter than 1 kpb from both DBs because they may produce
# sub-optimal alignments. If smaller sequences are still desired, adjust the -x
# option accordingly.
DBsplit -x1000 workdir/reference.dam
DBsplit -x1000 workdir/reads.db

cd workdir

# List suggestions for subsequent commands
dentist generate-dazzler-options

# Align the reference to itself
# NOTE: this may require parallelization on a cluster (see `HPC.daligner`).
daligner -I -l500 -e0.980100 reference.dam reference.dam

# Generate a repeat mask from the self-alignment
dentist mask reference.dam reference.reference.las dentist-self

# Align the reference to itself
# NOTE: this may require parallelization on a cluster (see `HPC.damapper`).
damapper -C -N -n.7 -e0.841500 -mdentist-self reference.dam reads.db

# Generate a repeat mask from the reads-alignment
# NOTE: adjust the read coverage to your data!
dentist mask -C25 reference.dam reads.db reference.reads.las dentist-reads

# Collect a set of candidates for gap filling (pileups)
dentist collect -m dentist-reads -m dentist-self \
                reference.dam reads.db reference.reads.las pileups.db

# Generate high-quality sequences from set of candidates.
# NOTE: this may require parallelization on a cluster. Use the `show-pile-ups`
#       sub-command to get the number of pileups and process them in batches
#       using the `--batch` option of the `process` sub-command. Subsequently, 
#       merge the results using the `merge-insertions` sub-command.
dentist process -m dentist-reads -m dentist-self \
                reference.dam reads.db reference.reads.las pileups.db \
                insertions.db

# Generate a gap-closed assembly.
dentist output reference.db insertions.db ../gap-closed.fasta
```


Maintainer
----------

Arne Ludwig &lt;<arne.ludwig@posteo.de>&gt;


Contributing
------------

Contributions are warmly welcome. Just create an [issue][gh-issues] or [pull request][gh-pr] on GitHub. If you submit a pull request please make sure that:

- the code compiles on Linux using the current release of [dmd][dmd-download],
- your code is covered with unit tests (if feasible) and
- `dub test` runs successfully.


[gh-issues]: https://github.com/a-ludi/dentist/issues
[gh-pr]: https://github.com/a-ludi/dentist/pulls
[dmd-download]: https://dlang.org/download.html#dmd


License
-------

This project is licensed under MIT License (see license in [LICENSE](./LICENSE).
