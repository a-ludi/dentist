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


[DUB]: https://code.dlang.org/download "Download DUB"

Usage
-----

Suppose we have the genome assembly `reference.fasta` that is to be updated and a set of reads `reads.fasta` with 25× coverage and a directory structure like this:

```
project
|-- input
|   |-- reference.fasta
|   `-- reads.fasta
|-- workdir
`-- result
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
dentist mask reference.dam reference.reference.las rep-self

# Align the reference to itself
# NOTE: this may require parallelization on a cluster (see `HPC.damapper`).
damapper -C -N -n.7 -e0.841500 reference.dam reads.db

# Generate a repeat mask from the reads-alignment
# NOTE: adjust the read coverage to your data!
dentist mask -C25 reference.dam reads.db reference.reads.las rep-reads

# Collect a set of candidates for gap filling (pileups)
dentist collect reference.dam reads.db reference.reads.las rep-reads pileups.db

# Generate high-quality sequences from set of candidates.
# NOTE: this may require parallelization on a cluster. Use the `show-pile-ups`
#       sub-command to get the number of pileups and process them in batches
#       using the `--batch` option of the `process` sub-command. Subsequently, 
#       merge the results using the `merge-insertions` sub-command.
dentist process reference.dam reads.db reference.reads.las rep-reads \\
        pileups.db rep-reads insertions.db

# Generate a gap-closed assembly.
dentist output reference.db insertions.db ../result/gap-closed.fasta
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
