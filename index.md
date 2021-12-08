DENTIST
=======

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat)](https://github.com/RichardLitt/standard-readme)
![GitHub](https://img.shields.io/github/license/a-ludi/dentist)
[![DUB](https://img.shields.io/dub/v/dentist)](https://code.dlang.org/packages/dentist)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/aludi/dentist?logo=docker&sort=semver)](https://hub.docker.com/repository/docker/aludi/dentist)
[![DOI 10.1101/2021.02.26.432990](https://img.shields.io/badge/DOI-(pre--print)_10.1101%2F2021.02.26.432990-informational)](https://doi.org/10.1101/2021.02.26.432990)

Long sequencing reads allow increasing contiguity and completeness of
fragmented, short-read based genome assemblies by closing assembly gaps,
ideally at high accuracy. DENTIST is a sensitive, highly-accurate and
automated pipeline method to close gaps in (short read) assemblies with long
reads.

**API documentation:** ([current][api-current], [v2.0.0][api-v2.0.0])

**First time here? Head over to [the example](#example) and make sure it works.**


[api-current]: ./api/current
[api-v2.0.0]: ./api/v2.0.0


Install
-------

### Use Conda (recommended)

Make sure [Conda][conda] is installed on your system. You can then use DENTIST like so:

```sh
# run the whole workflow on a cluster using Singularity
snakemake --configfile=snakemake.yml --use-conda -jall
snakemake --configfile=snakemake.yml --use-conda --profile=slurm
```

The last command is explained in more detail below in
[the usage section](#usage).


[conda]: https://docs.conda.io/projects/conda/


### Use a Singularity Container (recommended)

Make sure [Singularity][singularity] is installed on your system. You can then use the container like so:

```sh
# launch an interactive shell
singularity shell docker://aludi/dentist:stable

# execute a single command inside the container
singularity exec docker://aludi/dentist:stable dentist --version

# run the whole workflow on a cluster using Singularity
snakemake --configfile=snakemake.yml --use-singularity --profile=slurm
```

The last command is explained in more detail below in
[the usage section](#usage).


[singularity]: https://sylabs.io/guides/3.5/user-guide/index.html


### Use Pre-Built Binaries

Download the latest pre-built binaries from the [releases section][release]
and extract the contents. The pre-built binaries are stored in a subfolder
called `bin`. Here are the instructions for `v2.0.0`:

```sh
# download & extract pre-built binaries
wget https://github.com/a-ludi/dentist/releases/download/v2.0.0/dentist.v2.0.0.x86_64.tar.gz
tar -xzf dentist.v2.0.0.x86_64.tar.gz

# make binaries available to your shell
cd dentist.v2.0.0.x86_64
PATH="$PWD/bin:$PATH"

# check installation with
dentist -d
# Expected output:
# 
#daligner (part of `DALIGNER`; see https://github.com/thegenemyers/DALIGNER) [OK]
#damapper (part of `DAMAPPER`; see https://github.com/thegenemyers/DAMAPPER) [OK]
#DAScover (part of `DASCRUBBER`; see https://github.com/thegenemyers/DASCRUBBER) [OK]
#DASqv (part of `DASCRUBBER`; see https://github.com/thegenemyers/DASCRUBBER) [OK]
#DBdump (part of `DAZZ_DB`; see https://github.com/thegenemyers/DAZZ_DB) [OK]
#DBdust (part of `DAZZ_DB`; see https://github.com/thegenemyers/DAZZ_DB) [OK]
#DBrm (part of `DAZZ_DB`; see https://github.com/thegenemyers/DAZZ_DB) [OK]
#DBshow (part of `DAZZ_DB`; see https://github.com/thegenemyers/DAZZ_DB) [OK]
#DBsplit (part of `DAZZ_DB`; see https://github.com/thegenemyers/DAZZ_DB) [OK]
#fasta2DAM (part of `DAZZ_DB`; see https://github.com/thegenemyers/DAZZ_DB) [OK]
#fasta2DB (part of `DAZZ_DB`; see https://github.com/thegenemyers/DAZZ_DB) [OK]
#computeintrinsicqv (part of `daccord`; see https://gitlab.com/german.tischler/daccord) [OK]
#daccord (part of `daccord`; see https://gitlab.com/german.tischler/daccord) [OK]
```

The tarball additionally contains the Snakemake workflow, example config files
and this README. In short, everything you to run DENTIST.


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
troubles. These should be the same versions used in the
[Conda recipe](./conda/recipes/dentist-core/meta.yaml):

- [snakemake@5.32.1](https://snakemake.readthedocs.io/en/v5.32.1/getting_started/installation.html)
- [DENTIST@2.0.0](https://github.com/a-ludi/dentist/tree/v2.0.0)
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

Before you start producing wonderful scientific results, you should skip over
to the [example section](#example) and try to run the small example. This will
make sure your setup is working as expected.

### Quick execution with Snakemake (and Singularity or Conda)

> TL;DR
>
>     # edit dentist.yml and snakemake.yml
>     snakemake --configfile=snakemake.yml --use-conda --profile=slurm

Install [Snakemake][snakemake] version >=5.32.1 and copy these files into your
working directory:

```sh
cp -r -t .
    ./snakemake/dentist.yml \
    ./snakemake/Snakefile \
    ./snakemake/snakemake.yml \
    ./snakemake/envs \
    ./snakemake/scripts
```

Next edit `snakemake.yml` and `dentist.yml` to fit your needs and optionally
test your configuration with

    snakemake --configfile=snakemake.yml --use-conda --cores=1 -f -- validate_dentist_config

If no errors occurred the whole workflow can be executed using

    snakemake --configfile=snakemake.yml --use-conda --cores=all

For small genomes of a few 100 Mbp this should run on a regular workstation.
One may use Snakemake's `--cores` to run independent jobs in parallel. Larger
data sets may require a cluster in which case you can use Snakemake's
[cloud][snakemake-cloud] or [cluster][snakemake-cluster] facilities.


[snakemake]: https://snakemake.readthedocs.io/en/stable/index.html
[snakemake-cloud]: https://snakemake.readthedocs.io/en/stable/executing/cloud.html
[snakemake-cluster]: https://snakemake.readthedocs.io/en/stable/executing/cluster.html


#### Executing on a Cluster

To make execution on a cluster easy DENTIST comes with examples files to make
Snakemake use SLURM via DRMAA. Please read the [documentation of
Snakemake][snakemake-cluster] if this does not suit your needs. Another good
starting point is [the Snakemake-Profiles project][smp-project].

Start by copying these files to your working/home directory:

```sh
cp -r -t .
    ./snakemake/cluster.yml \
    ./snakemake/dentist.yml \
    ./snakemake/Snakefile \
    ./snakemake/snakemake.yml \
    ./snakemake/envs \
    ./snakemake/scripts
mkdir -p ~/.config/snakemake/slurm
# choose appropriate file from `snakemake/profile-slurm.*.yml`
cp ./snakemake/profile-slurm.submit-async.yml ~/.config/snakemake/slurm
```

Next [adjust the profile][snakemake-profiles] according to your cluster. This
should enable Snakemake to submit and track jobs on your cluster. You may use
the configuration values specified in `cluster.yml` to configure job names and
resource allocation for each step of the pipeline. Now, submit the workflow
to your cluster by

    snakemake --configfile=snakemake.yml --profile=slurm --use-conda

Note, parameters specified in the profile provide default values and can be
overridden by specifying different value on the CLI.


[smp-project]: https://github.com/snakemake-profiles/doc
[snakemake-profiles]: https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles


### Manual execution

Please inspect the Snakemake workflow to get all the details. It might be
useful to execute Snakemake with the `-p` switch which causes Snakemake to
print the shell commands. If you plan to write your own workflow management
for DENTIST please feel free to contact the maintainer!


Example
-------

Make sure you have [Snakemake][snakemake] 5.32.1 or later installed.

You can also use the convenient Singularity container to execute the rules.
Just make sure you have [Singularity][singularity] 3.5.x or later installed.

First of all download the test data and workflow and switch to the
`dentist-example` directory.

```sh
wget https://github.com/a-ludi/dentist-example/releases/download/v2.0.0-3/dentist-example.tar.gz
tar -xzf dentist-example.tar.gz
cd dentist-example
```


### Local Execution

Execute the entire workflow on your *local machine* using `all` cores:

```sh
# run the workflow
PATH="$PWD/bin:$PATH" snakemake --configfile=snakemake.yml --cores=all

# validate the files
md5sum -c checksum.md5
```

Execution takes approx. 7 minutes and a maximum of 1.7GB memory on my little
laptop with an Intel® Core™ i5-5200U CPU @ 2.20GHz.


### Execution with Conda

Execute the workflow inside a convenient Singularity image by adding `--use-conda` to the call to Snakemake:

```sh
# run the workflow
snakemake --configfile=snakemake.yml --use-conda --cores=all

# validate the files
md5sum -c checksum.md5
```

In more recent versions of Snakemake, you may need to also pass
`--conda-frontend=conda` unless you have [Mamba][mamba] installed. Mamba is a
faster alternative to Conda.


[mamba]: https://github.com/mamba-org/mamba


### Execution in Singularity Container

Execute the workflow inside a convenient Singularity image by adding `--use-singularity` to the call to Snakemake:

```sh
# run the workflow
snakemake --configfile=snakemake.yml --use-singularity --cores=all

# validate the files
md5sum -c checksum.md5
```


### Cluster Execution

Execute the workflow on a *SLURM cluster*:

```sh
mkdir -p "$HOME/.config/snakemake/slurm"
# select one of the profile-slurm.{drmaa,submit-async,submit-sync}.yml files
cp -v "profile-slurm.submit-async.yml" "$HOME/.config/snakemake/slurm/config.yaml"
# execute using the cluster profile
snakemake --configfile=snakemake.yml --use-conda --profile=slurm

# validate the files
md5sum -c checksum.md5
```

If you want to run with a different cluster manager or in the cloud, please
read the [advice below](#executing-on-a-cluster).


Configuration
-------------

DENTIST comprises a complex pipeline of with many options for tweaking. This
section points out some important parameters and their effect on the result or
performance.

The default parameters are rather **conservative**, i.e. they focus on
correctness of the result while not sacrificing too much sensitivity.

We also provide a **greedy** sample configuration
([`snakemake/dentist.greedy.yml`](./snakemake/dentist.greedy.yml)) which
focuses on sensitivity but may introduce more errors. _**Warning:** Use with
care! Always validate the closed gaps (e.g. manual inspection)._

In any case, the workflow creates an intermediate assembly
`workdir/{output_assembly}-preliminary.fasta` that contains all closed gaps,
i.e. before validation. It is accompanied by an AGP and BED file. You may
inspect these file for maximum sensitivity.


#### How to Choose DENTIST Parameters

While [the list of all commandline parameters][dentist-cli-summary] is a good
reference, it does not provide an overview of the important parameters.
Therefore, we provide this shorter list of important and influential
parameters. Please also consider adjusting the performance parameter in the
workflow configuration (`snakemake/snakemake.yml`).

- `--read-coverage`: This is the preferred way of providing values to
  `--max-coverage-reads`, `--max-improper-coverage-reads` and
  `--min-coverage-reads`. See below how their values are derived from
  `--read-coverage`.

    Ideally, the user provides the haploid read coverage which, can be
    inferred using a histogram of the alignment coverage across the assembly.
    Alternatively, the average raw read coverage can be used which is the
    number of base pairs in the reads divided by the number of base pairs
    in the assembly.

- `--ploidy`: Combined with `--read-coverage`, this parameters is the preferred
    way of providing `--min-coverage-reads`.

    We use the Wikipedia definition of ploidy, as "the number of complete sets of chromosomes in a cell" (https://en.wikipedia.org/wiki/Ploidy)

- `--max-coverage-reads`, `--max-improper-coverage-reads`: 
  These parameters are used to derive a repeat mask from the ref vs. reads
  alignment. If the coverage of (improper) alignments is larger than the given
  theshold it will be considered repetitive. If supplied, default values are derived from `--read-coverage` as follows:
  
    The maximum read coverage `C_max` is calculated from the global read
    coverage `C` (provided via --read-coverage) such that the probability of
    observing more than `C_max` alignments in a unique (non-repetitive) genomic
    region is very small (see [pre-print][dentist-bioarxiv], Methods section
    and Supplementary Table 2). In practice, this probability is approximated
    via

    ```
    C_max = floor(C / log(log(log(b * C + c) / log(a))))
    where
        a = 1.65
        b = 0.1650612
        c = 5.9354533
    ```

    To further increase the sensitivity, DENTIST searches for smaller
    repeat-induced local alignments. To this end, we define an alignment as
    proper if there are at most 100 bp (adjustable via
    --proper-alignment-allowance) of unaligned sequence on either end of the
    read. All other alignments, where only a smaller substring of the read
    aligns, are called improper. Improper alignments are often indicative of
    repetitive regions. Therefore, DENTIST considers genomic regions, where the
    number of improper read alignments is higher than a threshold to be
    repetitive. By default, this threshold equals half the global read coverage C.
    (see [pre-print][dentist-bioarxiv], Methods section). In practice, a smoothed
    version of `max(4, x/2)` is used to provide better performance for very low
    read coverage. The maximum improper read coverage `I_max` is computed as

    ```
    I_max = floor(a*x + exp(b*(c - x)))
    where
        a = 0.5
        b = 0.1875
        c = 8
    ```

- `--dust-{reads,ref}`, `--daligner-{consensus,reads-vs-reads,self}`,
  `--damapper-ref-vs-reads`, `--datander-ref`, `--daccord`:  
  These options allow passing parameters to the respective tools. They may have
  dramatic influence on the result. The default settings work well for PacBio
  CLR reads and should also work well with raw Nanopore data.

    In-depth discussion of each tool goes beyond the scope of this document,
  please refer to the respective documentations ([DBdust][DAZZ_DB],
  [daligner][daligner], [damapper][damapper], [datander][dascrubber],
  [daccord][daccord]).

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

- `--allow-single-reads`: May be used under careful consideration in
  combination with `--min-spanning-reads=1`. This is intended for one of the
  following scenarios:

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
  
    - `scaffoldGaps`: Closes only gaps that are marked by `N`s in the
      assembly. This is the default mode of operation. Use this if you do not
      want to alter the scaffolding of the assembly. See also
      `--existing-gap-bonus`.
    - `scaffolds`: Allows whole scaffolds to be joined in addition to the
      effects of `scaffoldGaps`. Use this if you have (many) scaffolds that
      are not yet full chromosome-scale.
    - `contigs`: Allows contigs to be rearranged freely. This is especially
      useful in _de novo_ assemblies **before** applying any other scaffolding
      methods as it increases the contiguity thus increasing the chance that
      large-scale scaffolding (e.g. Bionano or Hi-C) finds proper joins.

- `--min-coverage-reads`, `--min-spanning-reads`, `--region-context`:
  DENTIST validates closed gaps by mapping the reads to the gap-closed
  assembly. It requires for each gap and the base pairs down- and upstream
  (`--region-context`) are (1) covered by at least `--min-coverage-reads` reads
  at every position and (2) are spanned by at least `--min-spanning-reads`
  reads. Thus, increasing any of these numbers makes the *valid* gaps more
  robust but may reduce their number.

    If `--min-coverage-reads` is not provided, it will be derived from
    `--read-coverage` (see above) and `--ploidy`. Given (haploid) read coverage
    `C` and ploidy `p`, the minimum read coverage `C_min` is calculated as

    ```
      C_min = C / (2 * p)
    ```

    This corresponds to 50% of the long read coverage expected to be sequenced
    from a haploid locus (see [pre-print][dentist-bioarxiv], Methods section).


[dentist-cli-summary]: ./list-of-commandline-options.html


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

- `threads_per_process`: Sets the maximum number of threads/cores a single job
  may use. A single-threaded job will always allocate a single core but
  thread-parallel steps, e.g. the sequence alignments, will use up to
  `threads_per_process` if snakemake has been provided enough cores
  via `--cores`.
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


Troubleshooting
---------------


### Regular `ProtectedOutputException`

Snakemake has a [built-in facility to protect files][sm-protected-files] from
accidental overwrites. This is meant to avoid overwriting precious results
that took many CPU hours to produce. If executing a rule would overwrite a
protected file, Snakemake raises a `ProtectedOutputException`, e.g.:

```
ProtectedOutputException in line 1236 of /tmp/dentist-example/Snakefile:
Write-protected output files for rule collect:
workdir/pile-ups.db
  File "/usr/lib/python3.9/site-packages/snakemake/executors/__init__.py", line 136, in run_jobs
  File "/usr/lib/python3.9/site-packages/snakemake/executors/__init__.py", line 441, in run
  File "/usr/lib/python3.9/site-packages/snakemake/executors/__init__.py", line 230, in _run
  File "/usr/lib/python3.9/site-packages/snakemake/executors/__init__.py", line 155, in _run
```

Here `workdir/pile-ups.db` is the protected file that caused the error. If you
are sure of what you are doing, you can simply raise the protection by `chmod
-R +w ./workdir` and execute Snakemake again. Now, it will overwrite any files.


[sm-protected-files]: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#protected-and-temporary-files


### No internet connection on compute nodes

If you have no internet connection on your compute nodes or even the cluster
head node and want to use Singularity for execution, you will need to download
the container image manually and put it to a location accessible by all jobs.
Assume `/path/to/dir` is such a location on your cluster. Then download the
container image using

```sh
# IF internet connection on head node
singularity pull --dir /path/to/dir docker://aludi/dentist:stable

# ELSE (on local machine)
singularity pull docker://aludi/dentist:stable
# copy dentist_stable.sif to cluster
scp dentist_stable.sif cluster:/path/to/dir/dentist_stable.sif
```

When the image is in place you will need to adjust your configuration in
`snakemake.yml`:

```yaml
dentist_container: "/path/to/dir/dentist_stable.sif"
```

Now, you are ready for execution.

Note, if you want to use Conda without internet connection, you can just use the
pre-compiled binaries instead because they are just what Conda will install.
Be sure to adjust your `PATH` accordingly, e.g.:

```sh
PATH="$PWD/bin:$PATH" snakemake --configfile=snakemake.yml --profile=slurm
```


### Illegally formatted line from `DBshow -n`

This error message may appear in DENTIST's log files. It is a known bug that
will be fixed in a future release. In the meantime avoid FASTA headers that
contain a literal `" :: "`.


Citation
--------

> Arne Ludwig, Martin Pippel, Gene Myers, Michael Hiller. DENTIST – using long
> reads to close assembly gaps at high accuracy. __Submitted for peer review__.
> Pre-print at [https://doi.org/10.1101/2021.02.26.432990][dentist-bioarxiv]

[dentist-bioarxiv]: https://doi.org/10.1101/2021.02.26.432990 "Pre-print manuscript at bioRxiv"


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

This project is licensed under MIT License (see [LICENSE](./license.html)).
