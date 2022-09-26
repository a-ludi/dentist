DENTIST\: Mini Example
======================

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat)](https://github.com/RichardLitt/standard-readme)
![License](https://img.shields.io/github/license/a-ludi/dentist)
[![GitHub](https://img.shields.io/badge/GitHub-code-blue?logo=github)][dentist]

> A small example to test DENTIST's workflow

Quickly test [DENTIST][dentist] with this example workflow. It uses part of the
_D. melanogaster_ reference assembly (dm6) and simulated reads to demonstrate
the workflow. The full source code of DENTIST is available at
<https://github.com/a-ludi/dentist>.

If you experience issues, please search the issues at [DENTIST's
repository](https://github.com/a-ludi/dentist/issues) or create a new one if
you cannot find an answer to your problem.


Table of Contents
-----------------

- [Install](#install)
- [Usage](#usage)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Maintainer](#maintainer)
- [Contributing](#contributing)
- [License](#license)


Install
-------

Make sure you have [Snakemake][snakemake] 5.32.1 or later installed.

You can also use the convenient Conda package to execute the rules. Just make
sure you have [Mamba][mamba] installed.


Usage
-------

First of all download the test data and workflow and switch to the
`dentist-example` directory.

```sh
wget https://github.com/a-ludi/dentist/releases/download/v4.0.0/dentist-example.tar.gz
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

Make sure [Mamba][mamba] (a frontend for [Conda][conda]) is installed on your
system. Execute the workflow without explicit installation by adding
`--use-conda` to the call to Snakemake:

```sh
# run the workflow
snakemake --configfile=snakemake.yml --use-conda --cores=all

# validate the files
md5sum -c checksum.md5
```

*Note: If you do not have `mamba` installed, you may need to pass
`--conda-frontend=conda` to Snakemake.*


### Execution in Singularity Container (discouraged)

*Remark: the Singularity container may not work properly depending on your
system. (see issue #30)*

Execute the workflow inside a convenient Singularity image by adding
`--use-singularity` to the call to Snakemake:

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
cp -v "profile-slurm.sync.yml" "$HOME/.config/snakemake/slurm/config.yaml"
# execute using the cluster profile
snakemake --configfile=snakemake.yml --use-singularity --profile=slurm

# validate the files
md5sum -c checksum.md5
```

Find more advice on execution using a cluster in the
[online documentation][dentist-cluster].


[dentist-cluster]: https://a-ludi.github.io/dentist/#cluster-execution


Citation
--------

> Arne Ludwig, Martin Pippel, Gene Myers, Michael Hiller. DENTIST — using long
> reads for closing assembly gaps at high accuracy. _GigaScience_, Volume 11,
> 2022, giab100.
> [https://doi.org/10.1093/gigascience/giab100][dentist-gigascience]

[dentist-gigascience]: https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giab100/6514926 "Paper published at GigaScience"


Maintainer
----------

DENTIST is being developed by Arne Ludwig &lt;<ludwig@mpi-cbg.de>&gt; at
the Planck Institute of Molecular Cell Biology and Genetics, Dresden, Germany.


License
-------

This project is licensed under MIT License (see [LICENSE](./LICENSE)).


[dentist]: https://a-ludi.github.io/dentist/
[snakemake]: https://snakemake.readthedocs.io/en/v5.32.1/getting_started/installation.html
[singularity]: https://sylabs.io/guides/3.5/user-guide/quick_start.html
[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[mamba]: https://github.com/mamba-org/mamba
