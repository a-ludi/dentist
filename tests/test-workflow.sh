#!/bin/bash

# Unofficial Bash Strict Mode (http://redsymbol.net/articles/unofficial-bash-strict-mode/)
set -euo pipefail
IFS=$'\n\t'

SIMULATOR_ARGS=(
    -m25000
    -s12500
    -e.13
    -c20
    -r1724161952
)

REPO_ROOT="$(realpath "$(dirname "$0")/..")"
SNAKEFILE="$REPO_ROOT/snakemake/Snakefile"
ASSEMBLY_REF_DATA="\
H4sIAAAAAAAAA02XQZLkNgwE7/sd3/0I/4CBA3x24P9hIrOomdneboqSSKCAKoB/17///fXPnzp1
Zqaqa39mv8796a6qs6Px68ydOXtVfvr0feB+eeMO7tS9vsOzz937O9h/vn6y+o7Kl6dnN9nH9mVM
OWvSjrt3AVeZyqiGezV7/863pu4Fa645u247sffvE/fx/bpzvVOsO9w52jYTYzEAH7xYf7nYF/uO
ejdZs9e/NQs/RxvveF95m07WdDPevoac/jFxHy3NyZ5HvzHB+wnBGb5EdxZF99rf2sVYm791cZda
txbTGNi+3FkzccKJAfwX3qkXNez7HNm1yyjq5IlJ68Q8xNbSajIgbhO9dYcEughiPsljHpxA4I/h
30fWLJNigG40Sl/2usf/wXjjnDc1IzbsbJ24xOMmq9A0KddA6s/m96ZLEpBFCrtJ1bcOHvS8S57m
bjEvh0xKoh+EyK8gM2FIm6RjMMoQjuzq/ZCCLd+aDNpB9U/MT9inv03wd6H+ud8ExkQ2uMtA2IwT
MunLAlRh3d58J5YNLIPXj1WV7zUN7ah2rp3Rb+BeeyT2YDhPoAF5DuDCfBMtAGIOjpW7lAyaUPmT
MWebUGqX106GBiMxehSvTvB27/Nl2K5R/YIaOjfbXxdXXEi1iIliZuKpA7XhL8iZ9DD4S571wMg3
YH/KVijMLsTXAUZFcdQFsqEAA1Eq9XRYALrCv6Om4/rysxRYJGXvyLBm+db5IXYtnT8SLhtO0Fm0
mq3bXFKTlUDs3MXL9EOjA/QJSKBTHQmcCCreYBJmlfBIWvR80aZoRBdKHlZ3BJ8a8GVQNNxcV3TB
M8wmO3EA1PlJFnhb3siZVD5EUHVYhDql45XPu3BHWElhbG8qZcRe4kb0TOEdJIEJnKjUecpQqm0k
ZmMFdvPrDwqXCQVh+tewqBnHsjWpx/MFovrVpshudpdQQo0br+KY2iMeQIncgBlSPalsFMJK1Qx/
O4njDSWtLCTZgoiqLyk6rzzXV40V/KTVsXRAsYrcncRVLTuI8IFju01EZwUGZcUQxLVPGocax+YV
GpJ7Sm4qgHOQZ6KS52ugXFyjALLPKwnhK0mK6oerEynKtME7GeOyxHA0T3HQh3Lzthw8gFIdgIH8
s/kaOw4YPsH6lFvzNOpmltmWdZqNE/m3vmu1DD9pB9lF8CrNFU+BOzXGMrIbtikFHBr8k4gdJUyE
Ua8D65DEipCbupXcVOJpbx5XQShVXN9bMfhVM+fYuaiVbA5EdqQmnq3w686exE/aAVjadm9qjxw6
nVpC3VCPUOjuXyld9gDaB57U0jqW2lL9Jp3MSasepXu8ZE5pwMN00nLxtN2p6Y8KmWpR6nqQTmXL
fhBbhytKCeHhJ/mrINtUNCy2LqTamf9Z0tmjbliVKF/liqxd1lbT+HylFRMR3QqWr1E26iQBnb1c
7fD3RBxKCDHNSlQ57jA20aJAAJMkHfvtSo9q8U9FOLyflkOWBkkVsDoFI2coTwb9MFH4Qpm2QG/I
vhMKcE2qN11tJ2Ki/a6Pfa3g5KCQLpEUEgTb7g+TuOIdG78ogAe/CG1auHTLcsK+hgxNL0c/Y8uQ
J3XOqh1VDI7zMSCno3TOlcNfJN1GksLgueR9bITrfM3EeYr2PLPkhNovA+x0tYv39fWlZMtze0HC
BD7Hih7tVOBc8PVEKVGRPU9g8KJfdE5w8GySzrHUbvNy0hqXFlBIx8KV3uJXpX5nhfF0ld4h/tnu
8CCVNS/UvFY6p9qc5/PT9ed/P27fQgkQAAA="
GAP_SEQ_DATA="\
H4sIAAAAAAAAAyWMwQ0AMAgC/92S8GCB2z9FjVFQwCBPATGinSCFKJLgRnWHXevUHD1onaHcg311
sY1Uex/yrOV8YgAAAA=="
DENTIST_CONFIG_DATA="\
H4sIAAAAAAAAA6tWio9PSU1LLM0piY9XsqpWKkpNTNFNzi9LLUpMT1WyMjLQM9BRKsjJz0ypBPJq
a7kALygudjIAAAA="
WORKFLOW_CONFIG_DATA="\
H4sIAAAAAAAAA22QQU/DMAyF/4vPKWwILr0BJw4TCLghFLmJ24WFJIpdYEz976SUaCBxe8/2++zk
AJaCOBZtYujdAG0tnLxwDKCgH73Xb+idRXGl0koeSYELaRSG9gCZesoUDM1ZFDxFZnrt/L4RKpge
WbBwMqHlOvJt/ra07NOMuLu8vrq51Q+b+0eYFMRR6qJF6oovswOmxvjIZH9Ykzqeo23HyTuB9gma
j7MVPNdN/zTeY95ZlwuzKgU+DkupCC5etks8UdYpR0NczlorKDrhgEK6QzFbze6zPGR9ruC3v1gp
OH6j7nw0uzk/fQFgPnmIgwEAAA=="
# generate with: md5sum data/* data/.??* dentist.json snakemake.json checksums.output.md5 | bzip2 | base64
INPUT_CHECKSUMS_DATA="\
H4sIAAAAAAAAA33QOVIDQQwF0HxOMRfA9KZW6zhaCwMeKM8EcHu6ChMBDhQ+/S9xaiHFSrNgk1GE
CncuoNSKD+V1NT74kffdL/L6+XD18Ktv6ifjy4LWcyoFEtSSa62IBu5CTVRdgf7XwfvBC0qSpmEV
rSG1rmjmWZ1ihqP8Sj98P25UR9UZ1C31iEpCLKw1NR1Sc+o/0Vdn22+kNuAiCF6BB4EqG88LGKr0
LDVu5PRHXXnfl4SBKRcdLfqsXLvPm5G0QmNwpjv8ya5LSr1bt8Y6mQAMs5yiDS+cBCXf4Wf7WDJL
80ruHDiYxSywpw5Ec0limdy34zz/87y/bQvknAE9GpA654JiodaFQ8Mp9XXdN37xy5xv8AWelVCk
CQIAAA=="
# generate with: md5sum gap-closed.closed-gaps.bed gap-closed.fasta | gzip | base64
OUTPUT_CHECKSUMS_DATA="\
H4sIAAAAAAAAA1XLSwrDMAxF0XlWkQ0kyJblyMt51ieTQgvu/mmho44uHLiSOppIFiPtIJ9zFmlj
hJVrBPO+33gd9niu8POX4yvrnOGbsXJ3IwKnRIkKaFYoVUfzfv3NifXG9gHK8nuEcAAAAA=="


function _decompress()
{
    base64 -d | zcat
}


function _compress()
{
    gzip | base64
}


function _reference_assembly()
{
    _decompress <<<"$ASSEMBLY_REF_DATA"
}


function _gap_seq()
{
    _decompress <<<"$GAP_SEQ_DATA"
}


function _test_assembly()
{
    local GAP_SEQ
    GAP_SEQ="$(_gap_seq)"

    _reference_assembly \
    | sed "s/$GAP_SEQ/${GAP_SEQ//[acgt]/n}/"
}


function _reads()
{
    simulator "${SIMULATOR_ARGS[@]}" "$@"
}


function _dentist_config()
{
    _decompress <<<"$DENTIST_CONFIG_DATA"
}


function _workflow_config()
{
    _decompress <<<"$WORKFLOW_CONFIG_DATA"
}


function _input_checksums()
{
    _decompress <<<"$INPUT_CHECKSUMS_DATA"
}


function _output_checksums()
{
    _decompress <<<"$OUTPUT_CHECKSUMS_DATA"
}


function setup_workdir()
{
    EXECDIR="$PWD"
    WORKDIR="$(realpath "$(mktemp -d ".test-workflow-XXXXXX")")"
    trap clean_up exit
    (( ${DEBUG:-0} == 0 )) || trap save_workdir err

    cd "$WORKDIR"
}


function save_workdir()
{
    mv "$WORKDIR" "$EXECDIR/test-workflow-failed"

    echo "saved results under $EXECDIR/test-workflow-failed"
}


function clean_up()
{
    rm -rf "$WORKDIR"
}


function generate_files()
{
    mkdir data
    _reference_assembly \
    | tee data/assembly-reference.fasta \
    | fasta2DAM -i data/assembly-reference.dam
    _reads data/assembly-reference.dam > data/reads.fasta
    _test_assembly > data/assembly-test.fasta
    _dentist_config > dentist.json
    _workflow_config > snakemake.json
    _input_checksums > checksums.input.md5
    _output_checksums > checksums.output.md5

    md5sum -c checksums.input.md5 && echo "Input files are valid."
}


function run_workflow()
{
    snakemake --snakefile="$SNAKEFILE" --configfile=snakemake.json -j1 -p "$@"
}


function run_tests()
{
    md5sum -c checksums.output.md5 && echo "Output files are valid."
}


function main()
{
    setup_workdir
    generate_files
    run_workflow "$@"
    run_tests
}


main "$@"
