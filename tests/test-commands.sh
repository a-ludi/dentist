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
# generate with: md5sum data/* data/.??* dentist.json snakemake.json | gzip | base64
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


function dentist()
(
    [[ -v DENTIST_BINARY ]] \
    || DENTIST_BINARY="$(command which dentist)" \
    || DENTIST_BINARY="$REPO_ROOT/dentist"

    command "$DENTIST_BINARY" "$@"
)


function setup_workdir()
{
    EXECDIR="$PWD"
    WORKDIR="$(realpath "$(mktemp --tmpdir -d ".test-commands-XXXXXX")")"
    trap clean_up exit
    (( ${DEBUG:-0} == 0 )) || trap save_workdir err

    cd "$WORKDIR"
}


function save_workdir()
{
    mv "$WORKDIR" "$EXECDIR/test-commands-failed"

    echo "saved results under $EXECDIR/test-commands-failed"
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
(
    set -x
    mkdir -p workdir workdir/insertions
    fasta2DB workdir/reads.db data/reads.fasta && DBsplit -x20 workdir/reads.db
    fasta2DAM workdir/assembly-test.dam data/assembly-test.fasta && DBsplit -x20 workdir/assembly-test.dam
    ( cd workdir && datander '-T1' -s126 -l500 -e0.7 assembly-test.1 )
    TANmask -v -ntan workdir/assembly-test workdir/TAN.assembly-test.1.las
    Catrack -v workdir/assembly-test.dam tan
    DBdust workdir/assembly-test
    ( cd workdir && daligner -I '-T1' -l500 -e0.7 -mdust -mtan assembly-test.1 assembly-test.1 )
    LAmerge -v workdir/assembly-test.assembly-test.las workdir/assembly-test.1.assembly-test.1.las
    dentist mask --config=dentist.json  workdir/assembly-test.dam workdir/assembly-test.assembly-test.las dentist-self
    ( cd workdir && damapper -C '-T1' -e0.7 -mdust -mdentist-self -mtan assembly-test reads.1 )
    dentist propagate-mask --config=dentist.json  -m dentist-self workdir/assembly-test.dam workdir/reads.db workdir/assembly-test.reads.1.las dentist-self-1B
    dentist propagate-mask --config=dentist.json  -m dentist-self-1B workdir/reads.db workdir/assembly-test.dam workdir/reads.1.assembly-test.las dentist-self-H-1B
    dentist merge-masks --config=dentist.json  workdir/assembly-test.dam dentist-self-H dentist-self dentist-self-H-1B
    dentist propagate-mask --config=dentist.json  -m tan workdir/assembly-test.dam workdir/reads.db workdir/assembly-test.reads.1.las tan-1B
    dentist propagate-mask --config=dentist.json  -m tan-1B workdir/reads.db workdir/assembly-test.dam workdir/reads.1.assembly-test.las tan-H-1B
    dentist merge-masks --config=dentist.json  workdir/assembly-test.dam tan-H tan tan-H-1B
    LAmerge -v workdir/assembly-test.reads.las workdir/assembly-test.reads.1.las
    dentist mask --config=dentist.json  workdir/assembly-test.dam workdir/reads.db workdir/assembly-test.reads.las dentist-reads
    dentist propagate-mask --config=dentist.json  -m dentist-reads workdir/assembly-test.dam workdir/reads.db workdir/assembly-test.reads.1.las dentist-reads-1B
    dentist propagate-mask --config=dentist.json  -m dentist-reads-1B workdir/reads.db workdir/assembly-test.dam workdir/reads.1.assembly-test.las dentist-reads-H-1B
    dentist merge-masks --config=dentist.json  workdir/assembly-test.dam dentist-reads-H dentist-reads dentist-reads-H-1B
    dentist collect --config=dentist.json  --threads=1 --auxiliary-threads=1 --mask=dentist-self-H,tan-H,dentist-reads-H workdir/assembly-test.dam workdir/reads.db workdir/assembly-test.reads.las workdir/pile-ups.db
    dentist process --config=dentist.json  --threads=1 --auxiliary-threads=1 --mask=dentist-self-H,tan-H,dentist-reads-H --batch=0..2 workdir/assembly-test.dam workdir/reads.db workdir/pile-ups.db workdir/insertions/batch.0.db
    dentist merge-insertions workdir/insertions.db workdir/insertions/batch.0.db
    dentist output --config=dentist.json  --agp=workdir/gap-closed-preliminary.closed-gaps.agp --closed-gaps-bed=workdir/gap-closed-preliminary.closed-gaps.bed --revert=scaffolding,skip-gaps,skip-gaps-file,agp-dazzler --agp-dazzler workdir/assembly-test.dam workdir/reads.db workdir/insertions.db workdir/gap-closed-preliminary.fasta
    fasta2DAM workdir/gap-closed-preliminary.dam workdir/gap-closed-preliminary.fasta && DBsplit -x20 workdir/gap-closed-preliminary.dam
    dentist bed2mask --config=dentist.json  --data-comments --bed=workdir/gap-closed-preliminary.closed-gaps.bed workdir/gap-closed-preliminary.dam closed-gaps
    ( cd workdir && datander '-T1' -s126 -l500 -e0.7 gap-closed-preliminary.1 )
    TANmask -v -ntan workdir/gap-closed-preliminary workdir/TAN.gap-closed-preliminary.1.las
    Catrack -v workdir/gap-closed-preliminary.dam tan
    DBdust workdir/gap-closed-preliminary
    ( cd workdir && daligner -I '-T1' -l500 -e0.7 -mdust -mtan gap-closed-preliminary.1 gap-closed-preliminary.1 )
    LAmerge -v workdir/gap-closed-preliminary.gap-closed-preliminary.las workdir/gap-closed-preliminary.1.gap-closed-preliminary.1.las
    dentist mask --config=dentist.json  workdir/gap-closed-preliminary.dam workdir/gap-closed-preliminary.gap-closed-preliminary.las dentist-self
    ( cd workdir && damapper -C '-T1' -e0.7 -mdust -mdentist-self -mtan gap-closed-preliminary reads.1 )
    LAmerge -v workdir/gap-closed-preliminary.reads.las workdir/gap-closed-preliminary.reads.1.las
    LAsplit workdir/gap-closed-preliminary.@.reads.las 1 < workdir/gap-closed-preliminary.reads.las
    dentist validate-regions --config=dentist.json --threads=1 --weak-coverage-mask=1.dentist-weak-coverage workdir/gap-closed-preliminary.dam workdir/reads.db workdir/gap-closed-preliminary.1.reads.las closed-gaps > workdir/validation-report.1.json
    dentist merge-masks --config=dentist.json  workdir/gap-closed-preliminary.dam dentist-weak-coverage 1.dentist-weak-coverage
    cat workdir/validation-report.1.json > workdir/validation-report.json
    jq -r 'select(.isValid // false | not) | .contigIds | map(tostring) | join("-")' workdir/validation-report.json > workdir/skip-gaps.txt
    dentist output --config=dentist.json  --agp=gap-closed.agp --closed-gaps-bed=gap-closed.closed-gaps.bed --skip-gaps-file=workdir/skip-gaps.txt workdir/assembly-test.dam workdir/reads.db workdir/insertions.db gap-closed.fasta
)


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
