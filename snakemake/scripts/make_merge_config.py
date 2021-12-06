import json
import subprocess


def make_merge_config(config_file, insertions_batches):
    merge_config = {"merge-insertions": {"partitioned-insertions": insertions_batches}}

    with open(config_file, "w") as merge_config_file:
        json.dump(merge_config, merge_config_file)


def dentist_validate_file(config_file):
    subprocess.run(
        ["dentist", "validate-config", config_file],
        text=True,
        check=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )


make_merge_config(snakemake.output[0], snakemake.input.insertions_batches)
dentist_validate_file(snakemake.output[0])
