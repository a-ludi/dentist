import subprocess
import yaml


def full_validate_dentist_config(config_file):
    errors = list()
    config = None
    with open(config_file, "r") as cf:
        config = yaml.safe_load(cf)

    def prohibit_revert_in_default():
        if "__default__" in config:
            if "revert" in config["__default__"]:
                errors.append("highly discouraged use of `revert` in `__default__`")

    def __get_presence_of_flags(command, *inspect_flags):
        flags = dict()
        if "__default__" in config:
            for flag in inspect_flags:
                flags[flag] = flag in config["__default__"]
        if command in config:
            for flag in inspect_flags:
                flags[flag] = flag in config[command]
            if "revert" in config[command]:
                for flag in inspect_flags:
                    if flag in config[command]["revert"]:
                        flags[flag] = False
        return flags

    def check_presence_of_read_masking_threshold():
        flags = __get_presence_of_flags(
            "mask-repetitive-regions", "read-coverage", "max-coverage-reads"
        )

        if not (flags["read-coverage"] ^ flags["max-coverage-reads"]):
            errors.append(
                "must specify either --read-coverage or --max-coverage-reads for command `mask-repetitive-regions`"
            )

    def check_presence_of_validation_coverage_threshold():
        flags = __get_presence_of_flags(
            "validate-regions", "read-coverage", "ploidy", "min-coverage-reads"
        )

        if not (
            (flags["read-coverage"] and flags["ploidy"]) ^ flags["min-coverage-reads"]
        ):
            errors.append(
                "must specify either --read-coverage and --ploidy or --min-coverage-reads for command `validate-regions"
            )

    prohibit_revert_in_default()
    check_presence_of_read_masking_threshold()
    check_presence_of_validation_coverage_threshold()

    try:
        dentist_validate_file(config_file)
    except subprocess.CalledProcessError as e:
        errors.append(e.stderr)

    if len(errors) > 0:
        raise Exception("; ".join(errors))

    return config


def dentist_validate_file(config_file):
    subprocess.run(
        ["dentist", "validate-config", config_file],
        text=True,
        check=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )


full_validate_dentist_config(snakemake.input[0])
