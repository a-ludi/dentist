import re
import shlex
import subprocess

import yaml


def validate_config(config_file):
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
        subprocess.run(
            ["dentist", "validate-config", config_file],
            text=True,
            check=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        errors.append(e.stderr)

    if len(errors) > 0:
        raise Exception("; ".join(errors))

    return config


def generate_options_for(alignment_name, config_file, additional_flags=None):
    if not config_file.exists():
        return "false"

    generate_proc = subprocess.run(
        ["dentist", "generate", "--quiet", f"--config={config_file}"],
        text=True,
        check=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )
    commands = generate_proc.stdout.splitlines()

    if generate_proc.returncode != 0 or len(commands) == 0:
        error = generate_proc.stderr.decode()
        raise Exception(f"failed to get alignment commands: {error}")

    # iterate over pairs of lines
    for comment, command in zip(commands[::2], commands[1::2]):
        # find desired alignment type
        if alignment_name in comment:
            # clean up the command
            command = re.sub(r"<[^>]+>", "", command)
            command = command.rstrip()
            command = shlex.split(command)
            if not additional_flags is None:
                # add flags
                command.extend(additional_flags)
            return command

    raise Exception(
        "failed to get alignment command: unknown alignment_name `{alignment_name}`"
    )


def auxiliary_threads(threads):
    return max(1, threads // 4)


def main_threads(threads):
    return threads // auxiliary_threads(threads)
