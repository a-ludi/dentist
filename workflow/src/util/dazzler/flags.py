def ensure_threads_flag(flags, threads="{threads}"):
    return ensure_flags(flags, {"-T": threads})


def ensure_masks(flags, *masks):
    return ensure_flags(flags, {"-m": list(masks)})


def ensure_flags(flags, flags_with_values):
    def __keep_flag(flag):
        fname = flag[0:2]
        fvalue = flag[2:]

        if not fname in flags_with_values:
            return True
        elif fname in flags_with_values and not isinstance(
            flags_with_values[fname], list
        ):
            return False
        elif fname in flags_with_values and isinstance(flags_with_values[fname], list):
            return not fvalue in flags_with_values[fname]

    flags = list(f for f in flags if __keep_flag(f))

    for flag, value in flags_with_values.items():
        if value is True:
            flags.append(flag)
        elif isinstance(value, list):
            flags.extend((str(flag) + str(v) for v in value))
        elif value:
            flags.append(flag + str(value))

    return flags


def deduplicate_flags(flags):
    present_flags = set()
    deduplicated = list()
    for f in flags:
        if not f.startswith("-"):
            deduplicated.append(f)
        elif not f in present_flags:
            present_flags.add(f)
            deduplicated.append(f)

    return deduplicated


def assert_flag(
    flags, required_flag, message="required flag is missing: {missing_flags}"
):
    assert_flags(flags, [required_flag], message)


def assert_flags(
    flags, required_flags, message="required flag(s) are missing: {missing_flags}"
):
    present_flags = dict(((f[0:2], True) for f in flags))

    missing_flags = list()
    for required_flag in required_flags:
        if not present_flags.get(required_flag, False):
            missing_flags.append(required_flag)

    if len(missing_flags) > 0:
        e = AssertionError(message.format(missing_flags=", ".join(missing_flags)))
        e.missing_flags = missing_flags

        raise e
