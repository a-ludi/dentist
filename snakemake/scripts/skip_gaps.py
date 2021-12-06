import json
import subprocess


def generate_skip_gaps(validation_report, skip_gaps):
    validations = None
    with open(validation_report, "r") as validation_report_file:
        validations = multi_json_load(validation_report_file)

    with open(skip_gaps, "w") as skip_gaps_file:
        for validation in validations:
            if not validation.get("isValid", False):
                gap_spec = "-".join((str(cid) for cid in validation["contigIds"]))
                print(gap_spec, file=skip_gaps_file)


def multi_json_loads(s, **kwargs):
    """Deserialize ``s`` (a ``str``, ``bytes`` or ``bytearray`` instance
    containing whitespace-separated JSON documents) to a list of Python
    objects.
    """

    decoder = json.JSONDecoder(**kwargs)

    start = len(s) - len(s.lstrip())
    docs = list()
    while start < len(s):
        doc, pos = decoder.raw_decode(s, idx=start)
        docs.append(doc)
        lsep = len(s[pos:]) - len(s[pos:].lstrip())

        if lsep == 0:
            raise json.JSONDecodeError("Extra data", s, pos)

        start = pos + lsep

    return docs


def multi_json_load(fp, **kwargs):
    """Deserialize ``fp`` (a ``.read()``-supporting file-like object containing
    containing whitespace-separated JSON documents) to a list of Python
    objects.
    """

    return multi_json_loads(fp.read(), **kwargs)


generate_skip_gaps(snakemake.input[0], snakemake.output[0])
