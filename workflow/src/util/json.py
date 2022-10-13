import json


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
