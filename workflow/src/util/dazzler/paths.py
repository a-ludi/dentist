from pathlib import Path


def db_name(db_file):
    return Path(db_file).stem


def db_files(db):
    db = Path(db)
    assert db.suffix in {".db", ".dam"}

    suffixes = [".bps", ".idx"]
    if db.suffix == ".dam":
        suffixes.append(".hdr")
        suffixes.sort()

    hidden_file = db.with_name(f".{db.name}")

    return [db] + [hidden_file.with_suffix(suffix) for suffix in suffixes]


def homogenized_mask(mask):
    return mask + "-H"


def block_mask(mask, block):
    return "{}.{}".format(block, mask)


def block_masks(mask, db):
    return [block_mask(mask, b + 1) for b in range(get_num_blocks(db))]


def pseudo_block_mask(mask, block):
    return "{}-{}B".format(mask, block)


def pseudo_block_masks(mask, db, block=None):
    blocks = blocks_for(db, 1, block)

    return [pseudo_block_mask(mask, b) for b in blocks]


def mask_files(db, mask, *, block=None, pseudo_block=None):
    suffixes = ["anno", "data"]
    root = db.parent

    if pseudo_block is not None and block is None:
        block = str(pseudo_block)
        pseudo_block = True

    if not block is None:
        if pseudo_block:
            mask = pseudo_block_mask(mask, block)
        else:
            mask = block_mask(mask, block)

    def _mask_file(suffix):
        return db.with_name(f".{db.stem}.{mask}.{suffix}")

    return (_mask_file(suffix) for suffix in suffixes)


def alignment_file(db_a, db_b=None, block_a=None, block_b=None):
    if db_b is None:
        db_b = db_a
    db_a = db_name(db_a)
    db_b = db_name(db_b)
    parts = [db_a, block_a, db_b, block_b, "las"]

    return Path(".".join(str(part) for part in parts if part is not None))


def get_num_blocks(db):
    try:
        with db.open() as db_file:
            for line in db_file:
                if line.startswith("blocks ="):
                    return int(line.rpartition(" ")[2])

        raise EOFError("DB not split: {}".format(db))
    except FileNotFoundError:
        return 1


def get_blocks(db_or_num_blocks, batch_size=1, *, first=1):
    if isinstance(db_or_num_blocks, int):
        N = db_or_num_blocks
    else:
        N = get_num_blocks(db_or_num_blocks)

    return list(
        range(i, min(i + batch_size, N + 1)) for i in range(first, N + 1, batch_size)
    )


FULL_DB = object()


class block_range:
    def __init__(self, batch_size):
        self.batch_size = batch_size


def blocks_for(db, batch_size, block=None):
    def __dbblocks(db, batch_size):
        num_blocks = get_num_blocks(db)
        range_requested = isinstance(batch_size, block_range) or batch_size > 1
        if isinstance(batch_size, block_range):
            batch_size = batch_size.batch_size
        blocks = range(1, num_blocks + 1, batch_size)

        def __range(b):
            return "{}-{}".format(b, min(b + batch_size - 1, num_blocks))

        if range_requested:
            blocks = [__range(b) for b in blocks]

        return blocks

    def __block_or_range(block):
        if isinstance(block, int):
            return [block]
        elif isinstance(block, str):
            parts = block.split("-")

            if len(parts) == 1:
                return [int(parts[0])]
            elif len(parts) == 2:
                return range(int(parts[0]), int(parts[1]) + 1)
            else:
                raise ValueError("illegal value for block: " + block)
        else:
            return iter(block)

    if block is FULL_DB:
        return None
    elif block is None:
        return __dbblocks(db, batch_size)
    else:
        return __block_or_range(block)


def block_alignments(db_a, db_b=None, block_a=None, block_b=None):
    if db_b is None:
        db_b = db_a

    blocks_a = blocks_for(db_a, 1, block_a)
    blocks_b = blocks_for(db_b, 1, block_b)

    if not blocks_a is None and not blocks_b is None:
        return (alignment_file(db_a, db_b, i, j) for i in blocks_a for j in blocks_b)
    elif not blocks_b is None:
        return (alignment_file(db_a, db_b, block_b=j) for j in blocks_b)
    elif not blocks_a is None:
        return (alignment_file(db_a, db_b, block_a=i) for i in blocks_a)
    else:
        raise Exception("illegal block selection: both block_a and block_b are empty")
