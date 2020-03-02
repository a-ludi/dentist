#-----------------------------------------------------------------------------
# BEGIN functions
#-----------------------------------------------------------------------------

def rel_to_workdir(path):
    from os.path import relpath

    return relpath(path, config["workdir"])


def db_stub(db_file):
    from os.path import splitext

    return splitext(db_file)[0]


def db_name(db_file):
    from os.path import basename

    return basename(db_stub(db_file))


def prepend_ext(filename, ext):
    from os.path import splitext

    parts = splitext(filename)

    return parts[0] + ext + parts[1]


def fasta_to_workdb(fasta_file, ext):
    return join(config["workdir"], "{db}.{ext}".format(db=db_name(fasta_file), ext=ext))


def fasta2dazz_command(target_db):
    if target_db.endswith(".db"):
        return "fasta2DB"
    else:
        return "fasta2DAM"


def db_files(db):
    from os.path import dirname

    hidden_db_file_suffixes = [".bps", ".idx"]
    hidden_dam_file_suffixes = [".bps", ".hdr", ".idx"]
    root = dirname(db)
    suffixes = hidden_db_file_suffixes if db.endswith(".db") \
                else hidden_dam_file_suffixes

    def __hidden_file(suffix):
        if root:
            return "{}/.{}{}".format(root, db_name(db), suffix)
        else:
            return ".{}{}".format(db_name(db), suffix)

    return [db] + [__hidden_file(suffix)  for suffix in suffixes]


_fast2dam_checkpoint = dict()


def await_db_files(db):
    from os.path import splitext

    def __await_db_files(wildcards):
        if db in _fast2dam_checkpoint:
            _checkpoint = getattr(checkpoints, _fast2dam_checkpoint[db])

            return _checkpoint.get(workdir=workdir_).output
        else:
            raise Exception("unknown db: {}".format(db));

    return __await_db_files


def await_pile_ups():
    def __await_pile_ups(wildcards):
        return checkpoints.collect.get().output

    return __await_pile_ups


def alignment_file(db_a, db_b=None, block_a=None, block_b=None):
    if db_b is None:
        db_b = db_a
    db_a = db_name(db_a)
    db_b = db_name(db_b)

    filename = None
    if block_a is None and block_b is None:
        filename = "{}.{}.las".format(db_a, db_b)
    elif block_a is None:
        filename = "{}.{}.{}.las".format(db_a, db_b, block_b)
    elif block_b is None:
        filename = "{}.{}.{}.las".format(db_a, block_a, db_b)
    else:
        filename = "{}.{}.{}.{}.las".format(db_a, block_a, db_b, block_b)

    return join(config["workdir"], filename)


def make_flags(flags):
    try:
        from shlex import quote
    except ImportError:
        from pipes import quote

    if str(flags) == flags:
        return flags
    else:
        return " ".join(quote(flag)  for flag in flags)


def append_flags(flags, *new_flags):
    try:
        from shlex import quote
    except ImportError:
        from pipes import quote

    if len(flags) == 0:
        return make_flags(new_flags)
    else:
        return flags + " " + make_flags(new_flags)


def prepare_flags(flags):
    def secondary_expand(wildcards, input=None, output=None, threads=None, resources=None):
        return flags.format(wildcards=wildcards,
                            input=input,
                            output=output,
                            threads=threads,
                            resources=resources)

    return secondary_expand


def log_file(step):
    return join(config["logdir"], "{}.log".format(step))


def auxiliary_threads(wildcards, threads):
    if "auxiliary_threads" in config:
        return int(config["auxiliary_threads"])
    else:
        return max(1, threads // 4)


def main_threads(wildcards, threads):
    return threads // auxiliary_threads(wildcards, threads=threads)


def validate_dentist_config(config_file):
    subprocess.check_output(["dentist", "validate-config", config_file])


def generate_options_for(alignment_name, config_file):
    generate_template = "dentist generate --quiet --config={} | sed -nE '/^#.*{}/ {{ n; s/^([^<]+).*/\\1/p }}'"
    generate_cmd = generate_template.format(config_file, alignment_name)

    return subprocess.check_output(generate_cmd, shell=True).decode().rstrip()


def get_num_blocks(db):
    try:
        with open(db, 'r') as db_file:
            for line in db_file:
                if line.startswith("blocks ="):
                    return int(line.rpartition(" ")[2])

        raise EOFError("DB not split: {}".format(db))
    except FileNotFoundError:
        return 1


__num_contigs_regex = re.compile(r"^\+\s+R\s+(?P<num_contigs>\d+)\s*$", re.MULTILINE)


def get_num_contigs(db):
    try:
        dbdump = subprocess.check_output(["DBdump", db]).decode()
        match = __num_contigs_regex.match(dbdump)

        if not match:
            raise Exception("Could not read number of contigs in {}".format(db))

        return int(match.group("num_contigs"))
    except subprocess.CalledProcessError:
        return 1


def block_alignments(db_a, db_b=None, damapper=False):
    num_blocks_a = get_num_blocks(db_a) if not damapper else 0

    if db_b is None:
        db_b = db_a
        num_blocks_b = num_blocks_a
    else:
        num_blocks_b = get_num_blocks(db_b)

    blocks_a = range(1, num_blocks_a + 1)
    blocks_b = range(1, num_blocks_b + 1)

    if damapper:
        return (alignment_file(db_a, db_b, block_b=j)  for j in blocks_b)
    else:
        return (alignment_file(db_a, db_b, i, j)  for i in blocks_a for j in blocks_b)


def mask_files(db, mask):
    from os.path import dirname

    suffixes = ["anno", "data"]
    root = dirname(db)

    def __mask_files(suffix):
        if root:
            return "{}/.{}.{}.{}".format(root, db_name(db), mask, suffix)
        else:
            return ".{}.{}.{}".format(db_name(db), mask, suffix)

    return (__mask_files(suffix)  for suffix in suffixes)


def batch_range(wildcards, batch_size, num_elements):
    batch_id = int(wildcards.batch_id, base=10)
    from_index = batch_id * batch_size
    to_index = min(from_index + batch_size, num_elements)

    return "{}..{}".format(from_index, to_index)


def insertion_batch_range(wildcards):
    return batch_range(wildcards, config["batch_size"], get_num_pile_ups())


def get_num_pile_ups():
    from subprocess import CalledProcessError
    from os.path import exists

    try:
        info_cmd = ["dentist", "show-pile-ups", "-j", pile_ups]
        pile_ups_info = subprocess.check_output(info_cmd, stderr=subprocess.DEVNULL)

        return json.loads(pile_ups_info)["numPileUps"]
    except CalledProcessError as e:
        if exists(pile_ups):
            raise e

        return 1


def ceildiv(n, d):
    return (n + d - 1) // d


def insertions_batches():
    num_pile_ups = get_num_pile_ups()
    num_batches = ceildiv(num_pile_ups, config["batch_size"])
    num_digits = int(ceil(log10(num_batches)))
    batch_id_format = "{{:0{}d}}".format(num_digits)

    def get_insertions_batch(batch_id):
        workdir_ = config["workdir"]
        insertions_batch = config["workflow"]["insertions_batch"]
        batch_id = batch_id_format.format(batch_id)

        return join(workdir_, insertions_batch.format(batch_id=batch_id))

    return (get_insertions_batch(i)  for i in range(num_batches))

#-----------------------------------------------------------------------------
# END functions
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Variables for rules
#-----------------------------------------------------------------------------

# config shortcuts
inputs = config["inputs"]
outputs = config["outputs"]
workdir_ = config["workdir"]
workflow_ = config["workflow"]
max_threads = config["max_threads"]

# workflow files
reference_fasta = inputs["reference"]
reads_fasta = inputs["reads"]
output_assembly = outputs["output_assembly"]
reference = fasta_to_workdb(reference_fasta, "dam")
reads_db_type = "db" if inputs["reads_type"] == "PACBIO_SMRT" else "dam"
reads = fasta_to_workdb(reads_fasta, reads_db_type)
dentist_config_file = join(workdir_, workflow_["config"])
dentist_merge_config_file = prepend_ext(dentist_config_file, ".merge")
self_alignment = alignment_file(reference)
tandem_alignment = alignment_file("TAN", reference)
ref_vs_reads_alignment = alignment_file(reference, reads)
self_mask = workflow_["self_mask"]
tandem_mask = "tan"
reads_mask = workflow_["reads_mask"]
masks = [self_mask, tandem_mask, reads_mask]
pile_ups = join(workdir_, workflow_["pile_ups"])
insertions_batch = join(workdir_, workflow_["insertions_batch"])
insertions = join(workdir_, workflow_["insertions"])

# command-specific
from os import environ
dentist_flags = environ.get("DENTIST_FLAGS", "")
dbsplit_flags = make_flags(config["dbsplit_flags"])
dalign_flags = make_flags(config["dalign_flags"])
datander_flags = make_flags(config["datander_flags"])

if "TMPDIR" in environ:
    dalign_flags = append_flags(dalign_flags, "-P{}".format(environ["TMPDIR"]))
    datander_flags = append_flags(datander_flags, "-P{}".format(environ["TMPDIR"]))


localrules:
    extend_dentist_config,
    reference2dam


rule extend_dentist_config:
    output: dentist_config_file
    run:
        dentist_config = config["dentist_config"]

        if "__default__" not in dentist_config:
            dentist_config["__default__"] = dict()
        dentist_config["__default__"]["reference"] = reference
        dentist_config["__default__"]["reads"] = reads
        dentist_config["__default__"]["result"] = output_assembly
        dentist_config["__default__"]["ref-vs-reads-alignment"] = ref_vs_reads_alignment
        dentist_config["__default__"]["mask"] = masks
        dentist_config["__default__"]["pile-ups"] = pile_ups
        dentist_config["__default__"]["insertions"] = insertions

        if "mask-repetitive-regions" not in dentist_config:
            dentist_config["mask-repetitive-regions"] = dict()
        dentist_config["mask-repetitive-regions"]["reads"] = None

        with open(output[0], 'w') as ext_config_file:
            json.dump(dentist_config, ext_config_file)

        validate_dentist_config(output[0])


_fast2dam_checkpoint[reference] = "reference2dam"
checkpoint reference2dam:
    input: reference_fasta
    output: *db_files(reference)
    shell:
        "fasta2DAM {output[0]} {input} && DBsplit {dbsplit_flags} {output[0]}"
