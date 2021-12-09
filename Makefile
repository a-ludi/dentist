dentist_version = $(shell git describe)
arch = $(shell uname -m)
container = dentist_$(dentist_version).sif
container_def = singularity/dentist_$(dentist_version).def
dentist_env=dentist_$(firstword $(subst ., ,$(dentist_version)))

local_only = $(if $(findstring local,$(dentist_version)),$(1))

dentist_releases = $(shell git tag | sed -nE 's/v(.*)/\1/p' | tr '\n' ' ')
singularity_definitions = $(addprefix singularity/dentist_,$(addsuffix .def,$(dentist_releases)))

ddox = dub run -q ddox --
ddox_filter = \
    "--min-protection=Protected" \
    "--ex" "dentist.swinfo" "--ex" "dentist.modules" \
    "--ex" "dentist.common.binio._testdata"


dist_dir = dentist.$(dentist_version).$(arch)
dist_tarball = dentist.$(dentist_version).$(arch).tar.gz
binaries=$(addprefix $(dist_dir)/bin/,Catrack computeintrinsicqv daccord \
    daligner DAM2fasta damapper DAScover DASqv datander DB2fasta DBa2b DBb2a \
    DBdump DBdust DBmv DBrm DBshow DBsplit DBstats DBtrim DBwipe dentist \
    dumpLA fasta2DAM fasta2DB LAa2b LAb2a LAcat LAcheck LAdump LAmerge \
    lasfilteralignments LAshow LAsort LAsplit rangen simulator TANmask)
info_files = $(addprefix $(dist_dir)/,README.md CHANGELOG.md LICENSE)
workflow_files = $(addprefix $(dist_dir)/,\
    cluster.yml dentist.yml dentist.greedy.yml snakemake.yml \
    envs/$(dentist_env).yml Snakefile \
    $(patsubst snakemake/%,%,$(wildcard snakemake/profile-slurm.*.yml)) \
    $(patsubst snakemake/%,%,$(wildcard snakemake/scripts/*.py)))

dist_files = $(info_files) $(workflow_files) $(binaries)



.PHONY: all
all: singularity-image conda-env



.PHONY: singularity-image
singularity-image: $(container)

dentist_%.sif: singularity/dentist_%.def
	sudo singularity build $@ $<


.PHONY: conda-dentist-core
conda-dentist-core:
	conda-build $(call local_only,-c local) -c a_ludi -c bioconda \
	    $(call local_only,--clobber=conda/recipes/dentist-core/local-build.yaml) \
	    conda/recipes/dentist-core


.PHONY: conda-dentist
conda-dentist:
	conda-build $(call local_only,-c local) -c a_ludi -c bioconda -c conda-forge \
	    $(call local_only,--clobber=conda/recipes/dentist/local-build.yaml) \
	    conda/recipes/dentist


docs/list-of-commandline-options.md: $(binaries)
	$(dist_dir)/bin/dentist --list-options > $@


.PHONY: api-docs
api-docs: docs.json

docs.json:
	dub build --build=docs-json --config=testing
	$(ddox) filter $(ddox_filter) $@


.PHONY: dist
dist: $(dist_tarball)


.PHONY: prepare-dist
prepare-dist: $(dist_files)


.PHONY: clean-dist
clean-dist:
	rm -rf $(dist_dir)


.PHONY: binaries
binaries: $(binaries)


$(binaries) &:
	ROOT="$$(mktemp -d --tmpdir conda-dentist-core.XXXXXX)"; \
	PREFIX="$$ROOT/env"; \
	BINARIES=( $(addprefix "$$PREFIX/bin/",$(notdir $(binaries))) ); \
	trap 'rm -rf "$$ROOT"' exit; \
	conda create -y --copy -p "$$PREFIX" -c local -c a_ludi -c bioconda \
	    dentist-core==local \
	    jq==1.6 \
	    'python>=3,<4'; \
	echo 'Installing binaries into $(dist_dir)'; \
	install -Dt $(dist_dir)/bin "$${BINARIES[@]}"


$(dist_tarball): $(dist_files)
	@echo "Creating dist tarball ..."
	tar -czf $@ $^


$(dist_dir)/%: snakemake/%
	install -D $< $@


$(dist_dir)/%: %
	install -D $< $@
