dentist_version=$(shell git describe)
dentist_env=dentist_$(firstword $(subst ., ,$(dentist_version)))

distdir = dist
datadir = data
scriptdir = scripts
bindir = $(distdir)/bin

assembly_reference = $(datadir)/assembly-reference.fasta.gz
assembly_reference_dam = $(datadir)/assembly-reference.dam

simulated_reads = $(distdir)/reads.fasta
simulator_flags = -m25000 -s12500 -e.13 -c20 -r19339

gaps_bed = $(datadir)/gaps.bed
assembly_test = $(distdir)/reference.fasta
mk_test_assembly_script = $(scriptdir)/mk-test-assembly.awk

generated_inputs = $(assembly_test) $(simulated_reads)
checksum_inputs = checksum.inputs.md5
inputs_valid_marker = $(datadir)/.inputs-valid

workflow_config = $(distdir)/snakemake.yml
workflow_config_orig = ../snakemake/snakemake.yml
mk_workflow_config_script = $(scriptdir)/mk-workflow-config.awk

dentist_container=$(distdir)/dentist_$(dentist_version:v%=%).sif
binaries=$(addprefix $(bindir)/,Catrack computeintrinsicqv daccord daligner \
	DAM2fasta damapper DAScover DASqv datander DB2fasta DBa2b DBb2a DBdump \
	DBdust DBmv DBrm DBshow DBsplit DBstats DBtrim DBwipe dentist dumpLA \
	fasta2DAM fasta2DB LAa2b LAb2a LAcat LAcheck LAdump LAmerge \
	lasfilteralignments LAshow LAsort LAsplit rangen simulator TANmask)
runtime_environment = $(dentist_container) $(binaries)

workflow_results = $(shell tr '\n' ' ' < result-files.lst)
checksum_workflow = $(distdir)/checksum.md5

additional_dist_files = $(addprefix $(distdir)/,dentist.yml envs/$(dentist_env).yml README.md scripts Snakefile) \
	$(wildcard $(distdir)/profile-slurm.*.yml)

dist_tarball = dentist-example.tar.gz
dist_files = $(generated_inputs) $(workflow_config) $(checksum_workflow) \
    $(runtime_environment) $(additional_dist_files)


# Tools
simulator = simulator
snakemake = $(shell command -v snakemake)
singularity = $(shell command -v singularity)
conda = $(shell command -v mamba)

conda_prefix = $(shell mamba info --base | tail -n1)/envs/$(dentist_env)


# Recipe fragments
clean_workflow = [[ -$(noclean)- == -yes- ]] || $(MAKE) clean-workflow


.PHONY: prepare-dist
prepare-dist: $(dist_files)

.PHONY: input-data
input-data: $(generated_inputs) $(inputs_valid_marker)

$(assembly_reference_dam): $(assembly_reference)
	@echo "-- making DAM for reference assembly ..."
	mkdir -p $(@D)
	[[ ! -f $@ ]] || DBrm $@
	zcat -f $< | fasta2DAM -i $@ && DBsplit -x0 -a $@

$(simulated_reads): $(assembly_reference_dam)
	@echo "-- simulating reads ..."
	mkdir -p $(@D)
	$(simulator) $(simulator_flags) $< > $@

$(assembly_test): $(mk_test_assembly_script) $(gaps_bed) $(assembly_reference)
	@echo "-- building example assembly ..."
	mkdir -p $(@D)
	awk -f $(mk_test_assembly_script) $(gaps_bed) <(zcat -f $(assembly_reference)) > $@

$(inputs_valid_marker): $(checksum_inputs) $(generated_inputs)
	@echo "-- checking generated data ..."
	mkdir -p $(@D)
	md5sum -c $< && touch $@

$(workflow_config): $(workflow_config_orig) $(mk_workflow_config_script) $(dentist_container)
	awk -v 'dentist_container=$(notdir $(dentist_container))' \
	    -v 'dentist_env=envs/$(dentist_env).yml' \
	    -f $(mk_workflow_config_script) $< > $@

$(checksum_workflow): $(workflow_results) $(generated_inputs)
	cd $(distdir) && md5sum $(^:$(distdir)/%=%) > $(@:$(distdir)/%=%)

.PHONY: $(dentist_container)
$(dentist_container):
	$(MAKE) -C .. $(@F)
	ln -srf ../$(@F) $@

.PHONY: dist
dist: $(dist_tarball)

$(dist_tarball): $(dist_files)
	mkdir -p $(@D)
	tar --transform='s|^$(distdir)/|dentist-example/|' --dereference -czf $@ $^


.PHONY: clean-workflow
clean-workflow:
	rm -rf $(distdir)/.snakemake $(distdir)/gap-closed.{agp,fasta,closed-gaps.bed} $(distdir)/logs $(distdir)/workdir


# --- Tests ---


.PHONY: check-results
check-results:
	cd $(distdir) && ../check-results.sh


singularity_tests=.passed-singularity-dentist-dependency-check \
				  .passed-singularity-snakemake-syntax-check \
				  .passed-singularity-snakemake-workflow


.PHONY: test-singularity
test-singularity: $(singularity_tests)


.passed-singularity-dentist-dependency-check: $(dentist_container)
	$(singularity) exec $(dentist_container) dentist -d
	touch $@


.passed-singularity-snakemake-syntax-check: $(dentist_container) $(workflow_config)
	cd $(distdir) && $(snakemake) --configfile=snakemake.yml --use-singularity --config dentist_container=$(notdir $(dentist_container)) -nqj1
	touch $@


.passed-singularity-snakemake-workflow: $(generated_inputs) $(workflow_config) $(dentist_container)
	+$(clean_workflow)
	cd $(distdir) && $(snakemake) --configfile=snakemake.yml --use-singularity --config dentist_container=$(notdir $(dentist_container)) -jall
	$(MAKE) check-results
	touch $@


conda_tests=.passed-conda-dentist-dependency-check \
			.passed-conda-snakemake-syntax-check \
			.passed-conda-snakemake-workflow


.PHONY: test-conda
test-conda: $(conda_tests)


.PHONY: recreate-conda-env
recreate-conda-env:
	conda env remove -n $(dentist_env)
	$(MAKE) $(conda_prefix)


$(conda_prefix): $(distdir)/envs/$(dentist_env).yml
	conda env create -n $(dentist_env) --file=$(distdir)/envs/$(dentist_env).yml


.passed-conda-dentist-dependency-check: $(conda_prefix)
	$(conda) run -n $(dentist_env) dentist -d
	touch $@


.passed-conda-snakemake-syntax-check: $(conda_prefix) $(distdir)/envs/$(dentist_env).yml $(workflow_config)
	cd $(distdir) && $(snakemake) --use-conda --config 'dentist_env=envs/$(dentist_env).yml' --conda-frontend=$(notdir $(firstword $(conda))) --configfile=snakemake.yml -nqj1
	touch $@


.passed-conda-snakemake-workflow: $(generated_inputs) $(conda_prefix) $(distdir)/envs/$(dentist_env).yml $(workflow_config)
	+$(clean_workflow)
	cd $(distdir) && $(snakemake) --use-conda --config 'dentist_env=envs/$(dentist_env).yml' --conda-frontend=conda --configfile=snakemake.yml -jall
	$(MAKE) check-results
	touch $@


binaries_tests=.passed-binaries-path-correct-check \
               .passed-binaries-dentist-dependency-check \
               .passed-binaries-snakemake-syntax-check \
               .passed-binaries-snakemake-workflow


.PHONY: test-binaries
test-binaries: $(binaries_tests)


.PHONY: binaries
binaries: $(binaries)


$(binaries) &: $(conda_prefix) $(addprefix $(conda_prefix)/bin/,$(notdir $(binaries)))
	install -D -t $(bindir) $(filter-out $(conda_prefix),$^)


.passed-binaries-path-correct-check: export PATH := $(shell pwd)/$(bindir):$(PATH)
.passed-binaries-path-correct-check: $(addprefix .passed-binaries-path-correct-check-,$(notdir $(binaries)))
	touch $@


.passed-binaries-path-correct-check-%: export PATH := $(shell pwd)/$(bindir):$(PATH)
.passed-binaries-path-correct-check-%: $(bindir)/%
	[[ "$$(which $*)" == "$(shell pwd)/$(bindir)/$*" ]]
	touch $@


.passed-binaries-dentist-dependency-check: export PATH := $(shell pwd)/$(bindir):$(PATH)
.passed-binaries-dentist-dependency-check: $(binaries)
	dentist -d
	touch $@


.passed-binaries-snakemake-syntax-check: export PATH := $(shell pwd)/$(bindir):$(PATH)
.passed-binaries-snakemake-syntax-check: $(binaries) $(workflow_config)
	cd $(distdir) && $(snakemake) --configfile=snakemake.yml -nqj1
	touch $@


.passed-binaries-snakemake-workflow: export PATH := $(shell pwd)/$(bindir):$(PATH)
.passed-binaries-snakemake-workflow: $(generated_inputs) $(binaries) $(workflow_config)
	+$(clean_workflow)
	cd $(distdir) && $(snakemake) --configfile=snakemake.yml -jall
	$(MAKE) check-results
	touch $@
