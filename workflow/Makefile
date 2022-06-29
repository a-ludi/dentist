srcdir = src
source_files = $(shell find -L $(srcdir)/ -name '*.py')

workflow.pyz: $(source_files)
	cd $(srcdir) && zip $(abspath $@) $(^:src/%=%)
