srcdir = src
srcdir_dwe = external/dentist-workflow-engine/src
source_files = $(shell find -L $(srcdir)/ -name '*.py')
source_files_dwe = $(shell find -L $(srcdir_dwe)/ -name '*.py')

distdir = dist
dist_files = $(source_files:$(srcdir)/%=$(distdir)/%) $(source_files_dwe:$(srcdir_dwe)/%=$(distdir)/%)

shebang = \#!/usr/bin/env python3


workflow.pyz: $(dist_files)
	cd $(distdir) \
	&& zip $@ $(^:$(distdir)/%=%) \
	&& echo '$(shebang)' \
	   | cat - $@ > $(abspath $@) \
	&& chmod +x $(abspath $@)
	rm -f $(distdir)/$@


$(distdir)/%: $(srcdir)/%
	install -D $< $(<:$(srcdir)/%=$(distdir)/%)


$(distdir)/%: $(srcdir_dwe)/%
	install -D $< $(<:$(srcdir_dwe)/%=$(distdir)/%)
