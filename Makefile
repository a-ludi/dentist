dentist_readme = $(DENTIST)/README.md
dentist_license = $(DENTIST)/LICENSE
dentist_list_of_commandline_options = $(DENTIST)/docs/list-of-commandline-options.md

api_docs = $(wildcard docs.*.json)
api_versions = $(patsubst docs.%.json,%,$(api_docs))
api_htmls = $(addprefix api/,$(addsuffix /index.html,$(api_versions)))

site_files = index.md license.md list-of-commandline-options.md

ddox_generate_flags = --html-style=pretty

ddox = dub run -q ddox --
ddox_root = $(shell dub describe ddox | jq -r '.packages[] | select(.name == "ddox") | .path')

.PHONY: all
all: $(site_files) $(api_htmls)


index.md: $(dentist_readme)
license.md: $(dentist_license)
list-of-commandline-options.md: $(dentist_list_of_commandline_options)

%.md: scripts/make-%.sh scripts/include.rc
	bash scripts/make-$*.sh $(filter-out scripts/%,$^) > $@


.PHONY: api-docs
api-docs: $(api_htmls)


api/%/index.html: docs.%.json
	$(ddox) generate-html $(ddox_generate_flags) $< $(@D)
	cp -r -t $(@D) $(ddox_root)/public/*
