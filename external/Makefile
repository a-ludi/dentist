# Compiler configuration
CXXFLAGS=-O2 -std=c++11
LDLIBS=-lsdsl -ldivsufsort -ldivsufsort64

# Install paths
PREFIX=/usr/local
EXEC_PREFIX=$(PREFIX)
BINDIR=$(EXEC_PREFIX)/bin

# Binaries used by this Makefile
INSTALL_PROGRAM=install


all: fm-index

.PHONY: install install-strip clean
install: fm-index
	$(INSTALL_PROGRAM) $^ $(BINDIR)

install-strip:
	$(MAKE) INSTALL_PROGRAM='$(INSTALL_PROGRAM) -s' \
	        install

clean:
	$(RM) fm-index
