# quick n dirty makefile to build rust binary
# based off the output of `cargo build --release --verbose` and `cargo build --verbose`

SHELL := /bin/bash

PROJECT := linlexico
ROOT := $(shell realpath .)
COMMON_FLAGS = main.rs --edition=2018 --crate-type bin
#COMMON_FLAGS = simple.rs --edition=2018 --crate-type bin
RUSTC = rustc $(COMMON_FLAGS) --emit=dep-info,link

# release build
$(PROJECT):
	$(RUSTC) --crate-name $@ -C opt-level=3
-include $(PROJECT).d

# debug build
$(PROJECT)_dbg:
	$(RUSTC) --crate-name $@ -C debuginfo=2
-include $(PROJECT)_dbg.d

# Generate Documentation
# (Try opening doc/$(PROJECT)/index.html in a browser)
.PHONY: doc
doc:
	rustdoc $(COMMON_FLAGS) --crate-name $(PROJECT) --document-private-items

# Benchmark
# (Prefer hyperfine if it is present, but fall back to time)
.PHONY: bench
ifeq (, $(shell which hyperfine))
bench: $(PROJECT)
	time ./$(PROJECT)
else
bench: $(PROJECT)
	@# Hyperfine needs the command to be in quotes
	hyperfine './$(PROJECT)'
endif

# build everything
all: $(PROJECT) $(PROJECT)_dbg doc

# professor-proofing the makefile by adding aliases
.PHONY: docs build ball build_dbg debug dbg benchmark
docs: doc
build: $(PROJECT)
ball: $(PROJECT)
build_dbg: $(PROJECT)_dbg
debug: $(PROJECT)_dbg
dbg: $(PROJECT)_dbg
benchmark: bench


# Clean build dir
.PHONY: clean
clean:
	rm -rf $(PROJECT) $(PROJECT).d $(PROJECT)_dbg $(PROJECT)_dbg.d doc $(PROJECT)*.o
