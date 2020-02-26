# Pythia + ROOT + GSL Makefile

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Handle GZIP support.
GZIP_INC=
GZIP_FLAGS=
ifeq ($(GZIP_USE),true)
  GZIP_INC+= -DGZIPSUPPORT -I$(GZIP_INCLUDE)
  GZIP_FLAGS+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz
endif

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl 

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	@echo "Usage: make mainXX"

# The Makefile configuration.
Makefile.inc:
	$(error Error: PYTHIA must be configured, please run "./configure"\
                in the top PYTHIA directory)

# PYTHIA libraries.
$(PREFIX_LIB)/libpythia8.a :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

# *************************************************************************
# ROOT + GSL .
%: %.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(ROOT_USE),true)
	$(CXX) $< -o $@ -w -lgsl -lgslcblas -lm -I$(ROOT_INCLUDE) $(CXX_COMMON)\
	 `$(ROOTBIN)root-config --cflags`\
	 -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`
else
	@echo "Error: $@ requires ROOT"
endif
# *************************************************************************
