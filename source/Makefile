# TASMAN Makefile
#
# Examples:
#
#   make
#   make FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"
#   make FC=ifx FFLAGS="-O3"
#   make clean
#   make check

# Compiler-specific flags can be supplied by the user.
ifeq ($(origin FC), default)
  FC = gfortran
endif

ifeq ($(origin FFLAGS), undefined)
  ifneq ($(findstring gfortran,$(FC)),)
    FFLAGS = -w -O3 -ffp-contract=off
  else
    FFLAGS =
  endif
endif

# Installation directories.

TASMAN_DIR := $(abspath ..)
BIN_DIR    := $(TASMAN_DIR)/bin
TARGET     := $(BIN_DIR)/tasman

# Sibling TALYS and TEFAL installations used by the sample cases.

TALYS_DIR_CHECK := $(abspath $(TASMAN_DIR)/../talys)
TEFAL_DIR_CHECK := $(abspath $(TASMAN_DIR)/../tefal)

# TASMAN source files.
#
# Free-form source files must be named so that alphabetical order places
# module definitions, such as A0_tasman_mod.f90 and A1_error_handling_mod.f90,
# before their users.

SOURCES := $(sort $(wildcard *.f) $(wildcard *.f90))


.PHONY: all install clean check

all install: $(TARGET)

$(TARGET): $(SOURCES)
	@echo "Compiler: $(FC)"
	@echo "Flags:    $(FFLAGS)"
	@mkdir -p "$(BIN_DIR)"
	$(FC) $(FFLAGS) -o "$(TARGET)" $(SOURCES)
	@rm -f *.o *.mod
	@echo "Installed: $(TARGET)"

clean:
	rm -f *.o *.mod tasman "$(TARGET)"

check: all
	cd "$(TASMAN_DIR)/samples" && \
	  TASMAN_DIR="$(TASMAN_DIR)" \
	  TALYS_DIR="$(TALYS_DIR_CHECK)" \
	  TEFAL_DIR="$(TEFAL_DIR_CHECK)" \
	  PATH="$(TALYS_DIR_CHECK)/bin:$(TEFAL_DIR_CHECK)/bin:$$PATH" \
	  ./verify
