# File     : Makefile
# Purpose  : makefile for build subdir

# Default rule in case make is invoked in current dir

all:
	make -C ..

# Targets

LIB_TARGETS = libem

PKG_FILES = em.pc
INC_FILES = em_m.mod em.inc

ifeq ($(FRONTENDS),yes)
   BIN_TARGETS += em_test testffamp2
endif

# Source file lists & related variables
#
# xxx_SRCS = source files required to compile xxx
# xxx_USES = source files providing modules USE'd by xxx
# xxx_LIBS = library files required to link xxx

em_test_SRCS = em_test.fypp
em_test_USES = $(libem_SRCS)
em_test_LIBS = libem

testffamp2_SRCS = testffamp2.f ff_amp2.fypp
testffamp2_USES = $(libem_SRCS)
testffamp2_LIBS = libem

libem_SRCS = em_m.fypp em_mesa_m.fypp em_gyre_m.fypp freq_m.fypp \
             interp_m.fypp deriv_m.fypp utils_m.fypp

# Build flags

PKGS := gyre mesa lapack95

PKGCONF_FLAGS := --with-path=$(PKG_DIR)

ifeq ($(CRMATH),yes)
   PKGS += crmath
   MATH_SLOT = crmath
else
   MATH_SLOT = default
endif

PKGCONF_FLAGS += --define-variable=math_slot=$(MATH_SLOT)

FFLAGS := `pkgconf $(PKGCONF_FLAGS) --cflags $(PKGS)`
FPPFLAGS := `pkgconf $(PKGCONF_FLAGS) --cflags-only-I $(PKGS)`
LDFLAGS := `pkgconf $(PKGCONF_FLAGS) --libs $(PKGS)`

unexport MESASDK_MATH_SLOT

# Include standard rules

include Make.inc

# Target for linting

linter :
	fortitude check --file-extensions=f90

# Targets for source checking

check_src:
	@for src in $(ALL_SRCS); do \
             find $(SRC_DIRS) -path '*/forum/*' -prune -or -name $$src -print | grep -q  . || echo $$src MISSING; \
         done

check_src : SRC_DIRS ?= ../src

.PHONY: all linter check_src print_link
