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
   BIN_TARGETS += em_test
endif

# Source file lists & related variables
#
# xxx_SRCS = source files required to compile xxx
# xxx_USES = source files providing modules USE'd by xxx
# xxx_LIBS = library files required to link xxx

em_test_SRCS = em_test.fypp
em_test_USES = $(libem_SRCS)
em_test_LIBS = libem

libem_SRCS = em_m.fypp em_gyre_m.fypp freq_m.fypp util_m.fypp

# MESA stuff

MODPATH=${MESA_DIR}/include

MESA_LIB_DIR = $(MESA_DIR)/lib

LOAD_PGPLOT = `mesasdk_pgplot_link`
LOAD_LAPACK = `mesasdk_lapack_link`
LOAD_BLAS = `mesasdk_blas_link`
LOAD_HDF5 = `mesasdk_hdf5_link`
LOAD_CRMATH = `mesasdk_crmath_link`

LOAD_MATRIX = $(LOAD_CRMATH) -lmtx -lconst -lutils -lmath $(LOAD_LAPACK) $(LOAD_BLAS)
LOAD_MESA_NUMERICS = -linterp_2d -linterp_1d -lnum -lauto_diff -lforum $(LOAD_MATRIX) $(LOAD_HDF5)
LOAD_MESA_MICRO = -lnet -leos -lkap -lrates -lneu -lchem $(LOAD_MESA_NUMERICS)
LOAD_MESA_MACRO = -lionization -latm -lcolors -lturb $(LOAD_MESA_MICRO)
LOAD_MESA_STAR_SUPPORT = $(LOAD_MESA_MACRO)
LOAD_STAR_MODS = -lstar -lstar_data $(LOAD_MESA_STAR_SUPPORT)
STAR_LOAD_LIST = $(LOAD_STAR_MODS) $(LOAD_PGPLOT) $(LOAD_SE)
LOAD_MESA_STAR = -L$(MESA_LIB_DIR) $(STAR_LOAD_LIST)

# Build flags

PKGS := gyre lapack95

PKGCONF_FLAGS := --with-path=$(GYRE_DIR)/lib/pkgconfig

ifeq ($(CRMATH),yes)
   PKGS += crmath
   MATH_SLOT = crmath
else
   MATH_SLOT = default
endif

PKGCONF_FLAGS += --define-variable=math_slot=$(MATH_SLOT)

FFLAGS := `pkgconf $(PKGCONF_FLAGS) --cflags $(PKGS)` -I$(MESA_DIR)/include -I$(MESA_DIR)/star/make
FPPFLAGS := `pkgconf $(PKGCONF_FLAGS) --cflags-only-I $(PKGS)` -m os
LDFLAGS := `pkgconf $(PKGCONF_FLAGS) --libs $(PKGS)` $(LOAD_MESA_STAR)

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

.PHONY: all linter check_src
