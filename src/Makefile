# Variables

TARGETS=em_test em_test_warrick testffamp2

DEBUG=yes

# MESA stuff

MODPATH=${MESA_DIR}/include ${MESA_DIR}/star/make ${GYRE_DIR}/src/build

MESA_LIB_DIR = $(MESA_DIR)/lib

LOAD_PGPLOT = `mesasdk_pgplot_link`
LOAD_LAPACK = `mesasdk_lapack_link`
LOAD_BLAS = `mesasdk_blas_link`
LOAD_HDF5 = `mesasdk_hdf5_link`
LOAD_CRLIBM = -lcrlibm
LOAD_MATRIX = -lf2crlibm $(LOAD_CRLIBM) -lmtx -lconst -lutils $(LOAD_LAPACK) $(LOAD_BLAS)
LOAD_MESA_NUMERICS = -linterp_2d -linterp_1d -lnum $(LOAD_MATRIX) $(LOAD_HDF5)
LOAD_MESA_MICRO = -lnet -leos -lkap -lrates -lneu -lchem $(LOAD_MESA_NUMERICS)
LOAD_MESA_MACRO = -lionization -latm -lcolors $(LOAD_MESA_MICRO)
LOAD_MESA_STAR_SUPPORT = $(LOAD_MESA_MACRO)
LOAD_STAR_MODS = -lstar $(LOAD_MESA_STAR_SUPPORT)
STAR_LOAD_LIST = $(LOAD_STAR_MODS) $(LOAD_PGPLOT) $(LOAD_SE)
LOAD_MESA_STAR = -L$(MESA_LIB_DIR) $(STAR_LOAD_LIST)

LDFLAGS:=${LDFLAGS} ${LOAD_MESA_STAR} -lz

# Include other makefiles

SUBDIRS_GYRE=extern/core:common:par:diff:bvp:interp:model:matrix:search:context:ext:frontend:ad:rad:nad:rot:tar:tide:mode:grid:frontend:output:poly

SUBDIRS_MESA=

SUBDIRS=${GYRE_DIR}/src/$(subst :,:${GYRE_DIR}/src/,${SUBDIRS_GYRE}):${SUBDIRS_MESA}

include Make.inc

# Additional deps

testffamp2: ff_amp2.o

ff_amp2.o: em_lib.mod

