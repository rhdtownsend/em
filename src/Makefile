# Variables

TARGETS=em_test

DEBUG=yes

# MESA stuff

MODPATH=${MESA_DIR}/include ${GYRE_DIR}/src/build

MESA_LIB_DIR = $(MESA_DIR)/lib

LOAD_PGPLOT = `mesasdk_pgplot_link`
LOAD_LAPACK = `mesasdk_lapack_link`
LOAD_BLAS = `mesasdk_blas_link`
LOAD_HDF5 = `mesasdk_hdf5_link`
LOAD_CRLIBM = -lcrlibm
LOAD_MATRIX = -lf2crlibm $(LOAD_CRLIBM) -lmtx -lconst -lutils -lmesaklu $(LOAD_LAPACK) $(LOAD_BLAS)
LOAD_MESA_NUMERICS = -linterp_2d -linterp_1d -lnum $(LOAD_MATRIX) $(LOAD_HDF5)
LOAD_MESA_MICRO = -lnet -leos -lkap -lrates -lneu -lchem $(LOAD_MESA_NUMERICS)
LOAD_MESA_MACRO = -lionization -latm -lcolors -lmlt $(LOAD_MESA_MICRO)
LOAD_MESA_STAR_SUPPORT = $(LOAD_MESA_MACRO)
LOAD_STAR_MODS = -lstar $(LOAD_MESA_STAR_SUPPORT)
STAR_LOAD_LIST = $(LOAD_STAR_MODS) $(LOAD_PGPLOT) $(LOAD_SE)
LOAD_MESA_STAR = -L$(MESA_LIB_DIR) $(STAR_LOAD_LIST)

LDFLAGS:=${LDFLAGS} ${LOAD_MESA_STAR} -lz

# Include other makefiles

SUBDIRS=${GYRE_DIR}/src/extern/core:${GYRE_DIR}/src/extern/astro:${GYRE_DIR}/src/common:${GYRE_DIR}/src/diff:${GYRE_DIR}/src/bvp:${GYRE_DIR}/src/interp:${GYRE_DIR}/src/model:${GYRE_DIR}/src/matrix:${GYRE_DIR}/src/search:${GYRE_DIR}/src/ext:${GYRE_DIR}/src/frontend:${GYRE_DIR}/src/ad:${GYRE_DIR}/src/rad:${GYRE_DIR}/src/nad:${GYRE_DIR}/src/rot:${GYRE_DIR}/src/tar:${GYRE_DIR}/src/mode:${GYRE_DIR}/src/soln:${GYRE_DIR}/src/grid:${GYRE_DIR}/src/frontend:${GYRE_DIR}/src/output:${GYRE_DIR}/src/poly:${MESA_DIR}/star/job

include Make.inc
