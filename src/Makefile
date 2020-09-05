# Variables

TARGETS=em_test testffamp2

DEBUG=yes
FPE=no

export MESASDK_MATH_SLOT=crmath

# MESA stuff

MODPATH=${MESA_DIR}/include ${MESA_DIR}/star/make ${GYRE_DIR}/src/build

MESA_LIB_DIR = $(MESA_DIR)/lib

LOAD_PGPLOT = `mesasdk_pgplot_link`
LOAD_LAPACK = `mesasdk_lapack_link`
LOAD_BLAS = `mesasdk_blas_link`
LOAD_HDF5 = `mesasdk_hdf5_link`
LOAD_CRMATH = `mesasdk_crmath_link`
LOAD_MATRIX = $(LOAD_CRMATH) -lmtx -lconst -lutils -lmath $(LOAD_LAPACK) $(LOAD_BLAS) 
LOAD_MESA_NUMERICS = -linterp_2d -linterp_1d -lnum $(LOAD_MATRIX) $(LOAD_HDF5)
LOAD_MESA_MICRO = -lnet -leos -lkap -lrates -lneu -lchem $(LOAD_MESA_NUMERICS)
LOAD_MESA_MACRO = -lionization -latm -lcolors $(LOAD_MESA_MICRO)
LOAD_MESA_STAR_SUPPORT = $(LOAD_MESA_MACRO)
LOAD_STAR_MODS = -lstar -lstar_data $(LOAD_MESA_STAR_SUPPORT)
STAR_LOAD_LIST = $(LOAD_STAR_MODS) $(LOAD_PGPLOT) $(LOAD_SE)
LOAD_MESA_STAR = -L$(MESA_LIB_DIR) $(STAR_LOAD_LIST)

LDFLAGS:=${LDFLAGS} ${LOAD_MESA_STAR} -lz

# Include other makefiles

SUBDIRS_GYRE=extern/core:common:math:par:diff:bvp:interp:model:matrix:search:context:ext:frontend:ad:rad:nad:rot:tar:tide:mode:grid:frontend:output:poly
SUBDIRS_MESA=

SUBDIRS=${GYRE_DIR}/src/$(subst :,:${GYRE_DIR}/src/,${SUBDIRS_GYRE}):${SUBDIRS_MESA}

include ${GYRE_DIR}/src/build/Make.inc

# Additional rules

vpath %.mod ${MODPATH}

ifeq ($(CRMATH),yes)
gyre_math.fpp: gyre_math_crmath.fpp
	@ln -sf $< $@
else
gyre_math.fpp: gyre_math_intrinsic.fpp
	@ln -sf $< $@
endif

.INTERMEDIATE: gyre_math.fpp

.depend: gyre_math.fpp

# Additional deps

testffamp2: ff_amp2.o userff.o em_lib.o em_gyre.o em_freq.o em_util.o
	${F9XC} ${F9XFLAGS} -o ../../$@ $^ gyre_*.o core_*.o ${LDFLAGS}

ff_amp2.o: em_lib.mod

