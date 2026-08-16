F77 = gfortran
FFLAGS = -O2 -fall-intrinsics

# TSM - 06/04/2020 ###########################################################

FFINCL = -I$(MESA_DIR)/include -I$(MESA_DIR)/star/make -I$(GYRE_DIR)/src/build -finit-real=snan -fbacktrace -fmax-errors=25 -fcheck=all -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized -finline-limit=0 -ggdb -std=f2008 -fopenmp 

FFDEPS = ff_amp2.o em_lib.o em_gyre.o em_freq.o em_util.o gyre_*.o core_*.o

FFLIBS = -L$(MESA_DIR)/lib -lstar -lstar_data -lionization -latm -lcolors -lnet -leos -lkap -lrates -lneu -lchem -linterp_2d -linterp_1d -lnum -lmtx -lconst -lutils -lmath `mesasdk_lapack_link` `mesasdk_blas_link` `mesasdk_hdf5_link` `mesasdk_pgplot_link` -lz `lapack95_link` `crmath_link` `crlibm_link` `hdf5_link`

testffamp2: $(FFDEPS)
	$(F77) $(FFLAGS) $(FFINCL) -o ../../testffamp2 userff.f $(FFDEPS) $(FFLIBS)


# TSM - 08/15/2005 ###########################################################

MPIF77      = mpif77
MPIFFLAGS   = -fc=gfortran -O2
#MPIF77      = mpifort
#MPIFFLAGS   = -O2

mpi_pikaia.o: mpi_pikaia.f
	$(MPIF77) $(MPIFFLAGS) -c mpi_pikaia.f

pikaia_master.o: pikaia_master.f
	$(MPIF77) $(MPIFFLAGS) -c pikaia_master.f

mpi_fitness.o: mpi_fitness.f
	$(MPIF77) $(MPIFFLAGS) -c mpi_fitness.f

ff_slave.o: ff_slave.f
	$(MPIF77) $(MPIFFLAGS) -c ff_slave.f

pikaia_amp2: mpi_pikaia.o pikaia_master.o mpi_fitness.o ff_slave.o $(FFDEPS)
	$(MPIF77) $(MPIFFLAGS) $(FFINCL) -o ../../pikaia_amp2 mpi_pikaia.o \
	pikaia_master.o mpi_fitness.o ff_slave.o $(FFDEPS) $(FFLIBS)

pikaia_test: mpi_pikaia.o pikaia_master.o mpi_fitness.o ff_slave.o $(FFDEPS)
	$(MPIF77) $(MPIFFLAGS) $(FFINCL) -o ../../mesawork/pikaia_test mpi_pikaia.o \
	pikaia_master.o mpi_fitness.o ff_slave.o $(FFDEPS) $(FFLIBS)

