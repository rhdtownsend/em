F77 = gfortran
FFLAGS = -O2 -fall-intrinsics

# Example of how to compile/link testffamp2

PKG_DIR := ../../lib/pkgconfig

FFINCL := -fopenmp `pkgconf --with-path=$(PKG_DIR) --define-variable=math_slot=crmath --cflags em`

FFDEPS = ../../build/ff_amp2.o

FFLIBS = `pkgconf --with-path=$(PKG_DIR) --define-variable=math_slot=crmath --libs em crmath --static`

testffamp2: $(FFDEPS) ../frontend/testffamp2.f
	$(F77) $(FFLAGS) $(FFINCL) -o ../../testffamp2 $^ $(FFLIBS)

#

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

