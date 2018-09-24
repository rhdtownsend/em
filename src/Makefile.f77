F77 = gfortran
FFLAGS = -O2

# TSM - 11/17/2017 ###########################################################

FFINCL = -I$(MESA_DIR)/include -I$(GYRE_DIR)/src/build -I$(MESASDK_ROOT)/include  -finit-real=snan -ffpe-trap=invalid,overflow,zero -fbacktrace -fmax-errors=25 -O2 -march=native -fopenmp

FFDEPS = ff_amp2.o gyre_scan_par.o gyre_r_bvp.o gyre_sysmtx.o gyre_c_ext_func.o gyre_c_deriv.o gyre_nad_eqns.o gyre_deriv.o gyre_r_ext_func.o gyre_band_sysmtx.o core_linalg.o gyre_r_root.o gyre_r_diff.o gyre_c_interp.o gyre_bvp.o gyre_evol_model.o gyre_grid_spec.o gyre_c_block_sysmtx.o em_lib.o gyre_r_band_sysmtx.o gyre_c_sysmtx.o gyre_model_par.o gyre_r_sysmtx.o gyre_point.o gyre_grid.o gyre_c_trapz_diff.o gyre_rad_bound.o gyre_c_ext.o gyre_rad_eqns.o gyre_r_discrim_func.o gyre_r_magnus_diff.o gyre_r_eqns.o gyre_c_eqns.o gyre_nad_trans.o gyre_r_ext.o core_system.o gyre_rad_match.o gyre_r_state.o gyre_ad_trans.o gyre_discrim_func.o gyre_status.o gyre_mode.o gyre_c_diff.o gyre_ad_eqns.o gyre_interp.o gyre_freq.o gyre_diff.o gyre_dopp_rot.o gyre_c_discrim_func.o gyre_util.o gyre_ext.o gyre_c_root.o gyre_rad_trans.o em_util.o gyre_colloc_diff.o gyre_ext_func.o gyre_c_bound.o gyre_osc_par.o gyre_c_mirk_diff.o gyre_qad_eval.o gyre_cimplex.o gyre_wave.o gyre_atmos.o gyre_nad_bvp.o gyre_trapz_diff.o gyre_c_dopp_rot.o gyre_r_dopp_rot.o gyre_mode_par.o gyre_r_rot.o gyre_rot.o gyre_r_trapz_diff.o gyre_r_deriv.o gyre_c_search.o gyre_num_par.o core_memory.o gyre_context.o gyre_c_magnus_diff.o gyre_r_search.o gyre_r_colloc_diff.o gyre_nad_diff.o gyre_block_sysmtx.o gyre_rad_bvp.o gyre_state.o gyre_grid_util.o gyre_hom_model.o em_freq.o gyre_grid_factory.o gyre_rot_factory.o gyre_ad_bound.o gyre_mirk_diff.o gyre_nad_match.o gyre_ad_match.o gyre_rad_diff.o core_parallel.o gyre_constants.o gyre_twopt_model.o gyre_eqns.o gyre_c_band_sysmtx.o gyre_r_bound.o gyre_grid_weights.o gyre_r_interp.o core_order.o core_func.o gyre_ad_diff.o gyre_model.o gyre_diff_factory.o gyre_poly_model.o gyre_c_bvp.o gyre_bound.o em_gyre.o gyre_r_mirk_diff.o gyre_c_state.o gyre_search.o gyre_r_block_sysmtx.o gyre_ad_bvp.o gyre_nad_bound.o gyre_c_rot.o gyre_magnus_diff.o gyre_c_colloc_diff.o gyre_grid_par.o gyre_sysmtx_factory.o gyre_model_util.o gyre_linalg.o core_kinds.o gyre_root.o core_constants.o


FFLIBS = -L$(MESA_DIR)/lib -lstar -lionization -latm -lcolors -lnet -leos -lkap -lrates -lneu -lchem -lmlt -linterp_2d -linterp_1d -lnum -lf2crlibm -lcrlibm -lmtx -lconst -lutils `mesasdk_lapack_link` `mesasdk_blas_link` `mesasdk_hdf5_link` `mesasdk_pgplot_link`  -lz `mesasdk_lapack95_link`

testffamp2: $(FFDEPS)
	$(F77) $(FFLAGS) $(FFINCL) -o ../../mesawork/testffamp2 userff.f $(FFDEPS) $(FFLIBS)


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
	$(MPIF77) $(MPIFFLAGS) $(FFINCL) -o ../../mesawork/pikaia_amp2 mpi_pikaia.o \
	pikaia_master.o mpi_fitness.o ff_slave.o $(FFDEPS) $(FFLIBS)

pikaia_test: mpi_pikaia.o pikaia_master.o mpi_fitness.o ff_slave.o $(FFDEPS)
	$(MPIF77) $(MPIFFLAGS) $(FFINCL) -o ../../mesawork/pikaia_test mpi_pikaia.o \
	pikaia_master.o mpi_fitness.o ff_slave.o $(FFDEPS) $(FFLIBS)

