! Module   : em_gyre
! Purpose  : GYRE routines for Embedded MESA module
!
! Copyright 2016 Rich Townsend

$include 'core.inc'

module em_gyre

  ! Uses

  use core_system
  use core_constants

  use em_freq
  use em_util

  use star_lib
  use star_def
  use const_def, only: dp

  use gyre_ad_bvp
  use gyre_bvp
  use gyre_context
  use gyre_evol_model
  use gyre_ext
  use gyre_grid
  use gyre_grid_par
  use gyre_grid_factory
  use gyre_mode
  use gyre_mode_par
  use gyre_model
  use gyre_model_par
  use gyre_num_par
  use gyre_osc_par
  use gyre_search
  use gyre_scan_par
  use gyre_rad_bvp

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: run_gyre

contains

  subroutine run_gyre (id, l, fr_obs, fr_mod)

    integer, intent(in)       :: id
    integer, intent(in)       :: l
    type(freq_t), intent(in)  :: fr_obs
    type(freq_t), intent(out) :: fr_mod

    type(star_info), pointer      :: s
    integer                       :: ierr
    real(dp)                      :: freq_min
    real(dp)                      :: freq_max
    integer                       :: n_freq
    type(model_par_t)             :: ml_p
    type(mode_par_t)              :: md_p
    type(osc_par_t)               :: os_p
    type(num_par_t)               :: nm_p
    type(grid_par_t)              :: gr_p
    type(scan_par_t), allocatable :: sc_p(:)
    type(evol_model_t), pointer   :: ml
    type(context_t), pointer      :: cx
    type(mode_t), allocatable     :: md(:)
    integer                       :: n_md
    integer                       :: i_md

    ! Get the star pointer

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failure in star_ptr)

    ! Set frequency bounds

    freq_min = fr_obs%nu(1) - 0.5_dp*s%delta_nu
    freq_max = fr_obs%nu(fr_obs%n) + 0.5_dp*s%delta_nu
    n_freq = CEILING(5._dp*(freq_max - freq_min)/s%delta_nu)

    ! Set up parameters

    call set_gyre_pars(l, freq_min, freq_max, n_freq, 'UHZ', ml_p, md_p, os_p, nm_p, gr_p, sc_p)

    ! Make the model

    call make_gyre_model(id, ml_p, ml)

    ! Find modes

    call find_gyre_modes(ml, md_p, os_p, nm_p, gr_p, sc_p, cx, md)

    ! Store data

    n_md = SIZE(md)

    fr_mod = freq_t([(REAL(md(i_md)%freq('UHZ', 'INERTIAL')), i_md=1,n_md)], &
         [(1._dp, i_md=1,n_md)], &
         [(md(i_md)%E_norm(), i_md=1,n_md)], &
         md%n_pg)

    ! Free up the context and model

    deallocate(cx)
    deallocate(ml)

    ! Finish

    return

  end subroutine run_gyre

  !****

  subroutine make_gyre_model (id, ml_p, ml)

    integer, intent(in)           :: id
    type(model_par_t), intent(in) :: ml_p
    type(evol_model_t), pointer   :: ml

    real(dp), allocatable :: global_data(:)
    real(dp), allocatable :: point_data(:,:)
    integer               :: ierr
    integer               :: n
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: V_2(:)
    real(dp), allocatable :: As(:)
    real(dp), allocatable :: U(:)
    real(dp), allocatable :: c_1(:)
    real(dp), allocatable :: Omega_rot(:)

    ! Get the global and point data

    call star_get_pulse_data(id, 'GYRE', &
         add_center_point=.TRUE., &
         keep_surface_point=.TRUE., &
         add_atmosphere=.FALSE., &
         global_data=global_data, &
         point_data=point_data, &
         ierr=ierr)
    $ASSERT(ierr == 0,Failure in star_get_pulse_data)

    ! Make the model_t

    associate ( &
         M_star => global_data(1), &
         R_star => global_data(2), &
         L_star => global_data(3), &
         r => point_data(1,:), &
         M_r => point_data(2,:), &
         L => point_data(3,:), &
         P => point_data(4,:), &
         T => point_data(5,:), &
         rho => point_data(6,:), &
         N2 => point_data(8,:), &
         Gamma_1 => point_data(9,:))

      ! Calculate dimensionless structure data

      x = r/R_star

      n = SIZE(x)

      allocate(V_2(n))
      allocate(As(n))
      allocate(U(n))
      allocate(c_1(n))
      allocate(Omega_rot(n))

      where (x /= 0._dp)
         V_2 = G_GRAVITY*M_r*rho/(P*r*x**2)
         As = r**3*N2/(G_GRAVITY*M_r)
         U = 4._dp*PI*rho*r**3/M_r
         c_1 = (r/R_star)**3/(M_r/M_star)
      elsewhere
         V_2 = 4._dp*PI*G_GRAVITY*rho(1)**2*R_star**2/(3._dp*P(1))
         As = 0._dp
         U = 3._dp
         c_1 = 3._dp*(M_star/R_star**3)/(4._dp*PI*rho)
      end where

      Omega_rot = 0._dp

      allocate(ml, SOURCE=evol_model_t(x, M_star, R_star, L_star, ml_p))

      call ml%define(I_V_2, V_2)
      call ml%define(I_AS, As)
      call ml%define(I_U, U)
      call ml%define(I_C_1, c_1)

      call ml%define(I_GAMMA_1, Gamma_1)

      call ml%define(I_OMEGA_ROT, Omega_rot)

      call ml%commit()

    end associate

    ! Finish

  end subroutine make_gyre_model

  !****

  subroutine set_gyre_pars (l, freq_min, freq_max, n_freq, freq_units, ml_p, md_p, os_p, nm_p, gr_p, sc_p)

    integer, intent(in)                        :: l
    real(dp), intent(in)                       :: freq_min
    real(dp), intent(in)                       :: freq_max
    integer, intent(in)                        :: n_freq
    character(*), intent(in)                   :: freq_units
    type(model_par_t), intent(out)             :: ml_p
    type(mode_par_t), intent(out)              :: md_p
    type(osc_par_t), intent(out)               :: os_p
    type(num_par_t), intent(out)               :: nm_p
    type(grid_par_t), intent(out)              :: gr_p
    type(scan_par_t), allocatable, intent(out) :: sc_p(:)

    ! Set up parameters for a GYRE run

    ! Model

    ml_p = model_par_t(add_center=.FALSE., &
                       repair_As=.TRUE.)

    ! Mode

    md_p = mode_par_t(i=1, &
                      l=l, &
                      m=0)

    ! Oscillation

    os_p = osc_par_t(x_ref=1._dp, &
                     outer_bound='JCD', &
                     inertia_norm='RADIAL', &
                     eddington_approx=.TRUE.)

    ! Numerical

    nm_p = num_par_t(diff_scheme='COLLOC_GL4')

    ! Grid (the first setting is preferable, but can lead to huge
    ! models. So, go with the second one for now)

!    gr_p = grid_par_t(alpha_osc=10._dp, &
!                      alpha_exp=2._dp, &
!                      n_inner=5)

    gr_p = grid_par_t(n_inner=5)

    ! Scan

    allocate(sc_p(1))

    sc_p(1) = scan_par_t(freq_min=freq_min, &
                         freq_max=freq_max, &
                         n_freq=n_freq, &
                         freq_min_units=freq_units, &
                         freq_max_units=freq_units, &
                         freq_frame='INERTIAL', &
                         grid_type='LINEAR', &
                         grid_frame='INERTIAL')

    ! Finish

    return

  end subroutine set_gyre_pars

  !****

  subroutine find_gyre_modes (ml, md_p, os_p, nm_p, gr_p, sc_p, cx, md)

    type(evol_model_t), pointer, intent(in) :: ml
    type(mode_par_t), intent(in)            :: md_p
    type(osc_par_t), intent(in)             :: os_p
    type(num_par_t), intent(in)             :: nm_p
    type(grid_par_t), intent(in)            :: gr_p
    type(scan_par_t), intent(in)            :: sc_p(:)
    type(context_t), pointer, intent(out)   :: cx
    type(mode_t), allocatable, intent(out)  :: md(:)

    real(dp), allocatable       :: omega(:)
    type(grid_t)                :: gr
    real(dp)                    :: omega_min
    real(dp)                    :: omega_max
    class(r_bvp_t), allocatable :: bp
    integer                     :: n_md
    integer                     :: d_md

    ! Set up the context

    allocate(cx)
    cx = context_t(ml, gr_p, md_p, os_p)

    ! Set up the frequency array

    call build_scan(cx, md_p, os_p, sc_p, omega)

    ! Create the grid

    gr = grid_t(cx, omega, gr_p)

    if (gr%n_k > 25000) stop 'Grid exceeded hard-wired size limit of 25,000 points'

    ! Set frequency bounds and perform checks

    if (nm_p%restrict_roots) then
       omega_min = MINVAL(omega)
       omega_max = MAXVAL(omega)
    else
       omega_min = -HUGE(0._dp)
       omega_max = HUGE(0._dp)
    endif

    call check_scan(ml, gr, omega, md_p, os_p)

    ! Create the bvp_t

    if (md_p%l == 0 .AND. os_p%reduce_order) then
       allocate(bp, SOURCE=rad_bvp_t(cx, gr, md_p, nm_p, os_p))
    else
       allocate(bp, SOURCE=ad_bvp_t(cx, gr, md_p, nm_p, os_p))
    endif

    ! Find modes

    n_md = 0
    d_md = 16

    allocate(md(d_md))

    call scan_search(bp, omega, omega_min, omega_max, process_mode, nm_p)

    ! Resize the md array just to the modes found

    md = md(:n_md)

    ! Finish (this explicit deallocation shouldn't be necessary, but
    ! without it valgrind complains)

    deallocate(bp)

    return

  contains

    subroutine process_mode (md_new, n_iter, chi)

      type(mode_t), intent(in)  :: md_new
      integer, intent(in)       :: n_iter
      type(r_ext_t), intent(in) :: chi

      ! Process the mode

      if (md_new%n_pg < md_p%n_pg_min .OR. md_new%n_pg > md_p%n_pg_max) return

      ! Store it

      n_md = n_md + 1

      if (n_md > d_md) then
         d_md = 2*d_md
         call reallocate(md, [d_md])
      endif

      md(n_md) = md_new

      !call md(n_md)%prune()

      ! Finish

      return

    end subroutine process_mode

  end subroutine find_gyre_modes

end module em_gyre
