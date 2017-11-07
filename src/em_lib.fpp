
  ! Module   : em_lib
! Purpose  : main Embedded MESA module
!
! Copyright 2016 Rich Townsend

$include 'core.inc'

module em_lib

  ! Uses

  use core_system
  use core_constants

  use em_gyre
  use em_freq

  use star_lib
  use star_def
  use const_def, only: dp
  use atm_lib
  use eos_lib

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Module variables

  type(freq_t), save :: fr_obs_m(0:3) ! Observed frequencies
  type(freq_t), save :: fr_mod_m(0:3) ! Model frequencies
  type(freq_t), save :: fr_cor_m(0:3) ! Model frequencies

  character(64), save :: state_m     ! State-machine state
  real(dp), save      :: f_enter_m   ! Threshold for entering GYRE calcs
  real(dp), save      :: f_exit_m    ! Threshold for exiting GYRE calcs
  real(dp), save      :: y_tol_m     ! Frequency convergence tolerance
  real(dp), save      :: dt_frac_m   ! Timestep limit fraction

  real(dp), save :: Teff_m ! Effective temperature of best model
  real(dp), save :: logg_m ! Surface gravity of best model
  real(dp), save :: FeH_m  ! Metallicity of best model
  real(dp), save :: R_m    ! Radius of best model
  real(dp), save :: L_m    ! Luminosity of best model
  real(dp), save :: age_m  ! Age of best model
  real(dp), save :: a3     ! Coefficient of cubic term
  real(dp), save :: a1     ! Coefficient of inverse term

  ! Access specifiers

  private

  public :: dp
  public :: freq_t
  public :: init
  public :: create_star
  public :: evolve_star_to_zams
  public :: evolve_star_to_log_g
  public :: evolve_star_seismic
  public :: destroy_star
  public :: set_obs_freqs
  public :: clear_obs_freqs
  public :: get_mod_freqs
  public :: get_cor_freqs
  public :: get_mod_data
  public :: apply_cubic_correction
  public :: apply_combined_correction

contains

  subroutine init ()

    ! Initialize

    call clear_obs_freqs()

    ! Finish

    return

  end subroutine init

  !*****

  function create_star (M, Y, Z, alpha, f_ov) result (id)

    real(dp), intent(in) :: M
    real(dp), intent(in) :: Y
    real(dp), intent(in) :: Z
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: f_ov
    integer              :: id

    integer                  :: ierr
    type(star_info), pointer :: s
    real(dp)                 :: f0_ov_div_f_ov
    integer                  :: which_atm

    ! Create a star with the specified initial parameters
    ! mixing-length parameter

    ! Allocate the star
 
    id = alloc_star(ierr)
    $ASSERT(ierr == 0,Failure in alloc_star)

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failed in star_ptr)

    ! Set up job parameters (by passing an empty filename, no inlist
    ! will be read, but defaults are set)

    call read_star_job(s, '', ierr)
    $ASSERT(ierr == 0,Failed in read_star_job)

    ! Initialize starlib (okay to do extra calls on this)
    
    call starlib_init(s, ierr)
    $ASSERT(ierr == 0,Failed in starlib_init)

    ! Set up handles

    call star_set_kap_and_eos_handles(id, ierr)
    $ASSERT(ierr == 0,Failed in star_set_kap_and_eos_handles)

    ! Set up controls (by passing an empty filename, no inlist will be
    ! read, but defaults are set)

    call star_setup(id, '', ierr)
    $ASSERT(ierr == 0,Failed in star_steup)

    s%initial_mass = M
    s%initial_y = Y
    s%initial_Z = Z
    s%mixing_length_alpha = alpha
    s%max_age = 15.d9

    s%num_cells_for_smooth_brunt_B = 0
    s%add_double_points_to_pulse_data = .TRUE.
    s%threshold_grad_mu_for_double_point = 1d0
    s%max_number_of_double_points = 100

!    s%varcontrol_target = 1D-4
!    s%mesh_delta_coeff = 0.5D0

    ! (Note: this mirrors the code in star/astero/src/extras_support)

    f0_ov_div_f_ov = 1._dp

    s% overshoot_f_above_nonburn_core = f_ov
    s% overshoot_f_above_nonburn_shell = f_ov
    s% overshoot_f_below_nonburn_shell = f_ov
    s% overshoot_f_above_burn_h_core = f_ov
    s% overshoot_f_above_burn_h_shell = f_ov
    s% overshoot_f_below_burn_h_shell = f_ov
    s% overshoot_f_above_burn_he_core = f_ov
    s% overshoot_f_above_burn_he_shell = f_ov
    s% overshoot_f_below_burn_he_shell = f_ov
    s% overshoot_f_above_burn_z_core = f_ov
    s% overshoot_f_above_burn_z_shell = f_ov
    s% overshoot_f_below_burn_z_shell = f_ov
    s% overshoot_f0_above_nonburn_core = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_above_nonburn_shell = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_below_nonburn_shell = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_above_burn_h_core = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_above_burn_h_shell = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_below_burn_h_shell = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_above_burn_he_core = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_above_burn_he_shell = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_below_burn_he_shell = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_above_burn_z_core = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_above_burn_z_shell = f0_ov_div_f_ov*f_ov
    s% overshoot_f0_below_burn_z_shell = f0_ov_div_f_ov*f_ov
    
    ! Output controls

    s%write_profiles_flag = .FALSE.
    s%do_history_file = .FALSE.
    s%photo_interval = 0

    s%log_directory = '.'
    s%photo_directory = '.'

    s%terminal_interval = 0

    ! Set up procedure hooks

    s%extras_check_model => null()
    
    ! Set up the optical depth at the atmosphere base (this appears
    ! here because we don't use do_star_job_controls_before)

    which_atm = atm_option(s%which_atm_option, ierr)
    $ASSERT(ierr == 0,Failed in atm_option)
    
    s%tau_base = atm_tau_base(which_atm, ierr)
    $ASSERT(ierr == 0,Failed in atm_tau_base)

    ! Set up EOS blending

    call eos_set_HELM_OPAL_Zs(s%eos_handle, 1d0, 1d0, ierr)
    $ASSERT(ierr == 0,Failed in eos_set_HELM_OPAL_Zs)

    ! Create the star as a pre-MS object (cannot start with a ZAMS
    ! model, because we use non-standard helium abundances)

    call star_create_pre_ms_model(id, &
         s%job%pre_ms_T_c, &
         s%job%pre_ms_guess_rho_c, &
         s%job%pre_ms_d_log10_P, &
         s%job%pre_ms_logT_surf_limit, &
         s%job%pre_ms_logP_surf_limit, &
         s%job%initial_zfracs, &
         s%job%dump_missing_metals_into_heaviest, &
         s%job%change_net .OR. s%job%change_initial_net, &
         s%job%new_net_name, &
         s%job%pre_ms_relax_num_steps, ierr)
    $ASSERT(ierr == 0,Failed in star_create_pre_ms_model)

    ! Finish

    return

  end function create_star

  !****

  subroutine evolve_star (id)

    integer, intent(in) :: id

    logical, parameter :: DEBUG = .FALSE.

    type(star_info), pointer :: s
    integer                  :: ierr
    logical                  :: first_try
    logical                  :: just_did_backup
    integer                  :: result
    
    ! Run the evolve loop

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failure in star_ptr)

    evolve_loop: do

       first_try = .TRUE.
       just_did_backup = .FALSE.

       step_loop: do

          ! Take the step
            
          result = star_evolve_step(id, first_try, just_did_backup)

          ! Check the model
          
          if (result == keep_going) result = star_check_model(id)
          if (result == keep_going .AND. &
               ASSOCIATED(s%extras_check_model)) result = s%extras_check_model(id, 0)

          ! Pick the next timestep
          
          if (result == keep_going) result = star_pick_next_timestep(id)            
          if (result == keep_going) exit step_loop

          ! See what went wrong
               
          if (result == redo) then
             result = star_prepare_to_redo(id)
          end if
          if (result == retry) then
             result = star_prepare_to_retry(id)
          end if
          if (result == backup) then
             result = star_do1_backup(id)
             just_did_backup = .TRUE.
          else
             just_did_backup = .FALSE.
          end if
          if (result == terminate) then
             exit step_loop
          end if

          first_try = .FALSE.
               
       end do step_loop

       if (DEBUG) then
          write(OUTPUT_UNIT, *) 'Evolved:', s%star_age, s%time_step, s%why_TLim, result
       endif

       ! Once we get here, the only options are keep_going or terminate.
       ! redo, retry, or backup must be done inside the step_loop

       if (result == keep_going) then
          result = star_finish_step(id, 0, .FALSE., ierr)
          $ASSERT(ierr == 0,Failed in star_finish_step)
       endif

       if (result == terminate) then
          result = star_finish_step(id, 0, .FALSE., ierr)
          $ASSERT(ierr == 0,Failed in star_finish_step)
          exit evolve_loop
       endif

    end do evolve_loop

    ! Finish

    return

  end subroutine evolve_star

  !****
  
  subroutine evolve_star_to_zams (id)

    integer, intent(in) :: id

    type(star_info), pointer :: s
    integer                  :: ierr
    logical                  :: stp_save

    ! Evolve the star to the ZAMS

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failure in star_ptr)

    ! Set controls

    stp_save = s%stop_near_zams
    s%stop_near_zams = .TRUE.

    ! Evolve the star

    call evolve_star(id)

    ! Reset controls

    s%stop_near_zams = stp_save

    ! Finish

    return

  end subroutine evolve_star_to_zams

  !****

  subroutine evolve_star_to_log_g (id, log_g)

    integer, intent(in)  :: id
    real(dp), intent(in) :: log_g

    type(star_info), pointer :: s
    integer                  :: ierr
    real(dp)                 :: log_g_save
    
    ! Evolve the star until the surface gravity drops below log_g

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failure in star_ptr)

    ! Set controls

    log_g_save = s%log_g_lower_limit
    s%log_g_lower_limit = log_g

    ! Evolve

    call evolve_star(id)

    ! Reset controls

    s%log_g_lower_limit = log_g_save

    ! Finish

    return

  end subroutine evolve_star_to_log_g

  !****

  subroutine evolve_star_seismic (id, f_enter, f_exit, y_tol, dt_frac)

    integer, intent(in)            :: id
    real(dp), intent(in), optional :: f_enter
    real(dp), intent(in), optional :: f_exit
    real(dp), intent(in), optional :: y_tol
    real(dp), intent(in), optional :: dt_frac

    type(star_info), pointer                        :: s
    integer                                         :: ierr
    procedure(extras_check_model_seismic_), pointer :: ecm_save

    ! Set up parameters

    if (PRESENT(f_enter)) then
       f_enter_m = f_enter
    else
       f_enter_m = 0.1_dp
    endif

    if (PRESENT(f_exit)) then
       f_exit_m = f_exit
    else
       f_exit_m = 0.2_dp
    endif

    if (PRESENT(y_tol)) then
       y_tol_m = y_tol
    else
       y_tol_m = 1E-4_dp
    endif

    if (PRESENT(dt_frac)) then
       dt_frac_m = dt_frac
    else
       dt_frac_m = 0.1_dp
    endif

    ! Clear fit data

    call clear_mod_freqs_()
    call clear_mod_data_()
    
    ! Evolve the star until the seismic constraints are met

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failure in star_ptr)

    ! Set controls

    ecm_save => s%extras_check_model
    s%extras_check_model => extras_check_model_seismic_

    ! Initialize the state machine

    state_m = 'START'

    ! Evolve the star

    call evolve_star(id)

    ! Reset controls

    s%extras_check_model => ecm_save

    ! Finish

    return

  end subroutine evolve_star_seismic

  !****

  function extras_check_model_seismic_ (id, id_extra) result (action)

    integer, intent(in) :: id
    integer, intent(in) :: id_extra
    integer             :: action

    logical, parameter :: DEBUG = .FALSE.

    type(star_info), pointer :: s
    integer                  :: ierr
    logical, save            :: call_run_gyre = .FALSE.
    real(dp)                 :: Delta_nu_obs
    real(dp)                 :: Delta_nu_mod
    real(dp), save           :: Delta_nu_prev = 0._dp
    type(freq_t)             :: fr_mod(0:3)
    type(freq_t)             :: fr_cor(0:3)
    integer                  :: i_lo
    real(dp)                 :: nu_obs_lo
    real(dp)                 :: nu_mod_lo
    real(dp), save           :: nu_mod_prev
    integer                  :: n_pg_lo
    integer, save            :: n_pg_prev
    real(dp)                 :: y
    real(dp)                 :: y_prev = 0._dp
    real(dp)                 :: z
    real(dp), save           :: z_min = 0._dp
    real(dp), save           :: x_a = 0._dp
    real(dp), save           :: x_b = 0._dp
    real(dp), save           :: y_a = 0._dp
    real(dp), save           :: y_b = 0._dp
    real(dp), save           :: x_new = 0._dp
    real(dp), save           :: delta_t = 0._dp
    real(dp)                 :: dt_lim
    integer                  :: l

    ! Get the star pointer

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failure in star_ptr)

    ! Decide whether to run GYRE at all

    Delta_nu_obs = fr_obs_m(0)%Delta_nu()

    if (DEBUG) then
       write(OUTPUT_UNIT, *) 'Delta_nu_obs:',Delta_nu_obs
    end if

    if (.NOT. call_run_gyre) then
       call_run_gyre = ABS(Delta_nu_obs - s%delta_nu) <= f_enter_m*Delta_nu_obs
    endif
    
    if (.NOT. call_run_gyre) then
       action = keep_going
       return
    endif

    ! call clear_mod_freqs_
    ! call clear_cor_freqs_
    fr_mod = freq_t()
    fr_cor = freq_t()
    
    ! Run GYRE for radial modes

    call run_gyre(id, 0, fr_obs_m(0), fr_mod(0))

    ! Correct frequencies before bracketing

    call apply_cubic_correction(fr_mod, fr_obs_m, fr_cor)

    Delta_nu_mod = fr_cor(0)%Delta_nu()
    
    ! Find the best match to the lowest-frequency observed mode

    nu_obs_lo = fr_obs_m(0)%nu(1)

    i_lo = MINLOC(ABS(fr_cor(0)%nu - nu_obs_lo), DIM=1)

    nu_mod_lo = fr_cor(0)%nu(i_lo)
    n_pg_lo = fr_cor(0)%n_pg(i_lo)

    ! Update discriminants

    y = nu_mod_lo - nu_obs_lo
    z = ABS(Delta_nu_mod - Delta_nu_obs)

    write(OUTPUT_UNIT, *) 'y =', y, ';  z =', z, '; nu_mod_lo =', nu_mod_lo

    ! (Possibly) switch to a new state

    action = keep_going

    select case (state_m)

    case ('START')

       z_min = HUGE(0._dp)

       state_m = 'EVOLVE'

    case ('RESTART')

       state_m = 'EVOLVE'

    case ('EVOLVE')

       if (ABS(Delta_nu_mod - Delta_nu_obs) > f_exit_m*Delta_nu_obs) then

          state_m = 'FINISH'

       elseif (y*y_prev <= 0._dp .AND. n_pg_lo == n_pg_prev) then

          ! Set up a new bracket

          delta_t = s%dt

          x_a = 0._dp
          x_b = 1._dp

          y_a = y_prev
          y_b = y

          state_m = 'BRACKET'

       endif

    case ('BRACKET')

       if (n_pg_lo == n_pg_prev) then

          ! Update the bracket

          if (y_a*y <= 0._dp) then

             x_b = x_new

             y_b = y
             
          else

             x_a = x_new

             y_a = y

          endif
             
          if (DEBUG) then
             write(OUTPUT_UNIT, *) 'Bracket updated:', x_a, x_b, y_a, y_b
          endif

          ! Check for convergence

          if (ABS(y) < y_tol_m) then

             if (DEBUG) then
                write(OUTPUT_UNIT, *) 'Converged to solution:', y_a, y_b, y, z
             endif

             ! Store the model if this is the best so far

             if (z < z_min) then

                if (DEBUG) then
                   write(OUTPUT_UNIT, *) 'Best solution updated:', y_a, y_b, y, z
                endif

                call store_mod_data_(id)

                fr_mod_m(0) = fr_mod(0)

                ! Run GYRE for non-radial modes

                do l = 1, UBOUND(fr_obs_m, 1)
                   if (fr_obs_m(l)%n > 0) then
                      call run_gyre(id, l, fr_obs_m(l), fr_mod_m(l))
                   endif
                end do

                z_min = z

             endif

             ! Continue the evolution

             s%dt = delta_t*(1._dp - x_new)

             state_m = 'RESTART'

          endif

       else

          ! Abandon the bracket, because the fundamental mode wasn't found

          if (DEBUG) then
             write(OUTPUT_UNIT, *) 'Abandoning bracket (fundamental not found)'
          endif

          s%dt = delta_t

          state_m = 'RESTART'

       endif

    case default

       $ABORT(Invalid state)

    end select

    ! Prepare for the next step

    select case (state_m)

    case ('EVOLVE')

       ! Limit the timestep based on the change in Delta_nu (relative
       ! change in Delta_nu must not exceed dt_frac_m)

       dt_lim = dt_frac_m*ABS(Delta_nu_mod/(Delta_nu_mod - Delta_nu_prev))*s%dt

       if (s%dt > dt_lim) then
          s%dt = dt_lim
          action = redo
       endif

    case ('BRACKET')

       ! Set the timestep for bisection

       x_new = 0.5_dp*(x_a + x_b)

       s%dt = delta_t*x_new
       action = redo

    case ('FINISH')

       action = terminate

    end select

    ! Save the current state for the next call

    Delta_nu_prev = Delta_nu_mod

    nu_mod_prev = nu_mod_lo
    n_pg_prev = n_pg_lo

    y_prev = y

    ! Finish

    return
      
  end function extras_check_model_seismic_

  !****

  subroutine match_modes (fr_mod, fr_obs, fr_cor)
    ! match observed modes to model modes one by one
    ! currently naive: just finds nearest mode
    type(freq_t), intent(in)  :: fr_mod(0:3)
    type(freq_t), intent(in)  :: fr_obs(0:3)
    type(freq_t), intent(out) :: fr_cor(0:3)
    integer :: i, j, l
    do l = 0, 3
       if (fr_obs(l)%n < 1) cycle
       if (fr_mod(l)%n < 1) cycle

       fr_cor(l) = fr_obs(l)

       do i = 1, fr_obs(l)%n
          j = MINLOC(ABS(fr_mod(l)%nu - fr_obs(l)%nu(i)), DIM=1)
          fr_cor(l)%nu(i) = fr_mod(l)%nu(j)
          fr_cor(l)%E_norm(i) = fr_mod(l)%E_norm(j)
          fr_cor(l)%n_pg(i) = fr_mod(l)%n_pg(j)
       end do
    end do
  end subroutine match_modes

  !****
  
  subroutine apply_cubic_correction (fr_mod, fr_obs, fr_cor)
    type(freq_t), intent(in)  :: fr_mod(0:3)
    type(freq_t), intent(in)  :: fr_obs(0:3)
    type(freq_t), intent(out) :: fr_cor(0:3)
    real(dp) :: X, y, XtX, Xty
    integer :: i, j, l

    call match_modes(fr_mod, fr_obs, fr_cor)

    XtX = 0._dp
    Xty = 0._dp
    X = 0._dp
    y = 0._dp

    do l = 0, 3
       if (fr_obs(l)%n < 1) cycle
       if (fr_cor(l)%n < 1) cycle
       do i = 1, fr_obs(l)%n
          X = fr_cor(l)%nu(i)**3/fr_cor(l)%E_norm(i)/fr_obs(l)%dnu(i)
          y = (fr_obs(l)%nu(i)-fr_cor(l)%nu(i))/fr_obs(l)%dnu(i)

          XtX = XtX + X*X
          Xty = Xty + X*y
       end do
    end do

    a3 = Xty/XtX

    do l = 0, 3
       if (fr_obs(l)%n < 1) cycle
       if (fr_cor(l)%n < 1) cycle
       do i = 1, fr_obs(l)%n
          fr_cor(l)%nu(i) = fr_cor(l)%nu(i) &
               + a3*fr_cor(l)%nu(i)**3/fr_cor(l)%E_norm(i)
       end do
    end do

  end subroutine apply_cubic_correction
  
  !****

  subroutine apply_combined_correction (fr_mod, fr_obs, fr_cor)
    type(freq_t) :: fr_mod(0:3)
    type(freq_t) :: fr_obs(0:3)
    type(freq_t) :: fr_cor(0:3)
    real(dp) :: X(2), y, XtX(2,2), XtXi(2,2), Xty(2), detXtX
    integer :: i, j, l

    call match_modes(fr_mod, fr_obs, fr_cor)

    XtX = 0._dp
    Xty = 0._dp
    X = 0._dp
    y = 0._dp

    do l = 0, 3
       if (fr_obs(l)%n < 1) cycle
       if (fr_cor(l)%n < 1) cycle
       do i = 1, fr_obs(l)%n
          X(1) = fr_cor(l)%nu(i)**(-1)/fr_cor(l)%E_norm(i)/fr_obs(l)%dnu(i)
          X(2) = fr_cor(l)%nu(i)**3/fr_cor(l)%E_norm(i)/fr_obs(l)%dnu(i)
          y = (fr_obs(l)%nu(i)-fr_cor(l)%nu(i))/fr_obs(l)%dnu(i)

          XtX(1,1) = XtX(1,1) + X(1)*X(1)
          XtX(1,2) = XtX(1,2) + X(1)*X(2)
          XtX(2,2) = XtX(2,2) + X(2)*X(2)
          Xty(1) = Xty(1) + X(1)*y
          Xty(2) = Xty(2) + X(2)*y
       end do
    end do

    XtX(2,1) = XtX(1,2)

    XtXi(1,1) = XtX(2,2)
    XtXi(2,2) = XtX(1,1)
    XtXi(1,2) = -XtX(1,2)
    XtXi(2,1) = -XtX(2,1)

    detXtX = XtX(1,1)*XtX(2,2) - XtX(1,2)*XtX(2,1)
    XtXi = XtXi/detXtX

    a1 = XtXi(1,1)*Xty(1) + XtXi(1,2)*Xty(2)
    a3 = XtXi(2,1)*Xty(1) + XtXi(2,2)*Xty(2)

    do l = 0, 3
       if (fr_obs(l)%n < 1) cycle
       if (fr_cor(l)%n < 1) cycle
       do i = 1, fr_obs(l)%n
          fr_cor(l)%nu(i) = fr_cor(l)%nu(i) &
               + (a1*fr_cor(l)%nu(i)**(-1) + a3*fr_cor(l)%nu(i)**3) &
               /fr_cor(l)%E_norm(i)
       end do
    end do

  end subroutine apply_combined_correction
  
  !****

  subroutine destroy_star (id)

    integer, intent(in) :: id

    integer :: ierr

    ! Free the handle

    call free_star(id, ierr)
    $ASSERT(ierr == 0,Failure in free_star)

    ! Finish

    return

  end subroutine destroy_star

  !****

  subroutine set_obs_freqs (l, fr)

    integer, intent(in)      :: l
    type(freq_t), intent(in) :: fr

    integer :: n
    integer :: i

    ! Set observational frequency data for harmonic degree l

    $ASSERT(l >= LBOUND(fr_obs_m, 1) .AND. l <= UBOUND(fr_obs_m, 1),Invalid harmonic degree)

    fr_obs_m(l) = fr

    ! Finish

    return

  end subroutine set_obs_freqs

  !****

  subroutine clear_obs_freqs ()

    ! Clear observational frequency data

    fr_obs_m = freq_t()

    ! Finish

    return

  end subroutine clear_obs_freqs

  !****

  subroutine clear_cor_freqs_ ()

    ! Clear corrected frequency data

    fr_cor_m = freq_t()

    ! Finish

    return

  end subroutine clear_cor_freqs_

  !****

  function get_mod_freqs (l) result (fr)

    integer, intent(in) :: l
    type(freq_t)        :: fr

    ! Get model frequency data for harmonic degree l

    $ASSERT(l >= LBOUND(fr_obs_m, 1) .AND. l <= UBOUND(fr_obs_m, 1),Invalid harmonic degree)

    fr = fr_mod_m(l)

    ! Finish

    return

  end function get_mod_freqs

  !****

  function get_cor_freqs (l) result (fr)

    integer, intent(in) :: l
    type(freq_t)        :: fr

    ! Get model frequency data for harmonic degree l

    $ASSERT(l >= LBOUND(fr_obs_m, 1) .AND. l <= UBOUND(fr_obs_m, 1),Invalid harmonic degree)

    fr = fr_cor_m(l)

    ! Finish

    return

  end function get_cor_freqs

  !****

  subroutine clear_mod_freqs_ ()

    ! Clear model frequency data

    fr_mod_m = freq_t()

    ! Finish

    return
    
  end subroutine clear_mod_freqs_

  !****

  subroutine get_mod_data (Teff, logg, FeH, R, L, age)

    real(dp), intent(out) :: Teff
    real(dp), intent(out) :: logg
    real(dp), intent(out) :: FeH
    real(dp), intent(out) :: R
    real(dp), intent(out) :: L
    real(dp), intent(out) :: age
    
    ! Return best-fit model data

    Teff = Teff_m
    logg = logg_m
    FeH = FeH_m
    R = R_m
    L = L_m
    age = age_m

    ! Finish

    return

  end subroutine get_mod_data

  !****

  subroutine clear_mod_data_ ()

    ! Clear model data

    Teff_m = 0._dp
    logg_m = 0._dp
    FeH_m = 0._dp
    R_m = 0._dp
    L_m = 0._dp
    age_m = 0._dp

    ! Finish

    return

  end subroutine clear_mod_data_

  !****

  subroutine store_mod_data_ (id)

    integer, intent(in) :: id

    type(star_info), pointer :: s
    integer                  :: ierr

    ! Get the star pointer

    call star_ptr(id, s, ierr)
    $ASSERT(ierr == 0,Failure in star_ptr)

    ! Store model data

    Teff_m = s%Teff
    logg_m = s%log_surface_gravity
    FeH_m = log10((1.-s%surface_h1-s%surface_he3-s%surface_he4)/s%surface_h1)+1.64d0
    !FeH_m = log10((s%surface_z / s%surface_x) / 0.02293d0 )
    R_m = s%photosphere_r
    L_m = s%photosphere_L
    age_m = s%star_age

    ! Finish

    return

  end subroutine store_mod_data_

end module em_lib

  
