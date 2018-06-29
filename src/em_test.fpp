! Module   : em_test
! Purpose  : test frontend for Embedded MESA
!
! Copyright 2016-2017 Rich Townsend

$include 'core.inc'

program em_test

  ! Uses

  use em_lib
  
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  type(freq_t) :: fr_obs_0
  type(freq_t) :: fr_obs_1
  type(freq_t) :: fr_obs_2
  integer      :: id
  integer      :: t_code
  type(freq_t) :: fr_mod_0
  type(freq_t) :: fr_mod_1
  type(freq_t) :: fr_mod_2
  integer      :: i
  real(dp)     :: Teff
  real(dp)     :: logg
  real(dp)     :: FeH
  real(dp)     :: R
  real(dp)     :: L
  real(dp)     :: age

  ! Initialize EM

  call init_em()

  ! Define observational frequencies

  fr_obs_0 = freq_t([799.70d0,855.30d0,909.92d0,965.16d0,1021.81d0,1078.97d0,1135.32d0,1192.12d0,1250.12d0], &
                    [0.27d0,0.73d0,0.26d0,0.36d0,0.28d0,0.33d0,0.34d0,0.45d0,0.89d0], &
                    [1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0], &
                    [1,2,3,4,5,6,7,8,9])

  fr_obs_1 = freq_t([748.60d0,777.91d0,828.21d0,881.29d0,935.90d0,991.09d0,1047.79d0,1104.68d0,1161.27d0,1216.95d0], &
                    [0.23d0,0.24d0,0.42d0,0.29d0,0.23d0,0.22d0,0.24d0,0.22d0,0.33d0,0.53d0], &
                    [1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0], &
                    [1,2,3,4,5,6,7,8,9,10])

  fr_obs_2 = freq_t([794.55d0,905.31d0,961.47d0,1017.56d0,1075.01d0,1130.79d0,1187.55d0,1246.78d0], &
                    [0.52d0,0.35d0,0.49d0,0.27d0,0.27d0,0.61d0,0.32d0,0.84d0], &
                    [1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0], &
                    [1,2,3,4,5,6,7,8])
  
  call set_obs_freqs(0, fr_obs_0)
  call set_obs_freqs(1, fr_obs_1)
  call set_obs_freqs(2, fr_obs_2)
  
  ! Create a star

  id = create_star( &
       M=1.3120640689574825D+00, &
       Y=0.28D0, &
       Z=0.02D0, &
       alpha=1.6295047446620596D+00, &
       f_ov=1.4999999999999999D-02, &
       max_age=1D11)

  ! Evolve it to the ZAMS

  call evolve_star_to_zams(id, t_code)

  if (t_code == t_ok) then
     print *,'Evolve to ZAMS: OK'
  else
     print *,'Evolve to ZAMS: Failed, termination code =', t_code
     stop
  end if

  ! Evolve it until seismic constraints are met

  call evolve_star_seismic(id, t_code)

  if (t_code == t_ok) then
     print *,'Evolve to seismic: OK'
  elseif (t_code == t_max_age) then
     print *,'Evolve to seismic: Reached maximum age'
  else
     print *,'Evolve to seismic: Failed, termination code =', t_code
     stop
  endif

  ! Get model data

  call get_mod_data(Teff, logg, FeH, R, L, age)

  print *,'Teff:', Teff
  print *,'logg:', logg
  print *,'FeH:', FeH
  print *,'R:', R
  print *,'L:', L
  print *,'age:', age

  ! Get model frequencies

  fr_mod_0 = get_mod_freqs(0)
  fr_mod_1 = get_mod_freqs(1)
  fr_mod_2 = get_mod_freqs(2)

  print *,'l=0 results (n_pg, nu, E_norm):'
  call write_results(fr_mod_0)

  print *,'l=1 results (n_pg, nu, E_norm):'
  call write_results(fr_mod_1)

  print *,'l=2 results (n_pg, nu, E_norm):'
  call write_results(fr_mod_2)

  ! Destroy the star

  call destroy_star(id)

  ! Finalize EM

  call final_em()

  ! Finish

contains

  subroutine write_results (fr)

    type(freq_t), intent(in) :: fr

    integer :: i

    do i = 1, fr%n
       print *,'  ', fr%n_pg(i), fr%nu(i), fr%E_norm(i)
    end do

    ! Finish

    return

  end subroutine write_results

end program em_test

