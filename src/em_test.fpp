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

  integer      :: id
  type(freq_t) :: fr_obs(0:3)
  type(freq_t) :: fr_mod(0:3)
  type(freq_t) :: fr_cor(0:3)
  integer      :: i
  integer      :: j
  real(dp)     :: Teff
  real(dp)     :: logg
  real(dp)     :: FeH
  real(dp)     :: R
  real(dp)     :: L
  real(dp)     :: age

  real(dp), allocatable :: r_mod(:)
  real(dp), allocatable :: r_freq(:)
  real(dp), allocatable :: r_obs(:)
  real(dp), allocatable :: r_invcov(:,:)
  integer :: n_r, n_r010, n_r02, n_r13, ierr
  real(dp) :: chi2

  ! Define observational frequencies

  fr_obs = freq_t()
  ! Rich's original data
  ! fr_obs(0) = freq_t([799.70d0,855.30d0,909.92d0,965.16d0,1021.81d0,1078.97d0,1135.32d0,1192.12d0,1250.12d0], &
  !                   [0.27d0,0.73d0,0.26d0,0.36d0,0.28d0,0.33d0,0.34d0,0.45d0,0.89d0], &
  !                   [1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0], &
  !                   [1,2,3,4,5,6,7,8,9])

  ! fr_obs(1) = freq_t([748.60d0,777.91d0,828.21d0,881.29d0,935.90d0,991.09d0,1047.79d0,1104.68d0,1161.27d0,1216.95d0], &
  !                   [0.23d0,0.24d0,0.42d0,0.29d0,0.23d0,0.22d0,0.24d0,0.22d0,0.33d0,0.53d0], &
  !                   [1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0], &
  !                   [1,2,3,4,5,6,7,8,9,10])

  ! fr_obs(2) = freq_t([794.55d0,905.31d0,961.47d0,1017.56d0,1075.01d0,1130.79d0,1187.55d0,1246.78d0], &
  !                   [0.52d0,0.35d0,0.49d0,0.27d0,0.27d0,0.61d0,0.32d0,0.84d0], &
  !                   [1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0], &
  !                   [1,2,3,4,5,6,7,8])

  ! Travis' Sun-like data
  ! fr_obs(0) = freq_t([2093.53d0, 2228.95d0, 2363.03d0, 2496.29d0, 2629.84d0, 2764.32d0, 2899.11d0, &
  !                    3033.92d0, 3168.94d0, 3303.73d0, 3439.17d0, 3575.23d0, 3711.6d0, 3848.2d0, 3984.79d0], &
  !                   [0.14d0, 0.1d0, 0.07d0, 0.07d0, 0.06d0, 0.05d0, 0.04d0, 0.04d0, &
  !                    0.05d0, 0.06d0, 0.1d0, 0.15d0, 0.26d0, 0.5d0, 0.65d0], &
  !                   [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
  !                    1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
  !                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
  
  ! fr_obs(1) = freq_t([2156.91d0, 2292.05d0, 2425.82d0, 2559.31d0, 2693.56d0, 2828.36d0, 2963.45d0, &
  !                    3098.39d0, 3233.4d0, 3368.88d0, 3504.5d0, 3640.95d0, 3777.21d0, 3913.93d0, 4051.51d0], &
  !                   [0.08d0, 0.09d0, 0.09d0, 0.06d0, 0.06d0, 0.05d0, 0.05d0, 0.05d0, &
  !                    0.05d0, 0.08d0, 0.09d0, 0.16d0, 0.2d0, 0.36d0, 0.61d0], &
  !                   [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
  !                    1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
  !                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

  ! fr_obs(2) = freq_t([2082.35d0, 2217.71d0, 2352.32d0, 2486.11d0, 2619.87d0, 2754.72d0, 2889.66d0, &
  !                    3024.81d0, 3160.05d0, 3295.05d0, 3430.88d0, 3566.95d0, 3703.7d0, 3840.9d0, 3977.52d0], &
  !                   [0.55d0, 0.26d0, 0.14d0, 0.13d0, 0.08d0, 0.07d0, 0.06d0, 0.07d0, &
  !                    0.07d0, 0.09d0, 0.14d0, 0.24d0, 0.32d0, 0.64d0, 0.87d0], &
  !                   [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
  !                    1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
  !                   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

  ! fr_obs(3) = freq_t([2541.91d0, 2676.26d0, 2811.32d0, 2947.1d0, 3082.24d0, 3218.03d0, 3352.91d0, 3489.33d0, 3625.71d0], &
  !      [0.53d0, 0.41d0, 0.3d0, 0.16d0, 0.26d0, 0.16d0, 0.34d0, 0.43d0, 0.64d0], &
  !      [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
  !      [1, 2, 3, 4, 5, 6, 7, 8, 9])

  ! Warrick's dataset, reduced from Travis'
  fr_obs(0) = freq_t([2496.29d0, 2629.84d0, 2764.32d0, 2899.11d0, 3033.92d0, &
                      3168.94d0, 3303.73d0, 3439.17d0, 3575.23d0, 3711.6d0, 3848.2d0, 3984.79d0], &
                     [0.07d0, 0.06d0, 0.05d0, 0.04d0, 0.04d0, &
                      0.05d0, 0.06d0, 0.1d0, 0.15d0, 0.26d0, 0.5d0, 0.65d0], &
                     [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
                      1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

  fr_obs(1) = freq_t([2559.31d0, 2693.56d0, 2828.36d0, 2963.45d0, &
                      3098.39d0, 3233.4d0, 3368.88d0, 3504.5d0, 3640.95d0, 3777.21d0, 3913.93d0, 4051.51d0], &
                     [0.06d0, 0.06d0, 0.05d0, 0.05d0, 0.05d0, &
                      0.05d0, 0.08d0, 0.09d0, 0.16d0, 0.2d0, 0.36d0, 0.61d0], &
                     [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
                      1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

  fr_obs(2) = freq_t([2486.11d0, 2619.87d0, 2754.72d0, 2889.66d0, &
                      3024.81d0, 3160.05d0, 3295.05d0, 3430.88d0, 3566.95d0, 3703.7d0, 3840.9d0, 3977.52d0], &
                     [0.13d0, 0.08d0, 0.07d0, 0.06d0, 0.07d0, &
                      0.07d0, 0.09d0, 0.14d0, 0.24d0, 0.32d0, 0.64d0, 0.87d0], &
                     [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
                      1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
  
  ! Warrick's dataset, extra reduced from Travis'
  ! fr_obs(0) = freq_t([2764.32d0, 2899.11d0, 3033.92d0, 3168.94d0, 3303.73d0, &
  !                     3439.17d0, 3575.23d0, 3711.6d0, 3848.2d0, 3984.79d0], &
  !                    [0.05d0, 0.04d0, 0.04d0, 0.05d0, 0.06d0, &
  !                     0.1d0, 0.15d0, 0.26d0, 0.5d0, 0.65d0], &
  !                    [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
  !                     1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
  !                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

  ! fr_obs(1) = freq_t([2828.36d0, 2963.45d0, 3098.39d0, 3233.4d0, 3368.88d0, &
  !                     3504.5d0, 3640.95d0, 3777.21d0, 3913.93d0, 4051.51d0], &
  !                    [0.05d0, 0.05d0, 0.05d0, 0.05d0, 0.08d0, &
  !                     0.09d0, 0.16d0, 0.2d0, 0.36d0, 0.61d0], &
  !                    [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
  !                     1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
  !                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

  ! fr_obs(2) = freq_t([2754.72d0, 2889.66d0, 3024.81d0, 3160.05d0, 3295.05d0, &
  !                     3430.88d0, 3566.95d0, 3703.7d0, 3840.9d0, 3977.52d0], &
  !                    [0.07d0, 0.06d0, 0.07d0, 0.07d0, 0.09d0, &
  !                     0.14d0, 0.24d0, 0.32d0, 0.64d0, 0.87d0], &
  !                    [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
  !                     1.d0, 1.d0, 1.d0, 1.d0, 1.d0], &
  !                    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
  
  call set_obs_freqs(0, fr_obs(0))
  call set_obs_freqs(1, fr_obs(1))
  call set_obs_freqs(2, fr_obs(2))

  ! Read data and (inverse) covariance for ratios
  ! Produced externally with Python
  ! I should be better about IO units but let's first make something that works
  n_r010 = 20
  n_r02 = 11
  n_r13 = 0
  n_r = n_r010 + n_r02 + n_r13

  allocate(r_freq(n_r))
  allocate(r_obs(n_r))
  allocate(r_invcov(n_r, n_r))
  
  ierr = 0
  open(unit=10, file='ratio_means.txt', status='old', iostat=ierr)
  if (ierr /= 0) then
     write(*,*) 'failed to open ratio_means.txt'
     stop 1
  end if

  do i = 1, n_r
     read(10, *) r_freq(i), r_obs(i)
  end do

  close(10)

  open(unit=10, file='ratio_invcov.txt', status='old', iostat=ierr)
  if (ierr /= 0) then
     write(*,*) 'failed to open ratio_invcov.txt'
     stop 1
  end if

  do i = 1, n_r
     read(10, *) r_invcov(i,:)
  end do

  close(10)

  ! Done loading ratio data

  ! Create a star

  ! id = create_star( &
  !      M=1.3120640689574825D+00, &
  !      Y=0.28D0, &
  !      Z=0.02D0, &
  !      alpha=1.6295047446620596D+00, &
  !      f_ov=1.4999999999999999D-02)
  id = create_star( &
       M=1.00d0, &
       Z=0.017d0, &
       Y=0.265d0, &
       alpha=2.0d0, &
       f_ov=0.0d0)

  ! Evolve it to the ZAMS

  call evolve_star_to_zams(id)

  ! Evolve it until seismic constraints are met

  call evolve_star_seismic(id)

  ! Get model data

  call get_mod_data(Teff, logg, FeH, R, L, age)

  print *,'Teff:', Teff
  print *,'logg:', logg
  print *,'FeH:', FeH
  print *,'R:', R
  print *,'L:', L
  print *,'age:', age

  ! Get model frequencies

  fr_mod = freq_t()
  fr_mod(0) = get_mod_freqs(0)
  fr_mod(1) = get_mod_freqs(1)
  fr_mod(2) = get_mod_freqs(2)

  print *,'l=0 results (n_pg, nu, E_norm):'
  call write_results(fr_mod(0))

  print *,'l=1 results (n_pg, nu, E_norm):'
  call write_results(fr_mod(1))

  print *,'l=2 results (n_pg, nu, E_norm):'
  call write_results(fr_mod(2))

  ! Get corrected frequencies

  fr_cor = freq_t()

  print *, 'Cubic corrected frequencies'
  call apply_cubic_correction(fr_mod, fr_obs, fr_cor)

  print *,'l=0 results (n_pg, nu, E_norm):'
  call write_results(fr_cor(0))

  print *,'l=1 results (n_pg, nu, E_norm):'
  call write_results(fr_cor(1))

  print *,'l=2 results (n_pg, nu, E_norm):'
  call write_results(fr_cor(2))

  print *, 'Combined corrected frequencies'
  call apply_combined_correction(fr_mod, fr_obs, fr_cor)

  print *,'l=0 results (n_pg, nu, E_norm):'
  call write_results(fr_cor(0))

  print *,'l=1 results (n_pg, nu, E_norm):'
  call write_results(fr_cor(1))

  print *,'l=2 results (n_pg, nu, E_norm):'
  call write_results(fr_cor(2))

  print *,'complete results:'
  print *,'              fr_obs                         fr_err               fr_mod                  fr_mod_E_norm               fr_cor'
  do i = 0, 2
     do j = 1, fr_obs(i)%n
        print *, fr_obs(i)%nu(j), fr_obs(i)%dnu(j), fr_mod(i)%nu(j), &
             fr_mod(i)%E_norm(j), fr_cor(i)%nu(j)
     end do
  end do

  j = fr_obs(1)%n - 2
  allocate(r_mod(j))
  call get_r010(fr_obs(1)%nu(2:j+1), fr_mod(0), fr_mod(1), r_mod)

  write(*,*) 'r010 ratios'
  write(*,*) '          i               r_freq                           r010'
  do i = 1, size(r_mod)
     write(*,*) i, fr_obs(1)%nu(i+1), r_mod(i)
  end do
  deallocate(r_mod)

  j = fr_obs(2)%n - 1
  allocate(r_mod(j))
  call get_r02(fr_obs(0)%nu, fr_mod(0), fr_mod(1), fr_mod(2), r_mod)
  
  write(*,*) 'r02 ratios'
  write(*,*) '          i               r_freq                            r02'
  do i = 1, size(r_mod)
     write(*,*) i, fr_obs(2)%nu(i+1), r_mod(i)
  end do
  deallocate(r_mod)

  ! Compute chi^2 for ratios

  call chi2_ratios(n_r010, n_r02, n_r13, r_freq, r_obs, r_invcov, &
       fr_mod, chi2)
  write(*,*) "chi2_ratios =", chi2

  deallocate(r_obs, r_invcov)
  
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
