! Module   : em_freq
! Purpose  : frequency data & manipulation
!
! Copyright 2016 Rich Townsend

$include 'core.inc'

module em_freq

  ! Uses

  use core_kinds
  use core_order

  use em_util

  ! No implicit typing

  implicit none

  ! Derived-type definitions

  type :: freq_t
     real(dp), allocatable :: nu(:)
     real(dp), allocatable :: dnu(:)
     integer, allocatable  :: n_pg(:)
     integer               :: n = 0
   contains
     private
     procedure, public :: Delta_nu
  end type freq_t

  ! Interfaces

  interface freq_t
     module procedure freq_t_
  end interface freq_t

  ! Access specifiers

  private

  public :: freq_t

contains

  function freq_t_ (nu, dnu, n_pg) result (fr)

    real(dp), intent(in) :: nu(:)
    real(dp), intent(in) :: dnu(:)
    integer, intent(in)  :: n_pg(:)
    type(freq_t)         :: fr

    integer :: i(SIZE(nu))

    $CHECK_BOUNDS(SIZE(dnu),SIZE(nu))
    $CHECK_BOUNDS(SIZE(n_pg),SIZE(nu))

    ! Construct the freq_t

    i = sort_indices(nu)

    fr%nu = nu(i)
    fr%dnu = dnu(i)

    fr%n_pg = n_pg(i)

    fr%n = SIZE(nu)

    ! Finish

    return

  end function freq_t_

  !****

  function Delta_nu (fr)

    class(freq_t), intent(in)      :: fr
    real(dp)                       :: Delta_nu

    real(dp) :: a
    real(dp) :: b

    ! Determine Delta_nu using a linear fit

    call linear_fit(REAL(fr%n_pg, dp), fr%nu, fr%dnu, a, b)

    Delta_nu = b

    ! Finish

    return

  end function Delta_nu

end module em_freq

  
