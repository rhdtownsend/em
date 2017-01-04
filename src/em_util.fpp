! Module   : em_util
! Purpose  : utility routines for Embedded MESA module
!
! Copyright 2016 Rich Townsend

$include 'core.inc'

module em_util

  ! Uses

  use const_def, only: dp

  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Access specifiers

  private

  public :: linear_fit

contains

  subroutine linear_fit (x, y, dy, a, b)

    real(dp), intent(in)  :: x(:)
    real(dp), intent(in)  :: y(:)
    real(dp), intent(in)  :: dy(:)
    real(dp), intent(out) :: a
    real(dp), intent(out) :: b

    real(dp) :: S
    real(dp) :: S_x
    real(dp) :: S_y
    real(dp) :: S_xx
    real(dp) :: S_xy
    real(dp) :: Delta

    $CHECK_BOUNDS(SIZE(y),SIZE(x))
    $CHECK_BOUNDS(SIZE(dy),SIZE(x))

    ! Set up fitting sums

    S = SUM(1._dp/dy**2)

    S_x = SUM(x/dy**2)
    S_y = SUM(y/dy**2)

    S_xx = SUM(x**2/dy**2)
    S_xy = SUM(x*y/dy**2)

    ! Calculate the fit

    Delta = S*S_xx - S_x**2

    a = (S_xx*S_y - S_x*S_xy)/Delta
    b = (S*S_xy - S_x*S_y)/Delta

    ! Finish

    return

  end subroutine linear_fit

end module em_util
