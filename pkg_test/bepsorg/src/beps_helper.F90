

MODULE mo_helper
  use shr_kind_mod, only: r8 =>shr_kind_r8

  IMPLICIT NONE

  ! PARAMETERS FOR AUXILLIARY FUNCTIONS
  real, PARAMETER :: zmin = 1.e-18
!  REAL, PARAMETER :: eta = 0.99999            ! curvature parameter for mins/maxs

CONTAINS

  !*********************************************************
  !* FUNCTION errf
  !* the (cumulative) error function
  !* numerical recipes in Fortran 77, Chapter 6.2
  !*********************************************************

  real FUNCTION errf (x)
    real :: x, z, t
    z=ABS(x)
    t=1./(1.+0.5*z)
    errf=1.-0.5*t*EXP(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
         t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
         t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
    IF (x<0) errf=1.-errf
  END FUNCTION errf

  !*********************************************************
  !*  FUNCTION mins
  !*  smoothed minimum function
  !*********************************************************

!  REAL FUNCTION mins (x, y)
!    REAL :: x, y
!    REAL :: z
!    z = (x+y)**2 - 4.*eta*x*y
!    IF (z.GE.zmin) THEN
!       mins = (x + y - SQRT(z)) / (2.*eta)
!    ELSE
!       mins = 0.
!    ENDIF
!  END FUNCTION mins
  real FUNCTION mins (x, y, eta)
    real :: x, y, eta
    real :: z
    z = (x+y)**2 - 4.*eta*x*y
    z = max (z, zmin)
    mins = (x + y - SQRT(z)) / (2.*eta)
  END FUNCTION mins

  !*********************************************************
  !*  FUNCTION maxs
  !*  smoothed maximum function
  !*********************************************************

!  REAL FUNCTION maxs (x, y)
!    REAL :: x, y
!    REAL :: z
!    z = (x+y)**2 - 4.*eta*x*y
!    IF (z.GE.zmin) THEN
!       maxs = (x + y + SQRT(z)) / (2.*eta)
!    ELSE
!       maxs = 0.
!    ENDIF
!  END FUNCTION maxs
  real FUNCTION maxs (x, y, eta)
    real :: x, y, eta
    real :: z
    z = (x+y)**2 - 4.*eta*x*y
    z = max (z, zmin)
    maxs = (x + y + SQRT(z)) / (2.*eta)
  END FUNCTION maxs

  !*********************************************************
  !*  FUNCTION minx
  !*  minimum function with exponential transition
  !*********************************************************

  real FUNCTION minx (x, y, x0)
    real :: x, y, x0
    IF (x.LE.y+x0) THEN
       minx = x - x0*EXP((x-y)/x0-1.)
    ELSE
       minx = y
    ENDIF
  END FUNCTION minx

 real function fominef_ss(p1,p2,tune)
    real :: p1,p2
    real :: tune ! tune/2 >= min(p1,p2) - fominef
    fominef_ss = min(p1,p2)-0.5*tune*exp(-abs(p1-p2)/tune)
  end function fominef_ss

 real function fomaxef_ss(p1,p2,tune)
    real :: p1,p2
    real :: tune ! tune/2 >= fomaxef - max(p1,p2)
    fomaxef_ss = max(p1,p2)+0.5*tune*exp(-abs(p1-p2)/tune)
  end function fomaxef_ss

  !*********************************************************
  !*  FUNCTION maxx
  !*  maximum function with exponential transition
  !*********************************************************

  real FUNCTION maxx (x, y, x0)
    real,intent(in) :: x, y, x0 ! snb: to fix bug in lf95-optimisation (variable in function erroneously inherits intent(out) attribute from calling scope)
    IF (x.GE.y-x0) THEN
       maxx = x + x0*EXP(-(x-y)/x0-1.)
    ELSE
       maxx = y
    ENDIF
  END FUNCTION maxx

  !*********************************************************
  !*  FUNCTION mmin
  !*  maximum function with exponential transition
  !*********************************************************

  real FUNCTION mmin (x, y, x0) ! maximum error is x0/e for x==y-x0
    real,intent(in) :: x, y, x0 ! snb: to fix bug in lf95-optimisation (variable in function erroneously inherits intent(out) attribute from calling scope)
    IF (x.LT.y) THEN
       mmin = x + (y-x)*EXP((x-y)/x0)
    ELSE ! exact for (x.le.y):
       mmin = y
    ENDIF
  END FUNCTION mmin

  !*********************************************************
  !*  FUNCTION mmax
  !*  maximum function with exponential transition
  !*********************************************************

  real FUNCTION mmax (x, y, x0) ! maximum error is -x0/e for x==y+x0
    real,intent(in) :: x, y, x0 ! snb: to fix bug in lf95-optimisation (variable in function erroneously inherits intent(out) attribute from calling scope)
    IF (x.GT.y) THEN
       mmax = x + (y-x)*EXP((y-x)/x0)
    ELSE ! exact for (x.le.y):
       mmax = y
    ENDIF
  END FUNCTION mmax

  !*********************************************************
  !*  FUNCTION my_random
  !*  
  !*********************************************************


END MODULE mo_helper

