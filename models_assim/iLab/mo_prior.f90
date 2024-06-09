!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file mo_prior.f90
!> \brief defines the prior of the one-dimensional control vector
!>        together with the control vector uncertainty (in physical coordinates)
!>        and a (short) parameter name.
!>
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  January 2020
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mo_prior
  implicit none
  real(kind=8), allocatable :: x_pr(:)    !-- prior control vector in normalised coordinates
  real(kind=8), allocatable :: x_sigma(:) !-- uncertainty of prior (physical units)
  character(len=32), allocatable :: x_prname(:) !-- prior control vector component names
end module mo_prior
