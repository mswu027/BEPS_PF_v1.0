!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  !> \file mo_bepsfunc_ctl.f90
  !>
  !> DESCRIPTION: environment settings to control functional implemention of BEPS
  !               runtime behaviour.
  !
  !> \authors The Inversion Lab (Michael Vossbeck) 
  !
  !> \date  March 2020
  !                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mo_bepsfunc_ctl
  implicit none

  !=======================================
  !
  !         c o n s t a n t s
  !
  !--
  real(kind=8), parameter :: fillval_r8 = -99999._8 !-- fill value (partly used for
  !--
  ! MV/TK-20191119:
  ! - we will focus on variables SIF, Thetam, COS_flux and VOD as targets for
  !   the cost function for the sitelevel version then!
  ! - potential observation operators for atmospheric CO2 can be discussed at a later stage ...
  ! MSWU-20201026: confirmed to drop VOD from the target variables
  character(len=8), dimension(*), parameter :: beps_simvars = &
       (/&
       'SIF     ',&
       'Thetam  ',&
       'COS_flux'/)
  integer, parameter :: nsimvar = size(beps_simvars)

  !--
  logical, save :: enable_netcdf_out = .true.  !-- whether to write output to NetCDF files
  !--
  logical, save :: write_restart     = .true.  !-- whether to write restart files
  !-- overall number of simulated time-points
  integer, save :: ntp = -1
  integer, allocatable, save :: time_points(:,:) !-- yr,mn,dy,tod,caldy,yrdoys => per time-step
  !-- support for NetCDF output (time-axis)
  character(len=len('YYYY-MM-DDTHH:MM:SS')), save :: ref_date = ''
  real(kind=8), allocatable :: seconds_since_ref(:)   !-- seconds elapsed since reference date
end module mo_bepsfunc_ctl
