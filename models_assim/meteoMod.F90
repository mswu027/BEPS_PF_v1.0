module meteoMod
use shr_kind_mod, only: r8=>shr_kind_r8
use beps_con,only: rho_a
implicit none

real(r8),public :: density_air     ! kg/m3
real(r8),public :: cp_air          ! specific heat of air  J/kg/C
real(r8),public :: vpd             ! kPa
real(r8),public :: slope_vapor     ! slope  of vapor presure to temperature kPa/C
real(r8),public :: psy             ! psychrometer constant  kPa/C
real(r8),public :: e_saturate      ! saturate water vapor   kPa
real(r8),public :: e_actual        ! actual water vapor     kPa
real(r8),public :: sp_humidity     ! specific humidity g/g
real(r8),public :: latent_water    ! condensation    J/kg
real(r8),public :: latent_snow     ! sublimation/deposition    J/kg

contains

subroutine  meteo_pack(temp,rh)
implicit none
real(r8),intent(in)  :: temp
real(r8),intent(in)  :: rh

density_air  = rho_a         ! 0C 1.292
e_saturate   = 0.61078*exp(17.3*temp/(237.3+temp))
e_actual     = e_saturate*rh/100.
vpd          = e_saturate - e_actual
sp_humidity  = 0.622*e_actual/(101.35-0.378*e_actual)
cp_air       = 1004.65*(1.+0.84*sp_humidity)
slope_vapor  = 2503./(temp+237.3)**2*exp(17.27*temp/(temp+237.3))
psy          = 0.066
latent_water = (2.501 - 0.00237*temp)*1e6 
latent_snow  = 2.83*1e6

end subroutine


end module





