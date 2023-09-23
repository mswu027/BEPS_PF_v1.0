!!*******************************************
!! This module is used to get/put soilp data at single point
!! from global datasets
!! flag =0: get from ; flag=1: put into
!! Created by J.Wang
!!*******************************************

subroutine retrive_soilp(soilp,i,j,flag)
use bepstype
use beps_soilMod
implicit none

type(soil)            :: soilp
integer               :: i,j,flag
type(soils),pointer   :: p
  !iLab-added:function with optional argument requires an explicit interface!
  interface
     subroutine endrun(msg)
       implicit none
       character(len=*), intent(in), optional :: msg
     end subroutine endrun
  end interface

p  => soilstat

if(i > npoints .or. j > PFT) call endrun("out of the range of spatial points or PFTs")

if(flag == 0 ) then   ! getting data
    soilp%n_layer        = p%n_layer(i)
    soilp%Zp             = p%Zp(i,j)
    soilp%Zsp            = p%Zsp(i,j)
    soilp%r_rain_g       = p%r_rain_g(i,j)
    soilp%r_drainage     = p%r_drainage(i,j)
    soilp%r_root_decay   = p%r_root_decay(i,j)
    soilp%psi_min        = p%psi_min(i,j)
    soilp%alpha          = p%alpha(i,j)
    soilp%f_soilwater    = p%f_soilwater(i,j)

    soilp%d_soil(:)           = p%d_soil(i,:)
    soilp%f_root(:)           = p%f_root(i,:,j)
    soilp%dt(:)               = p%dt(i,:,j)
    soilp%thermal_cond(:)     = p%thermal_cond(i,:,j)
    soilp%theta_vfc(:)        = p%theta_vfc(i,:,j)
    soilp%theta_vwp(:)        = p%theta_vwp(i,:,j)
    soilp%fei(:)              = p%fei(i,:,j)
    soilp%Ksat(:)             = p%Ksat(i,:,j)
    soilp%psi_sat(:)          = p%psi_sat(i,:,j)
    soilp%b(:)                = p%b(i,:,j)
    soilp%density_soil(:)     = p%density_soil(i,:)
    soilp%f_org(:)            = p%f_org(i,:,j)
    soilp%ice_ratio(:)        = p%ice_ratio(i,:,j)
    soilp%thetam(:)           = p%thetam(i,:,j)
    soilp%thetam_prev(:)      = p%thetam_prev(i,:,j)
    soilp%temp_soil_p(:)      = p%temp_soil_p(i,:,j)
    soilp%temp_soil_c(:)      = p%temp_soil_c(i,:,j)
    soilp%f_ice(:)            = p%f_ice(i,:,j)
    soilp%psim(:)             = p%psim(i,:,j)
    soilp%thetab(:)           = p%thetab(i,:,j)
    soilp%psib(:)             = p%psib(i,:,j)
    soilp%r_waterflow(:)      = p%r_waterflow(i,:,j)
    soilp%km(:)               = p%km(i,:,j)
    soilp%kb(:)               = p%kb(i,:,j)
    soilp%KK(:)               = p%KK(i,:,j)
    soilp%Cs(:)               = p%Cs(i,:,j)
    soilp%lambda(:)           = p%lambda(i,:,j)
    soilp%Ett(:)              = p%Ett(i,:,j)
    soilp%G(:)                = p%G(i,:,j)

else if (flag ==1) then   !! storing data
    p%n_layer(i)         = soilp%n_layer
    p%Zp(i,j)            = soilp%Zp
    p%Zsp(i,j)           = soilp%Zsp
    p%r_rain_g(i,j)      = soilp%r_rain_g
    p%r_drainage(i,j)    = soilp%r_drainage
    p%r_root_decay(i,j)  = soilp%r_root_decay
    p%psi_min(i,j)       = soilp%psi_min
    p%alpha(i,j)         = soilp%alpha
    p%f_soilwater(i,j)   = soilp%f_soilwater

    p%d_soil(i,:)             = soilp%d_soil(:)
    p%f_root(i,:,j)           = soilp%f_root(:)
    p%dt(i,:,j)               = soilp%dt(:)
    p%thermal_cond(i,:,j)     = soilp%thermal_cond(:)
    p%theta_vfc(i,:,j)        = soilp%theta_vfc(:)
    p%theta_vwp(i,:,j)        = soilp%theta_vwp(:)
    p%fei(i,:,j)              = soilp%fei(:)
    p%Ksat(i,:,j)             = soilp%Ksat(:)
    p%psi_sat(i,:,j)          = soilp%psi_sat(:)
    p%b(i,:,j)                = soilp%b(:)
    p%density_soil(i,:)       = soilp%density_soil(:)
    p%f_org(i,:,j)            = soilp%f_org(:)
    p%ice_ratio(i,:,j)        = soilp%ice_ratio(:)
    p%thetam(i,:,j)           = soilp%thetam(:)
    p%thetam_prev(i,:,j)       = soilp%thetam_prev(:)
    p%temp_soil_p(i,:,j)      = soilp%temp_soil_p(:)
    p%temp_soil_c(i,:,j)      = soilp%temp_soil_c(:)
    p%f_ice(i,:,j)            = soilp%f_ice(:)
    p%psim(i,:,j)             = soilp%psim(:)
    p%thetab(i,:,j)           = soilp%thetab(:)
    p%psib(i,:,j)             = soilp%psib(:)
    p%r_waterflow(i,:,j)      = soilp%r_waterflow(:)
    p%km(i,:,j)               = soilp%km(:)
    p%kb(i,:,j)               = soilp%kb(:)
    p%KK(i,:,j)               = soilp%KK(:)
    p%Cs(i,:,j)               = soilp%Cs(:)
    p%lambda(i,:,j)           = soilp%lambda(:)
    p%Ett(i,:,j)              = soilp%Ett(:)
    p%G(i,:,j)                = soilp%G(:)
end if

end subroutine

