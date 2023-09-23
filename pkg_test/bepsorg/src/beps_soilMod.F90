

!*******************************************************
! Function: Soil module( Initialize/Update soil status,soil
!           water, soil thermal etc.)
! Created : Jun Wang
! Date    : 2016/12/5
! 2023/06/20
! Xing Xiuli: particle filter DA; Prof. Wu: simulating microwave VOD.
! Combination of PA and simulating VOD
!*******************************************************
module beps_soilMod
use shr_kind_mod, only: r8=>shr_kind_r8
use beps_par
use meteoMod
use beps_con
use mid_results
implicit none
integer,private       :: i           ! generic index

type,public:: soil
    integer   :: n_layer     !the number layers used in model

    real(r8)  :: Zp          ! depth of ponded water on the ground surface
    real(r8)  :: Zsp         ! snow depth
    real(r8)  :: r_rain_g    ! the rainfall rate on ground m/s
!    real(r8)  :: soil_r     ! not used
    real(r8)  :: r_drainage  ! units:??
    real(r8)  :: r_root_decay! decay_rate_of_root_distribution
    real(r8)  :: psi_min     ! for fw
    real(r8)  :: alpha       ! for fw
    real(r8)  :: f_soilwater
! following parameters are essential for the Soil-Plant-Atmosphere Continuum (SPAC)
    real(r8)  ::  Sp    ! initial plant water storage   @Xiuli Xing 20221113
    real(r8)  :: biomass_root     ! aboveground biomass @Xiuli Xing 20221113

!!! Properties belong to each soil horizon
    real(r8)  :: d_soil(0:MAX_LAYERS-1)
    real(r8)  :: f_root(0:MAX_LAYERS-1)   !root weight
    real(r8)  :: dt(0:MAX_LAYERS-1)       ! the weight calculated from soil_water_factor
   ! read from patameters
    real(r8)  :: thermal_cond(0:MAX_LAYERS-1)  ! thermal conductivity
    real(r8)  :: theta_vfc(0:MAX_LAYERS-1)     !  field capacity
    real(r8)  :: theta_vwp(0:MAX_LAYERS-1)     !  wiltng point
    real(r8)  :: fei(0:MAX_LAYERS-1)           ! porosity
    real(r8)  :: Ksat(0:MAX_LAYERS-1)          ! saturated hydraulic conductivity
    real(r8)  :: psi_sat(0:MAX_LAYERS-1)       ! water potential in sat
    real(r8)  :: b(0:MAX_LAYERS-1)             ! Cambell parameter b
    real(r8)  :: density_soil(0:MAX_LAYERS-1)  ! soil bulk density of layer
    real(r8)  :: f_org(0:MAX_LAYERS-1)         ! volume fraction of organic matter in layer (%)
   ! needed to save
    real(r8)  :: ice_ratio(0:MAX_LAYERS-1)     ! The ratio of ice of soil layer
    real(r8)  :: thetam(0:MAX_LAYERS-1)        ! soil water content in the layer
    real(r8)  :: thetam_prev(0:MAX_LAYERS-1)
    real(r8)  :: temp_soil_p(0:MAX_LAYERS-1)   ! soil temperature in this layer
    real(r8)  :: temp_soil_c(0:MAX_LAYERS-1)   !?

!!   ! derived variables
    real(r8)  :: f_ice(0:MAX_LAYERS-1)
    real(r8)  :: psim(0:MAX_LAYERS-1)          ! soil water suction in this layer
    real(r8)  :: thetab(0:MAX_LAYERS-1)        ! soil water content at the bottom of each layer
    real(r8)  :: psib(0:MAX_LAYERS-1)          ! soil water suction at the bottom this layer
    real(r8)  :: r_waterflow(0:MAX_LAYERS-1)   ! the liquid water flow rates at the soil layer interfaces
    real(r8)  :: km(0:MAX_LAYERS-1)
    real(r8)  :: kb(0:MAX_LAYERS-1)            ! the hydraulic conducitivity
    real(r8)  :: KK(0:MAX_LAYERS-1)            ! the averaged conductivity of two soil layer
    real(r8)  :: Cs(0:MAX_LAYERS-1)            !
    real(r8)  :: lambda(0:MAX_LAYERS-1)        ! thermal conductivity of each soil layer
    real(r8)  :: Ett(0:MAX_LAYERS-1)           ! ET in each layer
    real(r8)  :: G(0:MAX_LAYERS-1)             ! energy fluxes
end type soil

public :: Init_soil_parameters, &       ! initial
          Init_soil_status,     &
          UpdateHeatFlux,       &        !Soil thermal
          Update_Cs,            &
          UpdateSoilThermalConductivity, &
          SurfaceTemperature,   &
          UpdateSoilMoisture,   &        ! soil water
          Soil_water_uptake,    &
          soil_water_factor_v2, &
          Soil_evaporation

private ::  Init_soil_rootfraction, &
           Update_ice_ratio

contains
! Initialization process
 subroutine Init_soil_parameters(lc,lai_yr,stxt,f_r_decay,r_root_decay,p)
 implicit none
 real(r8) :: biomass,biomass_leaf_o,biomass_stem_o,biomass_root_o
 real(r8) :: biomass_leaf_u,biomass_stem_u,biomass_root_u,biomass_root
 integer,intent(in)  :: lc
 integer,intent(in)  :: stxt
 real(r8),intent(in) :: r_root_decay,lai_yr,f_r_decay
 type(soil)          :: p
 write(*,*) "lc=:", lc
 write(*,*) "stxt=:", stxt
 p%n_layer      = 5
 p%KK(0:4)=max(1.e-12,p%KK(0:4))         !constrain KK
 !write(*,*)  "p%KK = ",p%KK

 if(lc == 3 .or. lc == 4) then
    p%psi_min  = 10.0      ! for fw
    p%alpha    = 1.5
 else
    p%psi_min  = 33.0
    p%alpha    = 0.4
 end if

 p%d_soil(0:4)      = (/0.05,0.10,0.20,0.40,1.25/)    ! depth_layer
 p%r_root_decay     = f_r_decay*r_root_decay          !scaling root distribution, f_r_decay

 call Init_soil_rootfraction(p)

 p%density_soil(0:4) = (/1300.0,1500.0,1517.0,1517.0,1517.0/)
 p%f_org(0:4)        = (/5.,2.,1.,1.,0.3/)

 select case (stxt)
 case(1)   ! sand
    p%b(0:4)             = (/1.7,1.9,2.1,2.3,2.5/)
    p%Ksat(0:4)          = (/58.,52.,46.,35.,10./)*1e-6   ! saturated hydraulic conductivity
    p%fei(0:4)           = 0.437
    p%theta_vfc(0:4)     = 0.09  !field capacity
    p%theta_vwp(0:4)     = 0.03  !wilting point
    p%thermal_cond(0:4)  = 8.6 ! thermal conductivity
    p%psi_sat(0:4)       = (/0.07,0.08,0.09,0.10,0.12/)
 case(2)    !loamy sand
    p%b(0:4)             = (/2.1,2.3,2.5,2.7,2.9/)
    p%Ksat(0:4)          = (/17.,15.,14.,10.,3./)*1e-6
    p%fei(0:4)           = 0.437
    p%theta_vfc(0:4)     = 0.21
    p%theta_vwp(0:4)     = 0.06
    p%thermal_cond(0:4)  = 8.3
    p%psi_sat(0:4)       = (/0.09,0.10,0.11,0.12,0.14/)
 case(3)   ! sandy loam
    p%b(0:4)             = (/3.1,3.3,3.5,3.7,3.9/)
    p%Ksat(0:4)          = (/720.,648.,576.,432.,144./)*1e-8
    p%fei(0:4)           = 0.453
    p%theta_vfc(0:4)     = 0.21
    p%theta_vwp(0:4)     = 0.10
    p%thermal_cond(0:4)  = 8.0
    p%psi_sat(0:4)       = (/0.15,0.16,0.17,0.18,0.20/)
 case(4)    !loam
    p%b(0:4)             = (/4.5,4.7,4.9,5.1,5.3/)
    p%Ksat(0:4)          = (/370.,330.,296.,222.,74./)*1e-8
    p%fei(0:4)           = 0.463
    p%theta_vfc(0:4)     = 0.27
    p%theta_vwp(0:4)     = 0.12
    p%thermal_cond(0:4)  = 7.0
    p%psi_sat(0:4)       = (/0.11,0.12,0.13,0.14,0.16/)
 case(5)   !silty loam
    p%b(0:4)             = (/4.7,4.9,5.1,5.3,5.5/)
    p%Ksat(0:4)          = (/190.,170.,152.,114.,38./)*1e-8
    p%fei(0:4)           = 0.501
    p%theta_vfc(0:4)     = 0.33
    p%theta_vwp(0:4)     = 0.13
    p%thermal_cond(0:4)  = 6.3
    p%psi_sat(0:4)       = (/0.21,0.22,0.23,0.24,0.26/)
 case(6)   ! sandy caly loam
    p%b(0:4)             = (/4.0,4.2,4.4,4.6,4.8/)
    p%Ksat(0:4)          = (/12.,10.8,96.,72.,24./)*1e-7
    p%fei(0:4)           = 0.398
    p%theta_vfc(0:4)     = 0.26
    p%theta_vwp(0:4)     = 0.15
    p%thermal_cond(0:4)  = 7.0
    p%psi_sat(0:4)       = (/0.28,0.29,0.30,0.31,0.33/)
 case(7)    !clay loam
    p%b(0:4)             = (/5.2,5.4,5.6,5.8,6.0/)
    p%Ksat(0:4)          = (/64.,58.,51.,38.,13./)*1e-8
    p%fei(0:4)           = 0.464
    p%theta_vfc(0:4)     = 0.32
    p%theta_vwp(0:4)     = 0.20
    p%thermal_cond(0:4)  = (/5.8,5.8,5.7,5.8,5.8/)
    p%psi_sat(0:4)       = (/0.26,0.27,0.28,0.29,0.31/)
 case(8)  !silty clay loam
    p%b(0:4)             = (/6.6,6.8,7.0,7.2,7.4/)
    p%Ksat(0:4)          = (/42.,38.,34.,25.2,8.4/)*1e-8
    p%fei(0:4)           = 0.471
    p%theta_vfc(0:4)     = 0.37
    p%theta_vwp(0:4)     = 0.32
    p%thermal_cond(0:4)  = 4.2
    p%psi_sat(0:4)       = (/0.33,0.34,0.35,0.36,0.38/)
 case(9)  ! sandy clay
    p%b(0:4)             = (/6.,6.2,6.4,6.6,6.8/)
    p%Ksat(0:4)          = (/33.,30.,26.4,19.8,6.6/)*1e-8
    p%fei(0:4)           = 0.43
    p%theta_vfc(0:4)     = 0.34
    p%theta_vwp(0:4)     = 0.24
    p%thermal_cond(0:4)  = 6.3
    p%psi_sat(0:4)       = (/0.29,0.30,0.31,0.32,0.34/)
 case(10) !silty clay
    p%b(0:4)             = (/7.9,8.1,8.3,8.5,8.7/)
    p%Ksat(0:4)          = (/25.,22.5,20.,15.,5./)*1e-8
    p%fei(0:4)           = 0.479
    p%theta_vfc(0:4)     = 0.39
    p%theta_vwp(0:4)     = 0.25
    p%thermal_cond(0:4)  = 4.0
    p%psi_sat(0:4)       = (/0.34,0.35,0.36,0.37,0.39/)
 case(11) ! clay
    p%b(0:4)             = (/7.6,7.8,8.0,8.2,8.4/)
    p%Ksat(0:4)          = (/17.,15.3,13.6,10.2,3.4/)*1e-8
    p%fei(0:4)           = 0.475
    p%theta_vfc(0:4)     = 0.40
    p%theta_vwp(0:4)     = 0.27
    p%thermal_cond(0:4)  = 4.4
    p%psi_sat(0:4)       = (/0.37,0.38,0.39,0.40,0.42/)
 case default
    p%b(0:4)             = (/7.6,7.8,8.0,8.2,8.4/)
    p%Ksat(0:4)          = (/17.,15.3,13.6,10.2,3.4/)*1e-8
    p%fei(0:4)           = 0.475
    p%theta_vfc(0:4)     = 0.40
    p%theta_vwp(0:4)     = 0.27
    p%psi_sat(0:4)       = (/0.37,0.38,0.39,0.40,0.42/)
 end select

if(lc >=1 .and. lc <= 5) then
           !
        biomass=0.9097*lai_yr+0.125*lai_yr*lai_yr
        biomass_leaf_o=0.05*biomass    !
        biomass_stem_o=0.95*biomass    !
        biomass_root_o=0.454*biomass
        !
        biomass_leaf_u=0.3*biomass_leaf_o  !
        biomass_stem_u=0.02*biomass_stem_o     !
        biomass_root_u=0.05*biomass_root_o !
        p%biomass_root = biomass_root_o + biomass_root_u

 else if(lc ==6 .or. lc ==9) then
        !!
        biomass=1.545*lai_yr+0.183*lai_yr*lai_yr
        biomass_leaf_o=0.04*biomass    !
        biomass_stem_o=0.96*biomass    !
        biomass_root_o=1.432*biomass**0.639    !
        biomass_leaf_u=0.3*biomass_leaf_o  !
        biomass_stem_u=0.01*biomass_stem_o     !
        biomass_root_u=0.01*biomass_root_o  !
        p%biomass_root = biomass_root_o + biomass_root_u

else if (lc == 10) then
        biomass = 1.227*lai_yr+0.154*lai_yr*lai_yr
        biomass_leaf_o  = 0.045*biomass
        biomass_stem_o  = 0.95*biomass
        biomass_root_o  = (0.454*biomass+1.432*biomass**0.639)/2.
        biomass_leaf_u  = 0.3*biomass_leaf_o
        biomass_stem_u  = 0.015*biomass_stem_o
        biomass_root_u  = 0.03*biomass_root_o
        p%biomass_root = biomass_root_o + biomass_root_u

else if (lc ==13) then
        biomass=1.545*lai_yr+0.183*lai_yr*lai_yr
        biomass_leaf_o=0.1*biomass    !
        biomass_stem_o=0.90*biomass    !
        biomass_root_o=1.432*biomass**0.639    !
        biomass_leaf_u=0.3*biomass_leaf_o     !
        biomass_stem_u=0.01*biomass_stem_o    ! 
        biomass_root_u=0.01*biomass_root_o    !
        p%biomass_root = biomass_root_o + biomass_root_u

else if(lc == 14 .or. lc == 15 .or. lc ==25 .or. lc ==40 .or. lc==41) then
        biomass_leaf_o=0.05*lai_yr  ! 
        biomass_stem_o=0.0          !
        biomass_root_o=0.061*lai_yr    !
        biomass_leaf_u=0.0
        biomass_stem_u=0.0
        biomass_root_u=0.0
        p%biomass_root = biomass_root_o + biomass_root_u
else
        biomass_leaf_o=0.05*lai_yr  ! 
        biomass_stem_o=0.0          !
        biomass_root_o=0.061*lai_yr    !
        biomass_leaf_u=0.0
        biomass_stem_u=0.0
        biomass_root_u=0.0
        p%biomass_root = biomass_root_o + biomass_root_u

end if

 return
 end subroutine

 subroutine Init_soil_status(p,Tsoil,Tair,Ms,snowdepth)
 ! soil temperatures and moisutre for each layer
 ! ponded water,snow depth
 type(soil)  :: p
 real(r8),intent(in)   :: Tsoil
 real(r8),intent(in)   :: Tair
 real(r8),intent(in)   :: Ms      ! soil moisture
 real(r8),intent(in)   :: snowdepth

 real(r8)  :: d_t

 d_t     = Tsoil - Tair
 p%Zp    = 0.0    ! depth of ponded water on the surface
 p%Zsp   = snowdepth
 p%r_rain_g    = 0.0

 if(d_t >5.0)   d_t  = 5.0
 if(d_t <-5.0)  d_t  = -5.0

 p%temp_soil_c(0)   = Tair+0.4*d_t
 p%temp_soil_c(1)   = Tair+0.5*d_t
 p%temp_soil_c(2)   = Tair+d_t
 p%temp_soil_c(3)   = Tair+1.2*d_t
 p%temp_soil_c(4)   = Tair+1.4*d_t
 p%temp_soil_c(5)   = Tair+1.4*d_t

 p%temp_soil_p(0:5) = p%temp_soil_c(0:5)

 p%thetam(0)        = 0.8*Ms
 p%thetam(1)        = Ms
 p%thetam(2)        = 1.05*Ms
 p%thetam(3)        = 1.10*Ms
 p%thetam(4)        = 1.15*Ms
 p%thetam(5)        = 1.25*Ms
 p%thetam_prev(0:5) = p%thetam(0:5)

 do i = 0,p%n_layer     !-1
    if(p%temp_soil_c(i) < -1.0) then
       p%ice_ratio(i)   = 1.0
    else if( p%temp_soil_c(i) > 0) then
       p%ice_ratio(i)   = 0.
    else
       p%ice_ratio(i)   = (0 - p%temp_soil_c(i))/1.0
    end if
 end do
 return
 end subroutine

 subroutine Init_soil_rootfraction(p)
 !Function rewritten by LHE, Jan 31,2013
 !Fortran version by Jun Wang, 8/12/2016
 implicit none
 type(soil) :: p
 real(r8)   :: cum_depth(0:MAX_LAYERS-1)

 cum_depth(0)   = p%d_soil(0)
 p%f_root(0)    = 1-p%r_root_decay**(cum_depth(0)*100)

 do i  = 1,p%n_layer-2  ! change to n_layer-1 from n_layer-2, a bug@MOUSONG,20221105
   cum_depth(i)  = cum_depth(i-1)+p%d_soil(i)
   p%f_root(i)   = p%r_root_decay**(cum_depth(i-1)*100) - p%r_root_decay**(cum_depth(i)*100)
 end do

 p%f_root(p%n_layer-1) = p%r_root_decay**(cum_depth(p%n_layer-2)*100)

 return
 end subroutine

!****************************************************
!  soil thermal regime
! update the soil temperatures for each soil layer
!****************************************************
subroutine SurfaceTemperature(temp_air,rh_air,depth_snow,depth_water,capacity_heat_soil1,capacity_heat_soil0, &
                              Gheat_g,depth_soil1,density_snow,tempL_u,netRad_g,evapo_soil,evapo_water_g,&
                              evapo_snow_g,lambda_soil1,percent_snow_g,heat_flux_soil1,temp_ground_last,&
                              temp_soil1_last,temp_any0_last,temp_snow_last,temp_soil0_last,temp_snow1_last,&
                              temp_snow2_last,temp_ground,temp_any0,temp_snow,temp_soil0,temp_snow1,temp_snow2,&
                              heat_flux)
! This subroutine will simulate the surface temperature in each step, as well as heat flux for surface to soil layers
! the core idea is to separate the interface as different layers by depth of snow, then calculate the temperature
! gradient and at last calculate the heat flux from ground surface to soil
!
! original beps would use Xg_snow[kkk] at some places
implicit none
real(r8),intent(in) :: temp_air,rh_air,depth_snow,depth_water,capacity_heat_soil1,capacity_heat_soil0
real(r8),intent(in) :: Gheat_g ! aerodynamic conductance of heat at ground Gheat = 1/ra_g
real(r8),intent(in) :: depth_soil1,density_snow,tempL_u,netRad_g
real(r8),intent(in) :: evapo_soil,evapo_water_g,evapo_snow_g
real(r8),intent(in) :: lambda_soil1  ! thermal conductivity of first layer soil
real(r8),intent(in) :: percent_snow_g,heat_flux_soil1
real(r8),intent(in) :: temp_ground_last,temp_soil1_last,temp_any0_last,temp_snow_last,&
                       temp_soil0_last,temp_snow1_last,temp_snow2_last
! ground => ground surface; soil0=>temperature of soil surface right above the soil in last step,the part is not covered by snow; soil1=>temperature of first layer soil in last step
real(r8),intent(out) :: temp_ground  ! ground surface tem in current
real(r8),intent(out) :: temp_any0    ! temperature of any layer right abover the soil,could be a mixture of snow temperature and soil surface temperature
real(r8),intent(out) :: temp_snow
real(r8),intent(out) :: temp_soil0   ! temperature of soil surface right above the soil, the part not covered by snow
real(r8),intent(out) :: temp_snow1, temp_snow2   ! temperature of snow layer 2 and 3, used when depth_snow > 0.05m
real(r8),intent(out) :: heat_flux   ! heat_flux from ground to soil

real(r8)  :: Gg   ! radiation available for heating the ground
real(r8)  :: lambda_snow   ! thermal conductivity
real(r8)  :: heat_flux_soil,heat_flux_snow  ! heat flux through the soil and snow fraction on ground, separatively
real(r8)  :: heat_flux_snow1,heat_flux_snow2
real(r8)  :: ra_g   ! aerodynamic resistence of heat
real(r8)  :: ttt    ! temporary vars


call meteo_pack(temp_air,rh_air)
ra_g   = 1./Gheat_g
lambda_snow = 0.021+4.2*density_snow/10000.+2.2*(density_snow**3)*1e-9
! available energy on ground
Gg     = netRad_g - evapo_snow_g*latent_snow - (evapo_water_g+evapo_soil)*latent_water

!!case 1 snow depth < 2cm, snow temperature ,ground temperature, soil surface temperature are the same
if( depth_snow <= 0.02) then
    ttt  = capacity_heat_soil1*0.02/kstep
    temp_ground = (temp_ground_last*ttt*ra_g*depth_soil1+ &
                   Gg*ra_g*depth_soil1 + &
                   density_air*cp_air*temp_air*depth_soil1+&
                   ra_g*lambda_soil1*temp_soil1_last)
    temp_ground = temp_ground/(density_air*cp_air*depth_soil1+ra_g*lambda_soil1+ttt*ra_g*depth_soil1)
    temp_ground = max(temp_ground_last-25,temp_ground)
    temp_ground = min(temp_ground_last+25,temp_ground)

    temp_any0   = temp_ground
    temp_snow   = temp_any0
    temp_soil0  = temp_any0
    temp_snow1  = temp_any0
    temp_snow2  = temp_any0

    heat_flux   = 2*lambda_soil1*(temp_any0 - temp_soil1_last)/depth_soil1
    heat_flux   = min(100.,heat_flux)
    heat_flux   = max(-100.,heat_flux)

else if(depth_snow > 0.02 .and. depth_snow <=0.05) then
!! snow fraction on ground decide the snow temperature based on energy balance
!! soil fraction on ground decide the soil surface temperature based on energy balance
!! snow and soil fraction works in parallel to determine the ground surface temperature
     ttt = capacity_heat_soil1*0.02/kstep    ! for soil
     temp_soil0 = (temp_soil0_last*ttt*ra_g*depth_soil1+&
                   Gg*ra_g*depth_soil1 + &
                   density_air*cp_air*temp_air*depth_soil1 + &
                   2*ra_g*lambda_soil1*temp_soil1_last) /&
                   (density_air*cp_air*depth_soil1+2*ra_g*lambda_soil1+ttt*ra_g*depth_soil1)
     temp_soil0 = max(temp_air-25,temp_soil0)
     temp_soil0 = min(temp_air+25,temp_soil0)

     ttt  = cp_ice*density_snow*depth_snow/kstep
     temp_snow  = (temp_snow_last*ttt*ra_g*depth_snow+&
                   Gg*ra_g*depth_snow+&
                   density_air*cp_air*tempL_u*depth_snow+&
                   ra_g*lambda_snow*temp_any0_last)/ &
                  (density_air*cp_air*depth_snow + ra_g*lambda_snow+ttt*ra_g*depth_snow)
     temp_snow  = max(temp_air-25,temp_snow)
     temp_snow  = min(temp_air+25,temp_snow)

     ttt = (lambda_soil1*temp_soil1_last/depth_soil1 + &
            temp_snow*lambda_snow + &
            0.02*capacity_heat_soil1/kstep*temp_any0_last)/ &
           (lambda_soil1/depth_soil1+lambda_snow/depth_snow+0.02*capacity_heat_soil1/kstep)
     temp_any0  = temp_soil0*(1-percent_snow_g) + ttt*percent_snow_g
     heat_flux_snow  = lambda_snow/(depth_snow+0.5*depth_soil1)*((temp_snow)-temp_soil1_last)
     heat_flux_soil  = heat_flux_snow*(temp_any0-temp_soil1_last)/depth_soil1

     heat_flux  = heat_flux_snow*percent_snow_g+heat_flux_soil*(1-percent_snow_g) !!!!Wrong???
     heat_flux  = min(100.,heat_flux)
     heat_flux  = max(-100.,heat_flux)

  ! starting to melt
    if(temp_snow >zero .and. temp_snow_last <= zero .and. depth_snow > zero) temp_snow = 0
  ! frozen
    if(temp_snow <zero .and. temp_snow_last >= zero .and. depth_water > zero) temp_snow = 0

    temp_ground = temp_snow*percent_snow_g + temp_soil0*(1-percent_snow_g)
    temp_ground = max(temp_air - 25, temp_ground)
    temp_ground = min(temp_air + 25, temp_ground)

    temp_snow1  = temp_snow
    temp_snow2  = temp_snow

else if(depth_snow > 0.05) then
!! case 3
!! snow_cover on ground is 100%
!! teh first layer of snow is set as 2cm
!! second layer as 2cm, too
!! the depth of third snow layer is depth_snow-0.04
    ttt  = cp_ice*density_snow*0.02/kstep
    temp_snow  = (temp_snow_last*ttt*ra_g*0.04 + &
                  Gg*ra_g*0.02 + &
                  density_air*cp_air*temp_air*0.04 + &
                  ra_g*lambda_snow*temp_snow1_last)/&
                 (density_air*cp_air*0.04+ ra_g*lambda_snow+ttt*ra_g*0.04)
    temp_snow=max(temp_air-25,temp_snow)
    temp_snow=min(temp_air+25,temp_snow)

    heat_flux_snow=lambda_snow*(temp_snow-temp_snow1_last)/0.04    !why 0.04 here?
    heat_flux=heat_flux_snow
    heat_flux=min(100.,heat_flux)
    heat_flux=max(-100.,heat_flux)

    heat_flux_snow1 = lambda_snow*(temp_snow1_last-temp_snow2_last)/(depth_snow-0.02)
    temp_snow1 = temp_snow1_last+(heat_flux- heat_flux_snow1)/(cp_ice*density_snow*0.02 )*kstep
    heat_flux_snow2 = (temp_snow2_last-temp_any0_last)/(0.5*(depth_snow-0.04)/lambda_snow+0.02/lambda_soil1)
    temp_snow2 = temp_snow2_last+(heat_flux_snow1- heat_flux_snow2)/(cp_ice*density_snow*(depth_snow-0.04))*kstep
    temp_any0  = temp_any0_last+(heat_flux_snow2- heat_flux_soil1) / (capacity_heat_soil0 * 0.02)*kstep
    temp_soil0 = temp_any0

    if(temp_snow > zero .and. temp_snow_last <= zero .and. depth_snow > zero ) temp_snow = 0
    if(temp_snow < zero .and. temp_snow_last >= zero .and. depth_water> zero ) temp_snow = 0

    temp_ground = temp_snow
end if

return
end subroutine

 subroutine UpdateHeatFlux(p,Xg_snow,lambda_snow,Tsn0,Tair_annual_mean,period_in_seconds)
! Xg_snow,lambda_snow was not be used     wangjun (why?)
 implicit none
 type(soil)   :: p
 real(r8),intent(in)   :: Xg_snow
 real(r8),intent(in)   :: lambda_snow
 real(r8),intent(in)   :: Tsn0,Tair_annual_mean
 integer,intent(in)    :: period_in_seconds
 !-- iLab::no need to have variable 's' 'implicitl save',
 !         or may also switch to parameter?
 ! real(r8)              :: S = 0.    !what? Wangjun
 real(r8) :: S

 S = 0._r8

 do i = 1,p%n_layer
    if(i< p%n_layer) then
      p%G(i) = (p%temp_soil_p(i-1)-p%temp_soil_p(i))/(0.5*p%d_soil(i-1)/p%lambda(i-1)+0.5*p%d_soil(i)/p%lambda(i))
    else
      p%G(i) = p%lambda(i-1)*(p%temp_soil_p(i-1)-Tair_annual_mean)/(DEPTH_F+p%d_soil(i-1)*0.5)
    end if

    if(p%G(i) > 200) p%G(i) = 200
    if(p%G(i) <-200) p%G(i) = -200
 end do

 do i = 0,p%n_layer-1
    p%temp_soil_c(i) = p%temp_soil_p(i)+(p%G(i)-p%G(i+1)+S)/(p%Cs(i)*p%d_soil(i))*real(period_in_seconds)
    if(p%temp_soil_c(i) > 50.0) p%temp_soil_c(i)  = 50.
    if(p%temp_soil_c(i) < -50.0) p%temp_soil_c(i) = -50.

    ! avoid floating point overflows @ Xing Xiuli
    if(p%temp_soil_c(i) < 1.e-6 .and. p%temp_soil_c(i)>0.) p%temp_soil_c(i) = 0.
    if(p%temp_soil_c(i) > -1.e-6 .and. p%temp_soil_c(i)<0.) p%temp_soil_c(i) = 0.
 end do
! do i = 0,p%n_layer-1
!  write(*,*) 'DG0031: soil temperature diagnosis', p%temp_soil_c(i)
! end do

 call Update_ice_ratio(p)

 do i = 0,p%n_layer-1
    p%temp_soil_p(i) = p%temp_soil_c(i)
 end do

 return
 end subroutine

 subroutine Update_Cs(p)
 implicit none
 type(soil) :: p
 real(r8)   :: term1,term2,term3

 !Chen B. (2007) Ecological Modelling 209, 277-300  (equation 18)
 do i = 0,p%n_layer-1
    term1  = 2.*1.e3*p%density_soil(i)/2.65
    term2  = 1.e6*p%thetam(i)*(4.2*(1.-p%ice_ratio(i))+2.09*p%ice_ratio(i))
    term3  = 2.5*1.e6*p%f_org(i)
    p%Cs(i) = term1 + term2 + term3
 end do

 return
 end subroutine

 subroutine Update_ice_ratio(p)
 implicit none
 type(soil)  :: p
 real(r8)    :: Lf0 = 3.34*1.e5  ! latent heat of fusion at 0C
 real(r8)    :: tmp

 do i = 0,p%n_layer-1
   ! starting to frozen
   if(p%temp_soil_p(i) >= 0. .and. p%temp_soil_c(i) < 0. .and. p%ice_ratio(i) < 1.0 .and. p%thetam(i) >0.) then
   !! Add p%thetam(i) >0. by @J.Wang
     tmp = (0.-p%temp_soil_c(i))*p%Cs(i)*p%d_soil(i)
     p%ice_ratio(i) = p%ice_ratio(i) + tmp/Lf0/1000./(p%thetam(i)*p%d_soil(i))
     p%ice_ratio(i) = min(1.0,p%ice_ratio(i))
     p%temp_soil_c(i) = 0
   ! Be melting
   else if(p%temp_soil_p(i) <= 0 .and. p%temp_soil_c(i) > 0. .and. p%ice_ratio(i) > 0.) then
     tmp = (p%temp_soil_c(i) - 0.0)*p%Cs(i)*p%d_soil(i)
     p%ice_ratio(i)  = p%ice_ratio(i) - tmp/Lf0/1000./(p%thetam(i)*p%d_soil(i))
     p%ice_ratio(i)  = max(0.,p%ice_ratio(i))
     p%temp_soil_c(i) = 0
   end if

   p%ice_ratio(i)  = p%ice_ratio(i)*p%thetam_prev(i)/p%thetam(i)
   p%ice_ratio(i)  = min(1.,p%ice_ratio(i))
 end do

 return
 end subroutine

 subroutine UpdateSoilThermalConductivity(p)
!to calculate thermal conductivity of each soil layer, advances in water resources 26(2003), 79-93*/
 implicit none
 type(soil) :: p
 real(r8)   :: ki = 2.1     ! the thermal conductivity of ice
 real(r8)   :: kw = 0.61    ! the thermal conductivity of water
 real(r8)   :: tmp1,tmp2,tmp3,tmp4

 do i  = 0,p%n_layer-1
   tmp1 = p%thermal_cond(i)**(1-p%fei(i))  ! dry
   tmp2 = ki**(1.2*p%thetam(i)*p%ice_ratio(i)) !ice
   tmp3 = kw**(p%thetam(i)*(1-p%ice_ratio(i))) !water
   tmp4 = p%thetam(i)/p%fei(i)  !!Sr ??

   p%lambda(i)  = (tmp1*tmp2*tmp3-0.15)*tmp4+0.15  !eq. 8. LHE
   p%lambda(i)  = max(p%lambda(i),0.15)
 end do

 return
 end subroutine


!******************************************************
! Soil hydraulic state
!******************************************************
 subroutine UpdateSoilMoisture(p)  !remove the parameter 'kstep'@J.Wang
 ! Last revision: may 20, 2015, by LHE
 ! Given the current condition to calculate soil moisture after a period.
 ! Richards equation. Source: ET and rain
 ! kstep is definee in beps_par.F90: the total second in this step ( the period )
 ! kkk (outside of the fucntion) : step within an hour or half hour measurement.
 implicit none
 type(soil)::p
 real(r8) :: Infil, Infil_max  ! infiltration, and Maximum infiltration
 real(r8) :: this_step,total_t,max_Fb,kkstep
 real(r8) :: d1

 kkstep = 1.*kstep
 do i=0,p%n_layer
    p%thetam_prev(i) = p%thetam(i)  !save previous thetam
 end do

do i=0,p%n_layer
   if(p%temp_soil_c(i) >0.0) then
      p%f_ice(i) = 1.0   !! f_ice should be named as f_water
   else if(p%temp_soil_c(i) < -1.) then
      p%f_ice(i) = 0.1
   else
      p%f_ice(i) = 0.1+0.9*(p%temp_soil_c(i) + 1.0)
   end if
 end do

!write(*,*) p%f_ice(0),p%Zp,p%r_rain_g

 !! juweimin
 !! this part solve the upper boundary condition (Infiltration). LHE
 !! the maximum Infiltration. The Inf should be changing very fast during precipitation because thetam
 !! is changing. LHE
 Infil_max = p%f_ice(0)*p%Ksat(0)*(1.+(p%fei(0)-p%thetam_prev(0))/p%d_soil(0)*p%psi_sat(0)*p%b(0)/p%fei(0))
 Infil     = max(p%f_ice(0)*(p%Zp/kkstep+p%r_rain_g),0.)
 Infil     = min(Infil_max,Infil)
 Infil     = max(0.,Infil)

 p%Zp    = (p%Zp/kkstep + p%r_rain_g - Infil)*kkstep*p%r_drainage ! Ponded water after runoff. This one is related to runoff

 this_step = 0.
 total_t = 0.
 max_Fb = 0.

 do while (total_t < kkstep)
     do i = 0,p%n_layer-1
        p%km(i) = p%f_ice(i)*p%Ksat(i)*(p%thetam(i)/p%fei(i))**(2.*p%b(i)+3.)
     end do
     do i= 0,p%n_layer-1
        if(i < p%n_layer-1) then
           p%thetab(i)  = (p%thetam(i+1)/p%d_soil(i+1)+p%thetam(i)/p%d_soil(i))/(1./p%d_soil(i)+1./p%d_soil(i+1))
        else
           d1 = (p%thetam(i) - p%thetab(i-1))*2./p%d_soil(i)
           d1 = max(d1,0.)
           p%thetab(i)  = p%thetam(i) + d1*p%d_soil(i)/2.
           p%thetab(i)  = min(p%thetab(i),p%fei(i))
        end if
     end do

     do i = 0,p%n_layer-1
       if(i<p%n_layer-1) then  ! the unsaturated hydraulic conductivity at soil lower boundary
          p%Kb(i) = p%f_ice(i)*(p%Ksat(i)*p%d_soil(i)+p%Ksat(i+1)*p%d_soil(i+1))/(p%d_soil(i)+p%d_soil(i+1))* &
                    (p%thetab(i)/p%fei(i))**(2.*p%b(i)+3.)    !! Note: Kb(0) to Kb(n_layer-1) are not used in the model
       else    ! i= n_layer-1
          p%Kb(i) = 0.5*p%f_ice(i)*p%Ksat(i)*(p%thetab(i)/p%fei(i))**(2.*p%b(i)+3.)
       end if
     end do

     ! the unsaturated soil water retention
     do i=0,p%n_layer-1
        p%psim(i)  = p%psi_sat(i)*(p%thetam(i)/p%fei(i))**(-p%b(i))
        p%psim(i)  = max(p%psi_sat(i),p%psim(i))   !! I see no neccesity to use this line unless thetam > fei
     end do

     ! the unsaturated soil water retention @boundary LHE
     do  i = 0,p%n_layer-1
       p%psib(i) = p%psi_sat(i)*(p%thetab(i)/p%fei(i))**(-1.*p%b(i))
       p%psib(i) = max(p%psi_sat(i),p%psib(i))
     end do

     ! the unsaturated hydraulic conductivity of soil p%n_layer
     do i = 0,p%n_layer-1
       if(i<p%n_layer-1) then
          p%KK(i) = (p%km(i)*p%psim(i)+p%km(i+1)*p%psim(i+1))/ &
                    (p%psim(i)+p%psim(i+1))*(p%b(i)+p%b(i+1))/(p%b(i)+p%b(i+1)+6)   ! see seller's
       else
          p%KK(i) = (p%km(i)*p%psim(i)+p%Kb(i)*p%psib(i))/ &
                    (p%psim(i)+p%psib(i))*p%b(i)/(p%b(i)+3)
       end if
     end do

     ! Fb flow speed, Dancy's law LHE
     do i = 0,p%n_layer-1
        if(i<p%n_layer-1) then
           p%r_waterflow(i) = p%KK(i)*(2*(p%psim(i+1) - p%psim(i))/(p%d_soil(i)+p%d_soil(i+1))+1)   ! downwards positive , +1 accounts for gravitational drainage LHE
        else
           p%r_waterflow(i) = 0.  ! from Ju
        end if
     end do

     ! check the r_waterflow further
     do i = 0,p%n_layer-2
        p%r_waterflow(i) = min((p%fei(i+1)-p%thetam(i+1))*p%d_soil(i+1)/kkstep + p%Ett(i+1),&
                         p%r_waterflow(i))
       if(abs(p%r_waterflow(i)) > max_Fb) max_Fb  = abs(p%r_waterflow(i))    ! find max_Fb for all p%LAYERS
     end do

     if(max_Fb > 1.e-5) then
        this_step  = 1.   ! determine the sub_step according to order of Fb emirically
     else if(max_Fb > 1.e-6) then
        this_step  = 30.
     else
!       this_step  = 360.
        this_step  = kkstep    !!@J.Wang replace 360 by kstep
     end if

     total_t  = total_t+this_step
     if(total_t > kstep) this_step = this_step - (total_t - kkstep)

     do i  =0,p%n_layer-1
       if(i==0) then
           p%thetam(i) = p%thetam(i) + (Infil*this_step - p%r_waterflow(i)*this_step - p%Ett(i)*this_step)/p%d_soil(i)
       else
           p%thetam(i) = p%thetam(i) + (p%r_waterflow(i-1)-p%r_waterflow(i)-p%Ett(i))*this_step/p%d_soil(i)
       end if

       p%thetam(i) = max(p%theta_vwp(i),p%thetam(i))
       p%thetam(i) = min(p%fei(i),p%thetam(i))
       end do
 end do    ! end do while

 do i=0,p%n_layer-1
    p%ice_ratio(i) = p%ice_ratio(i)*p%thetam_prev(i)/p%thetam(i)
    p%ice_ratio(i) = min(1.0,p%ice_ratio(i))
 end do

 return
end subroutine

 subroutine Soil_water_uptake(lai,Hp,a,b,c,p,Trans_o,Trans_u,Evap_soil,vod,fei_leaf)
 implicit none
 type(soil) :: p
 real(r8)   :: lai,Hp,prlsp,vod,deltal_Sp,k_rs
 real(r8)   :: Trans_o,Trans_u,Evap_soil,a,b,c
 real(r8)   :: Source!,a,b,c
 real(r8)   :: z_depth(0:MAX_LAYERS-1)
 real(r8)   :: thetaox, pox, Sox, tWA, Ttrig, tWB, fei_c,p1,p2,ftheta,theta_Amin
 real(r8)   :: r_xylem,r_r,Lr,deltal_min,deltal_max,p_delta,f_deltal,p_excess
 real(r8)   :: qupt_sum,q1,q2,q3,ppsl,ppslh,fpmax,fei_leaf,fei_min,fei_th,f_feil,Eta !pmax,
 real(r8)   :: f_T(0:MAX_LAYERS-1),f_theta(0:MAX_LAYERS-1),rp(0:MAX_LAYERS-1)
 real(r8)   :: r_delta(0:MAX_LAYERS-1),rs(0:MAX_LAYERS-1),qupt(0:MAX_LAYERS-1)

!!! In this part, the transpiration for plants was calculated using the SPAC approach,
!!! we adopted the Darcy's law for calculating plant hydraulics as done in CoupModel,
!!! this will result in the simultion of plant water content or leaf potential, and
!!! will be linked to VOD data at daily scale. @MOUSONG, 20221106
!!! Below are parameters for the SPAC-based water uptake modeling used in CoupModel,
!!! units are converted to fit BEPS.

 !Hp = param(29)
 !f_deltal=0.
 k_rs=10.**(-3)                 !need to be calibrated
 !deltal_Sp=0.
 !vod=0.
 !fei_leaf=0.
 pox = 4.
 fei_c = 400./100.            ! cm water to m water
 p1 = 400*10.**(-3)/24./3600.     ! 1/d to /s and covert mm to m for transpiration
 p2 = 0.1/1000./24./3600.     ! kg/(m2 d) to m/s
 fei_min = 15000./100.        ! cm water to m water
 ppsl = 1./1000.              ! mm to m
 ppslh = 0.5/1000.            ! mm/m to m/m
 r_r = 1000.*24.*3600.        ! d/m to s/m
 r_xylem = 1.*24.*3600.       ! d/m to s/m
 p_delta = 0.5                ! m2
 deltal_max = 0.01            ! m
 deltal_min = 0.001           ! m
 tWA = 0.8
 tWB = 0.
 p_excess = 2.0/1000./24./3600.  ! mm/d to m/s
 theta_Amin = 5./100.            ! % to -
 Ttrig = 15.                      ! oC
 qupt_sum = 0.
 qupt(0:4) = (/0.,0.,0.,0.,0./)       ! m/s
 fei_th = 1000./100.                  ! leaf threshold suction, cm to m water
 !q1 = 0.
 !q2 = 0.
 !q3 = 0.
 prlsp = 0.0001                   ! specific root length, gC/m
 !Lr = 1000.*p%biomass_root/prlsp    ! root lenght,0.1 m/m2, calculated from root biomass of PFT
 Lr = 0.1
 !a = 0.3                          ! These three para. need further tuning, [0,50],[0,20],[0,50]
 !b = 0.64
 !c = 0.04

 Source  = Trans_o+Trans_u        ! g/s    transpiration calculated from energy balance by original BEPS

! implement the root distribution function to calcualte relative water uptake for each layer

 z_depth(0)   =   p%d_soil(0)
 p%f_root(0)  =   1-p%r_root_decay**(z_depth(0)*100.)

 do i  = 1,p%n_layer-1    ! change to n_layer-1 from n_layer-2, a bug@MOUSONG,20221105
   z_depth(i)  = z_depth(i-1)+p%d_soil(i)
   p%f_root(i)   = p%r_root_decay**(z_depth(i-1)*100.) - p%r_root_decay**(z_depth(i)*100.)
 end do

 p%f_root(p%n_layer-1) = p%r_root_decay**(z_depth(p%n_layer-2)*100.)

 !fpmax = ppsl*LAI
 fpmax = ppslh*lai*Hp    ! maximum water storage, estimated from LAI and canopy height, m

! Estimate water uptake for each layer using the Darcy's approach

 do i = 0,p%n_layer-1
    ! water stress function
    thetaox = p%fei(i) - theta_Amin
    Sox = (p%thetam(i)-thetaox)/(p%fei(i)-thetaox)
    Sox = min(1.,Sox)
    Sox = max(Sox,1.e-6)
    ftheta = 10**(-pox*Sox)
    f_theta(i) = min((fei_c/p%psim(i))**(p1*Source/rho_w+p2),ftheta)
    !write(*,*), 'i =', i+1
    !write(*,*), 'p%r_root_decay =', p%r_root_decay
    !write(*,*)  "f_theta(i) = ",f_theta(i)
    ! temperature stress function
    f_T(i) = 1 - exp(-tWA*(max(0.,p%temp_soil_p(i)-Ttrig))**tWB)
    !write(*,*)  "f_T(i) = ",f_T(i)
    ! plant resistance
    !write(*,*)  "p%f_root = ",p%f_root(i)
    !write(*,*)  "Lr = ",Lr
    rp(i) = (r_xylem*Hp/p%f_root(i) + r_r/(Lr*p%f_root(i)))*(1.0/f_T(i))*(1.0/f_theta(i))
    !write(*,*)  "rp(i) = ",rp(i)
    !write(*,*)  "p%f_root(i) = ", p%f_root(i)
    !write(*,*)  "p%d_soil(i) = ", p%d_soil(i)
    ! soil-root resistance
    r_delta(i) = Lr*p%f_root(i)/p%d_soil(i)
    !write(*,*)  "r_delta(i) = ",r_delta(i)
    !write(*,*)  "deltal_min = ",deltal_min
    !write(*,*)  "deltal_max = ",deltal_max
    !write(*,*)  "-p_delta = ",-p_delta
    !write(*,*)  "p%f_root(i) = ", p%f_root(i)
    !write(*,*)  "-p_delta*r_delta(i) = ", -p_delta*r_delta(i)
    !write(*,*)  "exp(-p_delta*r_delta(i)) = ", exp(-p_delta*r_delta(i))
    f_deltal = deltal_min + (deltal_max - deltal_min)*exp(-p_delta*r_delta(i))

    !write(*,*)  "p%KK(i) = ", p%KK(i)
    !write(*,*)  "p%f_root(i) = ", p%f_root(i)
    !write(*,*)  "*************** "
    rs(i) = k_rs*f_deltal/(p%KK(i)*p%f_root(i))
    !write(*,*)  "i = ",i
    !write(*,*)  "p%KK*p%f_root = ",p%KK(i)*p%f_root(i)
    !rs(i) = 1000.
    !write(*,*)  "f_deltal = ",f_deltal
    !write(*,*)  "p%f_root(i) = ",p%f_root(i)
    !write(*,*)  "rs(i) = ",rs(i)
    ! leaf water potential
    !write(*,*)  "fpmax = ",fpmax
    !write(*,*)  "p%Sp = ",p%Sp
    !write(*,*)  "fei_min = ",fei_min
    !write(*,*)  "Hp = ",Hp

    fei_leaf = (1-p%Sp/fpmax)*(fei_min + Hp) - Hp

    !write(*,*)  "fei_leaf = ",fei_leaf
    !write(*,*)  "psim = ", p%psim(i)
    !write(*,*)  "z_depth", z_depth(p%n_layer-1)
    ! water uptake as the minimum of three terms
    !q1 = p%f_root(i)*( p%psim(i) - fei_leaf - (Hp+z_depth(p%n_layer-1)) )/(rp(i) + rs(i))
    q1 = p%f_root(i)*(fei_leaf -p%psim(i) - (Hp+z_depth(p%n_layer-1)) )/(rp(i) + rs(i))
    q1=max(0.,q1)
    !write(* ,*)  "q1 = ",q1
    q2 = p%f_root(i)*Source/rho_w + p_excess*p%f_root(i)
    !write(*,*)  "q2 = ",q2
    q3 = fpmax*p%f_root(i) - p%Sp*p%f_root(i)
    !write(*,*)  "q3 = ",q3
    q3=max(0.,q3)
    qupt(i) = min(q1,q2)
    qupt(i) = min(qupt(i),q3)
    !qupt(i) =q2

    !write(*,*)  "qupt(i) = ",qupt(i)

 end do

 ! calculate vod based on lai and leaf water potential, based on Liu et al., 2021
 !"Global ecosystem-scale plant hydraulic traits retrieved using modelÃƒâ€šÃ‚Â¨Cdata fusion"

 vod = (a + b*lai)*(1 + c*fei_leaf/101.)

! Calculate the actual transpiration with relation to leaf water potential

 f_feil = min((fei_leaf - fei_min)/(fei_th - fei_min),1.)
 f_feil = max(0.,f_feil)
  !write(*,*)  "f_feil = ",f_feil
 Eta = f_feil * Source/rho_w
! write(*,*)  "Eta = ",Eta
! Update the plant water storage
 do i=0,p%n_layer-1
    qupt_sum = qupt_sum + qupt(i)
 end do
 !write(*,*)  "qupt_sum = ",qupt_sum

 deltal_Sp =  -(Eta - qupt_sum) * kstep*1.
 !if (deltal_Sp>0) then
 !   deltal_Sp=min(10.**(-3),deltal_Sp)
 !else
 !   deltal_Sp=min(-10.**(-3),deltal_Sp)
 !end if
 p%Sp=p%Sp+deltal_Sp
 !write(*,*)  "deltal_Sp = ",deltal_Sp

! Blow is the code for original BEPS
 ! for the top layer
 !p%Ett(0) = (Source/rho_w)*p%dt(0) + Evap_soil/rho_w
 p%Ett(0) = min((Source/rho_w)*p%dt(0),qupt(0)) + Evap_soil/rho_w
!  p%Ett(0) = 0.
 ! for each layer
 do i = 1,p%n_layer-1
!    p%Ett(i) = min(Source/rho_w*p%dt(i),qupt(i))
    p%Ett(i) = min(Source/rho_w*p%dt(i),qupt(i))
 end do

 return
 end subroutine

 subroutine soil_water_factor_v2(p)
 ! Compute soil water stress factor
 ! Last revision by LHE. May 22. 2015
 ! Rewritten by : Liming He, Jan 29,2013
 ! Modified by: Mustapha El Maayar - March 2008
 ! Written by : Weiming Jun
 ! Fortran version : J. Wang April 2017
 implicit none
 type(soil) :: p
 real(r8):: ft(0:MAX_LAYERS-1), fpsisr(0:MAX_LAYERS-1)
 real(r8):: dtt(0:MAX_LAYERS-1)
 real(r8):: t1 = -0.02,t2=2.0
 !--iLab::changed in order to avoid "implicit save" for these variables
 ! real(r8):: dtt_sum = 0.,fpsisr_sum = 0.
 real(r8) :: dtt_sum, fpsisr_sum
 dtt_sum = 0._8
 fpsisr_sum = 0._8

 !! change the rule for updating p%psim @MOUSONG.WU,2018.11

 do i = 0,p%n_layer-1
    p%psim(i) = p%psi_sat(i)*(p%thetam(i)/p%fei(i))**(-p%b(i))
    p%psim(i) = max(p%psi_sat(i),1.e-6)
 end do

 do i = 0,p%n_layer-1
    if(p%psim(i) > p%psi_min) then
      fpsisr(i) = 1./(1. + ((p%psim(i) - p%psi_min)/p%psi_min)**p%alpha) !eq. 10 only f(psim) and eq.11 in Ju 2006
    else
      fpsisr(i) = 1.
    end if

    if(p%temp_soil_p(i) > 0.) then
       ft(i) = (1.-exp(t1*p%temp_soil_p(i)**t2)) !eq. 12 in Ju 2006
    else
       ft(i) = 0.
    end if

    fpsisr(i) = fpsisr(i)*ft(i)
 end do

 if(FW_VERSION ==1) then
    do i=0,p%n_layer-1
      dtt(i) = p%f_root(i)*fpsisr(i) !eq. 14 in Ju 2006
    end do
 else
    do i = 0,p%n_layer-1
      dtt(i) = p%f_root(i)
    end do
 end if

 do i=0,p%n_layer-1
   dtt_sum = dtt_sum + dtt(i)
 end do

 if(dtt_sum < 1.e-6) then
   p%f_soilwater = 0.1
   do i = 0,p%n_layer-1
      p%dt(i) = 0.
   end do
 else
   do i = 0,p%n_layer-1
      p%dt(i) = dtt(i)/dtt_sum ! eq. 14, in  JU
      p%dt(i) = max(p%dt(i),1.e-6)
      !if(isnan(p%dt(i))) then
      !p%dt(i) = 0.
      !write(*,*) p%dt(i)
      !end if
   end do

   do i =0,p%n_layer-1
      fpsisr_sum = fpsisr_sum + fpsisr(i)*p%dt(i) ! eq. 12, in Chen 2012 GBC; eq 15 in JU
   end do

   p%f_soilwater  = max(0.1,fpsisr_sum)
 end if

 return
end subroutine


 subroutine Soil_evaporation(temp_air,temp_g,rh_air,netRad_g,Gheat_g,percent_snow_g,depth_water,&
                             depth_snow,mass_water_g,mass_snow_g,density_snow,swc_g,porosity_g,evapo_soil,&
                             evapo_water_g,evapo_snow_g)
 ! this module calculate evap from ground surface/top soil, and evaporation of snow and pond water on surface
 ! eidtted by XZ Luo, May 25,2015

 ! input
 ! air temperature,ground surface temperature,relative humidity of ground (BEPS take it as the air RH)
 ! percentage of snow cover, depth of water/snow, soil water content on first soil layer, porosity of first layer

 ! OUTPUT
 ! evap from soil suface
 ! depth of water and snow on ground after evaporation and sublimation
 implicit none
 real(r8),intent(in) :: temp_air,temp_g,rh_air
 real(r8),intent(in) :: netRad_g  ! net radiation on ground
 real(r8),intent(in) :: Gheat_g   ! aerodynamic conductantce of heat on ground surface
 real(r8),intent(inout) :: percent_snow_g
 real(r8),intent(inout) :: depth_water,depth_snow  ! depth of water and snow on ground after rainfall/snowfall stage1 befor evap
                                                   ! output after substacting evps
 real(r8),intent(inout) :: mass_water_g,mass_snow_g
 real(r8),intent(in)    :: density_snow
 real(r8),intent(in)    :: swc_g,porosity_g ! soil water content (from last step) and porosity on ground
 real(r8),intent(out)   :: evapo_soil,evapo_water_g,evapo_snow_g

 real(r8) :: density_air_g,cp_air_g,vpd_g,slope_vapor_g,psy_air_g
 real(r8) :: Gwater_g !conductance of water on soil surface
 real(r8) :: lantent_water,latent_snow
 real(r8) :: density_water
 real(r8) :: length_step
 call meteo_pack(temp_g,rh_air)
 density_air_g  = density_air
 cp_air_g       = cp_air
 vpd_g          = vpd
 slope_vapor_g  = slope_vapor
 psy_air_g      = psy

 latent_water   = (2.501-0.00237*temp_air)*1e6
 latent_snow    =  2.83*1e6
 density_water  = rho_w
 length_step    = kstep

 ! adjust the rs due to CO2 impacts in non-water-limited areas, according to Yang et al., 2019, Nature Climate Change
 ! rs = rs_300*(1+S_rs(CO2-300)), rs_300=55 s m-1, S_rs=0.09% ppm-1, Mousong.Wu@2019.04
 if(swc_g/porosity_g < 0.5) then
     Gwater_g       = 1./(4.0*exp(8.2-4.2*swc_g/porosity_g))
 else
     Gwater_g       = 1./(55*(1+0.09/100.0*(CO2_air-300.)))
 end if

 ! get the percentage of snow
 if(depth_snow > 0.02) then
     percent_snow_g = 1.
 else
     percent_snow_g = mass_snow_g/(0.025*density_snow)
 end if
 percent_snow_g      = max(percent_snow_g,0.)
 percent_snow_g      = min(percent_snow_g,1.)

 ! when there are pond water on ground, there is evaporation from the water
 if(depth_water > 0 .and. depth_snow ==0) then
    evapo_water_g = 1./latent_water*(slope_vapor_g*(netRad_g*0.8-0)+density_air_g*cp_air_g*vpd_g*Gheat_g)/&
                    (slope_vapor_g+psy_air_g*(1+Gheat_g/0.01))
 else
    evapo_water_g = 0
 end if
 evapo_water_g = max(0.,evapo_water_g)
 if(evapo_water_g>0) evapo_water_g = min(evapo_water_g,depth_water*density_water/length_step)

 depth_water  = depth_water - evapo_water_g/density_water*length_step
 depth_water  = max(0.,depth_water)
 mass_water_g = mass_water_g - evapo_water_g*length_step

 ! when there are snow on ground, there s ony evaporation from the snow
 if(depth_snow > 0) then
   evapo_snow_g = 1./latent_snow*(slope_vapor_g*(netRad_g*0.8-0)+density_air_g*cp_air_g*vpd_g*Gheat_g)/&
                  (slope_vapor_g+psy_air_g*(1+Gheat_g/0.01))*percent_snow_g
 else
   evapo_snow_g = 0
 end if
 evapo_snow_g  = max(0.,evapo_snow_g)
 if(evapo_snow_g > 0) evapo_snow_g  = min(evapo_snow_g,mass_snow_g/length_step)
 mass_snow_g   = mass_snow_g - evapo_snow_g*length_step
 mass_snow_g   = max(mass_snow_g,0.)

 if(mass_snow_g >0) then
   depth_snow  = depth_snow - evapo_snow_g/density_snow*length_step
 else
   depth_snow  = 0
 end if

 if(depth_water >0 .or. depth_snow >0) then
   evapo_soil = 0
 else
   evapo_soil = (1.-percent_snow_g)*1/latent_water*(slope_vapor_g*(netRad_g-0)+density_air_g*cp_air_g*vpd_g*Gheat_g)/ &
                (slope_vapor_g +psy_air_g*(1+Gheat_g/Gwater_g))
   evapo_soil = max(0.,evapo_soil)
 end if
 return
 end subroutine

end module
