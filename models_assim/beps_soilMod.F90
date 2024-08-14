!*******************************************************
! Function: Soil module( Initialize/Update soil status,soil
!           water, soil thermal etc.)
! Created : Jun Wang
! Date    : 2016/12/5
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
    real(r8)  :: f_feileaf
    real(r8)  :: Sp
    real(r8)  :: biomass_root
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

! parameters related to plant hydraulics
type,public:: Phydraulic
    real(r8) :: theta_Amin
    real(r8) :: pox
    real(r8) :: fei_c
    real(r8) :: spac_p1
    real(r8) :: spac_p2
    real(r8) :: tWA
    real(r8) :: tWB
    real(r8) :: Ttrig
    real(r8) :: r_xylem
    real(r8) :: r_r
    real(r8) :: Lr
    real(r8) :: deltal_min
    real(r8) :: deltal_max
    real(r8) :: p_delta
    real(r8) :: ppslh
    real(r8) :: fei_min
    real(r8) :: fei_th
    real(r8) :: p_excess
    ! currently I do not consider the effect of leaf temperature @ Lu HU 2024/03/19
    real(r8) :: Tleaf_H
    real(r8) :: Tleaf_L
    real(r8) :: Tleaf_O
    real(r8) :: fei50
    real(r8) :: plc_a
    real(r8) :: plc_min
end type Phydraulic

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
 subroutine Init_soil_parameters(lc,stxt,r_root_decay,p)
 implicit none
 integer,intent(in)  :: lc
 integer,intent(in)  :: stxt
 real(r8),intent(in) :: r_root_decay
 type(soil)          :: p

 p%n_layer      = 5

 if(lc == 3 .or. lc == 4) then
    p%psi_min  = 10.0      ! for fw
    p%alpha    = 1.5
 else
    p%psi_min  = 33.0
    p%alpha    = 0.4
 end if

 p%d_soil(0:4)      = (/0.05,0.10,0.20,0.40,1.25/)    ! depth_layer
 p%r_root_decay     = r_root_decay

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
           !/*  calculating aboveground biomass based on LAI  J. Liu 2002 */
        biomass=0.9097*lai_yr+0.125*lai_yr*lai_yr
        biomass_stem_u=0.01*biomass_stem_o     !/* stem C of understory */
        biomass_root_u=0.01*biomass_root_o  !/* root C of understory */
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
        biomass_leaf_o=0.1*biomass    !/* leaf C of overstory */
        biomass_stem_o=0.90*biomass    !/* stem C of overstory */
        biomass_root_o=1.432*biomass**0.639    !/* root C of overstory  Kurz 1996 */
        biomass_leaf_u=0.3*biomass_leaf_o     !/* leaf C of understory */
        biomass_stem_u=0.01*biomass_stem_o    ! /* stem C of understory */
        biomass_root_u=0.01*biomass_root_o    !/* root C of understory */
        p%biomass_root = biomass_root_o + biomass_root_u
else if(lc == 14 .or. lc == 15 .or. lc ==25 .or. lc ==40 .or. lc==41) then
        biomass_leaf_o=0.05*lai_yr  ! /* leaf C = lai/20  from W.Ju 05y11*/
        biomass_stem_o=0.0          !/* stem C */
        biomass_root_o=0.061*lai_yr    !/* root C = lai/20*0.55/0.45  from W.Ju 05y11*/
        biomass_leaf_u=0.0
        biomass_stem_u=0.0
        biomass_root_u=0.0
        p%biomass_root = biomass_root_o + biomass_root_u
else
        biomass_leaf_o=0.05*lai_yr  ! /* leaf C = lai/20  from W.Ju 05y11*/
        biomass_stem_o=0.0          !/* stem C */
        biomass_root_o=0.061*lai_yr    !/* root C = lai/20*0.55/0.45  from W.Ju 05y11*/
        biomass_leaf_u=0.0
        biomass_stem_u=0.0
        biomass_root_u=0.0
        p%biomass_root = biomass_root_o + biomass_root_u

end if

 return
 end subroutine

subroutine Init_planthydraulics_para(phydra)
    implicit none
    type(Phydraulic) :: phydra
    phydra%theta_Amin=5./100.
    phydra%pox=4.
    phydra%fei_c=400./100.
    phydra%spac_p1=0.3 ! checked for coup model manual
    phydra%spac_p2=0.1
    phydra%tWA=0.8    ! Mellander et al., 2006, Modelling the effect of low soil temperatures on transpiration by Scots pine
    phydra%tWB=0.4
    phydra%Ttrig=0.  ! [-2 2] Wu et al., 2012, The role of air and soil temperatures in the seasonality of photosynthesis and transpiration in a boreal Scots pine ecosystem
    phydra%r_xylem=1.
    phydra%r_r=1000.
    phydra%Lr=100.
    phydra%deltal_min=0.001
    phydra%deltal_max=0.01
    phydra%p_delta=0.5
    phydra%ppslh=0.5/1000.
    phydra%fei_min=15000./100.
    phydra%fei_th=1000./100.
    phydra%p_excess=2.0/1000./24./3600.
    phydra%Tleaf_H=35.  ! currently do not consider
    phydra%Tleaf_L=5.    ! currently do not consider
    phydra%Tleaf_O=20.   ! currently do not consider
    phydra%fei50 = -1.5  ! MPa, when stomatal conductance loss half
    phydra%plc_a= 3.0  ! shape parameter
    phydra%plc_min = 0.1 ! minimum plc when critical or minimum leaf water potential occurs

 return
 end subroutine

 subroutine Init_soil_status(p,Tsoil,Tair,Ms,snowdepth)
 ! soil temperatures and moisutre for each layer
 ! ponded water,snow depth
 implicit none
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
 p%Sp = 1.e-6
 
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
 p%f_feileaf = 1.0
! the calculation of ice_ratio can be improved based on SFCC
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

 do i  = 1,p%n_layer-2
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
    p%thetam(i) = p%thetam_prev(i)  !save previous thetam
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

 ! 2023/12/13 after checked and add to update thetam and psim
 do i=0,p%n_layer
    p%thetam_prev(i) = p%thetam(i)  !save previous thetam
 end do

 do i=0,p%n_layer
    p%psim_prev(i) = p%psim(i)  !save previous thetam
 end do

 return
end subroutine

 !subroutine Soil_water_uptake(p,Trans_o,Trans_u,Evap_soil)
 !implicit none
 !type(soil) :: p
 !real(r8)   :: Trans_o,Trans_u,Evap_soil
 !real(r8)   :: Source
!
! Source  = Trans_o+Trans_u
!
! ! for the top layer
! p%Ett(0) = (Source/rho_w)*p%dt(0) + Evap_soil/rho_w
!!  p%Ett(0) = 0.
! ! for each layer
! do i = 1,p%n_layer-1
!   p%Ett(i) = Source/rho_w*p%dt(i)
! end do
!
! return
! end subroutine

subroutine Soil_water_uptake(fw_flag,p,Trans_o,Trans_u,Evap_soil,lai,Hp,a,b,c,phydra,&
    vod,fei_leaf,qupt_sum,Eta)
 implicit none
 type(soil) :: p
 real(r8)   :: Trans_o,Trans_u,Evap_soil,lai,Hp
 real(r8)   :: a,b,c
 type(Phydraulic), intent(in) :: phydra
 real(r8), intent(out) :: qupt_sum,Eta,vod,fei_leaf
 real(r8)   :: Source
 real(r8)   :: qupt(0:MAX_LAYERS-1)
 integer,intent(in)   :: fw_flag

 ! initialize the output variable: vod, fei_leaf, qupt_sum, ETa
 vod = 0.0
 fei_leaf = 0.0
 qupt_sum = 0.0
 Eta = 0.0

 Source  = Trans_o+Trans_u ! kg/(s*m2)
! 2024/03/12, flag for choosing the form of water stress on stomatal conductance or Vcmax
! 0: fws = 1.0, do not consider the water stress from soil or leaf water potential
! 1: fws = f_soilwater, soil moisture stress on BWB slope
! 2: fws = 1.0, vcmax = vcmax*f_feileaf, leaf water potential stress on Vcmax
! 3: fws = f_feileaf, leaf water potential stress on BWB slope
! 4: fws = 1.0, but f_feileaf only works on Etp, Eta = Etp * f_feileaf
! 5: fws = 1.0, vcmax = vcmax*f_soilwater, soil moisture stress on Vcmax
 if (fw_flag==2 .or. fw_flag==3 .or. fw_flag==4 .or. fw_flag==6) then

    call plant_hydraulics(fw_flag,lai,Hp,a,b,c,p,phydra,Source,vod,fei_leaf,qupt,qupt_sum,Eta)

 !.......considering the plant hydraulics..............
 ! for the first layer
    p%Ett(0) = min((Source/rho_w)*p%dt(0),qupt(0)) + Evap_soil/rho_w
    ! for each layer
    do i = 1,p%n_layer-1
       p%Ett(i) = min(Source/rho_w*p%dt(i),qupt(i))
    end do

 else
    !.....original BEPS p%Ett calculation................
    p%Ett(0) = (Source/rho_w)*p%dt(0) + Evap_soil/rho_w
    do i = 1,p%n_layer-1
       p%Ett(i) = Source/rho_w*p%dt(i)
    end do

 end if
 !.....................................................

 return
 end subroutine

! 2024/01/07 modified
 subroutine plant_hydraulics(fw_flag,lai,Hp,a,b,c,p,phydra,Trans,vod,fei_leaf,qupt,qupt_sum,Eta)
 implicit none
 type(soil) :: p
 type(Phydraulic),intent(in) :: phydra
 real(r8)   :: lai,Hp,Trans
 real(r8)   :: a,b,c
 real(r8)   :: Etp,deltal_Sp,thetaox, Sox,ftheta,f_deltal,fpmax,f_feil
 real(r8)   :: z_depth(0:MAX_LAYERS-1)
 real(r8)   :: theta_Amin, pox, spac_p1,spac_p2,tWA, tWB, Ttrig
 real(r8)   :: r_xylem,r_r,Lr,deltal_min,deltal_max,p_delta,ppslh,p_excess
 real(r8)   :: fei_c,fei_min,fei_th
 real(r8)   :: q1,q2,q3
 real(r8)   :: f_T(0:MAX_LAYERS-1),f_theta(0:MAX_LAYERS-1),rp(0:MAX_LAYERS-1)
 real(r8)   :: r_delta(0:MAX_LAYERS-1),rs(0:MAX_LAYERS-1)
 real(r8)   :: kkstep,Etp_mm
 real(r8)   :: twa_trig
 real(r8), intent(out) :: qupt_sum,Eta,vod,fei_leaf
 real(r8), intent(out) :: qupt(0:MAX_LAYERS-1)
 real(r8)   :: fei50, plc_a, plc_min
 integer,intent(in)   :: fw_flag

!!! In this part, the transpiration for plants was calculated using the SPAC approach,
!!! we adopted the Darcy's law for calculating plant hydraulics as done in CoupModel,
!!! this will result in the simultion of plant water content or leaf potential, and
!!! will be linked to VOD data at daily scale. @MOUSONG, 20221106
!!! Below are parameters for the SPAC-based water uptake modeling used in CoupModel,
!!! units are converted to fit BEPS.

!............default value from Couple model...........................
 !pox = 4.
 !fei_c = 400./100.            ! cm water to m water
 !spac_p1 = 0.3               ! 1/d
 !spac_p2 = 0.1                ! kg/(m2 d) equivalent to mm/d water
 !fei_min = 15000./100.        ! cm water to m water
 !ppslh = 0.5/1000.            ! mm/m to m/m
 !r_r = 1000.                  ! d/m
 !r_xylem = 1.                 ! d/m
 !p_delta = 0.5                ! m2
 !deltal_max = 0.01            ! m
 !deltal_min = 0.001           ! m
 !tWA = 0.8
 !tWB = 0.
 !p_excess = 2.0/1000./24./3600.  ! mm/d to m/s
 !theta_Amin = 5./100.            ! % to decimal
 !Ttrig = 15.                     ! oC
 !fei_th = 1000./100.             ! leaf threshold suction, cm to m water
 !Lr =0.1                         !Lr = 0.1  !default value
 !................................................................................

 !................. input parameter values............................................

 pox = phydra%pox
 !write(*,*), 'pox=',pox
 fei_c = phydra%fei_c                      ! m water
 !write(*,*), 'fei_c=',fei_c
 spac_p1 = phydra%spac_p1                  ! 1/d
 !write(*,*), 'spac_p1=',spac_p1
 spac_p2 = phydra%spac_p2                  ! kg/(m2 d) equivalent to mm/d water
 !write(*,*), 'spac_p2=',spac_p2
 fei_min = phydra%fei_min                  ! m water
 !write(*,*), 'fei_min=',fei_min
 ppslh = phydra%ppslh                      ! m/m
 !write(*,*), 'ppslh=',ppslh
 r_r = phydra%r_r                          ! d/m
 !write(*,*), 'r_r=',r_r
 r_xylem = phydra%r_xylem                  ! d/m
 !write(*,*), 'r_xylem=',r_xylem
 p_delta = phydra%p_delta                  ! m2
 !write(*,*), 'p_delta=',p_delta
 deltal_max = phydra%deltal_max            ! m
 !write(*,*), 'deltal_max=',deltal_max
 deltal_min = phydra%deltal_min            ! m
 !write(*,*), 'deltal_min=',deltal_min
 tWA = phydra%tWA
 !write(*,*), 'tWA=',tWA
 tWB = phydra%tWB
 !write(*,*), 'tWB=',tWB
 p_excess = phydra%p_excess                ! m/s
 !write(*,*), 'p_excess=',p_excess
 theta_Amin = phydra%theta_Amin            ! % to decimal
 !write(*,*), 'theta_Amin=',theta_Amin
 Ttrig = phydra%Ttrig                      ! oC
 !write(*,*), 'Ttrig=',Ttrig
 fei_th = phydra%fei_th                    ! leaf threshold suction, m water
 !write(*,*), 'fei_th=',fei_th
 Lr =phydra%Lr                             ! m/m2
 !write(*,*), 'Lr=',Lr
 !Tleaf_H =phydra%Tleaf_H                   ! degree C
 !write(*,*), 'Tleaf_H=',Tleaf_H
 !Tleaf_L =phydra%Tleaf_L
 !write(*,*), 'Tleaf_L=',Tleaf_L
 !Tleaf_O =phydra%Tleaf_O
 !write(*,*), 'Tleaf_O=',Tleaf_O
 !fei50 = phydra%fei50
 !plc_a = phydra%plc_a
 !plc_min = phydra%plc_min
 !..................Initialization................................................
 fpmax = 0.
 thetaox = 0.
 Sox = 0.
 ftheta = 0.
 f_theta(0:4) = (/0.,0.,0.,0.,0./)
 f_T(0:4) = (/0.,0.,0.,0.,0./)
 rp(0:4) = (/0.,0.,0.,0.,0./)
 r_delta(0:4) = (/0.,0.,0.,0.,0./)
 f_deltal = 0.
 rs(0:4) = (/0.,0.,0.,0.,0./)
 fei_leaf = 0.
 q1 = 0.
 q2 = 0.
 q3 = 0.
 qupt(0:4) = (/0.,0.,0.,0.,0./)       ! m/s
 vod = 0.
 f_feil = 0.
 Eta = 0.
 qupt_sum = 0.
 deltal_Sp = 0.

!2024/03/10 if using logistic function for f_feileaf, fei_min could be calculated by blow Eqs
!fei_min=fei50*(-101.0)*((1.0/plc_min-1)**(1.0/plc_a)), plc_min is the minimum degree of stomata reducing
!.............................................................................
 !prlsp = 0.0001                   ! specific root length, gC/m
 !Lr = 1000.*p%biomass_root/prlsp    ! root lenght,0.1 m/m2, calculated from root biomass of PFT, need to be improved in the future

 !a = 0.3                          ! These three para. need further tuning, [0,50],[0,20],[0,50]
 !b = 0.64
 !c = 0.04

! 2024/03/12, flag for choosing the form of water stress on stomatal conductance or Vcmax
! 0: fws = 1.0, do not consider the water stress from soil or leaf water potential
! 1: fws = f_soilwater, soil moisture stress on BWB slope
! 2: fws = 1.0, vcmax = vcmax*f_feileaf,leaf water potential stress on Vcmax
! 3: fws = f_feileaf,leaf water potential stress on BWB slope
! 4: fws = 1.0, but f_feileaf only works on Etp, Eta = Etp * f_feileaf
! 5: fws = 1.0, vcmax = vcmax*f_soilwater, soil moisture stress on Vcmax
! 6: fws = 1.0, vcmax = vcmax*min(f_soilwater,f_feileaf)
 if (fw_flag ==2 .or. fw_flag ==3 .or. fw_flag==6) then
    Etp  = Trans/p%f_feileaf
 else
    Etp  = Trans
 end if
 !write(*,*), 'Source= ',Source

! implement the root distribution function to calcualte relative water uptake for each layer

 z_depth(0)   =   p%d_soil(0)
 p%f_root(0)  =   1-p%r_root_decay**(z_depth(0)*100.)

! f_root root fraction for each layer
 do i  = 1,p%n_layer-1    ! change to n_layer-1 from n_layer-2, a bug@MOUSONG,20221105
   z_depth(i)  = z_depth(i-1)+p%d_soil(i) !soil depth from boundary surface to layer
   p%f_root(i)   = p%r_root_decay**(z_depth(i-1)*100.) - p%r_root_decay**(z_depth(i)*100.)
 end do

 p%f_root(p%n_layer-1) = p%r_root_decay**(z_depth(p%n_layer-2)*100.)

 !fpmax = ppsl*LAI
 fpmax = ppslh*lai*Hp    ! maximum water storage, estimated from LAI and canopy height, m water
 !write(*,*), 'fpmax =', fpmax
 !write(*,*)  "p%Sp = ",p%Sp
 !write(*,*)  "Sp/fpmax = ",p%Sp/fpmax

! Estimate water uptake for each layer using the Darcy's approach
 kkstep=1.*kstep
 Etp_mm=Etp/rho_w*10.**3*kkstep ! kg/(s*m2) / (kg/m3) * s = mm water
!write(*,*), 'Source_mm= ',Source_mm
!write(*,*), 'p%psim(i) =', p%psim(0)
!write(*,*), 'p%psim_prev(i) =', p%psim_prev(0)

 do i = 0,p%n_layer-1
    !write(*,*), 'i =', i+1
    ! water stress function
    !write(*,*), 'p%fei(i) =', p%fei(i)
    thetaox = p%fei(i) - theta_Amin
    !write(*,*), 'thetaox =', thetaox
    !write(*,*), 'p%fei =', p%fei(i)
    !write(*,*), 'p%thetam_prev=', p%thetam_prev(i)
    Sox = (p%thetam_prev(i)-thetaox)/(p%fei(i)-thetaox)

    !write(*,*), 'Sox =', Sox
    Sox = min(1.,Sox)
    Sox = max(Sox,1.e-6)
    !write(*,*), 'Sox =', Sox
    ftheta = 10.**(-pox*Sox)
    ! 2024/03/01
    !ftheta = max(ftheta, 0.1)
    ftheta = max(ftheta, 0.01)
    !write(*,*), 'ftheta =', ftheta
    !write(*,*), 'p%psim_prev(i) =', p%psim_prev(i)
    !write(*,*), 'fei_c/p%psim =', fei_c/p%psim(i)
    !write(*,*), 'p1*Source_mm+p2 =', p1*Source_mm+p2
    f_theta(i) = min((fei_c/p%psim_prev(i))**(spac_p1*Etp_mm+spac_p2),ftheta)
    !write(*,*)  "f_theta(i) = ",f_theta(i)
    ! soil temperature stress function
    twa_trig = -tWA*(max(0.,p%temp_soil_p(i)-Ttrig))**tWB
    !twa_trig = min(twa_trig,-1.e-6) ! avoiding the case: exp(0)=1.0
    ! exp(-40.)=4.24e.-18., avoid floating-point overflow
    twa_trig = max(twa_trig,-40.)
    !f_T(i) = 1 - exp(-tWA*(max(0.,p%temp_soil_p(i)-Ttrig))**tWB)
    f_T(i) = 1 - exp(twa_trig)
    ! 2024/03/01
    !f_T(i) = max(f_T(i),0.1)
    f_T(i) = max(f_T(i),0.01)
    !write(*,*)  "f_T(i) = ",f_T(i)
    ! plant resistance
    !write(*,*)  "p%f_root = ",p%f_root(i)
    !write(*,*)  "Lr = ",Lr
    ! rp unit: d
    !write(*,*)  "rp1 = ",r_xylem*Hp/p%f_root(i)
    !write(*,*)  "rp2 = ",r_r/(Lr*p%f_root(i))
    ! add leaf node in plant hydraulics
    rp(i) = (r_xylem*Hp/p%f_root(i) + r_r/(Lr*p%f_root(i)))*(1.0/f_T(i))*(1.0/f_theta(i))
    !write(*,*)  "rp(i) = ",rp(i)
    rp(i) = rp(i)*24.*3600.      ! convert from day to second
    !write(*,*)  "rp(i)*24*3600 = ",rp(i)
    !write(*,*)  "p%f_root(i) = ", p%f_root(i)
    !write(*,*)  "p%d_soil(i) = ", p%d_soil(i)
    ! soil-root resistance
    r_delta(i) = Lr*p%f_root(i)/p%d_soil(i)   ! unit : m/m3

    !write(*,*)  "r_delta(i) = ",r_delta(i)
    !write(*,*)  "deltal_min = ",deltal_min
    !write(*,*)  "deltal_max = ",deltal_max
    !write(*,*)  "-p_delta = ",-p_delta
    !write(*,*)  "p%f_root(i) = ", p%f_root(i)
    !write(*,*)  "-p_delta*r_delta(i) = ", -p_delta*r_delta(i)
    !write(*,*)  "exp(-p_delta*r_delta(i)) = ", exp(-p_delta*r_delta(i))
    !rdelta_tmp=min(r_delta(i),10./p_delta)
    !rp_delta = -p_delta*rdelta_tmp
    !write(*,*)  "p%KK(i)*p%f_root(i) = ", p%KK(i)*p%f_root(i)
    !write(*,*)  "rp_delta = ", rp_delta
    !f_deltal = deltal_min + (deltal_max - deltal_min)*exp(-p_delta*r_delta(i)) ! unit: m
    ! exp(-40.)=4.24e.-18., avoid floating-point overflow
    f_deltal = deltal_min + (deltal_max - deltal_min)*exp(max(-p_delta*r_delta(i),-40.)) ! unit: m
    !write(*,*)  "f_deltal = ", f_deltal
    !write(*,*)  "p%KK(i) = ", p%KK(i)
    !write(*,*)  "p%f_root(i) = ", p%f_root(i)
    !write(*,*)  "*************** "
    rs(i) = f_deltal/(p%KK(i)*p%f_root(i))  ! unit: s
    !write(*,*)  "p%KK*p%f_root = ",p%KK(i)*p%f_root(i)
    !write(*,*)  "f_deltal = ",f_deltal
    !write(*,*)  "p%f_root(i) = ",p%f_root(i)
    !write(*,*)  "rs(i) = ",rs(i)
    ! leaf water potential
    fei_leaf = (1.0-p%Sp/fpmax)*(fei_min + Hp) - Hp
    ! constrain the fei_min <= fei_leaf <= fei_th
    fei_leaf = max(fei_leaf,fei_th)
    !write(*,*)  "fei_leaf = ",fei_leaf
    !write(*,*)  "psim = ", p%psim(i)
    !write(*,*)  "z_depth", z_depth(p%n_layer-1)
    q1 = p%f_root(i)*(fei_leaf -p%psim_prev(i) - (Hp+z_depth(p%n_layer-1)) )/(rp(i) + rs(i))
    q1=max(0.,q1)
    !write(* ,*)  "q1 = ",q1
    q2 = p%f_root(i)*Etp/rho_w + p_excess*p%f_root(i)
    q2 = max(0.,q2)
    !write(*,*)  "q2 = ",q2
    q3 = (fpmax*p%f_root(i) - p%Sp*p%f_root(i))/kkstep
    !write(*,*)  "q3 = ",q3
    q3=max(0.,q3)
    qupt(i) = min(q1,q2)
    !write(*,*)  "qupt(i) = ",qupt(i)
 end do
 !write(*,*)  "fei_leaf = ",fei_leaf
 ! calculate vod based on lai and leaf water potential, based on Liu et al., 2021
 !"Global ecosystem-scale plant hydraulic traits retrieved using modelÃƒâ€šÃ‚Â¨Cdata fusion"

 !fei_leaf m H2o to MPa, 1 m H2o = 0.0098 MPa or 1/101 MPa
 ! The relation between VOD and fei_leaf/LAI should be further explored
 vod = (a + b*lai)*(1.0 + c*fei_leaf/101.)


! leaf temperature stress factor
 !T_canopy=(Tc_o_sunlit*LAI_o_sunlit+Tc_o_shaded*LAI_o_shaded+Tc_u_sunlit*LAI_u_sunlit+ &
     !Tc_u_shaded*LAI_u_shaded)/(LAI_o_sunlit+LAI_o_shaded+LAI_u_sunlit+LAI_u_shaded)
 !fTc_p = (Tleaf_H-Tleaf_O)/(Tleaf_O-Tleaf_L)
 !fTleaf=max(T_canopy-Tleaf_L,1.e-6)/(Tleaf_O-Tleaf_L)*(max(Tleaf_H-T_canopy,1.e-6)/(Tleaf_H-Tleaf_O))**fTc_p
! note: currently, I do not consider the leaf temperature stress factor @ Lu Hu

! The calculation form of water stress factor could be linear or non-linear
 f_feil = min((fei_leaf - fei_min)/(fei_th - fei_min),1.)
 f_feil = max(0.001,f_feil)
 !f_feil = max(1.e-6,f_feil)
 !f_feil = max(0.001,f_feil)

 ! using the logistic function to calculate the water stress factor @ Lu Hu
 !f_feil = 1.0/(1.0+(fei_leaf*(-1.0/101.0)/fei50)**plc_a)
 !f_feil = max(1.e-6,f_feil)
 !write(*,*)  "f_feil = ",f_feil

 ! 2024/03/13 gs = g0 + m* hs*An/Cs, wherever f_feileaf works on slope(m) or An,
 ! gs changes and the transpiration will be influenced!
 !Eta = Source/rho_w * f_feil
 if (fw_flag == 4) then
     Eta = Trans/rho_w * f_feil

 else
     Eta = Trans/rho_w
 end if



! Update the plant water storage
 do i=0,p%n_layer-1
    qupt_sum = qupt_sum + qupt(i)
 end do
 !write(*,*)  "qupt_sum = ",qupt_sum

 !deltal_Sp =  -(Eta - qupt_sum) * kstep*1.
 deltal_Sp =  -(Eta - qupt_sum) * kkstep

 !if (deltal_Sp>0) then
 !   deltal_Sp=min(10.**(-3),deltal_Sp)
 !else
 !   deltal_Sp=min(-10.**(-3),deltal_Sp)
 !end if
 p%Sp=p%Sp+deltal_Sp
 ! avoid the case: Sp<0
 p%Sp=max(1.e-6,p%Sp)
!p%Sp=max(2.0/1000.0,p%Sp)
 !write(*,*)  "deltal_Sp = ",deltal_Sp
 ! 2024/01/07
 p%f_feileaf = f_feil

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
      fpsisr(i) = 1./(1. + ((p%psim(i) - p%psi_min)/p%psi_min)**p%alpha)
    else
      fpsisr(i) = 1.
    end if

    if(p%temp_soil_p(i) > 0.) then
       ft(i) = (1.-exp(t1*p%temp_soil_p(i)**t2))
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
      p%dt(i) = dtt(i)/dtt_sum
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
   p%f_soilwater = min(1.0,fpsisr_sum)
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
