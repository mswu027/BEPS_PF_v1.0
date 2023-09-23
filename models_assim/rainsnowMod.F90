module rainsnowMod


use shr_kind_mod,only:r8=>shr_kind_r8
use beps_con,only:cp_ice,latent_fusion,density_water=>rho_w
use beps_par
implicit none

public:: rainfall_stage1,&
         rainfall_stage2,&
         snowpack_stage1,&
         snowpack_stage2,&
         snowpack_stage3

contains
!!rainfall_stage1 happens before evporation of intercepted water from canopy (supply)
!!rainfall_stage2 happens after evaporation of intercepted water from canopy (demand)
!*output:
!percentage of canopy covered by rainfall,overstorey and understorey (provided to evaporation_canopy)
!mass of water available for evaporation on canopy in this step
!precipitation on ground
!optical output:intercepted mass of rainfall in this step
!*input:
!air temperature
!precipitation (m/s)
!remain of water on leaves from last step(kg/m2) per leaf area
!lead area index of overstorey and understorey,excluding stem
!length of this step(s),if 10min,then it is set as 600
!air temperature and humidity
!*
subroutine rainfall_stage1(temp_air,precipitation,mass_water_o_last,mass_water_u_last,lai_o,lai_u,clumping,&
                           mass_water_o,mass_water_u,percent_water_o,percent_water_u,precipitation_g)
implicit none
real(r8),intent(in):: temp_air,precipitation !(oC,m/s)  precipitation=liquidwater @J.Wang
real(r8),intent(in):: mass_water_o_last,mass_water_u_last  !remins of water from last step
real(r8),intent(in):: lai_o,lai_u,clumping
real(r8),intent(out) :: mass_water_o,mass_water_u  !mass of water on leaves kg/m2 per ground area
real(r8),intent(out) :: percent_water_o,percent_water_u
real(r8),intent(out) :: precipitation_g

real(r8) :: precipitation_o,precipitation_u
real(r8) :: massMax_water_o,massMax_water_u !Maximum mass of water could be intercepted per leaf area(kg/m2)
real(r8) :: massStep_water_o,massStep_water_u ! mass of water intercepted in this step per leaf area
real(r8) :: length_step 

length_step      = kstep
precipitation_g  = 0.

!! overstorey
precipitation_o  = precipitation
mass_water_o     = mass_water_o_last + precipitation_o*length_step*density_water*(1-exp(-lai_o*clumping))
massMax_water_o  = 0.1*lai_o

mass_water_o     = max(0.,mass_water_o)
mass_water_o     = min(massMax_water_o,mass_water_o)

massStep_water_o = mass_water_o - mass_water_o_last
massStep_water_o = max(0.,massStep_water_o)

percent_water_o  = min(1.,mass_water_o/massMax_water_o)

!! understorey
precipitation_u  = precipitation_o - massStep_water_o/density_water/length_step
mass_water_u     = mass_water_u_last+precipitation_u*length_step*density_water*(1-exp(-lai_u*clumping))
massMax_water_u  = 0.1*lai_u

mass_water_u     = max(0.,mass_water_u)
mass_water_u     = min(mass_water_u,massMax_water_u)

massStep_water_u = mass_water_u - mass_water_u_last
massStep_water_u = max(massStep_water_u,0.)

percent_water_u  = min(1.,mass_water_u/massMax_water_u)

!! ground
precipitation_g  = precipitation_u - massStep_water_u/density_water/length_step

return
end subroutine rainfall_stage1


!!this module will calculate the water remained on canopy surface after evaporation in this step
subroutine rainfall_stage2(evapo_water_o,evapo_water_u,mass_water_o,mass_water_u)
implicit none
real(r8),intent(in) :: evapo_water_o,evapo_water_u
real(r8),intent(inout):: mass_water_o,mass_water_u

mass_water_o  = mass_water_o - evapo_water_o*kstep
mass_water_o  = max(0.,mass_water_o)

mass_water_u  = mass_water_u - evapo_water_u*kstep
mass_water_u  = max(0.,mass_water_u)

return
end subroutine rainfall_stage2

!*****************************************************
!     Snow Pack
!*****************************************************
!! this module will calculate the percentage of canopy and ground covered by snow
!! and output albedo of snow (used in enverge balance) and density of snow in this step
!! by XZ Luo, May 25,2015
!! snowpack_stage1 happens before any consumption of snow in this step, after the snow fall (supply)
!! snowpack_stage2 happens after sublimation from ground and canopy (demand)
!! snowpack_stage3 happens after frozen and melt of snow pack (demand)
!*input:
! air temperature,preciipitation, depth of snow from last step
! density of snow from last step,mass of snow on canopy and ground  from last step
! length of step, LAI_o/u and albedo of snow from last step
!*output
! mass of snow on canopy and ground accululation of snowfall
! albedo of snow in this step
! density of snow in this step

subroutine snowpack_stage1(temp_air,precipitation,mass_snow_o_last,mass_snow_u_last,mass_snow_g_last,density_snow_last,&
                           area_snow_o_last,area_snow_u_last,&
                           mass_snow_o,mass_snow_u,mass_snow_g,lai_o,lai_u,clumping,area_snow_o,area_snow_u,&
                           percent_snow_o,percent_snow_u,percent_snow_g,density_snow,depth_snow,albedo_v_snow,albedo_n_snow)
implicit none
real(r8),intent(in) :: temp_air,precipitation ! here precipitation = solid water @J.Wang
real(r8),intent(in) :: mass_snow_o_last,mass_snow_u_last,mass_snow_g_last,density_snow_last,area_snow_o_last,area_snow_u_last
real(r8),intent(in) :: lai_o,lai_u,clumping
real(r8),intent(out):: mass_snow_o,mass_snow_u,mass_snow_g ! mass of intercepted snow on canopy and gournd 
real(r8),intent(out):: area_snow_o,area_snow_u
real(r8),intent(out):: percent_snow_o,percent_snow_u,percent_snow_g !percentage of snow cover on canopy and ground
real(r8),intent(inout) :: density_snow
real(r8),intent(inout) :: depth_snow,albedo_v_snow,albedo_n_snow

real(r8) :: massMax_snow_o, massMax_snow_u
real(r8) :: massStep_snow_o, massStep_snow_u
real(r8) :: areaMax_snow_o, areaMax_snow_u !Maximum area of snow at overstorey and understorey
real(r8) :: change_depth_snow !change of snow depth on ground
real(r8) :: density_water, density_new_snow !density of water, density of newly fallen snow
real(r8) :: snowrate, snowrate_o, snowrate_u, snowrate_g
real(r8) :: albedo_v_Newsnow, albedo_n_Newsnow ! albedo of newly fallen snow in visible and near infrared band

density_new_snow=67.9+51.3*exp(temp_air/2.6)
albedo_v_Newsnow=0.94
albedo_n_Newsnow=0.8
massMax_snow_o=0.1*lai_o
massMax_snow_u=0.1*lai_u
areaMax_snow_o=lai_o*0.01
areaMax_snow_u=lai_u*0.01

snowrate  = precipitation  !!@J.Wang
!write(*,*) 'snowrate as precp', snowrate
!-- iLab::snowrate_o is always used in the condition at the end of this routine,
!         so it must definitively be initialised.
!         We are applying the outcommented initialisation
!         as used in 'temp_air<0' branch (just below):
snowrate_o   = max(snowrate,1.e-6)         !! @MOUSONG.WU
if(temp_air < 0) then
  !! overstorey
!  snowrate_o   = max(snowrate,1.e-6)         !! @MOUSONG.WU
!  write(*,*) 'snowrate_o: ', snowrate_o
  mass_snow_o  = mass_snow_o_last + snowrate_o*kstep*density_new_snow*(1-exp(-lai_o*clumping))
!  write(*,*) 'mass_snow_o', mass_snow_o, massMax_snow_o
  percent_snow_o  = mass_snow_o/massMax_snow_o
  percent_snow_o  = max(0.,percent_snow_o )
  percent_snow_o  = min(1.,percent_snow_o )

  area_snow_o     = percent_snow_o*areaMax_snow_o
  massStep_snow_o = mass_snow_o - mass_snow_o_last

  !! understorey
  snowrate_u   = snowrate - massStep_snow_o/density_new_snow/kstep
  snowrate_u   = max(0.,snowrate_u)

  mass_snow_u  = mass_snow_u_last+snowrate_u*kstep*density_new_snow*(1-exp(-lai_u*clumping))
  percent_snow_u  = mass_snow_u/massMax_snow_u
  percent_snow_u  = max(0.,percent_snow_u)
  percent_snow_u  = min(1.,percent_snow_u)

  area_snow_u     = percent_snow_u*areaMax_snow_u
  massStep_snow_u = mass_snow_u - mass_snow_u_last

  !! ground
  snowrate_g  = snowrate_u - massStep_snow_u/density_new_snow/kstep
  snowrate_g  = max(0.,snowrate_g)
  change_depth_snow = snowrate_g*kstep
else
  !! overstorey
  mass_snow_o  = mass_snow_o_last
  percent_snow_o  = mass_snow_o/massMax_snow_o
  percent_snow_o  = max(0., percent_snow_o)
  percent_snow_o  = min(1., percent_snow_o)
  area_snow_o     = area_snow_o_last

  !! understorey
  mass_snow_u     = mass_snow_u_last
  percent_snow_u  = mass_snow_u/massMax_snow_u
  percent_snow_u  = max(0.,percent_snow_u)
  percent_snow_u  = min(1.,percent_snow_u)
  area_snow_u     = area_snow_u_last
 
  !!ground 
  change_depth_snow=0.
end if

  change_depth_snow=max(0.,change_depth_snow)
  mass_snow_g=mass_snow_g_last+change_depth_snow*density_new_snow
  mass_snow_g=max(0.,mass_snow_g)
!  write(*,*) 'change_depth_snow', change_depth_snow
!  write(*,*) 'density_snow', density_snow_last
  
  if (change_depth_snow > 0.) then
   density_snow=(density_snow_last*depth_snow+density_new_snow*change_depth_snow)/(depth_snow+change_depth_snow)
  else
   density_snow=(density_snow_last-250.)*exp(-0.001*kstep/3600.)+250.0   !!@J.Wang ???
!    density_snow= 250.
  end if

  if(mass_snow_g >0.) then
     depth_snow  = mass_snow_g/density_snow
  else 
     depth_snow  = 0.
  end if

  percent_snow_g = mass_snow_g/(0.05*density_snow)
  percent_snow_g = min(percent_snow_g,1.)

!! albedo of snow in this step 
 if(snowrate_o>0.) then
     albedo_v_snow   = (albedo_v_snow-0.70)*exp(-0.005*kstep/3600.)+0.7  !!@J.Wang cannot understand clearly
     albedo_n_snow   = (albedo_n_snow-0.42)*exp(-0.005*kstep/3600.)+0.42  
  else
     albedo_v_snow   = albedo_v_Newsnow
     albedo_n_snow   = albedo_n_Newsnow
  end if

  return
end subroutine snowpack_stage1

!this module will calculate the snow remained on canopy surface after evaporation in this step
subroutine snowpack_stage2(evapo_snow_o,evapo_snow_u,mass_snow_o,mass_snow_u)
implicit none
real(r8),intent(in)  :: evapo_snow_o,evapo_snow_u
real(r8),intent(inout) :: mass_snow_o,mass_snow_u

mass_snow_o  = max(0.,mass_snow_o - evapo_snow_o*kstep)
mass_snow_u  = max(0.,mass_snow_u - evapo_snow_u*kstep)
return
end subroutine snowpack_stage2

!This module simulates the process of snow melting and water frozen in this step
!*input
! depth of snow on ground after stage 1,
! air temperature
! ground surface temperature
! // output:
! the amount of the melted snow, frozen snow
!
subroutine snowpack_stage3(temp_air,temp_snow,temp_snow_last,density_snow,depth_snow,depth_water,mass_snow_g)
implicit none
real(r8),intent(in)  :: temp_air,temp_snow,temp_snow_last,density_snow
real(r8),intent(inout) :: depth_snow,depth_water,mass_snow_g

real(r8) :: depth_snow_sup, mass_snow_sup ! depth and mass of snow after stage1, and minus the amount of sublimation
real(r8) :: mass_snow_melt, mass_water_frozen
real(r8) :: melt_depth_snow, frozen_depth_snow
real(r8) :: melt_depth_water, frozen_depth_water

!! assumed sublimation happens before the melting and freezing 
depth_snow_sup  = depth_snow
mass_snow_sup   = mass_snow_g

!!simulate snow melt and freeze process
mass_snow_melt=0.
mass_water_frozen=0.
if(depth_snow_sup <=0.02) then
    if(temp_air >0. .and. depth_snow_sup >0.) then
       mass_snow_melt=temp_air*0.0075*kstep/3600*0.3 
       mass_snow_melt=min(mass_snow_sup, mass_snow_melt)
    else
       mass_snow_melt=0.
    end if
else if(depth_snow_sup > 0.02 .and. depth_snow_sup <=0.05) then
    if (temp_snow>0. .and. temp_snow_last<0. .and. mass_snow_sup > 0.) then
      mass_snow_melt=temp_snow*cp_ice*density_snow*depth_snow_sup/latent_fusion
      mass_snow_melt=min(mass_snow_sup, mass_snow_melt)
    else
      mass_snow_melt=0.
    end if

    if(temp_snow<=0. .and. temp_snow_last >0. .and. depth_water>0.) then
       mass_water_frozen = (0.-temp_snow)*cp_ice*density_snow*depth_snow_sup/latent_fusion
       mass_water_frozen = max(mass_water_frozen,depth_water*density_water)
    else
       mass_water_frozen=0.
    end if
else if(depth_snow_sup >0.05) then
    if(temp_snow >0. .and. temp_snow_last <=0. .and. mass_snow_sup >0.) then
      mass_snow_melt=temp_snow*cp_ice*density_snow*depth_snow_sup*0.02/latent_fusion
      mass_snow_melt=min(mass_snow_sup, mass_snow_melt)
    else
      mass_snow_melt=0.
    end if

    if(temp_snow <=0. .and. temp_snow_last >0. .and. depth_water>0.) then
       mass_water_frozen = (0.-temp_snow)*cp_ice*density_snow*depth_snow_sup*0.02/latent_fusion
       mass_water_frozen = max(mass_water_frozen,depth_water*density_water)
    else
       mass_water_frozen=0.
    end if
end if

!!change in mass of snow on ground
mass_snow_g = mass_snow_g - mass_snow_melt+mass_water_frozen
mass_snow_g = max(0.,mass_snow_g)

!!change of depth in snow
melt_depth_snow    = mass_snow_melt/density_snow
frozen_depth_snow  = mass_water_frozen/density_snow
depth_snow         = depth_snow_sup -melt_depth_snow+frozen_depth_snow
if (isnan(depth_snow)) depth_snow = 0.
!write(*,*) depth_snow
depth_snow         = max(0.,depth_snow)

!!channge of depth in water
melt_depth_water   = mass_snow_melt/density_water
frozen_depth_water = mass_water_frozen/density_water
depth_water        = depth_water+melt_depth_water-frozen_depth_water
depth_water        = max(0.,depth_water)

return
end subroutine 


end module
