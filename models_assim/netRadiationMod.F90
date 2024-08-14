! This module calculate net radiation at both canopy level and leaf level XZ luo may23 2015
! output: net radiation for canopy,over/under storey and ground
!         ............. on sunlit/shaded leaves of over/understorey
!
! inputs: global solar radiation,cosine value for solar zenith angle,albedo of leaves
!         albedo of snow,percentage of snow cover
!         leaf area index 
!         temperature of over/under storey and ground
!         temperature of air/rh

subroutine netRadiation(shortRad_df,shortRad_dir,CosZs,temp_o,temp_u,temp_g,lai_o,lai_u,lai_os,lai_us, &
                        lai_o_sunlit,lai_o_shaded,lai_u_sunlit,lai_u_shaded,clumping,temp_air,rh,      &
                        albedo_snow_v,albedo_snow_n,percentArea_snow_o,percentArea_snow_u,percent_snow_g, &
                        albedo_v_o,albedo_n_o,albedo_v_u,albedo_n_u,albedo_v_g,albedo_n_g, &
                        netRad_o,netRad_u,netRad_g,netRadLeaf_o_sunlit,netRadLeaf_o_shaded,netRadLeaf_u_sunlit, &
                        netRadLeaf_u_shaded,netShortRadLeaf_o_sunlit,netShortRadLeaf_o_shaded,&
                        netShortRadLeaf_u_sunlit, netShortRadLeaf_u_shaded,fapar)
use meteoMod
use shr_kind_mod, only: r8=>shr_kind_r8
use beps_con,only: zero,sb_constant=>sigma
implicit none

real(r8),intent(in) :: shortRad_df,shortRad_dir,CosZs,temp_o,temp_u,temp_g
real(r8),intent(in) :: lai_o,lai_u,lai_os,lai_us ! LAI of over/understorey with or without stem
real(r8),intent(in) :: lai_o_sunlit,lai_o_shaded,lai_u_sunlit,lai_u_shaded ! sunlit/shaded leaves with consideration of stem
real(r8),intent(in) :: clumping
real(r8),intent(in) :: temp_air,rh
real(r8),intent(in) :: albedo_snow_v,albedo_snow_n ! albedo of snow in this step
real(r8),intent(in) :: percentArea_snow_o,percentArea_snow_u,percent_snow_g
real(r8),intent(in) :: albedo_v_o,albedo_n_o,albedo_v_u,albedo_n_u,albedo_v_g,albedo_n_g 
real(r8),intent(out):: netRad_o,netRad_u,netRad_g ! net Radiation on over/understorey and ground
real(r8),intent(out):: netRadLeaf_o_sunlit,netRadLeaf_o_shaded,netRadLeaf_u_sunlit,netRadLeaf_u_shaded !leaf levels for ET
real(r8),intent(out):: netShortRadLeaf_o_sunlit,netShortRadLeaf_o_shaded,netShortRadLeaf_u_sunlit, netShortRadLeaf_u_shaded 
                       !! net shortwave radiation at leaf level for GPP. 
real(r8) :: fapar ! fraction of absorbed photosynthetically active radiation
real(r8) :: shortRad_global,netShortRad_o,netShortRad_u,netShortRad_g ! net short wave radiation
real(r8) :: netShortRad_o_dir,netShortRad_o_df,netShortRad_u_dir,netShortRad_u_df,netShortRad_g_dir,netShortRad_g_df
!real(r8) :: shortRad_dir,shortRad_df
real(r8) :: netLongRadLeaf_o_sunlit,netLongRadLeaf_o_shaded,netLongRadLeaf_u_sunlit,netLongRadLeaf_u_shaded
real(r8) :: netLongRad_o,netLongRad_u,netLongRad_g
real(r8) :: shortRadLeaf_o_dir,shortRadLeaf_u_dir,shortRadLeaf_o_df,shortRadLeaf_u_df
real(r8) :: albedo_o,albedo_u,albedo_g  !albedo of overstorey/understorey/groudn(considering snow)
real(r8) :: albedo_v_os,albedo_n_os,albedo_v_us,albedo_n_us,albedo_v_gs,albedo_n_gs ! albedo of three parts in visible and NIR band (considering snow)
real(r8) :: emissivity_air,emissivity_o,emissivity_u,emissivity_g  !emissivity of air,over/understorey, and ground
real(r8) :: longRad_air,longRad_o,longRad_u,longRad_g   ! longwave radiation emissted by different parts
real(r8) :: cosQ_o,cosQ_u ! indicators to describe leaf distribution angles in canopy. slightly related with LAI
real(r8) :: gap_o_dir,gap_u_dir,gap_o_df,gap_u_df !gap fraction of direct and diffuse radiation for over/unerstory (diffuse used for diffuse solar radiation and longwave radiation
real(r8) :: gap_os_dir,gap_us_dir,gap_os_df,gap_us_df  ! considering stem
!real(r8) :: fapar ! fraction of absorbed photosynthetically active radiation

!calculate albedo of canopy in this step
albedo_v_os  = albedo_v_o*(1.-percentArea_snow_o)+albedo_snow_v*percentArea_snow_o
albedo_n_os  = albedo_n_o*(1.-percentArea_snow_o)+albedo_snow_n*percentArea_snow_o
albedo_v_us  = albedo_v_u*(1.-percentArea_snow_u)+albedo_snow_v*percentArea_snow_u
albedo_n_us  = albedo_n_u*(1.-percentArea_snow_u)+albedo_snow_n*percentArea_snow_u

albedo_o     = 0.5*(albedo_v_os+albedo_n_os)
albedo_u     = 0.5*(albedo_v_us+albedo_n_us)

! calculate albedo of ground 
albedo_v_gs  = albedo_v_g*(1.-percent_snow_g)+albedo_snow_v*percent_snow_g
albedo_n_gs  = albedo_n_g*(1.-percent_snow_g)+albedo_snow_n*percent_snow_g
albedo_g     = 0.5*(albedo_v_gs+albedo_n_gs) 

! separate global solar radiation into df and dir  @orgin
!      Here we input df/dir directly               @J.Wang
shortRad_global = shortRad_df + shortRad_dir     ! @J.Wang

! fraction at each layer of canopy,df/dir, use LAI here
if (CosZs > zero) then
   if(-0.5*clumping*lai_o/CosZs < -10.) then
      gap_o_dir = 0.
   else
      gap_o_dir    = exp(-0.5*clumping*lai_o/CosZs)
   end if
   if(-0.5*clumping*lai_u/CosZs < -10.) then
      gap_u_dir = 0.
    else
      gap_u_dir    = exp(-0.5*clumping*lai_u/CosZs)
    end if
    if(-0.5*clumping*lai_os/CosZs < -10.) then
       gap_os_dir = 0.
    else
       gap_os_dir   = exp(-0.5*clumping*lai_os/CosZs)   !considering stem
    end if
    if(-0.5*clumping*lai_us/CosZs < -10.) then
       gap_us_dir = 0.
     else
       gap_us_dir   = exp(-0.5*clumping*lai_us/CosZs)
     end if
end if

cosQ_o       = 0.537+0.025*lai_o    !leaf distribution angles
cosQ_u       = 0.537+0.025*lai_u
gap_o_df     = exp(-0.5*clumping*lai_o/cosQ_o)
gap_u_df     = exp(-0.5*clumping*lai_u/cosQ_u)
gap_os_df    = exp(-0.5*clumping*lai_os/cosQ_o)
gap_us_df    = exp(-0.5*clumping*lai_us/cosQ_u)

!fpar_os = (1. - albedo_o) - (1. - albedo_u) * gap_os_dir
!fpar_us = (1. - albedo_u) - (1. - albedo_g) * gap_us_dir

!emissivity of each part
call meteo_pack(temp_air,rh)
emissivity_air = 1.-exp(-(e_actual*10.0)**((temp_air+273.15)/1200.0))
emissivity_air = min(1.,emissivity_air)
emissivity_air = max(0.7,emissivity_air)

emissivity_o   = 0.98
emissivity_u   = 0.98
emissivity_g   = 0.96

!net short direct radiation on canopy and ground
if(shortRad_global>zero .and. CosZs >zero) then  ! day time
    netShortRad_o_dir  = shortRad_dir*((1.-albedo_o)-(1.-albedo_u)*gap_o_dir)
    netShortRad_u_dir  = shortRad_dir*gap_o_dir*((1.-albedo_u)-(1.-albedo_g)*gap_u_dir)
    netShortRad_g_dir  = shortRad_dir*gap_o_dir*gap_u_dir*(1.-albedo_g)
else
    netShortRad_o_dir  = 0.
    netShortRad_u_dir  = 0
    netShortRad_g_dir  = 0
end if

!net short diffuse radiation on canopy and ground
if(shortRad_global>zero .and. CosZs >zero) then  ! day time
   netShortRad_o_df    = shortRad_df*((1.-albedo_o)-(1.-albedo_u)*gap_o_df)+&
                         0.21*clumping*shortRad_dir*(1.1-0.1*lai_o)*exp(-1.*CosZs)
   netShortRad_u_df    = shortRad_df*gap_o_df*((1.-albedo_u)-(1.-albedo_g)*gap_u_df)+ &
                         0.21*clumping*shortRad_dir*gap_o_dir*(1.1-0.1*lai_u)*exp(-1.*CosZs)
   netShortRad_g_df    = shortRad_df*gap_o_df*gap_u_df*(1.-albedo_g)
else
   netShortRad_o_df    = 0.
   netShortRad_u_df    = 0.
   netShortRad_g_df    = 0.
end if

!total net shortwave radiation at canopy level
netShortRad_o   = netShortRad_o_dir + netShortRad_o_df
netShortRad_u   = netShortRad_u_dir + netShortRad_u_df
netShortRad_g   = netShortRad_g_dir + netShortRad_g_df

!net longwave radiation on canopy and ground
longRad_air     = emissivity_air*sb_constant*(temp_air+273.15)**4
longRad_o       = emissivity_o*sb_constant*(temp_o+273.15)**4
longRad_u       = emissivity_u*sb_constant*(temp_u+273.15)**4  
longRad_g       = emissivity_g*sb_constant*(temp_g+273.15)**4

netLongRad_o    = (emissivity_o*(longRad_air+longRad_u*(1.-gap_u_df)+longRad_g*gap_u_df)-2.*longRad_o)*(1.-gap_o_df)+&
                  emissivity_o*(1.-emissivity_u)*(1.-gap_u_df)*(longRad_air*gap_o_df+longRad_o*(1.-gap_o_df))
netLongRad_u    = (emissivity_u*(longRad_air*gap_o_df+longRad_o*(1.-gap_o_df)+longRad_g)-2.*longRad_u)*(1.-gap_u_df)+&
                  (1.-emissivity_g)*((longRad_air*gap_o_df+longRad_o*(1.-gap_o_df))*gap_u_df+longRad_u*(1.-gap_u_df))+&
                  emissivity_u*(1.-emissivity_o)*(longRad_u*(1.-gap_u_df)+longRad_g*gap_u_df)*(1.-gap_o_df)
netLongRad_g    = emissivity_g*((longRad_air*gap_o_df+longRad_o*(1.-gap_o_df))*gap_u_df+longRad_u*(1.-gap_u_df)) - &
                  longRad_g+(1.-emissivity_u)*longRad_g*(1.-gap_u_df)

!total net radiation for overstorey/understorey/ground
netRad_o  = netShortRad_o+netLongRad_o
netRad_u  = netShortRad_u+netLongRad_u
netRad_g  = netShortRad_g+netLongRad_g

!leaf level net radiation updated way
! reference Chen2012 clumping index paper
if(shortRad_global>zero .and. CosZs > zero) then
   shortRadLeaf_o_dir = 0.5*shortRad_dir/CosZs
   shortRadLeaf_o_dir = min(shortRadLeaf_o_dir,0.7*1362.)
   shortRadLeaf_u_dir = shortRadLeaf_o_dir

   shortRadLeaf_o_df  = (shortRad_df-shortRad_df*gap_os_df)/lai_os+0.07*shortRad_dir*(1.1-0.1*lai_os)*exp(-CosZs)
   shortRadLeaf_u_df  = (shortRad_df*gap_o_df-shortRad_df*gap_o_df*gap_us_df)/lai_us+&
                        0.05*shortRad_dir*gap_o_dir*(1.1-0.1*lai_us)*exp(-CosZs)
else
   shortRadLeaf_o_dir = 0.
   shortRadLeaf_u_dir = 0.
   shortRadLeaf_o_df  = 0.
   shortRadLeaf_u_df  = 0.
end if

!overstorey sunlit leaves
if(lai_o_sunlit>0.) then
   netShortRadLeaf_o_sunlit = (shortRadLeaf_o_dir+shortRadLeaf_o_df)*(1.-albedo_o)
   netLongRadLeaf_o_sunlit  = netLongRad_o/lai_os   !leaf level net long
   netRadLeaf_o_sunlit      = netShortRadLeaf_o_sunlit + netLongRadLeaf_o_sunlit
else
   netShortRadLeaf_o_sunlit = (shortRadLeaf_o_dir+shortRadLeaf_o_df)*(1.-albedo_o)
   netLongRadLeaf_o_sunlit  = netLongRad_o
   netRadLeaf_o_sunlit      = netShortRadLeaf_o_sunlit+netLongRadLeaf_o_sunlit
end if

!overstorey shaded leaves
if(lai_o_shaded>0.) then
   netShortRadLeaf_o_shaded = shortRadLeaf_o_df*(1.-albedo_o)
   netLongRadLeaf_o_shaded  = netLongRad_o/lai_os
   netRadLeaf_o_shaded      = netShortRadLeaf_o_shaded + netLongRadLeaf_o_shaded
else
   netShortRadLeaf_o_shaded = shortRadLeaf_o_df*(1.-albedo_o)
   netLongRadLeaf_o_shaded  = netLongRad_o
   netRadLeaf_o_shaded      = netShortRadLeaf_o_shaded + netLongRadLeaf_o_shaded
end if

!understorey sunlit leaf
if(lai_u_sunlit>0.) then
  netShortRadLeaf_u_sunlit  = (shortRadLeaf_u_dir+shortRadLeaf_u_df)*(1.-albedo_u)
  netLongRadLeaf_u_sunlit   = netLongRad_u/lai_us
  netRadLeaf_u_sunlit       = netShortRadLeaf_u_sunlit + netLongRadLeaf_u_sunlit
else
  netShortRadLeaf_u_sunlit  = (shortRadLeaf_u_dir+shortRadLeaf_u_df)*(1.-albedo_u)
  netLongRadLeaf_u_sunlit   = netLongRad_u
  netRadLeaf_u_sunlit       = netShortRadLeaf_u_sunlit + netLongRadLeaf_u_sunlit
end if

!understorey shaded leaf
if(lai_u_shaded >0.) then
  netShortRadLeaf_u_shaded  = shortRadLeaf_u_df*(1.-albedo_u)
  netLongRadLeaf_u_shaded   = netLongRad_u/lai_us
  netRadLeaf_u_shaded       = netShortRadLeaf_u_shaded + netLongRadLeaf_u_shaded
else
  netShortRadLeaf_u_shaded  = shortRadLeaf_u_df*(1.-albedo_u)
  netLongRadLeaf_u_shaded   = netLongRad_u
  netRadLeaf_u_shaded       = netShortRadLeaf_u_shaded + netLongRadLeaf_u_shaded
end if

! calculate fapar according to Chen 1996, adopted by MousongWU 2024
fapar = 1. - albedo_o - (1. - albedo_g) * gap_o_dir

 end subroutine
