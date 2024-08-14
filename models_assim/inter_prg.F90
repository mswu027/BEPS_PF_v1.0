!! Last update: July,2015
!! Fortran version: 3/5/2017 @J.Wang

subroutine inter_prg(yr, mn, dy, tod, &
     lai,lai_input,lc,clumping,Vcmax0,VJ_slope,N_leaf,bh2o,mh2o,f_leaf,p_kc25,p_ko25,p_tau25,sif_alpha,sif_beta,&
     param,meteo,CosZs,var_o,var_n,soilp,mid_res,daylen)
use shr_kind_mod,only:r8=>shr_kind_r8
!--iLab::restrict use of beps_time_manager to required entities  
!--iLab-update::extended arguments to avoid time manager completely
! use beps_time_manager, only:is_end_curr_day, get_curr_date
use beps_soilMod
use bepstype
use mid_results
use beps_par
use AnGsMod
use rainsnowMod
implicit none

!--iLab::added date-elements as argument to avoid 'call get_curr_date' further below
integer, intent(in) :: yr,mn,dy,tod
integer,intent(in)   ::lc
real(r8),intent(in)  ::clumping,Vcmax0,VJ_slope,N_leaf,bh2o,mh2o,f_leaf,p_ko25,p_kc25,p_tau25,sif_alpha,sif_beta,CosZs,daylen   !!J.Wang
real(r8),intent(in)  ::param(0:49)
type(climatedata),intent(in) ::meteo
real(r8),intent(in)  ::var_o(0:40)
real(r8),intent(out) ::var_n(0:40)
real(r8),intent(inout)  :: lai
type(soil)           ::soilp
type(results)        ::mid_res
integer,intent(in)   :: lai_input
!real(r8),intent(out) :: gs_h2o,G_o_b
integer              :: num,kkk,i,j
integer,parameter    :: iter_max = 20
integer              :: landcover
real(r8)             :: lai_o,lai_u,stem_o,stem_u
real(r8)             :: d_soil(0:layer)
real(r8)             :: Zsp                       !the depth of snow on the surface
real(r8)             :: Zp,Zp1 = 0.,Zp2 = 0.      !depth of pounded water on the surface
real(r8)             :: height_wind_sp            !height of the Va measured for calculation of L
real(r8)             :: Qhc_o(0:iter_max),Qhc_u(0:iter_max),Qhg(0:iter_max)  !sensible heating
real(r8)             :: G(0:layer+1,0:iter_max)  !the heat flux into the canopy of over story --in W/m^2
real(r8)             :: Wcl_o(0:iter_max),Wcs_o(0:iter_max)  !the masses od rain and snow on the canopy
real(r8)             :: Xcl_o(0:iter_max),Xcs_o(0:iter_max)  !the fractoion of canopy coverd by liquid water and snow
real(r8)             :: Wcl_u(0:iter_max),Wcs_u(0:iter_max)
real(r8)             :: Xcl_u(0:iter_max),Xcs_u(0:iter_max)
real(r8)             :: r_rain_g(0:iter_max)   !the rainfall rate on ground surface m/s
real(r8)             :: rho_snow(0:iter_max)   !density of snow
real(r8)             :: alpha_v_sw(0:iter_max),alpha_n_sw(0:iter_max) !albedo of snow
real(r8)             :: Wg_snow(0:iter_max)    ! amount of snow on ground
real(r8)             :: Xg_snow(0:iter_max)    ! fraction of the ground surface overed by snow
real(r8)             :: Ac_snow_o(0:iter_max),Ac_snow_u(0:iter_max)    ! area of canopy covered by snow
real(r8)             :: Ts0(0:iter_max),Tsn0(0:iter_max),Tsm0(0:iter_max),Tsn1(0:iter_max),Tsn2(0:iter_max)  !surface temperature
real(r8)             :: Tc_u(0:iter_max)    !effective canopy temperature in K
real(r8)             :: Tm(0:layer+1,0:iter_max) !soil temperature at the bottom and middle of each layer
real(r8)             :: lambda_soil(0:layer+1)   !thermal conductivity of each soil layer
real(r8)             :: Cs(0:layer+1,0:iter_max) ! soil volumetric hear capacity of each soil layer J/m3/K
real(r8)             :: temp_air,temp_grd
real(r8)             :: rh_air,wind_sp,snow,rainfall   !%,m/s,m/s @J.Wang  make original precipitation into rainfall and snow
real(r8)             :: Eil_o(0:iter_max),EiS_o(0:iter_max) !the evaporation rate of intercepted water of overstory--in kg/m^2/s
real(r8)             :: Eil_u(0:iter_max),EiS_u(0:iter_max)
real(r8)             :: Trans_o(0:iter_max),Trans_u(0:iter_max)  !transpiration
real(r8)             :: Evap_soil(0:iter_max)
real(r8)             :: Evap_SW(0:iter_max)     !evaporation from water pond
real(r8)             :: Evap_SS(0:iter_max)     !evaporation from snow pack
real(r8)             :: lambda_snow(0:iter_max) !effective thermal conductivity of snow in m2/s
real(r8)             :: e_a10                   ! vapour partial pressure of water in kPa
real(r8)             :: Lv_liquid               ! the latent heat of vaporation from liquid at air temperature=Ta
real(r8)             :: Lv_solid = 2.83*1e6     ! the latent heat of vaporation from solid (snow/ice) at air temperature
real(r8)             :: Ks                      ! instantaneous total short wave radiation (Global radiation)
real(r8)             :: alpha_sat,alpha_dry
real(r8)             :: alpha_v_o,alpha_v_u    !visible albedo of overstory,  o--overstory, u--understory
real(r8)             :: alpha_n_o,alpha_n_u    !near_infrared albedo
real(r8)             :: alpha_g                ! all-wave ground surface albedo
real(r8)             :: alpha_v_g,alpha_n_g
real(r8)             :: Cp_ca                  ! specific heat of moist air above the canopy
real(r8)             :: ra_o,ra_u,ra_g         ! the aerodynamic resistance of overstory, understory and ground surface
real(r8)             :: q_ca                   ! actural canopy stomatal resistance  --in s/m
real(r8)             :: radiation_o, radiation_u, radiation_g
real(r8)             :: ip = 0.                !the cumulative infiltration at the time of ponding   --in m/s
real(r8)             :: Infil = 0.
real(r8)             :: zr  = 0.8
real(r8)             :: cpd = 1004.65
real(r8)             :: Cs_o_sunlit_old,Cs_o_shaded_old,Cs_u_sunlit_old,Cs_u_shaded_old
                        ! CO2 concentration on the surfaces of leaves
real(r8)             :: COSs_o_sunlit_old,COSs_o_shaded_old,COSs_u_sunlit_old,COSs_u_shaded_old
                        ! COS concentration on the surfaces of leaves
real(r8)             :: Tc_o_sunlit_old,Tc_o_shaded_old,Tc_u_sunlit_old,Tc_u_shaded_old
                        ! the effective canopy temperature in K
real(r8)             :: Gs_o_sunlit_new,Gs_o_shaded_new,Gs_u_sunlit_new,Gs_u_shaded_new
                        !stomatal conductance of the big leaf     for water
real(r8)             :: Gs_o_sunlit_old,Gs_o_shaded_old,Gs_u_sunlit_old,Gs_u_shaded_old
real(r8)             :: Ac_o_sunlit,Ac_o_shaded,Ac_u_sunlit,Ac_u_shaded  ! net photosynthesis rate
real(r8)             :: Cs_o_sunlit_new,Cs_o_shaded_new,Cs_u_sunlit_new,Cs_u_shaded_new
                        ! CO2 concentration on the surfaces of     leaves
real(r8)             :: COSs_o_sunlit_new,COSs_o_shaded_new,COSs_u_sunlit_new,COSs_u_shaded_new
                        ! COS concentration on the surfaces of     leaves
real(r8)             :: Ci_o_sunlit_new,Ci_o_shaded_new,Ci_u_sunlit_new,Ci_u_shaded_new
                        ! intercellular CO2 concentration pn th    e leaf
real(r8)             :: COSi_o_sunlit_new,COSi_o_shaded_new,COSi_u_sunlit_new,COSi_u_shaded_new
                        ! intercellular COS concentration pn the leaf
real(r8)             :: Ci_o_sunlit_old,Ci_o_shaded_old,Ci_u_sunlit_old,Ci_u_shaded_old
real(r8)             :: COSi_o_sunlit_old,COSi_o_shaded_old,COSi_u_sunlit_old,COSi_u_shaded_old
real(r8)             :: Cc_o_sunlit_new,Cc_o_shaded_new,Cc_u_sunlit_new,Cc_u_shaded_new
                        ! CO2 concentration in the chloroplast
real(r8)             :: COSc_o_sunlit_new,COSc_o_shaded_new,COSc_u_sunlit_new,COSc_u_shaded_new
                        ! COS concentration in the chloroplast
real(r8)             :: Tc_o_sunlit_new,Tc_o_shaded_new,Tc_u_sunlit_new,Tc_u_shaded_new
                        ! the effective canopy temperature in K
real(r8)             :: f_soilwater    ! an emperical parameter describin    g the relative availability of soil water for plants
real(r8)             :: Gw_o_sunlit,Gw_o_shaded,Gw_u_sunlit,Gw_u_shaded ! the total conductance for water from     the intercellular space of the leaves to the reference height above the canopy
real(r8)             :: Gc_o_sunlit,Gc_o_shaded,Gc_u_sunlit,Gc_u_shaded  ! the total conductance for CO2 from th    e intercellular space of the leaves to the reference height above the canopy
real(r8)             :: Gww_o_sunlit,Gww_o_shaded,Gww_u_sunlit,Gww_u_shaded ! the total conductance for water from     the surface of the leaves to the reference height above the canopy
real(r8)             :: Gh_o_sunlit,Gh_o_shaded,Gh_u_sunlit,Gh_u_shaded  !total conductance for heat transfer f    rom the leaf surface to the reference height above the canopy
real(r8)             :: psychrometer
real(r8)             :: R_o_sunlit,R_o_shaded,R_u_sunlit,R_u_shaded,fapar   !solar radiation absorbed by sunlit, s    haded leaves
real(r8)             :: Tco, Tcu,slope
real(r8)             :: H_o_sunlit,H_o_shaded         !sensible heat flux from leaves
real(r8)             :: LAI_o_sunlit,LAI_o_shaded,LAI_u_sunlit,LAI_u_shaded
real(r8)             :: LAIo_sunlit,LAIo_shaded,LAIu_sunlit,LAIu_shaded
real(r8)             :: radiation_o_sun, radiation_o_shaded,radiation_u_sun, radiation_u_shaded  ! net radiation of leaves
real(r8)             :: GPP_o_sunlit,GPP_o_shaded,GPP_u_sunlit,GPP_u_shaded
real(r8)             :: SIF_o_sunlit,SIF_o_shaded,SIF_u_sunlit,SIF_u_shaded  !canopy level
real(r8)             :: lSIF_o_sunlit,lSIF_o_shaded,lSIF_u_sunlit,lSIF_u_shaded  !leaf level
real(r8)             :: stSIF_o_sunlit,stSIF_o_shaded,stSIF_u_sunlit,stSIF_u_shaded !scat
real(r8)             :: COS_o_sunlit,COS_o_shaded,COS_u_sunlit,COS_u_shaded
real(r8)             :: lCOS_o_sunlit,lCOS_o_shaded,lCOS_u_sunlit,lCOS_u_shaded
real(r8)             :: VPS_air
real(r8)             :: gs_h2o
real(r8)             :: GH_o,G_o_a,G_o_b,G_u_a, G_u_b
real(r8)             :: canopyh_o,canopyh_u
real(r8)             :: VPD_air
real(r8)             :: mass_water_g
real(r8)             :: percentArea_snow_o, percentArea_snow_u
real(r8)             :: Gheat_g
real(r8)             :: b_h2o            !the intercept term in BWB model (mol H2O m-2 s-1)
real(r8)             :: m_h2o            ! the slope in BWB model
real(r8)             :: leleaf_o_sunlit,leleaf_o_shaded,leleaf_u_sunlit,leleaf_u_shaded !leaf latent heat flux (mol/m2/s)

!for the Vcmax-Nitrogen calculation
real(r8)             :: Kn = 0.3
real(r8)             :: G_theta = 0.5
!real(r8)             :: K,Vcmax0,Vcmax_sunlit,Vcmax_shaded,expr1,expr2,expr3
real(r8)             :: K,Vcmax_sunlit,Vcmax_shaded,expr1,expr2,expr3    !Vcmax0 as an input from outside
real(r8)             :: slope_Vcmax_N, leaf_N,Jmax_sunlit,Jmax_shaded
real(r8)             :: ffpa     ! for SIF simulation @JWang
real(r8)             :: temp_day ! for storing daily mean temperature
real(r8)             :: theta_day  ! for storing daily mean surface soil moisture
real(r8)             :: trans_day  ! for storing daily mean transpiration
real(r8)             :: cosa
real(r8)             :: cos_soil
!--iLab::introduced to avoid calling is_end_curr_day() from beps_time_manager
logical :: is_end_curr_day

is_end_curr_day = (tod == 0) !--iLab: taken from BEPS time manager
temp_day = 0.
theta_day = 0.
trans_day = 0.

psychrometer=0.066
alpha_sat   = param(24)
alpha_dry   = param(25)
canopyh_o   = param(29)  ! to be used for module aerodynamic_conductance
canopyh_u   = param(30)
height_wind_sp   = param(31)
!height_wind_sp   = 30.
!m_h2o       = param(33)    ! used for photosynthesis
!b_h2o       = param(34)
m_h2o       = mh2o
b_h2o       = bh2o

!-- iLab::g2_h2o is *only* set from other routines in case 'CosZs>0.',
!         so we *must* initialise it and have uncommented the initialiser
!         that was already present.
gs_h2o = 0.
!gs_h2o      = 0.
! Vcmax-Nitrogen calculations by G,Mo 2011.04
if(CosZs >0.) then
   K = G_theta*clumping/CosZs
!   Vcmax0  = param(36)    !an input from outside @J.Wang
   if(K >10.) then          !! adjust K range here to get rid of floating-point exceptions,@MOUSONG.WU
     expr1 = 1.
     expr2 = 1.
   else
     expr1   = 1.-exp(-K*lai)
     expr2   = 1.-exp(-lai*(Kn+K))
   end if
   expr3   = 1.-exp(-Kn*lai)

   if(expr1 >0.) then
      Vcmax_sunlit = Vcmax0*param(47)*param(46)*K*expr2/(Kn+K)/expr1
   else
      Vcmax_sunlit = Vcmax0
   end if

   if(K>0 .and. lai>expr1/K) then
      Vcmax_shaded = Vcmax0*param(47)*param(46)*(expr3/Kn-expr2/(Kn+K))/(lai-expr1/K)
   else
      Vcmax_shaded = Vcmax0
   end if
end if

!! LAI calculation module by B.Chen
lai_o  = lai
if(lai<0.1) lai_o = 0.1
landcover    = int(param(4))

! Calculate ffpa as a function of lai, this makes the ffpa vary with time, instead of being a constant as below,@MOUSONG WU, 2020-09-14
ffpa = 1. - exp(-0.45*lai)
ffpa = max(1.e-2,ffpa)
ffpa = min(1.0,ffpa)

! added for sif simulation @J. Wang
!select case (landcover)
!   case (1)    !conifer evergreen
!    ffpa = 0.6
!   case(2)      !conifer decidous
!    ffpa = 0.6
!   case(6)      !broadleaf decidous
!    ffpa = 0.8
!   case(9)      !broadleaf evergreen
!    ffpa = 0.8
!   case(10)     !mix
!    ffpa = 0.4
!   case(13)     !shrub
!    ffpa = 0.8
!   case(14)     ! grass
!    ffpa = 0.8
!   case(15)     ! crop
!    ffpa = 0.6
!   case(40)     ! c4 grass
!    ffpa = 0.6  
!   case(41)     ! C4 crop
!    ffpa = 0.6
!end select


!if(landcover == 25 .or. landcover ==40) then
if(landcover == 14 .or. landcover == 15 .or. landcover == 40 .or. landcover == 41) then   !14->grass 15->crop @JWang
   lai_u = 0.01
else
   lai_u = 1.18*exp(-0.99*lai_o)
end if

if(lai_u>lai_o) lai_u = 0.01

stem_o   = param(8)*0.2
stem_u   = param(9)*0.2

call lai_cal(stem_o,stem_u,landcover,CosZs,lai_o,clumping,lai_u,LAIo_sunlit,LAIo_shaded,LAIu_sunlit,LAIu_shaded,&
          LAI_o_sunlit,LAI_o_shaded,LAI_u_sunlit,LAI_u_shaded)    !Bing Chen


!-------initalization of this time step
Ks        = meteo%Srad
rh_air    = meteo%rh
wind_sp   = meteo%wind
rainfall  = meteo%rainfall  !m/s   liquid water
snow      = meteo%snow      !m/s   snow
temp_air  = meteo%temp

if(Ks <=0) then
   alpha_v_o   = 0.
   alpha_n_o   = 0.
   alpha_v_u   = 0.
   alpha_n_u   = 0.
else
   alpha_v_o   = param(22)
   alpha_n_o   = param(23)
   alpha_v_u   = param(22)
   alpha_n_u   = param(23)
end if

Qhc_o(0)       = var_o(1)

Ts0(0)         = var_o(3)
if((Ts0(0) - temp_air) >2.0)    Ts0(0) = temp_air+2.0
if((Ts0(0) - temp_air) <-2.0)   Ts0(0) = temp_air-2.0

Tsn0(0)        = var_o(4)
if((Tsn0(0) - temp_air)>2.0)    Tsn0(0)= temp_air+2.0
if((Tsn0(0) - temp_air)<-2.0)   Tsn0(0)= temp_air-2.0

Tsm0(0)         = var_o(5)
if((Tsm0(0)-temp_air)>2.)       Tsm0(0)= temp_air+2.0
if((Tsm0(0)-temp_air)<-2.)      Tsm0(0)= temp_air-2.0

Tsn1(0)         = var_o(6)
if((Tsn1(0)-temp_air)>2.0)      Tsn1(0)= temp_air+2.0
if((Tsn1(0)-temp_air)<-2.)      Tsn1(0)= temp_air-2.0

Tsn2(0)         = var_o(7)
if((Tsn2(0)-temp_air)>2.0)      Tsn2(0)= temp_air+2.0
if((Tsn2(0)-temp_air)<-2.)      Tsn2(0)= temp_air-2.0

Wcl_o(0)        = var_o(15)  !the mass of intercepted liquid water and snow, overstory
Wcs_o(0)        = var_o(16)

Wcl_u(0)        = var_o(18)
Wcs_u(0)        = var_o(19)
Wg_snow(0)      = var_o(20) !  fraction of ground surface covered by snow and snow mass

soilp%Zsp       = var_o(33)
soilp%Zp        = var_o(34)
soilp%r_rain_g  = var_o(35)

!-- iLab::*must* at least initialise arrays for complete number of iterations
!          (and not only the first element),
!         since elements are input/output(!) arguments to called routines
!         (snowpack_stage1,netRadiation)
Ac_snow_o(0:iter_max)    = var_o(36)
Ac_snow_u(0:iter_max)    = var_o(37)
rho_snow(0:iter_max)     = var_o(38)
alpha_v_sw(0:iter_max)   = var_o(39)
alpha_n_sw(0:iter_max)   = var_o(40)
! Ac_snow_o(0)    = var_o(36)
! Ac_snow_u(0)    = var_o(37)
! rho_snow(0)     = var_o(38)
! alpha_v_sw(0)   = var_o(39)
! alpha_n_sw(0)   = var_o(40)
Zsp             = soilp%Zsp
Zp              = soilp%Zp

if(Zp <0.001) Zp  = 0.
!if(Zp < 1.e-6) Zp = 0.
do i = 9,14
    soilp%temp_soil_p(i-9)   = var_o(i)
end do
do i = 21,26
    soilp%thetam_prev(i-21)  = var_o(i)
end do
do i = 27,32
    soilp%ice_ratio(i-27)    = var_o(i)
end do

! vcmax jmax module  by L. He
slope_Vcmax_N      = param(47)
!leaf_N             = param(46)
leaf_N             = N_leaf

call Vcmax_Jmax(lai_o,clumping,Vcmax0,VJ_slope,slope_Vcmax_N,leaf_N,CosZs,Vcmax_sunlit,Vcmax_shaded,Jmax_sunlit,Jmax_shaded)
! temperatures of overstorey and understorey canopies
Tc_o_sunlit_old=temp_air-0.5
Tc_o_shaded_old=temp_air-0.5
Tc_u_sunlit_old=temp_air-0.5
Tc_u_shaded_old=temp_air-0.5

do kkk = 1,kloop           !sub-time iteration @J.Wang
    ! Snow pack stage 1  by R.Luo
    call snowpack_stage1(temp_air,snow,Wcs_o(kkk-1),Wcs_u(kkk-1),Wg_snow(kkk-1),rho_snow(kkk-1),Ac_snow_o(kkk-1),&
                         Ac_snow_u(kkk-1),Wcs_o(kkk),Wcs_u(kkk),Wg_snow(kkk),&
                         lai_o,lai_u,clumping,Ac_snow_o(kkk),Ac_snow_u(kkk),Xcs_o(kkk),Xcs_u(kkk),Xg_snow(kkk),&
                         rho_snow(kkk),Zsp,alpha_v_sw(kkk),alpha_n_sw(kkk))
!    write(*,*) "DG01: Ac_snow_o(kkk) =",Ac_snow_o(kkk)

    ! rainfall stag 1
    call rainfall_stage1(temp_air,rainfall,Wcl_o(kkk-1),Wcl_u(kkk-1),lai_o,lai_u,clumping,Wcl_o(kkk),Wcl_u(kkk),&
                         Xcl_o(kkk),Xcl_u(kkk),r_rain_g(kkk))

    if(soilp%thetam_prev(1)<soilp%theta_vwp(1)*0.5) then
         alpha_g  = alpha_dry
    else
         alpha_g  = (soilp%thetam_prev(1)-soilp%theta_vwp(1)*0.5)/(soilp%fei(1)-soilp%theta_vwp(1)*0.5)*&
                    (alpha_sat-alpha_dry)+alpha_dry
    end if

    alpha_v_g     = 2./3.*alpha_g
    alpha_n_g     = 4./3.*alpha_g

    ! soil water factor module
    call soil_water_factor_v2(soilp)

    f_soilwater  = soilp%f_soilwater

    if(f_soilwater >1.0) f_soilwater = 1.0
    GH_o    = Qhc_o(kkk-1)        !used as the init. for module aerodynamic_conductance
    VPS_air = 0.61078*exp(17.3*temp_air/(237.3+temp_air))
    e_a10   = VPS_air*rh_air/100.
    VPD_air = VPS_air - e_a10   !water vapor deficit at the reference height @J.Wang maybe directly use meteo_pack

    q_ca    = 0.622 * e_a10 / (101.35- 0.378 * e_a10)           !g/g  no dimention
    Cp_ca   = cpd*(1.+0.84*q_ca)

    slope   = 2503.0 / (temp_air+ 237.3)**2 * exp(17.27 *temp_air/ (temp_air+ 237.3))

    Gs_o_sunlit_old=1./200.0
    Ci_o_sunlit_old=0.7*CO2_air
    Gs_o_shaded_old=1./200.0
    Ci_o_shaded_old=0.7*CO2_air
    Gs_u_sunlit_old=1./300.0
    Ci_u_sunlit_old=0.7*CO2_air
    Gs_u_shaded_old=1./300.0
    Ci_u_shaded_old=0.7*CO2_air
    COSi_o_sunlit_old = 0.7*COS_air
    COSi_o_shaded_old = 0.7*COS_air
    COSi_u_sunlit_old = 0.7*COS_air
    COSi_u_shaded_old = 0.7*COS_air
    percentArea_snow_o=Ac_snow_o(kkk)/lai_o/2.
    percentArea_snow_u=Ac_snow_u(kkk)/lai_u/2.

    temp_grd  = temp_air     !ground temperature substituted by air temperature

    num=0
    do while (.True.)
       num = num+1
       ! aerodynamic_conductance module by G.Mo
       call aerodynamic_conductance(canopyh_o,canopyh_u,height_wind_sp,clumping,temp_air,wind_sp,GH_o,lai_o+stem_o,&
                                    lai_u+stem_u,ra_o,ra_u,ra_g,G_o_a,G_o_b,G_u_a,G_u_b)
       Gh_o_sunlit = 1.0/(1.0/G_o_a+0.5/G_o_b)  !heat conductance of sunlit leaves of overstorey
       Gh_o_shaded = 1.0/(1.0/G_o_a+0.5/G_o_b)
       Gh_u_sunlit = 1.0/(1.0/G_u_a+0.5/G_u_b)
       Gh_u_shaded = 1.0/(1.0/G_u_a+0.5/G_u_b)

       Gww_o_sunlit=1.0/(1.0/G_o_a+1.0/G_o_b+100.)  ! conductance for intercepted water of sunlit leaves of overstorey
       Gww_o_shaded=1.0/(1.0/G_o_a+1.0/G_o_b+100.)
       Gww_u_sunlit=1.0/(1.0/G_u_a+1.0/G_u_b+100.)
       Gww_u_shaded=1.0/(1.0/G_u_a+1.0/G_u_b+100.)

       ! temperatures of overstorey and understorey canopies
       Tco=(Tc_o_sunlit_old*LAI_o_sunlit+Tc_o_shaded_old*LAI_o_shaded)/(LAI_o_sunlit+LAI_o_shaded)
       Tcu=(Tc_u_sunlit_old*LAI_u_sunlit+Tc_u_shaded_old*LAI_u_shaded)/(LAI_u_sunlit+LAI_u_shaded)

       ! net Radiation at canopy and leaf level module by X.Luo
       call netRadiation(meteo%S_dff,meteo%S_dir,CosZs,Tco,Tcu,temp_grd,lai_o,lai_u,lai_o+stem_o,lai_u+stem_u,&
                         LAI_o_sunlit,LAI_o_shaded,LAI_u_sunlit,LAI_u_shaded,clumping,temp_air,rh_air,alpha_v_sw(kkk),&
                         alpha_n_sw(kkk),percentArea_snow_o,percentArea_snow_u,Xg_snow(kkk),alpha_v_o,alpha_n_o,&
                         alpha_v_u,alpha_n_u,alpha_v_g,alpha_n_g,radiation_o,radiation_u,radiation_g,radiation_o_sun,&
                         radiation_o_shaded,radiation_u_sun,radiation_u_shaded,R_o_sunlit,R_o_shaded,R_u_sunlit,&
                         R_u_shaded,fapar)

       ! photosynthesis module by B. Chen
       Gw_o_sunlit  = 1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_sunlit_old) !conductance of sunlit leaves of overstorey for water
       Gw_o_shaded = 1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_shaded_old)
       Gw_u_sunlit = 1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_sunlit_old)
       Gw_u_shaded = 1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_shaded_old)

       leleaf_o_sunlit = Gw_o_sunlit*(VPD_air+slope*(Tc_o_sunlit_old-temp_air))*rho_a * Cp_ca/psychrometer
       leleaf_o_shaded = Gw_o_shaded*(VPD_air+slope*(Tc_o_shaded_old-temp_air))*rho_a * Cp_ca/psychrometer
       leleaf_u_sunlit = Gw_u_sunlit*(VPD_air+slope*(Tc_u_sunlit_old-temp_air))*rho_a * Cp_ca/psychrometer
       leleaf_u_shaded = Gw_u_shaded*(VPD_air+slope*(Tc_u_shaded_old-temp_air))*rho_a * Cp_ca/psychrometer

       if(CosZs>0.) then
          call photosynthesis(landcover,Tc_o_sunlit_old,f_leaf,p_kc25,p_ko25,p_tau25,R_o_sunlit,e_a10,&
                              G_o_b,Vcmax_sunlit,VJ_slope,f_soilwater,b_h2o,m_h2o,&
                              Ci_o_sunlit_old,temp_air,leleaf_o_sunlit,Gs_o_sunlit_new,gs_h2o,Ac_o_sunlit,&
                              Ci_o_sunlit_new,ffpa,sif_alpha,sif_beta,lSIF_o_sunlit,COSi_o_sunlit_old,&
                              COSi_o_sunlit_new,lCOS_o_sunlit)
          call photosynthesis(landcover,Tc_o_shaded_old,f_leaf,p_kc25,p_ko25,p_tau25,R_o_shaded,e_a10,&
                              G_o_b,Vcmax_shaded,VJ_slope,f_soilwater,b_h2o,m_h2o,&
                              Ci_o_shaded_old,temp_air,leleaf_o_shaded,Gs_o_shaded_new,gs_h2o,Ac_o_shaded,&
                              Ci_o_shaded_new,ffpa,sif_alpha,sif_beta,lSIF_o_shaded,COSi_o_shaded_old,&
                              COSi_o_shaded_new,lCOS_o_shaded)
          call photosynthesis(landcover,Tc_u_sunlit_old,f_leaf,p_kc25,p_ko25,p_tau25,R_u_sunlit,e_a10,&
                              G_u_b,Vcmax_sunlit,VJ_slope,f_soilwater,b_h2o,m_h2o,&
                              Ci_u_sunlit_old,temp_air,leleaf_u_sunlit,Gs_u_sunlit_new,gs_h2o,Ac_u_sunlit,&
                              Ci_u_sunlit_new,ffpa,sif_alpha,sif_beta,lSIF_u_sunlit,COSi_u_sunlit_old,&
                              COSi_u_sunlit_new,lCOS_u_sunlit)
          call photosynthesis(landcover,Tc_u_shaded_old,f_leaf,p_kc25,p_ko25,p_tau25,R_u_shaded,e_a10,&
                              G_u_b,Vcmax_shaded,VJ_slope,f_soilwater,b_h2o,m_h2o,&
                              Ci_u_shaded_old,temp_air,leleaf_u_shaded,Gs_u_shaded_new,gs_h2o,Ac_u_shaded,&
                              Ci_u_shaded_new,ffpa,sif_alpha,sif_beta,lSIF_u_shaded,COSi_u_shaded_old,&
                              COSi_u_shaded_new,lCOS_u_shaded)
        else
          Gs_o_sunlit_new=0.0001
          Ac_o_sunlit=0.0
          lSIF_o_sunlit=0.0
          lCOS_o_sunlit=0.0
          Ci_o_sunlit_new=CO2_air*0.7
          Cs_o_sunlit_new=CO2_air
          Cc_o_sunlit_new=CO2_air*0.7*0.8
          COSi_o_sunlit_new=COS_air*0.7
          COSs_o_sunlit_new=COS_air
          COSc_o_sunlit_new=COS_air*0.7*0.8

          Gs_o_shaded_new=0.0001
          Ac_o_shaded=0.0
          lSIF_o_shaded=0.0
          lCOS_o_shaded=0.0
          Ci_o_shaded_new=CO2_air*0.7
          Cs_o_shaded_new=CO2_air
          Cc_o_shaded_new=CO2_air*0.7*0.8
          COSi_o_shaded_new=COS_air*0.7
          COSs_o_shaded_new=COS_air
          COSc_o_shaded_new=COS_air*0.7*0.8

          Gs_u_sunlit_new=0.0001
          Ac_u_sunlit=0.0
          lSIF_u_sunlit=0.
          lCOS_u_sunlit=0.
          Ci_u_sunlit_new=CO2_air*0.7
          Cs_u_sunlit_new=CO2_air
          Cc_u_sunlit_new=CO2_air*0.7*0.8
          COSi_u_sunlit_new=COS_air*0.7
          COSs_u_sunlit_new=COS_air
          COSc_u_sunlit_new=COS_air*0.7*0.8

          Gs_u_shaded_new=0.0001
          Ac_u_shaded=0.0
          lSIF_u_shaded=0.
          lCOS_u_shaded=0.
          Ci_u_shaded_new=CO2_air*0.7
          Cs_u_shaded_new=CO2_air
          Cc_u_shaded_new=CO2_air*0.7*0.8
          COSi_u_shaded_new=COS_air*0.7
          COSs_u_shaded_new=COS_air
          COSc_u_shaded_new=COS_air*0.7*0.8

       end if
!       write(*,*) G_o_b,gs_h2o
       Ci_o_sunlit_old=Ci_o_sunlit_new
       Cs_o_sunlit_old=Cs_o_sunlit_new
       COSi_o_sunlit_old=COSi_o_sunlit_new
       COSs_o_sunlit_old=COSs_o_sunlit_new
       Gs_o_sunlit_old=Gs_o_sunlit_new
       Gw_o_sunlit=1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_sunlit_new) !conductance of sunlit leaves of overstorey for water
       Gc_o_sunlit=1.0/(1.0/G_o_a+1.4/G_o_b+1.6/Gs_o_sunlit_new) !conductance of sunlit leaves of overstorey for CO2

       Ci_o_shaded_old=Ci_o_shaded_new
       Cs_o_shaded_old=Cs_o_shaded_new
       COSi_o_shaded_old=COSi_o_shaded_new
       COSs_o_shaded_old=COSs_o_shaded_new
       Gs_o_shaded_old=Gs_o_shaded_new
       Gw_o_shaded=1.0/(1.0/G_o_a+1.0/G_o_b+1.0/Gs_o_shaded_new)
       Gc_o_shaded=1.0/(1.0/G_o_a+1.4/G_o_b+1.6/Gs_o_shaded_new)

       Ci_u_sunlit_old=Ci_o_sunlit_new
       Cs_u_sunlit_old=Cs_u_sunlit_new
       COSi_u_sunlit_old=COSi_o_sunlit_new
       COSs_u_sunlit_old=COSs_u_sunlit_new
       Gs_u_sunlit_old=Gs_o_sunlit_new
       Gw_u_sunlit=1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_sunlit_new)
       Gc_u_sunlit=1.0/(1.0/G_u_a+1.4/G_u_b+1.6/Gs_u_sunlit_new)

       Ci_u_shaded_old=Ci_u_shaded_new
       Cs_u_shaded_old=Cs_u_shaded_new
       COSi_u_shaded_old=COSi_u_shaded_new
       COSs_u_shaded_old=COSs_u_shaded_new
       Gs_u_shaded_old=Gs_u_shaded_new
       Gw_u_shaded=1.0/(1.0/G_u_a+1.0/G_u_b+1.0/Gs_u_shaded_new)
       Gc_u_shaded=1.0/(1.0/G_u_a+1.4/G_u_b+1.6/Gs_u_shaded_new)

       ! leaf tempratures module by L. He
       call Leaf_Temperatures(temp_air,slope,psychrometer,VPD_air,Cp_ca,Gw_o_sunlit,Gw_o_shaded,Gw_u_sunlit,Gw_u_shaded,&
                 Gww_o_sunlit,Gww_o_shaded,Gww_u_sunlit,Gww_u_shaded,Gh_o_sunlit,Gh_o_shaded,Gh_u_sunlit,Gh_u_shaded,&
                 Xcs_o(kkk),Xcl_o(kkk),Xcs_u(kkk),Xcl_u(kkk),radiation_o_sun,radiation_o_shaded,radiation_u_sun,&
                 radiation_u_shaded,Tc_o_sunlit_new,Tc_o_shaded_new,Tc_u_sunlit_new,Tc_u_shaded_new)

       H_o_sunlit =(Tc_o_sunlit_new-temp_air)*rho_a * Cp_ca*Gh_o_sunlit
       H_o_shaded =(Tc_o_shaded_new-temp_air)*rho_a * Cp_ca*Gh_o_shaded
       GH_o=H_o_sunlit*LAI_o_sunlit+H_o_shaded*LAI_o_shaded  !for next num aerodynamic conductance calculation

       if(abs(Tc_o_sunlit_new-Tc_o_sunlit_old)<0.02 .and. abs(Tc_o_shaded_new-Tc_o_shaded_old)<0.02 .and. &
          abs(Tc_u_sunlit_new-Tc_u_sunlit_old)<0.02 .and. abs(Tc_u_shaded_new-Tc_u_shaded_old)<0.02) then
          exit
       else
          if(num >22) then   !iteration does not converge.
             Tc_o_sunlit_old=temp_air
             Tc_o_shaded_old=temp_air
             Tc_u_sunlit_old=temp_air
             Tc_u_shaded_old=temp_air
             exit
           else
             Tc_o_sunlit_old=Tc_o_sunlit_new
             Tc_o_shaded_old=Tc_o_shaded_new
             Tc_u_sunlit_old=Tc_u_sunlit_new
             Tc_u_shaded_old=Tc_u_shaded_new
           end if
         end if
    end do     ! end do while
!    write(*,*) G_o_b,gs_h2o
    GPP_o_sunlit= Ac_o_sunlit*LAIo_sunlit
    GPP_o_shaded= Ac_o_shaded*LAIo_shaded
    GPP_u_sunlit= Ac_u_sunlit*LAIu_sunlit
    GPP_u_shaded= Ac_u_shaded*LAIu_shaded
    
    COS_o_sunlit= lCOS_o_sunlit*LAIo_sunlit
    COS_o_shaded= lCOS_o_shaded*LAIo_shaded
    COS_u_sunlit= lCOS_u_sunlit*LAIu_sunlit
    COS_u_shaded= lCOS_u_shaded*LAIu_shaded
  
!    stSIF_o_sunlit=lSIF_o_sunlit*0.5*clumping*(1.1-0.1*LAIo_sunlit)*exp(-CosZs)*LAIo_sunlit
!    stSIF_o_shaded=lSIF_o_shaded*0.5*clumping*(1.1-0.1*LAIo_shaded)*exp(-CosZs)*LAIo_shaded
!    stSIF_u_sunlit=lSIF_u_sunlit*0.5*clumping*(1.1-0.1*LAIu_sunlit)*exp(-CosZs)*LAIu_sunlit
!    stSIF_u_shaded=lSIF_u_shaded*0.5*clumping*(1.1-0.1*LAIu_shaded)*exp(-CosZs)*LAIu_shaded

!    SIF_o_sunlit=lSIF_o_sunlit*exp(-0.5*clumping*LAIo_sunlit/CosZs)*LAIo_sunlit+stSIF_o_sunlit
!    SIF_o_shaded=lSIF_o_shaded*exp(-0.5*clumping*LAIo_shaded/CosZs)*LAIo_shaded+stSIF_o_shaded
!    SIF_u_sunlit=lSIF_u_sunlit*exp(-0.5*clumping*LAIu_sunlit/CosZs)*LAIu_sunlit+stSIF_u_sunlit
!    SIF_u_shaded=lSIF_u_shaded*exp(-0.5*clumping*LAIu_shaded/CosZs)*LAIu_shaded+stSIF_u_shaded

    stSIF_o_sunlit=lSIF_o_sunlit*0.3*clumping*(1.1-0.1*LAIo_sunlit)*exp(-CosZs)
    stSIF_o_shaded=lSIF_o_shaded*0.3*clumping*(1.1-0.1*LAIo_shaded)*exp(-CosZs)
    stSIF_u_sunlit=lSIF_u_sunlit*0.3*clumping*(1.1-0.1*LAIu_sunlit)*exp(-CosZs)
    stSIF_u_shaded=lSIF_u_shaded*0.3*clumping*(1.1-0.1*LAIu_shaded)*exp(-CosZs)

!    stSIF_o_sunlit = 0.0   !!@JWang  for scattering correction,@MOUSONG,make sure this is correct???
!    stSIF_o_shaded = 0.0
!    stSIF_u_sunlit = 0.0
!    stSIF_u_shaded = 0.0

    SIF_o_sunlit=(lSIF_o_sunlit+stSIF_o_sunlit)*LAIo_sunlit*0.1
    SIF_o_shaded=(lSIF_o_shaded+stSIF_o_shaded)*LAIo_shaded*0.1
    SIF_u_sunlit=(lSIF_u_sunlit+stSIF_u_sunlit)*LAIu_sunlit*0.1
    SIF_u_shaded=(lSIF_u_shaded+stSIF_u_shaded)*LAIu_shaded*0.1

    !Transpiration module by X. Luo
    call transpiration(Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new,temp_air, rh_air,&
             Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded,LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,&
             Trans_o(kkk), Trans_u(kkk))

    ! Evaporation and sublimation from canopy by X. Luo
    call evaporation_canopy (Tc_o_sunlit_new, Tc_o_shaded_new, Tc_u_sunlit_new, Tc_u_shaded_new,temp_air, rh_air,&
            Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded,LAI_o_sunlit, LAI_o_shaded, LAI_u_sunlit, LAI_u_shaded,&
            Xcl_o(kkk),Xcl_u(kkk),Xcs_o(kkk),Xcs_u(kkk),Eil_o(kkk),Eil_u(kkk),EiS_o(kkk),EiS_u(kkk))

    ! Rainfall stage 2 by X. Luo
    call rainfall_stage2(Eil_o(kkk),Eil_u(kkk),Wcl_o(kkk),Wcl_u(kkk))
    ! Snow pack stage2
    call snowpack_stage2(EiS_o(kkk),EiS_u(kkk),Wcs_o(kkk),Wcs_u(kkk))

    ! Evaporation from soil module
    Gheat_g  = 1./ra_g
    mass_water_g  = rho_w*Zp
    call Soil_evaporation(temp_grd,Ts0(kkk-1),rh_air,radiation_g,Gheat_g,Xg_snow(kkk),Zp,Zsp,mass_water_g,&
           Wg_snow(kkk),rho_snow(kkk),soilp%thetam_prev(0),soilp%fei(0),Evap_soil(kkk),Evap_SW(kkk),Evap_SS(kkk))

    Zp        = mass_water_g/rho_w      ! update surface ponding after ponding evaporation calculation
    Zsp       = Wg_snow(kkk)/rho_snow(kkk)    ! update snow depth as well after snow evaporation calculation

    ! to be checked later:  why set these 4 to 0
    Eil_o(kkk) = 0.
    EiS_o(kkk) = 0.
    Eil_u(kkk) = 0.
    EiS_u(kkk) = 0.

    ! soil Thermal Conductivity module by L. He
    call UpdateSoilThermalConductivity(soilp)
    call Update_Cs(soilp)

    ! Surface temperature
    Cs(0,kkk)  = soilp%Cs(0)
    Cs(1,kkk)  = soilp%Cs(0)
    Tc_u(kkk)  = Tcu
    lambda_soil(1)  = soilp%lambda(0)
    d_soil(1)  = soilp%d_soil(0)
    Tm(1,kkk-1)= soilp%temp_soil_p(1)
    Tm(0,kkk-1)= soilp%temp_soil_p(0)
    G(1,kkk)   = soilp%G(0)

    call SurfaceTemperature(temp_air,rh_air,Zsp,Zp,Cs(1,kkk),Cs(0,kkk),Gheat_g,d_soil(1),rho_snow(kkk),Tc_u(kkk),&
                  radiation_g,Evap_soil(kkk),Evap_SW(kkk),Evap_SS(kkk),lambda_soil(1),Xg_snow(kkk),G(1,kkk),Ts0(kkk-1),&
                  Tm(1,kkk-1),Tm(0,kkk-1),Tsn0(kkk-1),Tsm0(kkk-1),Tsn1(kkk-1),Tsn2(kkk-1),Ts0(kkk),Tm(0,kkk),&
                  Tsn0(kkk),Tsm0(kkk),Tsn1(kkk),Tsn2(kkk),G(0,kkk))
    soilp%temp_soil_c(0)=Tm(0,kkk)

    ! Snow pack stage3 module
    call snowpack_stage3(temp_air, Tsn0(kkk), Tsn0(kkk-1),rho_snow(kkk),Zsp,Zp,Wg_snow(kkk))

    call SensibleHeat(Tc_o_sunlit_new,Tc_o_shaded_new,Tc_u_sunlit_new,Tc_u_shaded_new,Ts0(kkk),temp_air,rh_air,&
              Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded, Gheat_g,LAI_o_sunlit, LAI_o_shaded, &
              LAI_u_sunlit,LAI_u_shaded,Qhc_o(kkk),Qhc_u(kkk),Qhg(kkk))

    ! soil water module
    soilp%Zsp   = Zsp
    soilp%G(0)  = G(0,kkk)

    call UpdateHeatFlux(soilp,Xg_snow(kkk),lambda_snow(kkk),Tsn0(kkk),temp_air,kstep)
    call Soil_water_uptake(soilp,Trans_o(kkk),Trans_u(kkk),Evap_soil(kkk))
  
    soilp%r_rain_g  = r_rain_g(kkk)
    soilp%Zp        = Zp

    call UpdateSoilMoisture(soilp)
    Zp              = soilp%Zp
   
end do       !END kkk iteration
!    write(*,*) G_o_b,gs_h2o
    kkk  = kloop
    if(Tsn1(kkk) > 40.) Tsn2(kkk) = 40.   !True? or Tsn1 @J.Wang
    if(Tsn1(kkk) <-40.) Tsn2(kkk) = -40.
    if(Tsn2(kkk) > 40.) Tsn2(kkk) = 40.
    if(Tsn2(kkk) <-40.) Tsn2(kkk) = -40.

    var_n(1)    = Qhc_o(kkk)    ! SH
    var_n(3)    = Ts0(kkk)      ! The temperature of ground surface
    var_n(4)    = Tsn0(kkk)     !The temperature of ground surface
    var_n(5)    = Tsm0(kkk)
    var_n(6)    = Tsn1(kkk)
    var_n(7)    = Tsn2(kkk)

    do i=9,14
       var_n(i) = soilp%temp_soil_c(i-9)
    end do
    do i=21,26
       var_n(i) = soilp%thetam(i-21)
    end do
    do i=27,32
       var_n(i) = soilp%ice_ratio(i-27)
    end do
    var_n(15)   = Wcl_o(kkk)
    var_n(16)   = Wcs_o(kkk) !the mass of intercepted liquid water and snow, overstory
    var_n(18)   = Wcl_u(kkk)
    var_n(19)   = Wcs_u(kkk)
    var_n(20)   = Wg_snow(kkk) ! fraction of ground surface covered by snow and snow mass

    var_n(33)   = soilp%Zsp
    var_n(34)   = soilp%Zp
    var_n(35)   = soilp%r_rain_g
    var_n(36)   = Ac_snow_o(kkk)
    var_n(37)   = Ac_snow_u(kkk)
    var_n(38)   = rho_snow(kkk)
    var_n(39)   = alpha_v_sw(kkk)
    var_n(40)   = alpha_n_sw(kkk)

    Lv_liquid   = (2.501 - 0.00237 * temp_air) * 1000000.  !laten heat of water vaporization in j/kg

    !for output
    mid_res%Net_Rad   = radiation_o+radiation_u+radiation_g
    mid_res%LH        = Lv_liquid *(Trans_o(kkk)+Eil_o(kkk)+Trans_u(kkk)+Eil_u(kkk)+Evap_soil(kkk)+Evap_SW(kkk))+&
                        Lv_solid*(EiS_o(kkk)+EiS_u(kkk)+Evap_SS(kkk))
    mid_res%SH        = Qhc_o(kkk)+Qhc_u(kkk)+Qhg(kkk)
    mid_res%Trans     = (Trans_o(kkk)+Trans_u(kkk))/rho_w
    mid_res%Evap      = (Eil_o(kkk)+Eil_u(kkk)+Evap_soil(kkk)+Evap_SW(kkk))/rho_w+ &
                        (EiS_o(kkk)+EiS_u(kkk)+Evap_SS(kkk))/rho_snow(kkk)
    mid_res%gpp_o_sunlit = GPP_o_sunlit*12.*1.e-6*1.e-3    !J.Wang kg/m2/s// umol C/m2/s
    mid_res%gpp_u_sunlit = GPP_u_sunlit*12.*1.e-6*1.e-3   
    mid_res%gpp_o_shaded = GPP_o_shaded*12.*1.e-6*1.e-3
    mid_res%gpp_u_shaded = GPP_u_shaded*12.*1.e-6*1.e-3

    mid_res%GPP          = mid_res%gpp_o_sunlit+mid_res%gpp_u_sunlit+mid_res%gpp_o_shaded+mid_res%gpp_u_shaded
    mid_res%SIF          = SIF_o_sunlit+SIF_o_shaded+SIF_u_sunlit+SIF_u_shaded
    mid_res%thetam_surf  = soilp%thetam(0)
    mid_res%COS_plant    = COS_o_sunlit+COS_o_shaded+COS_u_sunlit+COS_u_shaded      ! pmol/m2/s

!    write(*,*) "thetam_surf = ", mid_res%thetam_surf
    
    if(lai_input < 0) then
       temp_day = temp_day + meteo%temp
       theta_day = theta_day + (soilp%theta_vfc(0) - soilp%theta_vwp(0))
       trans_day = trans_day + mid_res%Trans
       if(is_end_curr_day) then
          temp_day = temp_day/24.
          theta_day = theta_day/24.
          trans_day = trans_day/24.

    ! use the method in BETHY to calculate phenology, with a little modification, @Mousong.Wu,201905

          call beps_phenology(lc,daylen,temp_day,theta_day,trans_day,mid_res%lai_old)

          temp_day = 0.
          theta_day = 0.
          trans_day = 0.

          mid_res%lai_new = mid_res%lai_old
       else
          mid_res%lai_new = mid_res%lai_old
       end if
          lai = mid_res%lai_new 
          !mid_res%fAPAR = 1. - exp(-0.45*lai)    ! calculate fAPAR using the Lambert-Beer law, Benjamin Smith et al., 2008, &
                                                 ! & Forest Ecology and Management, @Mousong.Wu, 201905
!          write(*,*) 'lai = ', mid_res%lai_old
!          write(*,*) 'fAPAR = ', mid_res%fAPAR
    end if
    !lai = mid_res%lai_new 
    ! empirical relationship between fAPAR and LAI, @Mousong.Wu, 201905
    mid_res%fAPAR = 1. - exp(-0.45*lai)    ! calculate fAPAR using the Lambert-Beer law, Benjamin Smith et al., 2008, &
    ! implement the physical based fAPAR calculation, @Mousong.Wu, 2024
    !mid_res%fAPAR = fapar
   call cos_grnd(soilp,cos_soil)
   mid_res%COS_grnd = cos_soil             ! pmol/m2/s

    return
end subroutine


