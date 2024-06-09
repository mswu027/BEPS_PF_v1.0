module bepstype
  use shr_kind_mod,only:r8=>shr_kind_r8
  !--iLab::no need for use-statement
  ! use beps_par
  implicit none

  ! For restart
  !-- iLab::module field should have save attribute
  ! 'save' means values will not be changed
  real(r8), allocatable,public,save     ::  v2last(:,:,:)


  !*********************Climate Forcing datasets**********************
  type,public:: forc
     real(r8),pointer::  Temp(:)
     real(r8),pointer::  Tempmx(:)
     real(r8),pointer::  Tempmn(:)
     real(r8),pointer::  Wind(:)
#ifdef COUP_CSM
     real(r8),pointer::  Zref(:)   ! bottom atm level height
     real(r8),pointer::  Rain(:)   ! precipitation (liquid)
     real(r8),pointer::  Snow(:)   ! Snow rate
     real(r8),pointer::  Swndr(:)  ! SW: nir direct
     real(r8),pointer::  Swvdr(:)  ! SW: vis direct
     real(r8),pointer::  Swndf(:)  ! SW: nir diffuse
     real(r8),pointer::  Swvdf(:)  ! SW: vis diffuse
     real(r8),pointer::  Swdr(:)   ! SW direct radiation
     real(r8),pointer::  Swdf(:)   ! SW diffuse radiation
     real(r8),pointer::  Lwdn(:)   ! downward LW heat flux
     real(r8),pointer::  shum(:)   ! specific humidity(kg/kg)
     real(r8),pointer::  pres(:)   ! pressure (Pa)
#else
     real(r8),pointer::  Srad(:)   ! downward solar raditation
     real(r8),pointer::  Rh(:)     ! relative humidity(%)
     real(r8),pointer::  Rain(:)   ! liquid precipitation
     real(r8),pointer::  Snow(:)   ! Snow rate
     real(r8),pointer::  Swdr(:)   ! SW direct radiation
     real(r8),pointer::  Swdf(:)   ! SW diffuse radiation

     character(len=16) :: meteo_ref_yyyymmdd !--iLab::to consistenly handle tempooral settings
                                             !--      expected format yyyy-mm-dd
#endif
  end type forc
  ! clim is a derived type
  type(forc),save,target,public:: clim

  !***********************for CPL datasets****************************
#ifdef COUP_CSM
  type,public:: CPL
     real(r8),pointer:: rofliq(:)  ! lnd->rtm
     real(r8),pointer:: rofice(:)  ! lnd->rtm
     real(r8),pointer:: t_Rad(:)   ! RAD tem
     real(r8),pointer:: tref(:)    ! 2m temperature
     real(r8),pointer:: qref(:)    ! 2m specific humidity
     real(r8),pointer:: avsdr(:)   ! albedo: direct , visible
     real(r8),pointer:: anidr(:)   ! albedo: direct , near-ir
     real(r8),pointer:: avsdf(:)   ! albedo: diffuse, vis
     real(r8),pointer:: anidf(:)   ! abdedo: diffuse, nis
     real(r8),pointer:: snowh(:)   ! snow height
     real(r8),pointer:: u10(:)     ! 10m wind
     real(r8),pointer:: ddvel(:)   ! dry deposition velocity (optinal)
     real(r8),pointer:: fv(:)      ! friction velocity
     real(r8),pointer:: ram1(:)    ! aerodynamical resistance
     real(r8),pointer:: soilw(:)   ! volumetric soil water
     real(r8),pointer:: taux(:)    ! wind stress
     real(r8),pointer:: tauy(:)    !
     real(r8),pointer:: LH(:)      ! latent heat
     real(r8),pointer:: SH(:)      ! sensible heat
     real(r8),pointer:: lwup(:)    ! upward longwave heat flux
     real(r8),pointer:: evap(:)    ! evaporation
     real(r8),pointer:: swnet(:)   ! heat flux
     real(r8),pointer:: fco2(:)    ! co2 flux
     real(r8),pointer:: flxdst1(:) ! dust flux size bin 1
     real(r8),pointer:: flxdst2(:) !                    2
     real(r8),pointer:: flxdst3(:)
     real(r8),pointer:: flxdst4(:)
     real(r8),pointer:: flxvoc(:)
  end type CPL

  type(CPL),save,target,public :: lnd2atm
#endif

  !******************for boundary/yrdata/lai/cpools/
  type,public::surf
     integer,pointer::  lcno(:,:)     ! PFT types
     integer,pointer::  stext(:)      ! soil texture
     real(r8),pointer:: PCT_PFT(:,:) ! PFT fraction
     real(r8),pointer:: clumping(:)
     real(r8),pointer:: longitude(:)
     real(r8),pointer:: latitude(:)

     real(r8),pointer:: sdp(:)       ! snowdepth
     real(r8),pointer:: st(:)        !
     real(r8),pointer:: sw(:)        ! soil moisture

     real(r8),pointer:: laiyr(:,:)   ! for plant resp
     real(r8),pointer:: nppyr(:,:)   ! for soil resp

!!! soil carbon pools  units(g C/m2)
     real(r8),pointer:: ccd(:,:)
     real(r8),pointer:: cfmd(:,:)
     real(r8),pointer:: cfsd(:,:)
     real(r8),pointer:: cm(:,:)
     real(r8),pointer:: cp(:,:)
     real(r8),pointer:: cs(:,:)
     real(r8),pointer:: csm(:,:)
     real(r8),pointer:: csmd(:,:)
     real(r8),pointer:: cssd(:,:)
     real(r8),pointer:: lai(:,:)     ! for photosynthesis
     real(r8),pointer:: Vcmax(:,:)   ! for data assimilation
     ! 2024/03/29  add the canopy height
     real(r8),pointer:: HeightC(:)
     !     character,pointer:: name(:)
  end type surf

  type(surf),save,target,public:: bound   ! boundary conditions

  !******************for assimilation and parameter optimization/
  type,public::para
 ! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
     real(r8),pointer:: p_Vcmax(:,:)
     real(r8),pointer:: p_VJ_slope(:,:)
     real(r8),pointer:: p_VN_slope(:,:)
     real(r8),pointer:: p_b_h2o(:,:)
     real(r8),pointer:: p_m_h2o(:,:)
     real(r8),pointer:: p_f_leaf(:,:)
     real(r8),pointer:: p_kc25(:,:)
     real(r8),pointer:: p_ko25(:,:)
     real(r8),pointer:: p_tau25(:,:)

 ! parameters for SIF modelling
     real(r8),pointer:: p_sif_alpha(:,:)
     real(r8),pointer:: p_sif_beta(:,:)

! parameters related to plant respiration: plant_resp
     real(r8),pointer:: p_q10(:,:)
! parameters related to soil respiration: soil_resp
     real(r8),pointer:: p_f_resp(:,:)
! soil parameters
     real(r8),pointer:: p_f_decay(:,:)  ! scalar factor for root distribution
     real(r8),pointer:: p_Ksat_scalar(:,:)
     real(r8),pointer:: p_b_scalar(:,:)
     real(r8),pointer:: p_porosity_scalar(:,:)
     real(r8),pointer:: p_vfc_scalar(:,:)
     real(r8),pointer:: p_vwp_scalar(:,:)
     real(r8),pointer:: p_psisat_scalar(:,:)
     real(r8),pointer:: p_drainage_scalar(:,:)
 ! parameters for vod modelling
     real(r8),pointer:: p_vod_a(:,:)
     real(r8),pointer:: p_vod_b(:,:)
     real(r8),pointer:: p_vod_c(:,:)
! parameters for SPAC in soilwateruptake
     real(r8),pointer:: p_theta_Amin(:,:)
     real(r8),pointer:: p_pox(:,:)
     real(r8),pointer:: p_fei_c(:,:)
     real(r8),pointer:: p_spac_p1(:,:)
     real(r8),pointer:: p_spac_p2(:,:)
     real(r8),pointer:: p_tWA(:,:)
     real(r8),pointer:: p_tWB(:,:)
     real(r8),pointer:: p_Ttrig(:,:)
     real(r8),pointer:: p_r_xylem(:,:)
     real(r8),pointer:: p_r_r(:,:)
     real(r8),pointer:: p_Lr(:,:)
     real(r8),pointer:: p_deltal_min(:,:)
     real(r8),pointer:: p_deltal_max(:,:)
     real(r8),pointer:: p_p_delta(:,:)
     real(r8),pointer:: p_ppslh(:,:)
     real(r8),pointer:: p_fei_min(:,:)
     real(r8),pointer:: p_fei_th(:,:)
     real(r8),pointer:: p_p_excess(:,:)
     real(r8),pointer:: p_Tleaf_H(:,:)
     real(r8),pointer:: p_Tleaf_L(:,:)
     real(r8),pointer:: p_Tleaf_O(:,:)

 ! parameters' uncertainty
 ! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
     real(r8),pointer:: u_Vcmax(:,:)
     real(r8),pointer:: u_VJ_slope(:,:)
     real(r8),pointer:: u_VN_slope(:,:)
     real(r8),pointer:: u_b_h2o(:,:)
     real(r8),pointer:: u_m_h2o(:,:)
     real(r8),pointer:: u_f_leaf(:,:)
     real(r8),pointer:: u_kc25(:,:)
     real(r8),pointer:: u_ko25(:,:)
     real(r8),pointer:: u_tau25(:,:)

 ! parameters for SIF modelling
     real(r8),pointer:: u_sif_alpha(:,:)
     real(r8),pointer:: u_sif_beta(:,:)

! parameters related to plant respiration: plant_resp
     real(r8),pointer:: u_q10(:,:)
! parameters related to soil respiration: soil_resp
     real(r8),pointer:: u_f_resp(:,:)
! soil parameters
     real(r8),pointer:: u_f_decay(:,:)  ! scalar factor for root distribution
     real(r8),pointer:: u_Ksat_scalar(:,:)
     real(r8),pointer:: u_b_scalar(:,:)
     real(r8),pointer:: u_porosity_scalar(:,:)
     real(r8),pointer:: u_vfc_scalar(:,:)
     real(r8),pointer:: u_vwp_scalar(:,:)
     real(r8),pointer:: u_psisat_scalar(:,:)
     real(r8),pointer:: u_drainage_scalar(:,:)
 ! parameters for vod modelling
     real(r8),pointer:: u_vod_a(:,:)
     real(r8),pointer:: u_vod_b(:,:)
     real(r8),pointer:: u_vod_c(:,:)
! parameters for SPAC in soilwateruptake
     real(r8),pointer:: u_theta_Amin(:,:)
     real(r8),pointer:: u_pox(:,:)
     real(r8),pointer:: u_fei_c(:,:)
     real(r8),pointer:: u_spac_p1(:,:)
     real(r8),pointer:: u_spac_p2(:,:)
     real(r8),pointer:: u_tWA(:,:)
     real(r8),pointer:: u_tWB(:,:)
     real(r8),pointer:: u_Ttrig(:,:)
     real(r8),pointer:: u_r_xylem(:,:)
     real(r8),pointer:: u_r_r(:,:)
     real(r8),pointer:: u_Lr(:,:)
     real(r8),pointer:: u_deltal_min(:,:)
     real(r8),pointer:: u_deltal_max(:,:)
     real(r8),pointer:: u_p_delta(:,:)
     real(r8),pointer:: u_ppslh(:,:)
     real(r8),pointer:: u_fei_min(:,:)
     real(r8),pointer:: u_fei_th(:,:)
     real(r8),pointer:: u_p_excess(:,:)
     real(r8),pointer:: u_Tleaf_H(:,:)
     real(r8),pointer:: u_Tleaf_L(:,:)
     real(r8),pointer:: u_Tleaf_O(:,:)

  end type para

  type(para),save,target,public:: assim   ! optimization of parameters

  !******************for Particle Fliter/
  type,public::PF_para

     ! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
     real(r8),pointer:: Vcmax(:,:)
     real(r8),pointer:: VJ_slope(:,:)
     real(r8),pointer:: VN_slope(:,:)
     real(r8),pointer:: b_h2o(:,:)
     real(r8),pointer:: m_h2o(:,:)
     real(r8),pointer:: f_leaf(:,:)
     real(r8),pointer:: kc25(:,:)
     real(r8),pointer:: ko25(:,:)
     real(r8),pointer:: tau25(:,:)

 ! parameters for SIF modelling
     real(r8),pointer:: sif_alpha(:,:)
     real(r8),pointer:: sif_beta(:,:)

! parameters related to plant respiration: plant_resp
     real(r8),pointer:: q10(:,:)
! parameters related to soil respiration: soil_resp
     real(r8),pointer:: f_resp(:,:)
! soil parameters
     real(r8),pointer:: f_decay(:,:)  ! scalar factor for root distribution
     real(r8),pointer:: Ksat_scalar(:,:)
     real(r8),pointer:: b_scalar(:,:)
     real(r8),pointer:: porosity_scalar(:,:)
     real(r8),pointer:: vfc_scalar(:,:)
     real(r8),pointer:: vwp_scalar(:,:)
     real(r8),pointer:: psisat_scalar(:,:)
     real(r8),pointer:: drainage_scalar(:,:)
 ! parameters for vod modelling
     real(r8),pointer:: vod_a(:,:)
     real(r8),pointer:: vod_b(:,:)
     real(r8),pointer:: vod_c(:,:)
! parameters for SPAC in soilwateruptake
     real(r8),pointer:: theta_Amin(:,:)
     real(r8),pointer:: pox(:,:)
     real(r8),pointer:: fei_c(:,:)
     real(r8),pointer:: spac_p1(:,:)
     real(r8),pointer:: spac_p2(:,:)
     real(r8),pointer:: tWA(:,:)
     real(r8),pointer:: tWB(:,:)
     real(r8),pointer:: Ttrig(:,:)
     real(r8),pointer:: r_xylem(:,:)
     real(r8),pointer:: r_r(:,:)
     real(r8),pointer:: Lr(:,:)
     real(r8),pointer:: deltal_min(:,:)
     real(r8),pointer:: deltal_max(:,:)
     real(r8),pointer:: p_delta(:,:)
     real(r8),pointer:: ppslh(:,:)
     real(r8),pointer:: fei_min(:,:)
     real(r8),pointer:: fei_th(:,:)
     real(r8),pointer:: p_excess(:,:)
     real(r8),pointer:: Tleaf_H(:,:)
     real(r8),pointer:: Tleaf_L(:,:)
     real(r8),pointer:: Tleaf_O(:,:)

     real(r8),pointer:: pfweight(:,:)
     real(r8),pointer:: pfweightupdate(:,:)

  end type PF_para

  type(PF_para),save,target,public:: PF   ! PF parameters

  type,public::MC_para

     ! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
     real(r8),pointer:: Vcmax(:,:)
     real(r8),pointer:: VJ_slope(:,:)
     real(r8),pointer:: VN_slope(:,:)
     real(r8),pointer:: b_h2o(:,:)
     real(r8),pointer:: m_h2o(:,:)
     real(r8),pointer:: f_leaf(:,:)
     real(r8),pointer:: kc25(:,:)
     real(r8),pointer:: ko25(:,:)
     real(r8),pointer:: tau25(:,:)

 ! parameters for SIF modelling
     real(r8),pointer:: sif_alpha(:,:)
     real(r8),pointer:: sif_beta(:,:)

! parameters related to plant respiration: plant_resp
     real(r8),pointer:: q10(:,:)
! parameters related to soil respiration: soil_resp
     real(r8),pointer:: f_resp(:,:)
! soil parameters
     real(r8),pointer:: f_decay(:,:)  ! scalar factor for root distribution
     real(r8),pointer:: Ksat_scalar(:,:)
     real(r8),pointer:: b_scalar(:,:)
     real(r8),pointer:: porosity_scalar(:,:)
     real(r8),pointer:: vfc_scalar(:,:)
     real(r8),pointer:: vwp_scalar(:,:)
     real(r8),pointer:: psisat_scalar(:,:)
     real(r8),pointer:: drainage_scalar(:,:)
 ! parameters for vod modelling
     real(r8),pointer:: vod_a(:,:)
     real(r8),pointer:: vod_b(:,:)
     real(r8),pointer:: vod_c(:,:)
! parameters for SPAC in soilwateruptake
     real(r8),pointer:: theta_Amin(:,:)
     real(r8),pointer:: pox(:,:)
     real(r8),pointer:: fei_c(:,:)
     real(r8),pointer:: spac_p1(:,:)
     real(r8),pointer:: spac_p2(:,:)
     real(r8),pointer:: tWA(:,:)
     real(r8),pointer:: tWB(:,:)
     real(r8),pointer:: Ttrig(:,:)
     real(r8),pointer:: r_xylem(:,:)
     real(r8),pointer:: r_r(:,:)
     real(r8),pointer:: Lr(:,:)
     real(r8),pointer:: deltal_min(:,:)
     real(r8),pointer:: deltal_max(:,:)
     real(r8),pointer:: p_delta(:,:)
     real(r8),pointer:: ppslh(:,:)
     real(r8),pointer:: fei_min(:,:)
     real(r8),pointer:: fei_th(:,:)
     real(r8),pointer:: p_excess(:,:)
     real(r8),pointer:: Tleaf_H(:,:)
     real(r8),pointer:: Tleaf_L(:,:)
     real(r8),pointer:: Tleaf_O(:,:)

     ! performance metrics for SM
     real(r8),pointer:: SM_R2(:,:)
     real(r8),pointer:: SM_R(:,:)
     real(r8),pointer:: SM_RMSE(:,:)
     real(r8),pointer:: SM_ME(:,:)
     real(r8),pointer:: SM_AIC(:,:)
     real(r8),pointer:: SM_BIC(:,:)

     ! performance metrics for leaf water potential
     real(r8),pointer:: LWP_R2(:,:)
     real(r8),pointer:: LWP_R(:,:)
     real(r8),pointer:: LWP_RMSE(:,:)
     real(r8),pointer:: LWP_ME(:,:)
     real(r8),pointer:: LWP_AIC(:,:)
     real(r8),pointer:: LWP_BIC(:,:)

  end type MC_para
! Derived Types, MC
  type(MC_para),save,target,public:: MC

  !******************for Particle Fliter observation/
  type,public::PF_obs0
    ! real(r8),pointer:: obs_GPP(:)
    real(r8),pointer:: obs_var(:) ! using obs_var represents generalized target

  end type PF_obs0

  type(PF_obs0),save,target,public:: PF_obs

  type,public::MC_obs0
    ! real(r8),pointer:: obs_GPP(:)
    !real(r8),pointer:: obs_VOD(:) ! changed to VOD from obs_GPP @hulu 2023/06/30
    real(r8),pointer:: obs_SM(:)
    real(r8),pointer:: obs_LWP(:)

  end type MC_obs0

  type(MC_obs0),save,target,public:: MC_obs

  !******************for Particle Fliter resample/
  type,public::PF_resample0

    real(r8),pointer:: outparticles(:,:)
    real(r8),pointer:: resample_weight(:)
    real(r8),pointer:: resample_weight_update(:)

  end type PF_resample0

  type(PF_resample0),save,target,public:: PF_resample

  !*********************************************************************************************

  !------------------Soil Status----------
  !*******************************************************************
  type,public  :: soils
     integer,pointer   ::  n_layer(:)
     real(r8),pointer  ::  Zp(:,:)
     real(r8),pointer  ::  Zsp(:,:)
     real(r8),pointer  ::  r_rain_g(:,:)
     real(r8),pointer  ::  r_drainage(:,:)
     real(r8),pointer  ::  r_root_decay(:,:)
     real(r8),pointer  ::  psi_min(:,:)
     real(r8),pointer  ::  alpha(:,:)
     real(r8),pointer  ::  f_soilwater(:,:)

     real(r8),pointer  ::  Sp(:,:)  !! 2023/06/30


     real(r8),pointer  ::  d_soil(:,:)
     real(r8),pointer  ::  f_root(:,:,:)
     real(r8),pointer  ::  dt(:,:,:)
     real(r8),pointer  ::  thermal_cond(:,:,:)
     real(r8),pointer  ::  theta_vfc(:,:,:)
     real(r8),pointer  ::  theta_vwp(:,:,:)
     real(r8),pointer  ::  fei(:,:,:)
     real(r8),pointer  ::  Ksat(:,:,:)
     real(r8),pointer  ::  psi_sat(:,:,:)
     real(r8),pointer  ::  b(:,:,:)
     real(r8),pointer  ::  density_soil(:,:)   !!(npoints,0:Max_LayerS-1)
     real(r8),pointer  ::  f_org(:,:,:)
     real(r8),pointer  ::  ice_ratio(:,:,:)
     real(r8),pointer  ::  thetam(:,:,:)
     real(r8),pointer  ::  thetam_prev(:,:,:)
     real(r8),pointer  ::  temp_soil_p(:,:,:)
     real(r8),pointer  ::  temp_soil_c(:,:,:)
     real(r8),pointer  ::  f_ice(:,:,:)
     real(r8),pointer  ::  psim(:,:,:)
     real(r8),pointer  ::  thetab(:,:,:)
     real(r8),pointer  ::  psib(:,:,:)
     real(r8),pointer  ::  r_waterflow(:,:,:)
     real(r8),pointer  ::  km(:,:,:)
     real(r8),pointer  ::  kb(:,:,:)
     real(r8),pointer  ::  KK(:,:,:)
     real(r8),pointer  ::  Cs(:,:,:)
     real(r8),pointer  ::  lambda(:,:,:)
     real(r8),pointer  ::  Ett(:,:,:)
     real(r8),pointer  ::  G(:,:,:)
  end type soils

  type(soils),target,save,public :: soilstat

  !*************************************************for C4 xiuli*******************************
  !*************************************************param_gdd/
  type,public::param_gdd

    real(r8),pointer:: tt_veg(:,:)  ! thermal requirement of stage 1 of crop development (degree days).
    real(r8),pointer:: tt_rep(:,:)  ! thermal requirement of stage 2 of crop development (degree days).
    real(r8),pointer:: phot_type(:,:)  ! photoperiod yes=1 not=0
    real(r8),pointer:: emer_doy(:,:)  ! emergence doy
    real(r8),pointer:: har_doy(:,:)  ! harvest doy

  end type param_gdd

  type(param_gdd),target,save,public :: pgdd

  !*************************************************for C4 xiuli*******************************

  !*********************************************************
  !-----------Interest Variables For output-----------------
  !*********************************************************
  type,public :: res
     real(r8),pointer::  GPPpft(:,:)
     real(r8),pointer::  SIFpft(:,:)
     real(r8),pointer::  SIFpft_sat(:,:)   !accord with the satellite data
     real(r8),pointer::  NPPpft(:,:)
     real(r8),pointer::  NEPpft(:,:)
     real(r8),pointer::  SHpft(:,:)
     real(r8),pointer::  LHpft(:,:)
     real(r8),pointer::  Transpft(:,:)
     real(r8),pointer::  Evappft(:,:)
     real(r8),pointer::  Net_Radpft(:,:)
     real(r8),pointer::  GPP(:)
     real(r8),pointer::  SIF(:)
     real(r8),pointer::  SIF_sat(:)
     real(r8),pointer::  NPP(:)
     real(r8),pointer::  NEP(:)
     real(r8),pointer::  LAIpft(:,:)
     real(r8),pointer::  LAI(:)
     real(r8),pointer::  SH(:)
     real(r8),pointer::  LH(:)
     real(r8),pointer::  Trans(:)
     real(r8),pointer::  Evap(:)
     real(r8),pointer::  Net_Rad(:)
     real(r8),pointer::  Thetampft(:,:)
     real(r8),pointer::  Thetam(:) ! first soil layer
     real(r8),pointer::  fAPARpft(:,:)
     real(r8),pointer::  fAPAR(:)
     real(r8),pointer::  VODpft(:,:)
     real(r8),pointer::  VOD(:)
     real(r8),pointer::  COS_fluxpft(:,:)
     real(r8),pointer::  COS_flux(:)
     real(r8),pointer::  NPP_yr_acc(:,:)
     real(r8),pointer::  PWSpft(:,:)
     real(r8),pointer::  PWS(:)
     ! Eta =transpiration from SPAC
     real(r8),pointer::  Quptpft(:,:)
     real(r8),pointer::  Qupt(:)
     real(r8),pointer::  ETapft(:,:)
     real(r8),pointer::  ETa(:)
     real(r8),pointer::  fei_leaf(:)
     real(r8),pointer::  fei_leafpft(:,:)
     ! test for output different soil moisture and soil water potential
     real(r8),pointer::  Thetam2_pft(:,:)
     real(r8),pointer::  Thetam_layer2(:)
     real(r8),pointer::  Thetam3_pft(:,:)
     real(r8),pointer::  Thetam_layer3(:)
     real(r8),pointer::  Thetam4_pft(:,:)
     real(r8),pointer::  Thetam_layer4(:)
     real(r8),pointer::  Thetam5_pft(:,:)
     real(r8),pointer::  Thetam_layer5(:)
     real(r8),pointer::  SWP1_pft(:,:)
     real(r8),pointer::  SWP2_pft(:,:)
     real(r8),pointer::  SWP3_pft(:,:)
     real(r8),pointer::  SWP4_pft(:,:)
     real(r8),pointer::  SWP5_pft(:,:)
     real(r8),pointer::  SWP_layer1(:)
     real(r8),pointer::  SWP_layer2(:)
     real(r8),pointer::  SWP_layer3(:)
     real(r8),pointer::  SWP_layer4(:)
     real(r8),pointer::  SWP_layer5(:)
     real(r8),pointer::  TS_layer1(:)
     real(r8),pointer::  TS_layer2(:)
     real(r8),pointer::  TS_layer3(:)
     real(r8),pointer::  TS_layer4(:)
     real(r8),pointer::  TS_layer5(:)
     real(r8),pointer::  PondWater(:)
     real(r8),pointer::  Rain_g(:)
     real(r8),pointer::  TS1_pft(:,:)
     real(r8),pointer::  TS2_pft(:,:)
     real(r8),pointer::  TS3_pft(:,:)
     real(r8),pointer::  TS4_pft(:,:)
     real(r8),pointer::  TS5_pft(:,:)
     real(r8),pointer::  PondWater_pft(:,:)
     real(r8),pointer::  Rain_g_pft(:,:)
     real(r8),pointer::  f_soilwater(:)
     real(r8),pointer::  f_feileaf(:)
     real(r8),pointer::  f_soilwater_pft(:,:)
     real(r8),pointer::  f_feileaf_pft(:,:)
     real(r8),pointer::  f_Tleaf(:)
     real(r8),pointer::  f_Tleaf_pft(:,:)
     real(r8),pointer::  LHa(:) ! actual latent heat when f_feileaf works on Vcmax
     real(r8),pointer::  LHapft(:,:)
     ! add sunlit and shaded T and gpp
     real(r8),pointer::  GPP_o_sunlit_pft(:,:)
     real(r8),pointer::  GPP_o_shaded_pft(:,:)
     real(r8),pointer::  GPP_u_sunlit_pft(:,:)
     real(r8),pointer::  GPP_u_shaded_pft(:,:)

     real(r8),pointer::  GPP_o_sunlit(:)
     real(r8),pointer::  GPP_o_shaded(:)
     real(r8),pointer::  GPP_u_sunlit(:)
     real(r8),pointer::  GPP_u_shaded(:)

     real(r8),pointer::  TR_o_sunlit_pft(:,:)
     real(r8),pointer::  TR_o_shaded_pft(:,:)
     real(r8),pointer::  TR_u_sunlit_pft(:,:)
     real(r8),pointer::  TR_u_shaded_pft(:,:)

     real(r8),pointer::  TR_o_sunlit(:)
     real(r8),pointer::  TR_o_shaded(:)
     real(r8),pointer::  TR_u_sunlit(:)
     real(r8),pointer::  TR_u_shaded(:)

     real(r8),pointer::  Eil(:)
     real(r8),pointer::  Evap_soil(:)
     real(r8),pointer::  Evap_SW(:)
     real(r8),pointer::  EiS(:)
     real(r8),pointer::  Evap_SS(:)

     real(r8),pointer::  Eil_pft(:,:)
     real(r8),pointer::  Evap_soil_pft(:,:)
     real(r8),pointer::  Evap_SW_pft(:,:)
     real(r8),pointer::  EiS_pft(:,:)
     real(r8),pointer::  Evap_SS_pft(:,:)

  end type res

  type(res),save,target,public:: output

  type,public :: mc_res

     real(r8),pointer::  SMpft(:,:)
     real(r8),pointer::  SM(:,:)
     real(r8),pointer::  obsSM(:,:)
     ! tempSM for output SM in MC analysis
     real(r8),pointer::  tempSM(:)
     ! for leaf water potential analysis
     real(r8),pointer::  LWPpft(:,:)
     real(r8),pointer::  LWP(:,:)
     real(r8),pointer::  obsLWP(:,:)

  end type mc_res

  type(mc_res),save,target,public:: mc_output

end module bepstype
