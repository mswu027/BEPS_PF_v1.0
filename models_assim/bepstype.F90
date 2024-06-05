module bepstype
  use shr_kind_mod,only:r8=>shr_kind_r8
  !--iLab::no need for use-statement
  ! use beps_par
  implicit none

  ! For restart
  !-- iLab::module field should have save attribute
  real(r8), allocatable,public,save     ::  v2last(:,:,:)


  !*********************Climate Forcing datasets**********************
  type,public:: forc
     real(r8),pointer::  Temp(:)
     real(r8),pointer::  Tempmx(:)
     real(r8),pointer::  Tempmn(:)
     real(r8),pointer::  Wind(:)
#ifdef COUP_CSM
     real(r8),pointer::  Zref(:)   !bottom atm level height
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
     real(r8),pointer::  Srad(:)   !downward solar raditation
     real(r8),pointer::  Rh(:)     !relative humidity(%)
     real(r8),pointer::  Rain(:)   ! liquid precipitation
     real(r8),pointer::  Snow(:)   ! Snow rate
     real(r8),pointer::  Swdr(:)   ! SW direct radiation
     real(r8),pointer::  Swdf(:)   ! SW diffuse radiation
     character(len=10) :: meteo_ref_yyyymmdd = '' !--iLab::format 'yyyy-mm-dd'
                                                  !        met forcing time units is
                                                  !        "seconds since yyyy-mm-dd"
                                                  !        and we need to keep track on that
     integer :: meteo_ref_time(4)                 !-- iLab::reference time year/mon/day/hour
#endif
  end type forc

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
     !     real(r8),pointer:: p_Vcmax(:)
     !     real(r8),pointer:: p_q10(:)
     !     real(r8),pointer:: p_drainage(:)
     !     real(r8),pointer:: p_beta(:)
     !     real(r8),pointer:: p_Ksat(:)
     !     real(r8),pointer:: p_b(:)
     !     real(r8),pointer:: u_Vcmax(:)
     !     real(r8),pointer:: u_q10(:)
     !     real(r8),pointer:: u_drainage(:)
     !     real(r8),pointer:: u_beta(:)
     !     real(r8),pointer:: u_Ksat(:)
     !     real(r8),pointer:: u_b(:)

     real(r8),pointer:: lai(:,:)     ! for photosynthesis
     real(r8),pointer:: Vcmax(:,:)   ! for data assimilation
     !     character,pointer:: name(:)
     character(len=10) :: lai_ref_yyyymmdd = ''  !-- iLab::format 'yyyy-mm-dd'
                                                 !         expected LAI time unit is
                                                 !         "days since yyyy-mm-dd"
                                                 !         and we need to keep track on that
     integer :: lai_ref_time(4)                  !-- iLab::reference time year/mon/day/hour
  end type surf

  type(surf),save,target,public:: bound   ! boundary conditions

  !******************for assimilation and parameter optimization/
  type,public::para
     real(r8),pointer:: p_Vcmax(:)
     real(r8),pointer:: p_VJ_slope(:)
     real(r8),pointer:: p_q10(:)
     real(r8),pointer:: p_sif_alpha(:)
     real(r8),pointer:: p_sif_beta(:)
     real(r8),pointer:: p_taweff(:)
     real(r8),pointer:: p_D0(:)
     real(r8),pointer:: p_Ksat_scalar(:)
     real(r8),pointer:: p_b_scalar(:)
     real(r8),pointer:: p_f_leaf
     real(r8),pointer:: p_kc25
     real(r8),pointer:: p_ko25
     real(r8),pointer:: p_tau25
 !    real(r8),pointer:: p_f_lr
     real(r8),pointer:: p_agb2vod

     real(r8),pointer:: u_Vcmax(:)
     real(r8),pointer:: u_VJ_slope(:)
     real(r8),pointer:: u_q10(:)
     real(r8),pointer:: u_sif_alpha(:)
     real(r8),pointer:: u_sif_beta(:)
     real(r8),pointer:: u_taweff(:)
     real(r8),pointer:: u_D0(:)
     real(r8),pointer:: u_Ksat_scalar(:)
     real(r8),pointer:: u_b_scalar(:)
     real(r8),pointer:: u_f_leaf
     real(r8),pointer:: u_kc25
     real(r8),pointer:: u_ko25
     real(r8),pointer:: u_tau25
!     real(r8),pointer:: u_f_lr
     real(r8),pointer:: u_agb2vod

  end type para

  type(para),save,target,public:: assim   ! optimization of parameters


  !*******************************************************************
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
     real(r8),pointer::  Thetam(:)
     real(r8),pointer::  fAPARpft(:,:)
     real(r8),pointer::  fAPAR(:)
     real(r8),pointer::  VODpft(:,:)
     real(r8),pointer::  VOD(:)
     real(r8),pointer::  COS_fluxpft(:,:)
     real(r8),pointer::  COS_flux(:)
     real(r8),pointer::  NPP_yr_acc(:,:)
  end type res

  type(res),save,target,public:: output

end module bepstype
