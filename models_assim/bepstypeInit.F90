! This module is used to initialize bepstype variables
! Editted by J.Wang
! Date: 10May2017

module bepstypeInit
  use shr_kind_mod,only: r8=>shr_kind_r8
  !--iLab::added further entities (as module bepstype does no longer use beps_par)
  use beps_par,only:npoints,PFT,max_layers,texture,nparameters, mc_length,PF_np
  use bepstype
  use controlInput_mod
  implicit none
  !--iLab::can avoid 'save' here since no variables are declared,
  !        use-associated entities should already have the attribute
  ! save

  public  :: Initbepstype
#ifdef COUP_CSM
  private :: initlnd2atm
#endif
  private :: initatm2lnd
  private :: InitSurf
  private :: InitOuput
  private :: InitMCOuput
  private :: InitAssim
  private :: InitPF
  private :: InitMC
  public  :: InitPF_obs
  public  :: InitMC_obs
  private :: InitPF_resample
  private :: Initparam_gdd
contains

  subroutine Initbepstype()
    implicit none

    allocate(v2last(npoints,0:40,PFT))
    v2last  = 0.

    call initatm2lnd()
#ifdef COUP_CSM
    call initlnd2atm()
#endif

    call InitSurf()

    call InitSoilstat()

    call InitOuput()

    call InitMCOuput()

    call InitAssim()

    call InitPF()

    call InitMC()

    call InitPF_obs()

    call InitMC_obs()

    call InitPF_resample()

    call Initparam_gdd()

    return
  end subroutine Initbepstype


  subroutine initatm2lnd()
    implicit none
    !--iLab::avoid pointer
    ! type(forc),pointer ::p
    ! p=>clim

    allocate(clim%Temp(npoints))
    allocate(clim%Tempmx(npoints))
    allocate(clim%Tempmn(npoints))
    allocate(clim%Wind(npoints))

#ifdef COUP_CSM
    allocate(clim%Zref(npoints))
    allocate(clim%Rain(npoints))
    allocate(clim%Snow(npoints))
    allocate(clim%Swndr(npoints))
    allocate(clim%Swvdr(npoints))
    allocate(clim%Swndf(npoints))
    allocate(clim%Swvdf(npoints))
    allocate(clim%Swdr(npoints))
    allocate(clim%Swdf(npoints))
    allocate(clim%Lwdn(npoints))
    allocate(clim%shum(npoints))
    allocate(clim%pres(npoints))
#else
    allocate(clim%Srad(npoints))
    allocate(clim%Rh(npoints))
    allocate(clim%Rain(npoints))
    allocate(clim%Snow(npoints))
    allocate(clim%Swdr(npoints))
    allocate(clim%Swdf(npoints))
#endif

    clim%Temp(:)      = 0.
    clim%Tempmx(:)      = 0.
    clim%Tempmn(:)      = 0.
    clim%Wind(:)      = 0.
#ifdef COUP_CSM
    clim%Zref(:)     = 0.
    clim%Rain(:)     = 0.
    clim%Snow(:)     = 0.
    clim%Swndr(:)    = 0.
    clim%Swvdr(:)    = 0.
    clim%Swndf(:)    = 0.
    clim%Swvdf(:)    = 0.
    clim%Swdr(:)     = 0.
    clim%Swdf(:)     = 0.
    clim%Lwdn(:)     = 0.
    clim%shum(:)     = 0.
    clim%pres(:)     = 0.
#else
    clim%Srad(:)     = 0.
    clim%Rh(:)       = 0.
    clim%Rain(:)     = 0.
    clim%Snow(:)     = 0.
    clim%Swdr(:)     = 0.
    clim%Swdf(:)     = 0.
#endif

    return
  end subroutine initatm2lnd

#ifdef COUP_CSM

  subroutine initlnd2atm()
    implicit none

    !--iLab::avoid pointer
    ! type(CPL),pointer :: p
    ! p=>lnd2atm

    allocate(lnd2atm%rofliq(npoints))
    allocate(lnd2atm%rofice(npoints))
    allocate(lnd2atm%t_Rad(npoints))
    allocate(lnd2atm%tref(npoints))
    allocate(lnd2atm%qref(npoints))
    allocate(lnd2atm%avsdr(npoints))
    allocate(lnd2atm%anidr(npoints))
    allocate(lnd2atm%avsdf(npoints))
    allocate(lnd2atm%anidf(npoints))
    allocate(lnd2atm%snowh(npoints))
    allocate(lnd2atm%u10(npoints))
    allocate(lnd2atm%ddvel(npoints))
    allocate(lnd2atm%fv(npoints))
    allocate(lnd2atm%ram1(npoints))
    allocate(lnd2atm%soilw(npoints))
    allocate(lnd2atm%taux(npoints))
    allocate(lnd2atm%tauy(npoints))
    allocate(lnd2atm%LH(npoints))
    allocate(lnd2atm%SH(npoints))
    allocate(lnd2atm%lwup(npoints))
    allocate(lnd2atm%evap(npoints))
    allocate(lnd2atm%swnet(npoints))
    allocate(lnd2atm%fco2(npoints))
    allocate(lnd2atm%flxdst1(npoints))
    allocate(lnd2atm%flxdst2(npoints))
    allocate(lnd2atm%flxdst3(npoints))
    allocate(lnd2atm%flxdst4(npoints))
    allocate(lnd2atm%flxvoc(npoints))

    lnd2atm%rofliq(:)     = 0.
    lnd2atm%rofice(:)     = 0.
    lnd2atm%t_Rad(:)      = 0.
    lnd2atm%tref(:)       = 0.
    lnd2atm%qref(:)       = 0.
    lnd2atm%avsdr(:)      = 0.
    lnd2atm%anidr(:)      = 0.
    lnd2atm%avsdf(:)      = 0.
    lnd2atm%anidf(:)      = 0.
    lnd2atm%snowh(:)      = 0.
    lnd2atm%u10(:)        = 0.
    lnd2atm%ddvel(:)      = 0.
    lnd2atm%fv(:)         = 0.
    lnd2atm%ram1(:)       = 0.
    lnd2atm%soilw(:)      = 0.
    lnd2atm%taux(:)       = 0.
    lnd2atm%tauy(:)       = 0.
    lnd2atm%LH(:)         = 0.
    lnd2atm%SH(:)         = 0.
    lnd2atm%lwup(:)       = 0.
    lnd2atm%evap(:)       = 0.
    lnd2atm%swnet(:)      = 0.
    lnd2atm%fco2(:)       = 0.
    lnd2atm%flxdst1(:)    = 0.
    lnd2atm%flxdst2(:)    = 0.
    lnd2atm%flxdst3(:)    = 0.
    lnd2atm%flxdst4(:)    = 0.
    lnd2atm%flxvoc(:)     = 0.

  end subroutine initlnd2atm

#endif

  subroutine InitSurf()
    implicit none
    !--iLab::avoid pointer
    ! type(surf),pointer::p
    ! p=>bound

    allocate(bound%lcno(npoints,PFT))
    allocate(bound%stext(npoints))
    allocate(bound%PCT_PFT(npoints,PFT))
    allocate(bound%clumping(npoints))
    allocate(bound%longitude(npoints))
    allocate(bound%latitude(npoints))

    allocate(bound%sdp(npoints))
    allocate(bound%st(npoints))
    allocate(bound%sw(npoints))

    allocate(bound%laiyr(npoints,PFT))
    allocate(bound%nppyr(npoints,PFT))
    allocate(bound%ccd(npoints,PFT))
    allocate(bound%cfmd(npoints,PFT))
    allocate(bound%cfsd(npoints,PFT))
    allocate(bound%cm(npoints,PFT))
    allocate(bound%cp(npoints,PFT))
    allocate(bound%cs(npoints,PFT))
    allocate(bound%csm(npoints,PFT))
    allocate(bound%csmd(npoints,PFT))
    allocate(bound%cssd(npoints,PFT))
    allocate(bound%lai(npoints,PFT))
    allocate(bound%Vcmax(npoints,PFT))
    ! 2024/03/29
    allocate(bound%HeightC(npoints))

    bound%lcno      = 0
    bound%stext     = 0
    bound%PCT_PFT   = 0
    bound%clumping  = 0.
    bound%longitude = 0.
    bound%latitude  = 0.

    bound%sdp       = 0.
    bound%st        = 0.
    bound%sw        = 0.

    bound%laiyr     = 0.
    bound%nppyr     = 0.
    bound%ccd       = 0.
    bound%cfmd      = 0.
    bound%cfsd      = 0.
    bound%cm        = 0.
    bound%cp        = 0.
    bound%cs        = 0.
    bound%csm       = 0.
    bound%csmd      = 0.
    bound%cssd      = 0.

    bound%lai       = 0.
    bound%Vcmax     = 0.
    !2024/03/29
    bound%HeightC   = 0.

  end subroutine InitSurf

  subroutine InitAssim()
    implicit none
    !--iLab::avoid pointer
    ! type(para),pointer  ::p
    ! p => assim

! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
    allocate(assim%p_Vcmax(PFT,npoints))
    allocate(assim%p_VJ_slope(PFT,npoints))
    allocate(assim%p_VN_slope(PFT,npoints))
    allocate(assim%p_b_h2o(PFT,npoints))
    allocate(assim%p_m_h2o(PFT,npoints))
    allocate(assim%p_f_leaf(PFT,npoints))
    allocate(assim%p_kc25(PFT,npoints))
    allocate(assim%p_ko25(PFT,npoints))
    allocate(assim%p_tau25(PFT,npoints))
 ! parameters for SIF modelling
    allocate(assim%p_sif_alpha(PFT,npoints))
    allocate(assim%p_sif_beta(PFT,npoints))
! parameters related to plant respiration: plant_resp
    allocate(assim%p_q10(PFT,npoints))
! parameters related to soil respiration: soil_resp
    allocate(assim%p_f_resp(PFT,npoints))
! soil parameters
    allocate(assim%p_f_decay(PFT,npoints))
    allocate(assim%p_Ksat_scalar(texture,npoints))
    allocate(assim%p_b_scalar(texture,npoints))
    allocate(assim%p_porosity_scalar(texture,npoints))
    allocate(assim%p_vfc_scalar(texture,npoints))
    allocate(assim%p_vwp_scalar(texture,npoints))
    allocate(assim%p_psisat_scalar(texture,npoints))
    allocate(assim%p_drainage_scalar(PFT,npoints))
 ! parameters for vod modelling
    allocate(assim%p_vod_a(PFT,npoints))
    allocate(assim%p_vod_b(PFT,npoints))
    allocate(assim%p_vod_c(PFT,npoints))
! parameters for SPAC in soilwateruptake
    allocate(assim%p_theta_Amin(PFT,npoints))
    allocate(assim%p_pox(PFT,npoints))
    allocate(assim%p_fei_c(PFT,npoints))
    allocate(assim%p_spac_p1(PFT,npoints))
    allocate(assim%p_spac_p2(PFT,npoints))
    allocate(assim%p_tWA(PFT,npoints))
    allocate(assim%p_tWB(PFT,npoints))
    allocate(assim%p_Ttrig(PFT,npoints))
    allocate(assim%p_r_xylem(PFT,npoints))
    allocate(assim%p_r_r(PFT,npoints))
    allocate(assim%p_Lr(PFT,npoints))
    allocate(assim%p_deltal_min(PFT,npoints))
    allocate(assim%p_deltal_max(PFT,npoints))
    allocate(assim%p_p_delta(PFT,npoints))
    allocate(assim%p_ppslh(PFT,npoints))
    allocate(assim%p_fei_min(PFT,npoints))
    allocate(assim%p_fei_th(PFT,npoints))
    allocate(assim%p_p_excess(PFT,npoints))
    allocate(assim%p_Tleaf_H(PFT,npoints))
    allocate(assim%p_Tleaf_L(PFT,npoints))
    allocate(assim%p_Tleaf_O(PFT,npoints))

! parameters' uncertainty
! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
    allocate(assim%u_Vcmax(PFT,npoints))
    allocate(assim%u_VJ_slope(PFT,npoints))
    allocate(assim%u_VN_slope(PFT,npoints))
    allocate(assim%u_b_h2o(PFT,npoints))
    allocate(assim%u_m_h2o(PFT,npoints))
    allocate(assim%u_f_leaf(PFT,npoints))
    allocate(assim%u_kc25(PFT,npoints))
    allocate(assim%u_ko25(PFT,npoints))
    allocate(assim%u_tau25(PFT,npoints))
 ! parameters for SIF modelling
    allocate(assim%u_sif_alpha(PFT,npoints))
    allocate(assim%u_sif_beta(PFT,npoints))
! parameters related to plant respiration: plant_resp
    allocate(assim%u_q10(PFT,npoints))
! parameters related to soil respiration: soil_resp
    allocate(assim%u_f_resp(PFT,npoints))
! soil parameters
    allocate(assim%u_f_decay(PFT,npoints))
    allocate(assim%u_Ksat_scalar(texture,npoints))
    allocate(assim%u_b_scalar(texture,npoints))
    allocate(assim%u_porosity_scalar(texture,npoints))
    allocate(assim%u_vfc_scalar(texture,npoints))
    allocate(assim%u_vwp_scalar(texture,npoints))
    allocate(assim%u_psisat_scalar(texture,npoints))
    allocate(assim%u_drainage_scalar(PFT,npoints))
 ! parameters for vod modelling
    allocate(assim%u_vod_a(PFT,npoints))
    allocate(assim%u_vod_b(PFT,npoints))
    allocate(assim%u_vod_c(PFT,npoints))
! parameters for SPAC in soilwateruptake
    allocate(assim%u_theta_Amin(PFT,npoints))
    allocate(assim%u_pox(PFT,npoints))
    allocate(assim%u_fei_c(PFT,npoints))
    allocate(assim%u_spac_p1(PFT,npoints))
    allocate(assim%u_spac_p2(PFT,npoints))
    allocate(assim%u_tWA(PFT,npoints))
    allocate(assim%u_tWB(PFT,npoints))
    allocate(assim%u_Ttrig(PFT,npoints))
    allocate(assim%u_r_xylem(PFT,npoints))
    allocate(assim%u_r_r(PFT,npoints))
    allocate(assim%u_Lr(PFT,npoints))
    allocate(assim%u_deltal_min(PFT,npoints))
    allocate(assim%u_deltal_max(PFT,npoints))
    allocate(assim%u_p_delta(PFT,npoints))
    allocate(assim%u_ppslh(PFT,npoints))
    allocate(assim%u_fei_min(PFT,npoints))
    allocate(assim%u_fei_th(PFT,npoints))
    allocate(assim%u_p_excess(PFT,npoints))
    allocate(assim%u_Tleaf_H(PFT,npoints))
    allocate(assim%u_Tleaf_L(PFT,npoints))
    allocate(assim%u_Tleaf_O(PFT,npoints))

    assim%p_Vcmax   = 0.
    assim%p_VJ_slope   = 0.
    assim%p_VN_slope   = 0.
    assim%p_b_h2o=0.
    assim%p_m_h2o=0.
    assim%p_f_leaf   = 0.
    assim%p_kc25   = 0.
    assim%p_ko25   = 0.
    assim%p_tau25   = 0.

    assim%p_sif_alpha   = 0.
    assim%p_sif_beta   = 0.

    assim%p_q10     = 0.
    assim%p_f_resp   = 0.

    assim%p_f_decay   = 0.
    assim%p_Ksat_scalar   = 0.
    assim%p_b_scalar  = 0.
    assim%p_porosity_scalar  = 0.
    assim%p_vfc_scalar  = 0.
    assim%p_vwp_scalar  = 0.
    assim%p_psisat_scalar  = 0.
    assim%p_drainage_scalar  = 0.

    assim%p_vod_a   = 0.
    assim%p_vod_b   = 0.
    assim%p_vod_c   = 0.

    assim%p_theta_Amin = 0.
    assim%p_pox = 0.
    assim%p_fei_c = 0.
    assim%p_spac_p1 = 0.
    assim%p_spac_p2 = 0.
    assim%p_tWA = 0.
    assim%p_tWB = 0.
    assim%p_Ttrig = 0.
    assim%p_r_xylem = 0.
    assim%p_r_r = 0.
    assim%p_Lr = 0.
    assim%p_deltal_min = 0.
    assim%p_deltal_max = 0.
    assim%p_p_delta = 0.
    assim%p_ppslh = 0.
    assim%p_fei_min = 0.
    assim%p_fei_th = 0.
    assim%p_p_excess = 0.
    assim%p_Tleaf_H = 0.
    assim%p_Tleaf_L = 0.
    assim%p_Tleaf_O = 0.

    assim%u_Vcmax   = 0.
    assim%u_VJ_slope   = 0.
    assim%u_VN_slope   = 0.
    assim%u_b_h2o=0.
    assim%u_m_h2o=0.
    assim%u_f_leaf   = 0.
    assim%u_kc25   = 0.
    assim%u_ko25   = 0.
    assim%u_tau25   = 0.

    assim%u_sif_alpha   = 0.
    assim%u_sif_beta   = 0.

    assim%u_q10     = 0.
    assim%u_f_resp   = 0.

    assim%u_f_decay   = 0.
    assim%u_Ksat_scalar   = 0.
    assim%u_b_scalar  = 0.
    assim%u_porosity_scalar  = 0.
    assim%u_vfc_scalar  = 0.
    assim%u_vwp_scalar  = 0.
    assim%u_psisat_scalar  = 0.
    assim%u_drainage_scalar  = 0.

    assim%u_vod_a   = 0.
    assim%u_vod_b   = 0.
    assim%u_vod_c   = 0.

    assim%u_theta_Amin = 0.
    assim%u_pox = 0.
    assim%u_fei_c = 0.
    assim%u_spac_p1 = 0.
    assim%u_spac_p2 = 0.
    assim%u_tWA = 0.
    assim%u_tWB = 0.
    assim%u_Ttrig = 0.
    assim%u_r_xylem = 0.
    assim%u_r_r = 0.
    assim%u_Lr = 0.
    assim%u_deltal_min = 0.
    assim%u_deltal_max = 0.
    assim%u_p_delta = 0.
    assim%u_ppslh = 0.
    assim%u_fei_min = 0.
    assim%u_fei_th = 0.
    assim%u_p_excess = 0.
    assim%u_Tleaf_H = 0.
    assim%u_Tleaf_L = 0.
    assim%u_Tleaf_O = 0.

  end subroutine InitAssim

  subroutine InitPF()
    implicit none
    !--iLab::avoid pointer
    ! type(para),pointer  ::p
    ! p => assim

    allocate(PF%Vcmax(parloop,npoints))
    allocate(PF%VJ_slope(parloop,npoints))
    allocate(PF%VN_slope(parloop,npoints))
    allocate(PF%b_h2o(parloop,npoints))
    allocate(PF%m_h2o(parloop,npoints))
    allocate(PF%f_leaf(parloop,npoints))
    allocate(PF%kc25(parloop,npoints))
    allocate(PF%ko25(parloop,npoints))
    allocate(PF%tau25(parloop,npoints))

    allocate(PF%sif_alpha(parloop,npoints))
    allocate(PF%sif_beta(parloop,npoints))

    allocate(PF%q10(parloop,npoints))
    allocate(PF%f_resp(parloop,npoints))

    allocate(PF%f_decay(parloop,npoints))
    allocate(PF%Ksat_scalar(parloop,npoints))
    allocate(PF%b_scalar(parloop,npoints))
    allocate(PF%porosity_scalar(parloop,npoints))
    allocate(PF%vfc_scalar(parloop,npoints))
    allocate(PF%vwp_scalar(parloop,npoints))
    allocate(PF%psisat_scalar(parloop,npoints))
    allocate(PF%drainage_scalar(parloop,npoints))

    allocate(PF%vod_a(parloop,npoints))
    allocate(PF%vod_b(parloop,npoints))
    allocate(PF%vod_c(parloop,npoints))

    allocate(PF%theta_Amin(parloop,npoints))
    allocate(PF%pox(parloop,npoints))
    allocate(PF%fei_c(parloop,npoints))
    allocate(PF%spac_p1(parloop,npoints))
    allocate(PF%spac_p2(parloop,npoints))
    allocate(PF%tWA(parloop,npoints))
    allocate(PF%tWB(parloop,npoints))
    allocate(PF%Ttrig(parloop,npoints))
    allocate(PF%r_xylem(parloop,npoints))
    allocate(PF%r_r(parloop,npoints))
    allocate(PF%Lr(parloop,npoints))
    allocate(PF%deltal_min(parloop,npoints))
    allocate(PF%deltal_max(parloop,npoints))
    allocate(PF%p_delta(parloop,npoints))
    allocate(PF%ppslh(parloop,npoints))
    allocate(PF%fei_min(parloop,npoints))
    allocate(PF%fei_th(parloop,npoints))
    allocate(PF%p_excess(parloop,npoints))
    allocate(PF%Tleaf_H(parloop,npoints))
    allocate(PF%Tleaf_L(parloop,npoints))
    allocate(PF%Tleaf_O(parloop,npoints))

    allocate(PF%pfweight(parloop,npoints))
    allocate(PF%pfweightupdate(parloop,npoints))

    PF%Vcmax   = 0.
    PF%VJ_slope   = 0.
    PF%VN_slope = 0.
    PF%b_h2o=0.
    PF%m_h2o=0.
    PF%f_leaf   = 0.
    PF%kc25   = 0.
    PF%ko25   = 0.
    PF%tau25   = 0.

    PF%sif_alpha   = 0.
    PF%sif_beta   = 0.

    PF%q10     = 0.
    PF%f_resp = 0.

    PF%f_decay = 0.
    PF%Ksat_scalar   = 0.
    PF%b_scalar  = 0.
    PF%porosity_scalar = 0.
    PF%vfc_scalar = 0.
    PF%vwp_scalar = 0.
    PF%psisat_scalar = 0.
    PF%drainage_scalar = 0.

    PF%vod_a = 0.
    PF%vod_b = 0.
    PF%vod_c = 0.

    PF%theta_Amin = 0.
    PF%pox = 0.
    PF%fei_c = 0.
    PF%spac_p1 = 0.
    PF%spac_p2 = 0.
    PF%tWA = 0.
    PF%tWB = 0.
    PF%Ttrig = 0.
    PF%r_xylem = 0.
    PF%r_r = 0.
    PF%Lr = 0.
    PF%deltal_min = 0.
    PF%deltal_max = 0.
    PF%p_delta = 0.
    PF%ppslh = 0.
    PF%fei_min = 0.
    PF%fei_th = 0.
    PF%p_excess = 0.
    PF%Tleaf_H = 0.
    PF%Tleaf_L = 0.
    PF%Tleaf_O = 0.

    PF%pfweight   = 0.
    PF%pfweightupdate   = 0.

  end subroutine InitPF

  subroutine InitMC()
    implicit none
    !--iLab::avoid pointer
    ! type(para),pointer  ::p
    ! p => assim
   allocate(MC%Vcmax(nparameters,npoints))
    allocate(MC%VJ_slope(nparameters,npoints))
    allocate(MC%VN_slope(nparameters,npoints))
    allocate(MC%b_h2o(nparameters,npoints))
    allocate(MC%m_h2o(nparameters,npoints))
    allocate(MC%f_leaf(nparameters,npoints))
    allocate(MC%kc25(nparameters,npoints))
    allocate(MC%ko25(nparameters,npoints))
    allocate(MC%tau25(nparameters,npoints))

    allocate(MC%sif_alpha(nparameters,npoints))
    allocate(MC%sif_beta(nparameters,npoints))

    allocate(MC%q10(nparameters,npoints))
    allocate(MC%f_resp(nparameters,npoints))

    allocate(MC%f_decay(nparameters,npoints))
    allocate(MC%Ksat_scalar(nparameters,npoints))
    allocate(MC%b_scalar(nparameters,npoints))
    allocate(MC%porosity_scalar(nparameters,npoints))
    allocate(MC%vfc_scalar(nparameters,npoints))
    allocate(MC%vwp_scalar(nparameters,npoints))
    allocate(MC%psisat_scalar(nparameters,npoints))
    allocate(MC%drainage_scalar(nparameters,npoints))

    allocate(MC%vod_a(nparameters,npoints))
    allocate(MC%vod_b(nparameters,npoints))
    allocate(MC%vod_c(nparameters,npoints))

    allocate(MC%theta_Amin(nparameters,npoints))
    allocate(MC%pox(nparameters,npoints))
    allocate(MC%fei_c(nparameters,npoints))
    allocate(MC%spac_p1(nparameters,npoints))
    allocate(MC%spac_p2(nparameters,npoints))
    allocate(MC%tWA(nparameters,npoints))
    allocate(MC%tWB(nparameters,npoints))
    allocate(MC%Ttrig(nparameters,npoints))
    allocate(MC%r_xylem(nparameters,npoints))
    allocate(MC%r_r(nparameters,npoints))
    allocate(MC%Lr(nparameters,npoints))
    allocate(MC%deltal_min(nparameters,npoints))
    allocate(MC%deltal_max(nparameters,npoints))
    allocate(MC%p_delta(nparameters,npoints))
    allocate(MC%ppslh(nparameters,npoints))
    allocate(MC%fei_min(nparameters,npoints))
    allocate(MC%fei_th(nparameters,npoints))
    allocate(MC%p_excess(nparameters,npoints))
    allocate(MC%Tleaf_H(nparameters,npoints))
    allocate(MC%Tleaf_L(nparameters,npoints))
    allocate(MC%Tleaf_O(nparameters,npoints))

    allocate(MC%SM_R(nparameters,npoints))
    allocate(MC%SM_R2(nparameters,npoints))
    allocate(MC%SM_RMSE(nparameters,npoints))
    allocate(MC%SM_ME(nparameters,npoints))
    allocate(MC%SM_AIC(nparameters,npoints))
    allocate(MC%SM_BIC(nparameters,npoints))

    allocate(MC%LWP_R(nparameters,npoints))
    allocate(MC%LWP_R2(nparameters,npoints))
    allocate(MC%LWP_RMSE(nparameters,npoints))
    allocate(MC%LWP_ME(nparameters,npoints))
    allocate(MC%LWP_AIC(nparameters,npoints))
    allocate(MC%LWP_BIC(nparameters,npoints))

    MC%Vcmax   = 0.
    MC%VJ_slope   = 0.
    MC%VN_slope = 0.
    MC%b_h2o=0.
    MC%m_h2o=0.
    MC%f_leaf   = 0.
    MC%kc25   = 0.
    MC%ko25   = 0.
    MC%tau25   = 0.

    MC%sif_alpha   = 0.
    MC%sif_beta   = 0.

    MC%q10     = 0.
    MC%f_resp = 0.

    MC%f_decay = 0.
    MC%Ksat_scalar   = 0.
    MC%b_scalar  = 0.
    MC%porosity_scalar = 0.
    MC%vfc_scalar = 0.
    MC%vwp_scalar = 0.
    MC%psisat_scalar = 0.
    MC%drainage_scalar = 0.

    MC%vod_a = 0.
    MC%vod_b = 0.
    MC%vod_c = 0.

    MC%theta_Amin = 0.
    MC%pox = 0.
    MC%fei_c = 0.
    MC%spac_p1 = 0.
    MC%spac_p2 = 0.
    MC%tWA = 0.
    MC%tWB = 0.
    MC%Ttrig = 0.
    MC%r_xylem = 0.
    MC%r_r = 0.
    MC%Lr = 0.
    MC%deltal_min = 0.
    MC%deltal_max = 0.
    MC%p_delta = 0.
    MC%ppslh = 0.
    MC%fei_min = 0.
    MC%fei_th = 0.
    MC%p_excess = 0.
    MC%Tleaf_H = 0.
    MC%Tleaf_L = 0.
    MC%Tleaf_O = 0.

    ! performance metrics
    MC%SM_R = 0.
    MC%SM_R2 = 0.
    MC%SM_RMSE = 0.
    MC%SM_ME = 0.
    MC%SM_AIC = 0.
    MC%SM_BIC = 0.

    MC%LWP_R = 0.
    MC%LWP_R2 = 0.
    MC%LWP_RMSE = 0.
    MC%LWP_ME = 0.
    MC%LWP_AIC = 0.
    MC%LWP_BIC = 0.

  end subroutine InitMC

  subroutine InitPF_obs()
    implicit none

    ! allocate(PF_obs%obs_GPP(npoints))
    allocate(PF_obs%obs_var(npoints)) ! changed to VOD from obs_GPP

     ! PF_obs%obs_GPP   = 0.
     PF_obs%obs_var   = 0.

  end subroutine InitPF_obs

subroutine InitMC_obs()
    implicit none

    !allocate(MC_obs%obs_VOD(npoints)) ! changed to VOD from obs_GPP
    allocate(MC_obs%obs_SM(npoints))
    allocate(MC_obs%obs_LWP(npoints))

     MC_obs%obs_SM   = 0.
     MC_obs%obs_LWP = 0.

  end subroutine InitMC_obs

  subroutine InitPF_resample()
    implicit none
    ! only consider 5 soil hydraulic parameters
    allocate(PF_resample%outparticles(parloop,PF_np+1)) ! 45 pamameters and one for pf_weight
    allocate(PF_resample%resample_weight(parloop))
    allocate(PF_resample%resample_weight_update(parloop))

    PF_resample%outparticles   = 0.
    PF_resample%resample_weight=0.
    PF_resample%resample_weight_update=0.

  end subroutine InitPF_resample

  subroutine InitSoilstat()
    implicit none
    !--iLab::avoid pointer
    ! type(soils),pointer  :: p
    ! p => soilstat

    allocate(soilstat%n_layer(npoints))
    allocate(soilstat%Zp(npoints,PFT))
    allocate(soilstat%Zsp(npoints,PFT))
    allocate(soilstat%r_rain_g(npoints,PFT))
    allocate(soilstat%r_drainage(npoints,PFT))
    allocate(soilstat%r_root_decay(npoints,PFT))
    allocate(soilstat%psi_min(npoints,PFT))
    allocate(soilstat%alpha(npoints,PFT))
    allocate(soilstat%f_soilwater(npoints,PFT))
    allocate(soilstat%d_soil(npoints,0:max_layers-1))
    allocate(soilstat%f_root(npoints,0:max_layers-1,PFT))
    allocate(soilstat%dt(npoints,0:max_layers-1,PFT))
    allocate(soilstat%thermal_cond(npoints,0:max_layers-1,PFT))
    allocate(soilstat%theta_vfc(npoints,0:max_layers-1,PFT))
    allocate(soilstat%theta_vwp(npoints,0:max_layers-1,PFT))
    allocate(soilstat%fei(npoints,0:max_layers-1,PFT))
    allocate(soilstat%Ksat(npoints,0:max_layers-1,PFT))
    allocate(soilstat%psi_sat(npoints,0:max_layers-1,PFT))
    allocate(soilstat%b(npoints,0:max_layers-1,PFT))
    allocate(soilstat%density_soil(npoints,0:max_layers-1))
    allocate(soilstat%f_org(npoints,0:max_layers-1,PFT))
    allocate(soilstat%ice_ratio(npoints,0:max_layers-1,PFT))
    allocate(soilstat%thetam(npoints,0:max_layers-1,PFT))
    allocate(soilstat%thetam_prev(npoints,0:max_layers-1,PFT))
    allocate(soilstat%temp_soil_p(npoints,0:max_layers-1,PFT))
    allocate(soilstat%temp_soil_c(npoints,0:max_layers-1,PFT))
    allocate(soilstat%f_ice(npoints,0:max_layers-1,PFT))
    allocate(soilstat%psim(npoints,0:max_layers-1,PFT))
    allocate(soilstat%thetab(npoints,0:max_layers-1,PFT))
    allocate(soilstat%psib(npoints,0:max_layers-1,PFT))
    allocate(soilstat%r_waterflow(npoints,0:max_layers-1,PFT))
    allocate(soilstat%km(npoints,0:max_layers-1,PFT))
    allocate(soilstat%kb(npoints,0:max_layers-1,PFT))
    allocate(soilstat%KK(npoints,0:max_layers-1,PFT))
    allocate(soilstat%Cs(npoints,0:max_layers-1,PFT))
    allocate(soilstat%lambda(npoints,0:max_layers-1,PFT))
    allocate(soilstat%Ett(npoints,0:max_layers-1,PFT))
    allocate(soilstat%G(npoints,0:max_layers-1,PFT))
    ! 2023/06/30
    allocate(soilstat%Sp(npoints,PFT))

    soilstat%n_layer(:)           = 0
    soilstat%Zp(:,:)              = 0.
    soilstat%Zsp(:,:)             = 0.
    soilstat%r_rain_g(:,:)        = 0.
    soilstat%r_drainage(:,:)      = 0.
    soilstat%r_root_decay(:,:)    = 0.
    soilstat%psi_min(:,:)         = 0.
    soilstat%alpha(:,:)           = 0.
    soilstat%f_soilwater(:,:)     = 0.
    ! 2023/06/30
    soilstat%Sp(:,:)              = 0.

    soilstat%d_soil(:,:)          = 0.
    soilstat%f_root(:,:,:)        = 0.
    soilstat%dt(:,:,:)            = 0.
    soilstat%thermal_cond(:,:,:)  = 0.
    soilstat%theta_vfc(:,:,:)     = 0.
    soilstat%theta_vwp(:,:,:)     = 0.
    soilstat%fei(:,:,:)           = 0.
    soilstat%Ksat(:,:,:)          = 0.
    soilstat%psi_sat(:,:,:)       = 0.
    soilstat%b(:,:,:)             = 0.
    soilstat%density_soil(:,:)    = 0.
    soilstat%f_org(:,:,:)         = 0.
    soilstat%ice_ratio(:,:,:)     = 0.
    soilstat%thetam(:,:,:)        = 0.
    soilstat%thetam_prev(:,:,:)   = 0.
    soilstat%temp_soil_p(:,:,:)   = 0.
    soilstat%temp_soil_c(:,:,:)   = 0.
    soilstat%f_ice(:,:,:)         = 0.
    soilstat%psim(:,:,:)          = 0.
    soilstat%thetab(:,:,:)        = 0.
    soilstat%psib(:,:,:)           = 0.
    soilstat%r_waterflow(:,:,:)   = 0.
    soilstat%km(:,:,:)            = 0.
    soilstat%kb(:,:,:)            = 0.
    soilstat%KK(:,:,:)            = 0.
    soilstat%Cs(:,:,:)            = 0.
    soilstat%lambda(:,:,:)        = 0.
    soilstat%Ett(:,:,:)           = 0.
    soilstat%G(:,:,:)             = 0.

  end subroutine InitSoilstat
  !******************************************for C4 crop Xiuli**************************

  subroutine Initparam_gdd()
    implicit none
    allocate(pgdd%tt_veg(npoints,PFT))
    allocate(pgdd%tt_rep(npoints,PFT))
    allocate(pgdd%phot_type(npoints,PFT))
    allocate(pgdd%emer_doy(npoints,PFT))
    allocate(pgdd%har_doy(npoints,PFT))

    pgdd%tt_veg(:,:)      = 0.
    pgdd%tt_rep(:,:)     = 0.
    pgdd%phot_type(:,:)   = 0.
    pgdd%emer_doy(:,:)  = 0.
    pgdd%har_doy(:,:) = 0.

  end subroutine Initparam_gdd

  !******************************************for C4 crop Xiuli**************************
  subroutine InitOuput()
    implicit none
    !--iLab::avoid pointer
    ! type(res),pointer::p
    ! p=>output

    allocate(output%GPPpft(npoints,PFT))
    allocate(output%SIFpft(npoints,PFT))
    allocate(output%SIFpft_sat(npoints,PFT))
    allocate(output%NPPpft(npoints,PFT))
    allocate(output%NEPpft(npoints,PFT))
    allocate(output%SHpft(npoints,PFT))
    allocate(output%LHpft(npoints,PFT))
    allocate(output%Transpft(npoints,PFT))
    allocate(output%Evappft(npoints,PFT))
    allocate(output%Net_Radpft(npoints,PFT))
    allocate(output%GPP(npoints))
    allocate(output%SIF(npoints))
    allocate(output%SIF_sat(npoints))
    allocate(output%NPP(npoints))
    allocate(output%NEP(npoints))
    allocate(output%LAIpft(npoints,PFT))
    allocate(output%LAI(npoints))
    allocate(output%SH(npoints))
    allocate(output%LH(npoints))
    allocate(output%Trans(npoints))
    allocate(output%Evap(npoints))
    allocate(output%Net_Rad(npoints))
    allocate(output%Thetampft(npoints,PFT))
    allocate(output%Thetam(npoints))
    allocate(output%fAPARpft(npoints,PFT))
    allocate(output%fAPAR(npoints))
    allocate(output%VODpft(npoints,PFT))  !!! why?
    allocate(output%VOD(npoints))
    allocate(output%COS_fluxpft(npoints,PFT))
    allocate(output%COS_flux(npoints))
    allocate(output%NPP_yr_acc(npoints,PFT))

    ! 2023/06/30
    allocate(output%PWS(npoints))
    allocate(output%ETa(npoints))
    allocate(output%Qupt(npoints))
    allocate(output%PWSpft(npoints,pft))
    allocate(output%ETapft(npoints,pft))
    allocate(output%Quptpft(npoints,pft))
    allocate(output%fei_leaf(npoints))
    allocate(output%fei_leafpft(npoints,PFT))
    ! 2023/10/30
    allocate(output%Thetam2_pft(npoints,PFT))
    allocate(output%Thetam3_pft(npoints,PFT))
    allocate(output%Thetam4_pft(npoints,PFT))
    allocate(output%Thetam5_pft(npoints,PFT))
    allocate(output%Thetam_layer2(npoints))
    allocate(output%Thetam_layer3(npoints))
    allocate(output%Thetam_layer4(npoints))
    allocate(output%Thetam_layer5(npoints))

    allocate(output%SWP1_pft(npoints,PFT))
    allocate(output%SWP2_pft(npoints,PFT))
    allocate(output%SWP3_pft(npoints,PFT))
    allocate(output%SWP4_pft(npoints,PFT))
    allocate(output%SWP5_pft(npoints,PFT))
    allocate(output%SWP_layer1(npoints))
    allocate(output%SWP_layer2(npoints))
    allocate(output%SWP_layer3(npoints))
    allocate(output%SWP_layer4(npoints))
    allocate(output%SWP_layer5(npoints))
    allocate(output%TS_layer1(npoints))
    allocate(output%TS_layer2(npoints))
    allocate(output%TS_layer3(npoints))
    allocate(output%TS_layer4(npoints))
    allocate(output%TS_layer5(npoints))
    allocate(output%PondWater(npoints))
    allocate(output%Rain_g(npoints))
    allocate(output%TS1_pft(npoints,PFT))
    allocate(output%TS2_pft(npoints,PFT))
    allocate(output%TS3_pft(npoints,PFT))
    allocate(output%TS4_pft(npoints,PFT))
    allocate(output%TS5_pft(npoints,PFT))
    allocate(output%PondWater_pft(npoints,PFT))
    allocate(output%Rain_g_pft(npoints,PFT))

    allocate(output%f_soilwater(npoints))
    allocate(output%f_feileaf(npoints))
    allocate(output%f_soilwater_pft(npoints,PFT))
    allocate(output%f_feileaf_pft(npoints,PFT))
    allocate(output%f_Tleaf(npoints))
    allocate(output%f_Tleaf_pft(npoints,PFT))
    ! 2024/03/12
    allocate(output%LHa(npoints))
    allocate(output%LHapft(npoints,PFT))

    ! sunlit and shaded gpp and tr
    allocate(output%GPP_o_sunlit_pft(npoints,PFT))
    allocate(output%GPP_o_shaded_pft(npoints,PFT))
    allocate(output%GPP_u_sunlit_pft(npoints,PFT))
    allocate(output%GPP_u_shaded_pft(npoints,PFT))
    allocate(output%GPP_o_sunlit(npoints))
    allocate(output%GPP_o_shaded(npoints))
    allocate(output%GPP_u_sunlit(npoints))
    allocate(output%GPP_u_shaded(npoints))

    allocate(output%TR_o_sunlit_pft(npoints,PFT))
    allocate(output%TR_o_shaded_pft(npoints,PFT))
    allocate(output%TR_u_sunlit_pft(npoints,PFT))
    allocate(output%TR_u_shaded_pft(npoints,PFT))
    allocate(output%TR_o_sunlit(npoints))
    allocate(output%TR_o_shaded(npoints))
    allocate(output%TR_u_sunlit(npoints))
    allocate(output%TR_u_shaded(npoints))

    allocate(output%Eil(npoints))
    allocate(output%Evap_soil(npoints))
    allocate(output%Evap_SW(npoints))
    allocate(output%EiS(npoints))
    allocate(output%Evap_SS(npoints))

    allocate(output%Eil_pft(npoints,PFT))
    allocate(output%Evap_soil_pft(npoints,PFT))
    allocate(output%Evap_SW_pft(npoints,PFT))
    allocate(output%EiS_pft(npoints,PFT))
    allocate(output%Evap_SS_pft(npoints,PFT))

    output%Eil(:) = 0.
    output%Evap_soil(:) = 0.
    output%Evap_SW(:) = 0.
    output%EiS(:) = 0.
    output%Evap_SS(:) = 0.

    output%Eil_pft(:,:) = 0.
    output%Evap_soil_pft(:,:) = 0.
    output%Evap_SW_pft(:,:) = 0.
    output%EiS_pft(:,:) = 0.
    output%Evap_SS_pft(:,:) = 0.

    output%GPP_o_sunlit_pft(:,:)=0.
    output%GPP_o_shaded_pft(:,:)=0.
    output%GPP_u_sunlit_pft(:,:)=0.
    output%GPP_u_shaded_pft(:,:)=0.
    output%TR_o_sunlit_pft(:,:)=0.
    output%TR_o_shaded_pft(:,:)=0.
    output%TR_u_sunlit_pft(:,:)=0.
    output%TR_u_shaded_pft(:,:)=0.

    output%GPP_o_sunlit(:)=0.
    output%GPP_o_shaded(:)=0.
    output%GPP_u_sunlit(:)=0.
    output%GPP_u_shaded(:)=0.
    output%TR_o_sunlit(:)=0.
    output%TR_o_shaded(:)=0.
    output%TR_u_sunlit(:)=0.
    output%TR_u_shaded(:)=0.


    output%Thetam2_pft(:,:)=0.
    output%Thetam3_pft(:,:)=0.
    output%Thetam4_pft(:,:)=0.
    output%Thetam5_pft(:,:)=0.
    output%Thetam_layer2(:)=0.
    output%Thetam_layer3(:)=0.
    output%Thetam_layer4(:)=0.
    output%Thetam_layer5(:)=0.
    output%SWP1_pft(:,:)=0.
    output%SWP2_pft(:,:)=0.
    output%SWP3_pft(:,:)=0.
    output%SWP4_pft(:,:)=0.
    output%SWP5_pft(:,:)=0.
    output%SWP_layer1(:)=0.
    output%SWP_layer2(:)=0.
    output%SWP_layer3(:)=0.
    output%SWP_layer4(:)=0.
    output%SWP_layer5(:)=0.
    output%Qupt(:)=0.
    output%Quptpft(:,:)=0.
    output%TS_layer1(:)=0.
    output%TS_layer2(:)=0.
    output%TS_layer3(:)=0.
    output%TS_layer4(:)=0.
    output%TS_layer5(:)=0.
    output%PondWater(:)=0.
    output%Rain_g(:)=0.
    output%TS1_pft(:,:)=0.
    output%TS2_pft(:,:)=0.
    output%TS3_pft(:,:)=0.
    output%TS4_pft(:,:)=0.
    output%TS5_pft(:,:)=0.
    output%PondWater_pft(:,:)=0.
    output%Rain_g_pft(:,:)=0.
    !!!
    output%GPPpft(:,:)  = 0.
    output%SIFpft(:,:)  = 0.
    output%SIFpft_sat(:,:)  = 0.
    output%NPPpft(:,:)  = 0.
    output%NEPpft(:,:)  = 0.
    output%SHpft(:,:)   = 0.
    output%LHpft(:,:)   = 0.
    output%Transpft(:,:)= 0.
    output%Evappft(:,:) = 0.
    output%Net_Radpft(:,:) = 0.
    output%LAIpft(:,:)  = 0.
    output%Thetampft(:,:)  = 0.
    output%fAPARpft(:,:)  = 0.
    output%VODpft(:,:)  = 0.
    output%COS_fluxpft(:,:)  = 0.
    output%NPP_yr_acc(:,:)  = 0.

    ! 2023/06/30
    output%PWSpft(:,:)      = 0.
    output%ETapft(:,:)      = 0.

    output%GPP(:)       = 0.
    output%SIF(:)       = 0.
    output%SIF_sat(:)   = 0.
    output%NPP(:)       = 0.
    output%NEP(:)       = 0.
    output%LAI(:)       = 0.
    output%SH(:)        = 0.
    output%LH(:)        = 0.
    output%Trans(:)     = 0.
    output%Evap(:)      = 0.
    output%Net_Rad(:)   = 0.
    output%Thetam(:)   = 0.
    output%fAPAR(:)    = 0.
    output%COS_flux(:) = 0.
    output%VOD(:)      = 0.

    ! 2023/06/30
    output%PWS(:)      = 0.
    output%ETa(:)      = 0.
    output%fei_leaf(:)      = 0.
    output%fei_leafpft(:,:)      = 0.

    output%f_soilwater(:) = 0.
    output%f_feileaf(:) = 0.
    output%f_Tleaf(:) = 0.
    output%f_soilwater_pft(:,:) = 0.
    output%f_feileaf_pft(:,:) = 0.
    output%f_Tleaf_pft(:,:) = 0.

    ! 2024/03/12
    output%LHa(:)        = 0.
    output%LHapft(:,:)        = 0.

  end subroutine InitOuput

  subroutine InitMCOuput()
    implicit none
    !--iLab::avoid pointer
    ! type(res),pointer::p
    ! p=>output

    allocate(mc_output%SMpft(npoints,PFT))
    allocate(mc_output%SM(mc_length,npoints))
    allocate(mc_output%obsSM(mc_length,npoints))
    allocate(mc_output%tempSM(npoints))

    allocate(mc_output%LWPpft(npoints,PFT))
    allocate(mc_output%LWP(mc_length,npoints))
    allocate(mc_output%obsLWP(mc_length,npoints))

    mc_output%SM(:,:)      = 0.
    mc_output%SMpft(:,:)      = 0.
    mc_output%obsSM(:,:)      = 0.
    mc_output%tempSM(:)  = 0.

    mc_output%LWP(:,:)      = 0.
    mc_output%LWPpft(:,:)      = 0.
    mc_output%obsLWP(:,:)      = 0.

  end subroutine InitMCOuput


end module bepstypeInit
