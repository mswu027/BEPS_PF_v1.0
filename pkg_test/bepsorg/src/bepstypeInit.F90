

! This module is used to initialize bepstype variables
! Editted by J.Wang
! Date: 10May2017

module bepstypeInit
  use shr_kind_mod,only: r8=>shr_kind_r8
  !--iLab::added further entities (as module bepstype does no longer use beps_par)
  use beps_par,only:npoints,PFT,max_layers,texture,nparameters
  use bepstype
  use controlInput_mod
  implicit none
  !--iLab::can avoid 'save' here since no variables are declared,
  !        use-associated entities should already have the attribute
  ! save

  public  :: Initbepstype



  private :: initatm2lnd
  private :: InitSurf
  private :: InitOuput
  private :: InitAssim
  private :: InitPF
  public  :: InitPF_obs
  private :: InitPF_resample
  private :: Initparam_gdd
contains

  subroutine Initbepstype()
    implicit none

    allocate(v2last(npoints,0:40,PFT))
    v2last  = 0.

    call initatm2lnd()




    call InitSurf()

    call InitSoilstat()

    call InitOuput()

    call InitAssim()

    call InitPF()

    call InitPF_obs()

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

    allocate(clim%Srad(npoints))
    allocate(clim%Rh(npoints))
    allocate(clim%Rain(npoints))
    allocate(clim%Snow(npoints))
    allocate(clim%Swdr(npoints))
    allocate(clim%Swdf(npoints))

    clim%Temp(:)      = 0.
    clim%Tempmx(:)      = 0.
    clim%Tempmn(:)      = 0.
    clim%Wind(:)      = 0.
    clim%Srad(:)     = 0.
    clim%Rh(:)       = 0.
    clim%Rain(:)     = 0.
    clim%Snow(:)     = 0.
    clim%Swdr(:)     = 0.
    clim%Swdf(:)     = 0.

    return
  end subroutine initatm2lnd


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
    !allocate(bound%p_Vcmax(PFT))
    !allocate(bound%p_q10(PFT))
    !allocate(bound%p_drainage(PFT))
    !allocate(bound%p_beta(PFT))
    !allocate(bound%p_Ksat(texture))
    !allocate(bound%p_b(texture))
    !allocate(bound%u_Vcmax(PFT))
    !allocate(bound%u_q10(PFT))
    !allocate(bound%u_drainage(PFT))
    !allocate(bound%u_beta(PFT))
    !allocate(bound%u_Ksat(texture))
    !allocate(bound%u_b(texture))

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

    !bound%p_Vcmax   = 0.
    !bound%p_q10     = 0.
    !bound%p_drainage   = 0.
    !bound%p_beta   = 0.
    !bound%p_Ksat   = 0.
    !bound%p_b   = 0.
    !bound%u_Vcmax   = 0.
    !bound%u_q10   = 0.
    !bound%u_drainage   = 0.
    !bound%u_beta   = 0.
    !bound%u_Ksat   = 0.
    !bound%u_b   = 0.

  end subroutine InitSurf

  subroutine InitAssim()
    implicit none
    !--iLab::avoid pointer
    ! type(para),pointer  ::p
    ! p => assim

    allocate(assim%p_Vcmax(PFT,npoints))
    allocate(assim%p_q10(PFT,npoints))
    allocate(assim%p_VJ_slope(PFT,npoints))
    allocate(assim%p_N_leaf(PFT,npoints))
    allocate(assim%p_r_decay(PFT,npoints))
    allocate(assim%p_b_h2o(PFT,npoints))
    allocate(assim%p_sif_alpha(PFT,npoints))
    allocate(assim%p_sif_beta(PFT,npoints))
    allocate(assim%p_taweff(PFT,npoints))
    allocate(assim%p_D0(PFT,npoints))
    allocate(assim%p_Ksat_scalar(texture,npoints))
    allocate(assim%p_b_scalar(texture,npoints))
    allocate(assim%p_m_h2o(PFT,npoints))
    allocate(assim%p_f_leaf(PFT,npoints))
    allocate(assim%p_kc25(PFT,npoints))
    allocate(assim%p_ko25(PFT,npoints))
    allocate(assim%p_tau25(PFT,npoints))
    !allocate(assim%p_f_lr)

    ! allocate(assim%p_agb2vod(PFT,npoints)) not using the agb to vod code

    ! 2023/06/30
    allocate(assim%p_f_resp(PFT,npoints))
    allocate(assim%p_VN_slope(PFT,npoints))
    allocate(assim%p_f_decay(PFT,npoints))
    allocate(assim%p_bwb(PFT,npoints))
    allocate(assim%p_a(PFT,npoints))
    allocate(assim%p_b(PFT,npoints))
    allocate(assim%p_c(PFT,npoints))

    allocate(assim%u_Vcmax(PFT,npoints))
    allocate(assim%u_q10(PFT,npoints))
    allocate(assim%u_VJ_slope(PFT,npoints))
    allocate(assim%u_N_leaf(PFT,npoints))
    allocate(assim%u_r_decay(PFT,npoints))
    allocate(assim%u_b_h2o(PFT,npoints))
    allocate(assim%u_sif_alpha(PFT,npoints))
    allocate(assim%u_sif_beta(PFT,npoints))
    allocate(assim%u_taweff(PFT,npoints))
    allocate(assim%u_D0(PFT,npoints))
    allocate(assim%u_Ksat_scalar(texture,npoints))
    allocate(assim%u_b_scalar(texture,npoints))
    allocate(assim%u_m_h2o(PFT,npoints))
    allocate(assim%u_f_leaf(PFT,npoints))
    allocate(assim%u_kc25(PFT,npoints))
    allocate(assim%u_ko25(PFT,npoints))
    allocate(assim%u_tau25(PFT,npoints))
    !allocate(assim%u_f_lr)

    ! allocate(assim%u_agb2vod(PFT,npoints))

    ! 2023/06/30
    allocate(assim%u_f_resp(PFT,npoints))
    allocate(assim%u_VN_slope(PFT,npoints))
    allocate(assim%u_f_decay(PFT,npoints))
    allocate(assim%u_bwb(PFT,npoints))
    allocate(assim%u_a(PFT,npoints))
    allocate(assim%u_b(PFT,npoints))
    allocate(assim%u_c(PFT,npoints))


    assim%p_Vcmax   = 0.
    assim%p_q10     = 0.
    assim%p_VJ_slope   = 0.
    assim%p_N_leaf=0.
    assim%p_r_decay=0.
    assim%p_b_h2o=0.
    assim%p_sif_alpha   = 0.
    assim%p_sif_beta   = 0.
    assim%p_taweff   = 0.
    assim%p_D0   = 0.
    assim%p_Ksat_scalar   = 0.
    assim%p_b_scalar  = 0.
    assim%p_m_h2o=0.
    assim%p_f_leaf   = 0.
    assim%p_kc25   = 0.
    assim%p_ko25   = 0.
    assim%p_tau25   = 0.

    !assim%p_f_lr   = 0.
    ! assim%p_agb2vod   = 0.

    ! 2023/06/30
    assim%p_f_resp   = 0.
    assim%p_VN_slope   = 0.
    assim%p_f_decay   = 0.
    assim%p_bwb   = 0.
    assim%p_a   = 0.
    assim%p_b   = 0.
    assim%p_c   = 0.

    assim%u_Vcmax   = 0.
    assim%u_q10     = 0.
    assim%u_VJ_slope   = 0.
    assim%u_N_leaf=0.
    assim%u_r_decay=0.
    assim%u_b_h2o=0.
    assim%u_sif_alpha   = 0.
    assim%u_sif_beta   = 0.
    assim%u_taweff   = 0.
    assim%u_D0   = 0.
    assim%u_Ksat_scalar   = 0.
    assim%u_b_scalar  = 0.
    assim%u_m_h2o=0.
    assim%u_f_leaf   = 0.
    assim%u_kc25   = 0.
    assim%u_ko25   = 0.
    assim%u_tau25   = 0.

    !assim%u_f_lr   = 0.
    ! assim%u_agb2vod   = 0.

	! 2023/06/30
    assim%u_f_resp   = 0.
    assim%u_VN_slope   = 0.
    assim%u_f_decay   = 0.
    assim%u_bwb   = 0.
    assim%u_a   = 0.
    assim%u_b   = 0.
    assim%u_c   = 0.


  end subroutine InitAssim

  subroutine InitPF()
    implicit none
    !--iLab::avoid pointer
    ! type(para),pointer  ::p
    ! p => assim

    allocate(PF%Vcmax(parloop,npoints))
    allocate(PF%q10(parloop,npoints))
    allocate(PF%VJ_slope(parloop,npoints))
    allocate(PF%N_leaf(parloop,npoints))
    allocate(PF%r_decay(parloop,npoints))
    allocate(PF%b_h2o(parloop,npoints))
    allocate(PF%sif_alpha(parloop,npoints))
    allocate(PF%sif_beta(parloop,npoints))
    allocate(PF%taweff(parloop,npoints))
    allocate(PF%D0(parloop,npoints))
    allocate(PF%Ksat_scalar(parloop,npoints))
    allocate(PF%b_scalar(parloop,npoints))
    allocate(PF%m_h2o(parloop,npoints))
    allocate(PF%f_leaf(parloop,npoints))
    allocate(PF%kc25(parloop,npoints))
    allocate(PF%ko25(parloop,npoints))
    allocate(PF%tau25(parloop,npoints))

    !allocate(assim%p_f_lr)
    ! allocate(PF%agb2vod(parloop,npoints))

    allocate(PF%pfweight(parloop,npoints))
    allocate(PF%pfweightupdate(parloop,npoints))

	! 2023/06/30
    allocate(PF%f_resp(parloop,npoints))
    allocate(PF%VN_slope(parloop,npoints))
    allocate(PF%f_decay(parloop,npoints))
    allocate(PF%bwb(parloop,npoints))
    allocate(PF%a(parloop,npoints))
    allocate(PF%b(parloop,npoints))
    allocate(PF%c(parloop,npoints))

    PF%Vcmax   = 0.
    PF%q10     = 0.
    PF%VJ_slope   = 0.
    PF%N_leaf=0.
    PF%r_decay=0.
    PF%b_h2o=0.
    PF%sif_alpha   = 0.
    PF%sif_beta   = 0.
    PF%taweff   = 0.
    PF%D0   = 0.
    PF%Ksat_scalar   = 0.
    PF%b_scalar  = 0.
    PF%m_h2o=0.
    PF%f_leaf   = 0.
    PF%kc25   = 0.
    PF%ko25   = 0.
    PF%tau25   = 0.

    !assim%p_f_lr   = 0.
    !PF%agb2vod   = 0.

    PF%pfweight   = 0.
    PF%pfweightupdate   = 0.

    ! 2023/06/30
    PF%f_resp = 0.
    PF%VN_slope = 0.
    PF%f_decay = 0.
    PF%bwb = 0.
    PF%a = 0.
    PF%b = 0.
    PF%c = 0.

  end subroutine InitPF

  subroutine InitPF_obs()
    implicit none

    ! allocate(PF_obs%obs_GPP(npoints))
    allocate(PF_obs%obs_VOD(npoints)) ! changed to VOD from obs_GPP

     ! PF_obs%obs_GPP   = 0.
     PF_obs%obs_VOD   = 0.

  end subroutine InitPF_obs

  subroutine InitPF_resample()
    implicit none

    allocate(PF_resample%outparticles(parloop,20)) ! 19 pamameters and one for pf_weight
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
    allocate(output%PWSpft(npoints,pft))
    allocate(output%ETapft(npoints,pft))
    allocate(output%fei_leaf(npoints))
    allocate(output%fei_leafpft(npoints,PFT))

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

  end subroutine InitOuput

end module bepstypeInit
