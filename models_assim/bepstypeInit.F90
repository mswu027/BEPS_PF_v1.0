! This module is used to initialize bepstype variables
! Editted by J.Wang
! Date: 10May2017

module bepstypeInit
  use shr_kind_mod,only: r8=>shr_kind_r8
  !--iLab::added further entities (as module bepstype does no longer use beps_par)
  use beps_par,only:npoints,PFT,max_layers,texture
  use bepstype
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
  private :: InitAssim
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

    call InitAssim()

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

    allocate(assim%p_Vcmax(PFT))
    allocate(assim%p_q10(PFT))
    allocate(assim%p_VJ_slope(PFT))
    allocate(assim%p_sif_alpha(PFT))
    allocate(assim%p_sif_beta(PFT))
    allocate(assim%p_taweff(PFT))
    allocate(assim%p_D0(PFT))
    allocate(assim%p_Ksat_scalar(texture))
    allocate(assim%p_b_scalar(texture))
    allocate(assim%p_f_leaf)
    allocate(assim%p_kc25)
    allocate(assim%p_ko25)
    allocate(assim%p_tau25)
    !allocate(assim%p_f_lr)
    allocate(assim%p_agb2vod)

    allocate(assim%u_Vcmax(PFT))
    allocate(assim%u_q10(PFT))
    allocate(assim%u_VJ_slope(PFT))
    allocate(assim%u_sif_alpha(PFT))
    allocate(assim%u_sif_beta(PFT))
    allocate(assim%u_taweff(PFT))
    allocate(assim%u_D0(PFT))
    allocate(assim%u_Ksat_scalar(texture))
    allocate(assim%u_b_scalar(texture))
    allocate(assim%u_f_leaf)
    allocate(assim%u_kc25)
    allocate(assim%u_ko25)
    allocate(assim%u_tau25)
    !allocate(assim%u_f_lr)
    allocate(assim%u_agb2vod)

    assim%p_Vcmax   = 0.
    assim%p_q10     = 0.
    assim%p_VJ_slope   = 0.
    assim%p_sif_alpha   = 0.
    assim%p_sif_beta   = 0.
    assim%p_taweff   = 0.
    assim%p_D0   = 0.
    assim%p_Ksat_scalar   = 0.
    assim%p_b_scalar  = 0.
    assim%p_f_leaf   = 0.
    assim%p_kc25   = 0.
    assim%p_ko25   = 0.
    assim%p_tau25   = 0.
    !assim%p_f_lr   = 0.
    assim%p_agb2vod   = 0.
    
    assim%u_Vcmax   = 0.
    assim%u_q10     = 0.
    assim%u_VJ_slope   = 0.
    assim%u_sif_alpha   = 0.
    assim%u_sif_beta   = 0.
    assim%u_taweff   = 0.
    assim%u_D0   = 0.
    assim%u_Ksat_scalar   = 0.
    assim%u_b_scalar  = 0.
    assim%u_f_leaf   = 0.
    assim%u_kc25   = 0.
    assim%u_ko25   = 0.
    assim%u_tau25   = 0.
    !assim%u_f_lr   = 0.
    assim%u_agb2vod   = 0.

  end subroutine InitAssim

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

    soilstat%n_layer(:)           = 0
    soilstat%Zp(:,:)              = 0.
    soilstat%Zsp(:,:)             = 0.
    soilstat%r_rain_g(:,:)        = 0.
    soilstat%r_drainage(:,:)      = 0.
    soilstat%r_root_decay(:,:)    = 0.
    soilstat%psi_min(:,:)         = 0.
    soilstat%alpha(:,:)           = 0.
    soilstat%f_soilwater(:,:)     = 0.

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
    allocate(output%VODpft(npoints,PFT))
    allocate(output%VOD(npoints))
    allocate(output%COS_fluxpft(npoints,PFT))
    allocate(output%COS_flux(npoints))
    allocate(output%NPP_yr_acc(npoints,PFT))

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
  end subroutine InitOuput

end module bepstypeInit
