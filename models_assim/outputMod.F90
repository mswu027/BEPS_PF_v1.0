module outputMod
  use shr_kind_mod,only:r8=>shr_kind_r8
  use bepstype
  use controlInput_mod, only:beps_out_dir, beps_rst_dir, &
       nhtfrq, nstpd, nscale, nlat, nlon, &
       check
  !--iLab::avoid beps_time_manager, all temporal information now passed as actual arguments
  ! use beps_time_manager
  use beps_par
  implicit none

  !! Here add variables for output by USERS
  real(r8),allocatable :: NEP9(:)
  real(r8),allocatable :: GPP9(:)
  real(r8),allocatable :: SIF9(:)
  real(r8),allocatable :: SIF9_sat(:)
  real(r8),allocatable :: NPP9(:)
  real(r8),allocatable :: GPPpft9(:,:)
  real(r8),allocatable :: SIFpft9(:,:)
  real(r8),allocatable :: LHpft9(:,:)
  real(r8),allocatable :: SHpft9(:,:)
  real(r8),allocatable :: Transpft9(:,:)
  real(r8),allocatable :: Evappft9(:,:)
  real(r8),allocatable :: Thetampft9(:,:)
  real(r8),allocatable :: SIFpft9_sat(:,:)
  real(r8),allocatable :: fAPARpft9(:,:)
  real(r8),allocatable :: VODpft9(:,:)
  real(r8),allocatable :: COS_fluxpft9(:,:)
  real(r8),allocatable :: temp9(:)
  real(r8),allocatable :: Wind9(:)
  real(r8),allocatable :: Rh9(:)
  real(r8),allocatable :: Rain9(:)
  real(r8),allocatable :: Snow9(:)
  real(r8),allocatable :: Swdr9(:)
  real(r8),allocatable :: Swdf9(:)
  real(r8),allocatable :: lai9(:)
  real(r8),allocatable :: LH9(:)
  real(r8),allocatable :: SH9(:)
  real(r8),allocatable :: Trans9(:)
  real(r8),allocatable :: Evap9(:)
  real(r8),allocatable :: Thetam9(:)
  real(r8),allocatable :: fAPAR9(:)
  real(r8),allocatable :: VOD9(:)
  real(r8),allocatable :: COS_flux9(:)
  ! 2023/07/19
  real(r8),allocatable :: fei_leaf9(:)
  real(r8),allocatable :: ETa9(:)
  real(r8),allocatable :: PWS9(:)
  real(r8),allocatable :: fei_leafpft9(:,:)
  real(r8),allocatable :: laipft9(:,:)
  real(r8),allocatable :: Qupt9(:)

  !2023/10/30
  real(r8),allocatable :: Thetam2_9(:)
  real(r8),allocatable :: Thetam3_9(:)
  real(r8),allocatable :: Thetam4_9(:)
  real(r8),allocatable :: Thetam5_9(:)
  real(r8),allocatable :: SWP1_9(:)
  real(r8),allocatable :: SWP2_9(:)
  real(r8),allocatable :: SWP3_9(:)
  real(r8),allocatable :: SWP4_9(:)
  real(r8),allocatable :: SWP5_9(:)

  real(r8),allocatable :: TS1_9(:)
  real(r8),allocatable :: TS2_9(:)
  real(r8),allocatable :: TS3_9(:)
  real(r8),allocatable :: TS4_9(:)
  real(r8),allocatable :: TS5_9(:)
  real(r8),allocatable :: PondWater_9(:)
  real(r8),allocatable :: Rain_g_9(:)
  real(r8),allocatable :: f_soilwater9(:)
  real(r8),allocatable :: f_feileaf9(:)
  !real(r8),allocatable :: f_Tleaf9(:)

  real(r8),allocatable :: LHa9(:)
  ! 2024/03/30
  real(r8),allocatable :: GPP_o_sunlit9(:)
  real(r8),allocatable :: GPP_o_shaded9(:)
  real(r8),allocatable :: GPP_u_sunlit9(:)
  real(r8),allocatable :: GPP_u_shaded9(:)

  real(r8),allocatable :: TR_o_sunlit9(:)
  real(r8),allocatable :: TR_o_shaded9(:)
  real(r8),allocatable :: TR_u_sunlit9(:)
  real(r8),allocatable :: TR_u_shaded9(:)

  real(r8),allocatable :: Eil9(:)
  real(r8),allocatable :: Evap_soil9(:)
  real(r8),allocatable :: Evap_SW9(:)
  real(r8),allocatable :: EiS9(:)
  real(r8),allocatable :: Evap_SS9(:)


  integer :: nst     = 0    ! for counting the simulation steps for monthly output

contains

  subroutine Init_output()
    implicit none

    allocate(NEP9(npoints))
    allocate(GPP9(npoints))
    allocate(NPP9(npoints))
    allocate(GPPpft9(npoints,PFT))
    allocate(SIFpft9(npoints,PFT))
    allocate(LHpft9(npoints,PFT))
    allocate(SHpft9(npoints,PFT))
    allocate(Transpft9(npoints,PFT))
    allocate(Evappft9(npoints,PFT))
    allocate(Thetampft9(npoints,PFT))
    allocate(SIFpft9_sat(npoints,PFT))
    allocate(SIF9(npoints))
    allocate(SIF9_sat(npoints))
    allocate(fAPARpft9(npoints,PFT))
    allocate(VODpft9(npoints,PFT))
    allocate(COS_fluxpft9(npoints,PFT))
    allocate(temp9(npoints))
    allocate(Wind9(npoints))
    allocate(Rh9(npoints))
    allocate(Rain9(npoints))
    allocate(Snow9(npoints))
    allocate(Swdr9(npoints))
    allocate(Swdf9(npoints))
    allocate(lai9(npoints))
    allocate(LH9(npoints))
    allocate(SH9(npoints))
    allocate(Trans9(npoints))
    allocate(Evap9(npoints))
    allocate(Thetam9(npoints))
    allocate(fAPAR9(npoints))
    allocate(VOD9(npoints))
    allocate(COS_flux9(npoints))
    ! 2023/07/19
    allocate(fei_leaf9(npoints))
    allocate(fei_leafpft9(npoints,PFT))
    allocate(ETa9(npoints))
    allocate(PWS9(npoints))
    allocate(laipft9(npoints,PFT))
    allocate(Qupt9(npoints))

    allocate(Thetam2_9(npoints))
    allocate(Thetam3_9(npoints))
    allocate(Thetam4_9(npoints))
    allocate(Thetam5_9(npoints))
    allocate(SWP1_9(npoints))
    allocate(SWP2_9(npoints))
    allocate(SWP3_9(npoints))
    allocate(SWP4_9(npoints))
    allocate(SWP5_9(npoints))

    allocate(TS1_9(npoints))
    allocate(TS2_9(npoints))
    allocate(TS3_9(npoints))
    allocate(TS4_9(npoints))
    allocate(TS5_9(npoints))
    allocate(PondWater_9(npoints))
    allocate(Rain_g_9(npoints))

    allocate(f_soilwater9(npoints))
    allocate(f_feileaf9(npoints))
    !allocate(f_Tleaf9(npoints))

    ! 2024/03/12
    allocate(LHa9(npoints))

    ! 2024/30/30
    allocate(GPP_o_sunlit9(npoints))
    allocate(GPP_o_shaded9(npoints))
    allocate(GPP_u_sunlit9(npoints))
    allocate(GPP_u_shaded9(npoints))

    allocate(TR_o_sunlit9(npoints))
    allocate(TR_o_shaded9(npoints))
    allocate(TR_u_sunlit9(npoints))
    allocate(TR_u_shaded9(npoints))

    allocate(Eil9(npoints))
    allocate(Evap_soil9(npoints))
    allocate(Evap_SW9(npoints))
    allocate(EiS9(npoints))
    allocate(Evap_SS9(npoints))

    Eil9(:) = 0.
    Evap_soil9(:) = 0.
    Evap_SW9(:) = 0.
    EiS9(:) = 0.
    Evap_SS9(:) = 0.

    GPP_o_sunlit9(:)=0.
    GPP_o_shaded9(:)=0.
    GPP_u_sunlit9(:)=0.
    GPP_u_shaded9(:)=0.

    TR_o_sunlit9(:)=0.
    TR_o_shaded9(:)=0.
    TR_u_sunlit9(:)=0.
    TR_u_shaded9(:)=0.

    Thetam2_9(:)=0.
    Thetam3_9(:)=0.
    Thetam4_9(:)=0.
    Thetam5_9(:)=0.
    SWP1_9(:)=0.
    SWP2_9(:)=0.
    SWP3_9(:)=0.
    SWP4_9(:)=0.
    SWP5_9(:)=0.

    TS1_9(:)=0.
    TS2_9(:)=0.
    TS3_9(:)=0.
    TS4_9(:)=0.
    TS5_9(:)=0.
    PondWater_9(:)=0.
    Rain_g_9(:)=0.

    NEP9(:)   = 0.0
    GPP9(:)   = 0.
    NPP9(:)   = 0.
    SIF9(:)   = 0.
    SIF9_sat(:)  = 0.
    GPPpft9(:,:)  = 0.
    SIFpft9(:,:)  = 0.
    LHpft9(:,:)  = 0.
    SHpft9(:,:)  = 0.
    Transpft9(:,:)  = 0.
    Evappft9(:,:)  = 0.
    Thetampft9(:,:)  = 0.
    SIFpft9_sat(:,:)  = 0.
    fAPARpft9(:,:) = 0.
    VODpft9(:,:)  = 0.
    COS_fluxpft9(:,:) = 0.
    temp9(:)      = 0.
    Wind9(:)      = 0.
    Rh9(:)        = 0.
    Rain9(:)      = 0.
    Snow9(:)      = 0.
    Swdr9(:)      = 0.
    Swdf9(:)      = 0.
    lai9(:)       = 0.
    LH9(:)       = 0.
    SH9(:)       = 0.
    Trans9(:)       = 0.
    Evap9(:)       = 0.
    Thetam9(:)       = 0.
    fAPAR9(:)     = 0.
    VOD9(:)      = 0.
    COS_flux9(:)    = 0.
    fei_leaf9(:) = 0. ! 2023/07/19
    fei_leafpft9(:,:) = 0.
    laipft9(:,:) = 0.
    ETa9(:) = 0.
    Qupt9(:) = 0.
    PWS9(:) = 0.
    f_soilwater9(:)=0.
    f_feileaf9(:)=0.
    !f_Tleaf9(:)=0.
    !2024/03/12
    LHa9(:) = 0.

  end subroutine Init_output

  !! average variables according to user's definition
  subroutine av_output(yr, mon, day, tod, kount, is_end_curr_month, ref_date, secs_since_ref)
    implicit none
    !-- iLab::turned yr,mon,day,tod to arguments and added the further arguments
    integer, intent(in) :: yr,mon,day,tod
    integer, intent(in) :: kount
    logical, intent(in) :: is_end_curr_month
    character(len=*), intent(in) :: ref_date
    real(r8), intent(in) :: secs_since_ref
    integer    :: ii,iii,p
    type(res),pointer   :: pp
    type(forc),pointer  :: ff
    type(surf),pointer  :: ss

    pp => output
    ff => clim
    ss => bound

    NEP9   = NEP9 + pp%NEP   !! accumulate
    GPP9   = GPP9 + pp%GPP
    SIF9   = SIF9 + pp%SIF
    SIF9_sat = SIF9_sat + pp%SIF_sat
    NPP9   = NPP9 + pp%NPP

    GPPpft9   = GPPpft9 + pp%GPPpft
    SIFpft9   = SIFpft9 + pp%SIFpft
    LHpft9   = LHpft9 + pp%LHpft
    SHpft9   = SHpft9 + pp%SHpft
    Transpft9   = Transpft9 + pp%Transpft
    Evappft9   = Evappft9 + pp%Evappft
    Thetampft9   = Thetampft9 + pp%Thetampft
    SIFpft9_sat = SIFpft9_sat + pp%SIFpft_sat
    fAPARpft9  =  fAPARpft9 + pp%fAPARpft
    VODpft9   = VODpft9 + pp%VODpft
    COS_fluxpft9 = COS_fluxpft9 + pp%COS_fluxpft
    temp9     = temp9+ff%Temp
    Wind9     = Wind9+ff%Wind
    Rh9       = Rh9+ff%Rh
    Rain9     = Rain9+ff%Rain
    Snow9     = Snow9+ff%Snow
    Swdr9     = Swdr9+ff%Swdr
    Swdf9     = Swdf9+ff%Swdf
    lai9      = lai9+pp%LAI
    LH9      = LH9+pp%LH
    SH9      = SH9+pp%SH
    Trans9      = Trans9+pp%Trans
    Evap9      = Evap9+pp%Evap
    Thetam9      = Thetam9+pp%Thetam
    fAPAR9     = fAPAR9+pp%fAPAR
    VOD9      = VOD9+pp%VOD
    COS_flux9 = COS_flux9+pp%COS_flux
    ! 2023/07/19
    fei_leaf9 = fei_leaf9+pp%fei_leaf
    fei_leafpft9 = fei_leafpft9+pp%fei_leafpft
    laipft9 = laipft9+pp%LAIpft
    ETa9 = ETa9+pp%ETa
    Qupt9 = Qupt9+pp%Qupt
    PWS9 = PWS9+pp%PWS

    ! 2023/10/30
    Thetam2_9=Thetam2_9+pp%Thetam_layer2
    Thetam3_9=Thetam3_9+pp%Thetam_layer3
    Thetam4_9=Thetam4_9+pp%Thetam_layer4
    Thetam5_9=Thetam5_9+pp%Thetam_layer5
    SWP1_9=SWP1_9+pp%SWP_layer1
    SWP2_9=SWP2_9+pp%SWP_layer2
    SWP3_9=SWP3_9+pp%SWP_layer3
    SWP4_9=SWP4_9+pp%SWP_layer4
    SWP5_9=SWP5_9+pp%SWP_layer5

    TS1_9=TS1_9+pp%TS_layer1
    TS2_9=TS2_9+pp%TS_layer2
    TS3_9=TS3_9+pp%TS_layer3
    TS4_9=TS4_9+pp%TS_layer4
    TS5_9=TS5_9+pp%TS_layer5
    PondWater_9=PondWater_9+pp%PondWater
    Rain_g_9=Rain_g_9+pp%Rain_g

    f_soilwater9=f_soilwater9+pp%f_soilwater
    f_feileaf9 = f_feileaf9+pp%f_feileaf
    !f_Tleaf9 = f_Tleaf9+pp%f_Tleaf

    LHa9 = LHa9 + pp%LHa

    ! 2024/03/30
    GPP_o_sunlit9 = GPP_o_sunlit9 + pp%GPP_o_sunlit
    GPP_o_shaded9 = GPP_o_shaded9 + pp%GPP_o_shaded
    GPP_u_sunlit9 = GPP_u_sunlit9 + pp%GPP_u_sunlit
    GPP_u_shaded9 = GPP_u_shaded9 + pp%GPP_u_shaded

    TR_o_sunlit9 = TR_o_sunlit9 + pp%TR_o_sunlit
    TR_o_shaded9 = TR_o_shaded9 + pp%TR_o_shaded
    TR_u_sunlit9 = TR_u_sunlit9 + pp%TR_u_sunlit
    TR_u_shaded9 = TR_u_shaded9 + pp%TR_u_shaded

    Eil9 = Eil9 + pp%Eil
    Evap_soil9 = Evap_soil9 + pp%Evap_soil
    Evap_SW9 = Evap_SW9 + pp%Evap_SW
    EiS9 = EiS9 + pp%EiS
    Evap_SS9 = Evap_SS9 + pp%Evap_SS


    !! currently I did not include the satellite SIF when nhtfrq < 0 @J.Wang
    if(nhtfrq < 0) then
       ! kount  = get_nstep()

       if(mod(kount,nstpd) ==0) then
          NEP9   = NEP9/nstpd    !! average
          GPP9   = GPP9/nstpd
          SIF9   = SIF9/nstpd
          NPP9   = NPP9/nstpd

          GPPpft9  = GPPpft9/nstpd
          SIFpft9  = SIFpft9/nstpd
          LHpft9  = LHpft9/nstpd
          SHpft9  = SHpft9/nstpd
          Transpft9  = Transpft9/nstpd
          Evappft9  = Evappft9/nstpd
          Thetampft9  = Thetampft9/nstpd
          fAPARpft9 = fAPARpft9/nstpd
          VODpft9 = VODpft9/nstpd
          COS_fluxpft9 = COS_fluxpft9/nstpd

          if(nhtfrq == -24) then
             SIF9_sat       = SIF9_sat
             SIFpft9_sat    = SIFpft9_sat
          else
             SIF9_sat = 0.
             SIFpft9_sat = 0.
          end if

          temp9    = temp9/nstpd
          Wind9    = Wind9/nstpd
          Rh9      = Rh9/nstpd
          Rain9    = Rain9/nstpd
          Snow9    = Snow9/nstpd
          Swdr9    = Swdr9/nstpd
          Swdf9    = Swdf9/nstpd
          lai9     = lai9/nstpd
          LH9     = LH9/nstpd
          SH9     = SH9/nstpd
          Trans9     = Trans9/nstpd
          Evap9     = Evap9/nstpd
          Thetam9     = Thetam9/nstpd
          fAPAR9    = fAPAR9/nstpd
          VOD9      = VOD9/nstpd
          COS_flux9  = COS_flux9/nstpd
          fei_leaf9 = fei_leaf9/nstpd ! 2023/07/19
          fei_leafpft9 = fei_leafpft9/nstpd
          laipft9 = laipft9/nstpd
          ETa9 = ETa9/nstpd
          Qupt9 = Qupt9/nstpd
          PWS9 = PWS9/nstpd

          ! 2023/10/30
          Thetam2_9=Thetam2_9/nstpd
          Thetam3_9=Thetam3_9/nstpd
          Thetam4_9=Thetam4_9/nstpd
          Thetam5_9=Thetam5_9/nstpd
          SWP1_9=SWP1_9/nstpd
          SWP2_9=SWP2_9/nstpd
          SWP3_9=SWP3_9/nstpd
          SWP4_9=SWP4_9/nstpd
          SWP5_9=SWP5_9/nstpd

          TS1_9=TS1_9/nstpd
          TS2_9=TS2_9/nstpd
          TS3_9=TS3_9/nstpd
          TS4_9=TS4_9/nstpd
          TS5_9=TS5_9/nstpd
          PondWater_9=PondWater_9/nstpd
          Rain_g_9=Rain_g_9/nstpd

          f_soilwater9=f_soilwater9/nstpd
          f_feileaf9=f_feileaf9/nstpd
          !f_Tleaf9=f_Tleaf9/nstpd

          LHa9 = LHa9/nstpd

          GPP_o_sunlit9 = GPP_o_sunlit9/nstpd
          GPP_o_shaded9 = GPP_o_shaded9/nstpd
          GPP_u_sunlit9 = GPP_u_sunlit9/nstpd
          GPP_u_shaded9 = GPP_u_shaded9/nstpd

          TR_o_sunlit9 = TR_o_sunlit9/nstpd
          TR_o_shaded9 = TR_o_shaded9/nstpd
          TR_u_sunlit9 = TR_u_sunlit9/nstpd
          TR_u_shaded9 = TR_u_shaded9/nstpd

          Eil9 = Eil9/nstpd
          Evap_soil9 = Evap_soil9/nstpd
          Evap_SW9 = Evap_SW9/nstpd
          EiS9 = EiS9/nstpd
          Evap_SS9 = Evap_SS9/nstpd


          if (nscale == 0) then
             call write_output_global(yr, mon, day, tod)
          else
             call write_output_site(yr, mon, day, tod, ref_date, secs_since_ref)
          end if

          NEP9   = 0.
          GPP9   = 0.
          SIF9   = 0.
          NPP9   = 0.
          SIF9_sat = 0.

          GPPpft9  = 0.
          SIFpft9  = 0.
          LHpft9  = 0.
          SHpft9  = 0.
          Transpft9  = 0.
          Evappft9  = 0.
          Thetampft9  = 0.
          SIFpft9_sat = 0.
          fAPARpft9 = 0.
          VODpft9  = 0.
          COS_fluxpft9 = 0.
          temp9    = 0.
          Wind9    = 0.
          Rh9      = 0.
          Rain9    = 0.
          Snow9    = 0.
          Swdr9    = 0.
          Swdf9    = 0.
          lai9     = 0.
          LH9     = 0.
          SH9     = 0.
          Trans9     = 0.
          Evap9     = 0.
          Thetam9     = 0.
          fAPAR9   = 0.
          VOD9     = 0.
          COS_flux9  = 0.
          fei_leaf9 = 0. ! 2023/07/19
          fei_leafpft9 = 0.
          laipft9 = 0.
          ETa9=0.
          Qupt9=0.
          PWS9 = 0.
         ! 2023/10/30
          Thetam2_9=0.
          Thetam3_9=0.
          Thetam4_9=0.
          Thetam5_9=0.
          SWP1_9=0.
          SWP2_9=0.
          SWP3_9=0.
          SWP4_9=0.
          SWP5_9=0.

          TS1_9=0.
          TS2_9=0.
          TS3_9=0.
          TS4_9=0.
          TS5_9=0.
          PondWater_9=0.
          Rain_g_9=0.

          f_soilwater9=0.
          f_feileaf9=0.
          !f_Tleaf9=0.

          LHa9 = 0.

          GPP_o_sunlit9 = 0.
          GPP_o_shaded9 = 0.
          GPP_u_sunlit9 = 0.
          GPP_u_shaded9 = 0.

          TR_o_sunlit9 = 0.
          TR_o_shaded9 = 0.
          TR_u_sunlit9 = 0.
          TR_u_shaded9 = 0.

          Eil9 = 0.
          Evap_soil9 = 0.
          Evap_SW9 = 0.
          EiS9 = 0.
          Evap_SS9 = 0.


       end if
    else if(nhtfrq ==0) then   !!monthly output
       nst      = nst +1
       if(is_end_curr_month) then
          NEP9   = NEP9/nst    !! average
          GPP9   = GPP9/nst
          SIF9   = SIF9/nst
          NPP9   = NPP9/nst

          GPPpft9  = GPPpft9/nst     !kg/m2/s
          SIFpft9  = SIFpft9/nst
          LHpft9  = LHpft9/nst
          SHpft9  = SHpft9/nst
          Transpft9  = Transpft9/nst
          Evappft9  = Evappft9/nst
          Thetampft9  = Thetampft9/nst
          fAPARpft9   = fAPARpft9/nst
          VODpft9     = VODpft9/nst
          COS_fluxpft9 = COS_fluxpft9/nst
          temp9    = temp9/nst
          Wind9    = Wind9/nst
          Rh9      = Rh9/nst
          Rain9    = Rain9/nst
          Snow9    = Snow9/nst
          Swdr9    = Swdr9/nst
          Swdf9    = Swdf9/nst
          lai9     = lai9/nst
          LH9     = LH9/nst
          SH9     = SH9/nst
          Trans9     = Trans9/nst
          Evap9     = Evap9/nst
          Thetam9     = Thetam9/nst
          fAPAR9      = fAPAR9/nst
          VOD9        = VOD9/nst
          COS_flux9   = COS_flux9/nst
          fei_leaf9 = fei_leaf9/nst ! 2023/07/19
          fei_leafpft9 = fei_leafpft9/nst
          laipft9 = laipft9/nst

          ETa9 = ETa9/nst
          Qupt9 = Qupt9/nst
          PWS9 = PWS9/nst

          !--iLab::yr,mon,day,tod now provided as arguments
          ! call get_prev_date(yr, mon, day, tod)
          !!              write(*,*) "write out data on ",yr,mon,day
          SIFpft9_sat   = SIFpft9_sat/day
          SIF9_sat      = SIF9_sat/day

          ! 2023/10/30
          Thetam2_9=Thetam2_9/nst
          Thetam3_9=Thetam3_9/nst
          Thetam4_9=Thetam4_9/nst
          Thetam5_9=Thetam5_9/nst
          SWP1_9=SWP1_9/nst
          SWP2_9=SWP2_9/nst
          SWP3_9=SWP3_9/nst
          SWP4_9=SWP4_9/nst
          SWP5_9=SWP5_9/nst

          TS1_9=TS1_9/nst
          TS2_9=TS2_9/nst
          TS3_9=TS3_9/nst
          TS4_9=TS4_9/nst
          TS5_9=TS5_9/nst
          PondWater_9=PondWater_9/nst
          Rain_g_9=Rain_g_9/nst

          f_soilwater9=f_soilwater9/nst
          f_feileaf9=f_feileaf9/nst
          !f_Tleaf9=f_Tleaf9/nst

          LHa9 = LHa9/nst

          GPP_o_sunlit9 = GPP_o_sunlit9/nst
          GPP_o_shaded9 = GPP_o_shaded9/nst
          GPP_u_sunlit9 = GPP_u_sunlit9/nst
          GPP_u_shaded9 = GPP_u_shaded9/nst

          TR_o_sunlit9 = TR_o_sunlit9/nst
          TR_o_shaded9 = TR_o_shaded9/nst
          TR_u_sunlit9 = TR_u_sunlit9/nst
          TR_u_shaded9 = TR_u_shaded9/nst

          Eil9 = Eil9/nst
          Evap_soil9 = Evap_soil9/nst
          Evap_SW9 = Evap_SW9/nst
          EiS9 = EiS9/nst
          Evap_SS9 = Evap_SS9/nst

          if (nscale == 0) then
             call write_output_global(yr,mon,day,tod)
          else
             call write_output_site(yr, mon, day, tod, ref_date, secs_since_ref)
          end if

          NEP9     = 0.
          GPP9     = 0.
          SIF9     = 0.
          SIF9_sat = 0.
          NPP9     = 0.

          GPPpft9  = 0.
          SIFpft9  = 0.
          LHpft9  = 0.
          SHpft9  = 0.
          Transpft9  = 0.
          Evappft9  = 0.
          Thetampft9  = 0.
          SIFpft9_sat  = 0.
          fAPARpft9 = 0.
          VODpft9  = 0.
          COS_fluxpft9  = 0.
          temp9    = 0.
          Wind9    = 0.
          Rh9      = 0.
          Rain9    = 0.
          Snow9    = 0.
          Swdr9    = 0.
          Swdf9    = 0.
          lai9     = 0.
          LH9     = 0.
          SH9     = 0.
          Trans9     = 0.
          Evap9     = 0.
          Thetam9     = 0.
          fAPAR9   = 0.
          VOD9     = 0.
          COS_flux9 = 0.
          nst      = 0
          fei_leaf9 = 0. ! 2023/07/19
          fei_leafpft9 = 0.
          laipft9 = 0.
          ETa9 = 0.
          Qupt9 = 0.
          PWS9 = 0.

        ! 2023/10/30
          Thetam2_9=0.
          Thetam3_9=0.
          Thetam4_9=0.
          Thetam5_9=0.
          SWP1_9=0.
          SWP2_9=0.
          SWP3_9=0.
          SWP4_9=0.
          SWP5_9=0.

          TS1_9=0.
          TS2_9=0.
          TS3_9=0.
          TS4_9=0.
          TS5_9=0.
          PondWater_9=0.
          Rain_g_9=0.

          f_soilwater9=0.
          f_feileaf9=0.
          !f_Tleaf9=0.

          LHa9 = 0.

          GPP_o_sunlit9 = 0.
          GPP_o_shaded9 = 0.
          GPP_u_sunlit9 = 0.
          GPP_u_shaded9 = 0.

          TR_o_sunlit9 = 0.
          TR_o_shaded9 = 0.
          TR_u_sunlit9 = 0.
          TR_u_shaded9 = 0.

          Eil9 = 0.
          Evap_soil9 = 0.
          Evap_SW9 = 0.
          EiS9 = 0.
          Evap_SS9 = 0.

       end if
    end if

    !deallocate(NEP9)
    !deallocate(GPP9)
    !deallocate(NPP9)
    !deallocate(GPPpft9)
    !deallocate(SIFpft9)
    !deallocate(Thetampft9)
    !deallocate(SIFpft9_sat)
    !deallocate(SIF9)
    !deallocate(SIF9_sat)
    !deallocate(temp9)
    !deallocate(Wind9)
    !deallocate(Rh9)
    !deallocate(Rain9)
    !deallocate(Snow9)
    !deallocate(Swdr9)
    !deallocate(Swdf9)
    !deallocate(lai9)
    !deallocate(Thetam9)

  end subroutine av_output
! write global need to be checked! @Lu HU
  subroutine write_output_global(yy,mm,dd,tod)
    use netcdf
    implicit none
    !--iLab::yy,mm,dd,tod turned to arguments
    integer, intent(in) :: yy,mm,dd,tod
    real(r8),dimension(nlp)        :: NEP1,GPP1,SIF1,SIF_sat1,NPP1,temp1,Wind1,Rh1,Rain1,Snow1,Swdr1,Swdf1,lai1,LH1,SH1, &
         Trans1,Evap1,Thetam1,fAPAR1,VOD1,COS_flux1, fei_leaf1  ! add fei_leaf 2023/07/20
    real(r8),dimension(nlon*nlat)  :: NEP2,GPP2,SIF2,SIF_sat2,NPP2,temp2,Wind2,Rh2,Rain2,Snow2,Swdr2,Swdf2,lai2,LH2,SH2, &
         Trans2,Evap2,Thetam2,fAPAR2,VOD2,COS_flux2, fei_leaf2
    real(r8),dimension(nlon,nlat)  :: NEP3,GPP3,SIF3,SIF_sat3,NPP3,temp3,Wind3,Rh3,Rain3,Snow3,Swdr3,Swdf3,lai3,LH3,SH3, &
         Trans3,Evap3,Thetam3,fAPAR3,VOD3,COS_flux3, fei_leaf3
    real(r8),dimension(nlp,PFT)           :: GPPpft1,LHpft1,SHpft1,Transpft1,Evappft1,Thetampft1, &
         fAPARpft1,VODpft1,COS_fluxpft1, fei_leafpft1, laipft1
    real(r8),dimension(nlon*nlat,PFT)     :: GPPpft2,LHpft2,SHpft2,Transpft2,Evappft2,Thetampft2,&
         fAPARpft2,VODpft2,COS_fluxpft2, fei_leafpft2, laipft2
    real(r8),dimension(nlon,nlat,PFT)     :: GPPpft3,LHpft3,SHpft3,Transpft3,Evappft3,Thetampft3, &
         fAPARpft3,VODpft3,COS_fluxpft3,fei_leafpft3, laipft3

    real(r8) :: lon(nlon),lat(nlat)

    integer   :: ierr
    integer   :: ncid,dimid_lon,dimid_lat,dimid_time,dimid_PFT,varid
    character(len=8)    :: datestr
    character(len=255)  :: fln1,fln2,name,unit
    integer   :: nt,status
    integer   :: i

    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_gatherv(NEP9(1),npoints,mpi_real8,NEP1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(GPP9(1),npoints,mpi_real8,GPP1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(SIF9(1),npoints,mpi_real8,SIF1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(SIF9_sat(1),npoints,mpi_real8,SIF_sat1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(NPP9(1),npoints,mpi_real8,NPP1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(fAPAR9(1),npoints,mpi_real8,fAPAR1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(VOD9(1),npoints,mpi_real8,VOD1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(COS_flux9(1),npoints,mpi_real8,COS_flux1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(temp9(1),npoints,mpi_real8,temp1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Wind9(1),npoints,mpi_real8,Wind1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Rh9(1),npoints,mpi_real8,Rh1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Rain9(1),npoints,mpi_real8,Rain1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Snow9(1),npoints,mpi_real8,Snow1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Swdr9(1),npoints,mpi_real8,Swdr1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Swdf9(1),npoints,mpi_real8,Swdf1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(lai9(1),npoints,mpi_real8,lai1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(LH9(1),npoints,mpi_real8,LH1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(SH9(1),npoints,mpi_real8,SH1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Trans9(1),npoints,mpi_real8,Trans1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Evap9(1),npoints,mpi_real8,Evap1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Thetam9(1),npoints,mpi_real8,Thetam1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)


    !do i = 1,PFT
    !   call mpi_gatherv(GPPpft9(1,i),npoints,mpi_real8,GPPpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(Evappft9(1,i),npoints,mpi_real8,Evappft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(Thetampft9(1,i),npoints,mpi_real8,Thetampft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(fAPARpft9(1,i),npoints,mpi_real8,fAPARpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(VODpft9(1,i),npoints,mpi_real8,VODpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(COS_fluxpft9(1,i),npoints,mpi_real8,COS_fluxpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !end do
    !call mpi_barrier(mpi_comm_world,ierr)
    NEP1 = NEP9
    GPP1 = GPP9
    SIF1 = SIF9
    SIF_sat1 = SIF9_sat
    NPP1 = NPP9
    fAPAR1 = fAPAR9
    VOD1 = VOD9
    COS_flux1 = COS_flux9
    temp1 = temp9
    Wind1 = Wind9
    Rh1 = Rh9
    Rain1 = Rain9
    Snow1 = Snow9
    Swdr1 = Swdr9
    Swdf1 = Swdf9
    lai1 = lai9
    LH1 = LH9
    SH1 = SH9
    Trans1 = Trans9
    Evap1 = Evap9
    Thetam1 = Thetam9
    GPPpft1 = GPPpft9
    Evappft1 = Evappft9
    Thetampft1 = Thetampft9
    fAPARpft1 = fAPARpft9
    VODpft1 = VODpft9
    COS_fluxpft1 = COS_fluxpft9
    fei_leaf1 = fei_leaf9 ! 2023/07/20
    fei_leafpft1 = fei_leafpft9
    laipft1 =laipft9
    !if(myid ==0) then
    NEP2   = 0.
    GPP2   = 0.
    NPP2   = 0.

    GPPpft2    = 0.
    LHpft2    = 0.
    SHpft2    = 0.
    Transpft2    = 0.
    Evappft2    = 0.
    Thetampft2    = 0.
    fAPARpft2   = 0.
    VODpft2     = 0.
    COS_fluxpft2 = 0.
    temp2      = 0.
    Wind2      = 0.
    Rh2        = 0.
    Rain2      = 0.
    Snow2      = 0.
    Swdr2      = 0.
    Swdf2      = 0.
    lai2       = 0.
    LH2       = 0.
    SH2       = 0.
    Trans2       = 0.
    Evap2       = 0.
    Thetam2    = 0.
    SIF2       = 0.
    SIF_sat2   = 0.
    fAPAR2    = 0.
    VOD2      = 0.
    COS_flux2 = 0.
    fei_leaf2 = 0. ! 2023/07/20
    fei_leafpft2 = 0.
    laipft2 = 0.

    SIF2(mapping)  = SIF1
    SIF_sat2(mapping) = SIF_sat1
    NEP2(mapping)  = NEP1
    GPP2(mapping)  = GPP1
    NPP2(mapping)  = NPP1

    GPPpft2(mapping,:)   = GPPpft1
    LHpft2(mapping,:)   = LHpft1
    SHpft2(mapping,:)   = SHpft1
    Transpft2(mapping,:)   = Transpft1
    Evappft2(mapping,:)   = Evappft1
    Thetampft2(mapping,:)   = Thetampft1
    fAPARpft2(mapping,:)   = fAPARpft1
    VODpft2(mapping,:)     = VODpft1
    COS_fluxpft2(mapping,:) = COS_fluxpft1
    temp2(mapping)       = temp1
    Wind2(mapping)       = Wind1
    Rh2(mapping)         = Rh1
    Rain2(mapping)       = Rain1
    Snow2(mapping)       = Snow1
    Swdr2(mapping)       = Swdr1
    Swdf2(mapping)       = Swdf1
    lai2(mapping)        = lai1
    LH2(mapping)        = LH1
    SH2(mapping)        = SH1
    Trans2(mapping)        = Trans1
    Evap2(mapping)        = Evap1
    Thetam2(mapping)     = Thetam1
    fAPAR2(mapping)     = fAPAR1
    VOD2(mapping)       = VOD1
    COS_flux2(mapping)  = COS_flux1
    ! 2023/07/20
    fei_leaf2(mapping) = fei_leaf1
    fei_leafpft2(mapping,:) = fei_leafpft1
    laipft2(mapping,:) = laipft1

    NEP3    = reshape(NEP2,(/nlon,nlat/))
    GPP3    = reshape(GPP2,(/nlon,nlat/))
    NPP3    = reshape(NPP2,(/nlon,nlat/))
    SIF3    = reshape(SIF2,(/nlon,nlat/))
    SIF_sat3 = reshape(SIF_sat2,(/nlon,nlat/))

    GPPpft3 = reshape(GPPpft2,(/nlon,nlat,PFT/))
    LHpft3 = reshape(LHpft2,(/nlon,nlat,PFT/))
    SHpft3 = reshape(SHpft2,(/nlon,nlat,PFT/))
    Transpft3 = reshape(Transpft2,(/nlon,nlat,PFT/))
    Evappft3 = reshape(Evappft2,(/nlon,nlat,PFT/))
    Thetampft3 = reshape(Thetampft2,(/nlon,nlat,PFT/))
    fAPARpft3 = reshape(fAPARpft2,(/nlon,nlat,PFT/))
    VODpft3 = reshape(VODpft2,(/nlon,nlat,PFT/))
    COS_fluxpft3 = reshape(COS_fluxpft2,(/nlon,nlat,PFT/))

    temp3   = reshape(temp2,(/nlon,nlat/))
    Wind3   = reshape(Wind2,(/nlon,nlat/))
    Rh3     = reshape(Rh2,(/nlon,nlat/))
    Rain3   = reshape(Rain2,(/nlon,nlat/))
    Snow3   = reshape(Snow2,(/nlon,nlat/))
    Swdr3   = reshape(Swdr2,(/nlon,nlat/))
    Swdf3   = reshape(Swdf2,(/nlon,nlat/))
    lai3    = reshape(lai2,(/nlon,nlat/))
    LH3    = reshape(LH2,(/nlon,nlat/))
    SH3    = reshape(SH2,(/nlon,nlat/))
    Trans3    = reshape(Trans2,(/nlon,nlat/))
    Evap3    = reshape(Evap2,(/nlon,nlat/))
    Thetam3    = reshape(Thetam2,(/nlon,nlat/))
    fAPAR3    = reshape(fAPAR2,(/nlon,nlat/))
    VOD3    = reshape(VOD2,(/nlon,nlat/))
    COS_flux3    = reshape(COS_flux2,(/nlon,nlat/))
    ! 2023/07/20
    fei_leaf3 = reshape(fei_leaf2, (/nlon,nlat/))
    fei_leafpft3 = reshape(fei_leafpft2, (/nlon,nlat,PFT/))
    laipft3 = reshape(laipft2, (/nlon,nlat,PFT/))
    ! LWP is usually negative 2023/10/23
    ! 1 m H2o = 0.0098 MPa or 1/101 MPa
    fei_leaf3 = fei_leaf3/101. * (-1.)
    fei_leafpft3 = fei_leafpft3/101. * (-1.)
    !--iLab::yy,mm,dd,tod are arguments now
    ! call get_prev_date(yy,mm,dd,tod)
    if(nhtfrq <0) then
       write(datestr,"(i8)") yy*10000+mm*100+dd
       nt   = (tod/3600+1)/(-nhtfrq)
    else if(nhtfrq ==0) then
       write(datestr,"(i6)") yy*100+mm
       nt   = 1
    end if

    write(*,*) "Writing out simulation file now!"
    fln1  = trim(beps_out_dir)//"beps_global_"//trim(datestr)//".nc"
    status =  nf90_open(fln1,nf90_write,ncid)
    if(status .ne. nf90_noerr) then

       call check(nf90_create(fln1,nf90_share,ncid))
       call check(nf90_def_dim(ncid,"lon",nlon,dimid_lon))
       call check(nf90_def_dim(ncid,"lat",nlat,dimid_lat))
       call check(nf90_def_dim(ncid,"PFT",PFT,dimid_PFT))
       call check(nf90_def_dim(ncid,"time",nf90_unlimited,dimid_time))

       call check(nf90_def_var(ncid,"time",nf90_double,(/dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"calendar","Gregorian"))

       call check(nf90_def_var(ncid,"lon",nf90_double,(/dimid_lon/),varid))
       call check(nf90_put_att(ncid,varid,"units","degree_east"))
       call check(nf90_put_att(ncid,varid,"long_name","longitude"))
       call check(nf90_put_att(ncid,varid,"axis","X"))

       call check(nf90_def_var(ncid,"lat",nf90_double,(/dimid_lat/),varid))
       call check(nf90_put_att(ncid,varid,"units","degree_north"))
       call check(nf90_put_att(ncid,varid,"long_name","latitude"))
       call check(nf90_put_att(ncid,varid,"axis","Y"))

       call check(nf90_def_var(ncid,"NEP",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Net Ecosystem Productivity"))

       call check(nf90_def_var(ncid,"GPP",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity"))

       call check(nf90_def_var(ncid,"VOD",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","-"))
       call check(nf90_put_att(ncid,varid,"name","Vegetation Optical Depth"))

       call check(nf90_def_var(ncid,"fAPAR",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","-"))
       call check(nf90_put_att(ncid,varid,"name","Fraction of Absorbed Photosynthetically Active Radiation"))

       call check(nf90_def_var(ncid,"COS_flux",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","mol/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","COS flux for soil and plant"))

       !   call check(nf90_def_var(ncid,"SIF_sat",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       !   call check(nf90_put_att(ncid,varid,"units","mW/m2/nm/sr"))
       !   call check(nf90_put_att(ncid,varid,"name","solar-induced SIF over the OCO2 pass-time"))

       call check(nf90_def_var(ncid,"SIF",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","mW/m2/nm/sr"))
       call check(nf90_put_att(ncid,varid,"name","solar-induced SIF"))

       call check(nf90_def_var(ncid,"GPPpft",nf90_double,(/dimid_lon,dimid_lat,dimid_PFT,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity"))

       call check(nf90_put_att(ncid,NF90_GLOBAL,"model","Beps runs"))
       call check(nf90_put_att(ncid,NF90_GLOBAL,"institution","Nanjing University"))

       call check(nf90_def_var(ncid,"Thetam",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"name","Surface soil moisture"))

       call check(nf90_def_var(ncid,"LH",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"name","Latent heat flux"))

       call check(nf90_def_var(ncid,"SH",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"name","Sensible heat flux"))

       call check(nf90_def_var(ncid,"Trans",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Transpiration"))

       call check(nf90_def_var(ncid,"Evap",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Evaporation"))

       call check(nf90_def_var(ncid,"Evappft",nf90_double,(/dimid_lon,dimid_lat,dimid_PFT,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Evaporation"))

       call check(nf90_def_var(ncid,"Thetampft",nf90_double,(/dimid_lon,dimid_lat,dimid_PFT,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"name","Surface soil moisture"))
       ! 2023/07/20
       call check(nf90_def_var(ncid,"fei_leaf",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","MPa"))
       call check(nf90_put_att(ncid,varid,"name","Leaf water potential"))

       call check(nf90_def_var(ncid,"LAI",nf90_double,(/dimid_lon,dimid_lat,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m2/m2"))
       call check(nf90_put_att(ncid,varid,"name","Leaf area index"))

       call check(nf90_def_var(ncid,"fei_leafpft",nf90_double,(/dimid_lon,dimid_lat,dimid_PFT,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","MPa"))
       call check(nf90_put_att(ncid,varid,"name","Leaf water potential"))

       call check(nf90_enddef(ncid))
    end if
    !! For temporary output  , Should be improved later @J.Wang
    do i = 1,nlat
       lat(i) = -89.5+i-1
    end do
    call check(nf90_inq_varid(ncid,"lat",varid))
    call check(nf90_put_var(ncid,varid,lat))

    do i = 1,nlon
       lon(i)  = 0.5+i-1.
    end do
    call check(nf90_inq_varid(ncid,"lon",varid))
    call check(nf90_put_var(ncid,varid,lon))

    !    call check(nf90_inq_varid(ncid,"time",varid))
    !    call check(nf90_put_var(ncid,varid,nt,start=(/nt/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"NEP",varid))
    call check(nf90_put_var(ncid,varid,NEP3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"GPP",varid))
    call check(nf90_put_var(ncid,varid,GPP3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"VOD",varid))
    call check(nf90_put_var(ncid,varid,VOD3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"fAPAR",varid))
    call check(nf90_put_var(ncid,varid,fAPAR3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"COS_flux",varid))
    call check(nf90_put_var(ncid,varid,COS_flux3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    !    call check(nf90_inq_varid(ncid,"SIF_sat",varid))
    !    call check(nf90_put_var(ncid,varid,SIF_sat3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"SIF",varid))
    call check(nf90_put_var(ncid,varid,SIF3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"GPPpft",varid))
    call check(nf90_put_var(ncid,varid,GPPpft3,start=(/1,1,1,nt/),count=(/nlon,nlat,PFT,1/)))
    call check(nf90_inq_varid(ncid,"Thetam",varid))
    call check(nf90_put_var(ncid,varid,Thetam3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"LH",varid))
    call check(nf90_put_var(ncid,varid,LH3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"SH",varid))
    call check(nf90_put_var(ncid,varid,SH3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"Trans",varid))
    call check(nf90_put_var(ncid,varid,Trans3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"Evap",varid))
    call check(nf90_put_var(ncid,varid,Evap3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"Thetampft",varid))
    call check(nf90_put_var(ncid,varid,Thetampft3,start=(/1,1,1,nt/),count=(/nlon,nlat,PFT,1/)))
    call check(nf90_inq_varid(ncid,"Evappft",varid))
    call check(nf90_put_var(ncid,varid,Evappft3,start=(/1,1,1,nt/),count=(/nlon,nlat,PFT,1/)))
    ! 2023/07/20
    call check(nf90_inq_varid(ncid,"fei_leaf",varid))
    call check(nf90_put_var(ncid,varid,fei_leaf3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"LAI",varid))
    call check(nf90_put_var(ncid,varid,lai3,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_inq_varid(ncid,"fei_leafpft",varid))
    call check(nf90_put_var(ncid,varid,fei_leafpft3,start=(/1,1,1,nt/),count=(/nlon,nlat,PFT,1/)))
    call check(nf90_close(ncid))
    !end if
    !call mpi_barrier(mpi_comm_world,ierr)

  end subroutine write_output_global


  subroutine write_output_site(yy,mm,dd,tod,ref_date,secs_since_ref)
   !subroutine write_output_site(yy,mm,dd,tod,ref_date,secs_since_ref,p)
    use netcdf
    implicit none
    !--iLab::yy,mm,dd,tod turned into arguments
    integer, intent(in) :: yy,mm,dd,tod
    !integer, intent(in) :: yy,mm,dd,tod,p
    character(len=*), intent(in) :: ref_date
    real(r8), intent(in) :: secs_since_ref
    character(len=*), parameter :: sub = 'write_output_site'
    real(r8),dimension(nlp)        :: NEP1,GPP1,SIF1,SIF_sat1,NPP1,temp1,Wind1,Rh1,Rain1,Snow1,Swdr1,Swdf1,lai1,LH1,SH1, &
         Trans1,Evap1,Thetam1,fAPAR1,VOD1,COS_flux1, fei_leaf1, ETa1,Qupt1, PWS1
    real(r8),dimension(nlp,PFT)    :: GPPpft1,LHpft1,SHpft1,Transpft1,Evappft1,Thetampft1,fAPARpft1,VODpft1,COS_fluxpft1, &
         fei_leafpft1,laipft1
     ! 2023/10/30
    real(r8),dimension(nlp)  :: Thetam2,Thetam3,Thetam4,Thetam5, SWP1,SWP2,SWP3,SWP4,SWP5
    real(r8),dimension(nlp)  :: TS1,TS2,TS3,TS4,TS5
    real(r8),dimension(nlp)  :: pond,rain_g
    ! 2024/01/08
    real(r8),dimension(nlp)  :: f_soilwater1, f_feileaf1, f_Tleaf1
    real(r8),dimension(nlp)  :: LHa1
    ! 2024/03/30
    real(r8),dimension(nlp)  :: GPP_o_sunlit1,GPP_o_shaded1,GPP_u_sunlit1,GPP_u_shaded1
    real(r8),dimension(nlp)  :: TR_o_sunlit1,TR_o_shaded1,TR_u_sunlit1,TR_u_shaded1
    real(r8),dimension(nlp)  :: Eil1,Evap_soil1,Evap_SW1,EiS1,Evap_SS1

    integer   :: ierr
    integer   :: ncid,dimid_site,dimid_time,dimid_PFT,varid
    integer   :: nsite(nlp)
    !character(len=8)    :: datestr,ppp
    character(len=8)    :: datestr
    character(len=255)  :: fln1,fln2,name,unit
    integer   :: nt,status
    integer   :: i
    !-- iLab::reduce amount of terminal output (to be reactivated on purpose)
    logical :: ldebug = .False.
    !-- iLab::added for consistent initialisation of NetCDF variables
    real(r8), parameter :: fill_value = -99999._r8

    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_gatherv(NEP9(1),npoints,mpi_real8,NEP1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(GPP9(1),npoints,mpi_real8,GPP1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(SIF9(1),npoints,mpi_real8,SIF1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(SIF9_sat(1),npoints,mpi_real8,SIF_sat1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(NPP9(1),npoints,mpi_real8,NPP1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(fAPAR9(1),npoints,mpi_real8,fAPAR1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(VOD9(1),npoints,mpi_real8,VOD1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(COS_flux9(1),npoints,mpi_real8,COS_flux1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(temp9(1),npoints,mpi_real8,temp1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Wind9(1),npoints,mpi_real8,Wind1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Rh9(1),npoints,mpi_real8,Rh1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Rain9(1),npoints,mpi_real8,Rain1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Snow9(1),npoints,mpi_real8,Snow1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Swdr9(1),npoints,mpi_real8,Swdr1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Swdf9(1),npoints,mpi_real8,Swdf1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(lai9(1),npoints,mpi_real8,lai1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(LH9(1),npoints,mpi_real8,LH1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(SH9(1),npoints,mpi_real8,SH1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Trans9(1),npoints,mpi_real8,Trans1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Evap9(1),npoints,mpi_real8,Evap1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_gatherv(Thetam9(1),npoints,mpi_real8,Thetam1(1),dp,sp,mpi_real8,0,mpi_comm_world,ierr)


    !do i = 1,PFT
    !   call mpi_gatherv(GPPpft9(1,i),npoints,mpi_real8,GPPpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(Evappft9(1,i),npoints,mpi_real8,Evappft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(Thetampft9(1,i),npoints,mpi_real8,Thetampft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(fAPARpft9(1,i),npoints,mpi_real8,fAPARpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(VODpft9(1,i),npoints,mpi_real8,VODpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_gatherv(COS_fluxpft9(1,i),npoints,mpi_real8,COS_fluxpft1(1,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
    !end do
    !call mpi_barrier(mpi_comm_world,ierr)
    NEP1 = NEP9
    GPP1 = GPP9
    SIF1 = SIF9
    SIF_sat1 = SIF9_sat
    NPP1 = NPP9
    fAPAR1 = fAPAR9
    VOD1 = VOD9
    COS_flux1 = COS_flux9
    temp1 = temp9
    Wind1 = Wind9
    Rh1 = Rh9
    Rain1 = Rain9
    Snow1 = Snow9
    Swdr1 = Swdr9
    Swdf1 = Swdf9
    lai1 = lai9
    LH1 = LH9
    SH1 = SH9
    Trans1 = Trans9
    Evap1 = Evap9
    Thetam1 = Thetam9
    GPPpft1 = GPPpft9
    Evappft1 = Evappft9
    Thetampft1 = Thetampft9
    fAPARpft1 = fAPARpft9
    VODpft1 = VODpft9
    COS_fluxpft1 = COS_fluxpft9
    ! 2023/07/20
    fei_leaf1 = fei_leaf9
    fei_leafpft1 = fei_leafpft9
    laipft1 = laipft9
    ! LWP is usually negative
    ! 1 m H2o = 0.0098 MPa or 1/101 MPa
    fei_leaf1 = fei_leaf1/101. * (-1.)
    fei_leafpft1 = fei_leafpft1/101. * (-1.)
    ETa1 = ETa9
    Qupt1 = Qupt9
    PWS1 = PWS9
    !if(myid ==0) then
    ! 2023/10/30
    Thetam2=Thetam2_9
    Thetam3=Thetam3_9
    Thetam4=Thetam4_9
    Thetam5=Thetam5_9
    SWP1=SWP1_9/101.*(-1.)
    SWP2=SWP2_9/101.*(-1.)
    SWP3=SWP3_9/101.*(-1.)
    SWP4=SWP4_9/101.*(-1.)
    SWP5=SWP5_9/101.*(-1.)

    TS1=TS1_9
    TS2=TS2_9
    TS3=TS3_9
    TS4=TS4_9
    TS5=TS5_9
    !pond=PondWater_9
    !rain_g=Rain_g_9
    f_soilwater1 = f_soilwater9
    f_feileaf1 = f_feileaf9
    !f_Tleaf1 = f_Tleaf9

    ! 2024/03/12
    LHa1 = LHa9

    !2024/03/30
    GPP_o_sunlit1 = GPP_o_sunlit9
    GPP_o_shaded1 = GPP_o_shaded9
    GPP_u_sunlit1 = GPP_u_sunlit9
    GPP_u_shaded1 = GPP_u_shaded9

    TR_o_sunlit1 = TR_o_sunlit9
    TR_o_shaded1 = TR_o_shaded9
    TR_u_sunlit1 = TR_u_sunlit9
    TR_u_shaded1 = TR_u_shaded9

    Eil1 = Eil9
    Evap_soil1 = Evap_soil9
    Evap_SW1 = Evap_SW9
    EiS1 = EiS9
    Evap_SS1 = Evap_SS9
    !write (*,*),'EiS1 =',EiS1

    !--iLab::yy,mm,dd,tod are arguments now
    ! call get_prev_date(yy,mm,dd,tod)
    if(nhtfrq <0) then
       !datestr="US_Ne1"
       datestr='US-MOz' ! site name as output NC file
       !write(datestr,"(i8)") yy
       !write(datestr,"(i8)") yy*10000+mm*100+dd
       nt   = (secs_since_ref/3600+1)/(-nhtfrq)
       !nt   = (tod/3600+1)/(-nhtfrq)
       ! !-- iLab::seconds elapsed since reference time (added for time-variable output)
       ! call timemgr_diff_secs(yy_ref*10000+mm_ref*100+dd_ref, tod_ref, yy*10000+mm*100+dd, tod,&
       !      secs_since_ref(1))
    else if(nhtfrq ==0) then
       write(datestr,"(i6)") yy*100+mm
       nt   = 1
    end if

    !-- iLab::make logging output depend on flag
    if(ldebug) then
       write(*,*) "Writing out simulation file now!"
    endif
    !write(ppp,"(i8)") p
    !fln1  = trim(beps_out_dir)//"beps_site_"//trim(adjustl(ppp))//"_"//trim(datestr)//".nc"
    fln1  = trim(beps_out_dir)//"beps_site_"//trim(datestr)//".nc"
    status =  nf90_open(fln1,nf90_write,ncid)
    if(status .ne. nf90_noerr) then

       call check(nf90_create(fln1,nf90_share,ncid))
       call check(nf90_def_dim(ncid,"nsite",nlp,dimid_site))
       call check(nf90_def_dim(ncid,"PFT",PFT,dimid_PFT))
       call check(nf90_def_dim(ncid,"time",nf90_unlimited,dimid_time))

       call check(nf90_def_var(ncid,"time",nf90_double,(/dimid_time/),varid))
       call check(nf90_put_att(ncid, varid, "long_name", "time"))
       !-- iLab::added for proper time-variable (hourly output only)
       if( nhtfrq<0 ) then
           call check(nf90_put_att(ncid, varid, "units", "seconds since "//ref_date))
          !call check(nf90_put_att(ncid, varid, "units", "hours since "//ref_date))
       endif
       call check(nf90_put_att(ncid,varid,"calendar","Gregorian"))

       !-- iLab::added NetCDF output initialisation with prescribed _FillValue
       !         *UPDATE*: disabled for now, since maybe not compatible on
       !                   platform which Mousong is using.
       call check(nf90_def_var(ncid,"nsite",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","-"))
       call check(nf90_put_att(ncid,varid,"long_name","site_number"))
       call check(nf90_put_att(ncid,varid,"axis","X"))

       call check(nf90_def_var(ncid,"Net_Rad",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"name","Net Radiation"))

       !call check(nf90_def_var(ncid,"NEP",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","Net Ecosystem Productivity"))

       call check(nf90_def_var(ncid,"GPP",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity"))

       !call check(nf90_def_var(ncid,"VOD",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"units","-"))
       !call check(nf90_put_att(ncid,varid,"name","Vegetation Optical Depth"))

       !call check(nf90_def_var(ncid,"fAPAR",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"units","-"))
       !call check(nf90_put_att(ncid,varid,"name","Fraction of Absorbed Photosynthetically Active Radiation"))

       !call check(nf90_def_var(ncid,"COS_flux",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"units","pmol/m2/s"))
       !call check(nf90_put_att(ncid,varid,"name","COS flux for soil and plant"))

       !   call check(nf90_def_var(ncid,"SIF_sat",nf90_double,(/dimid_site,dimid_time/),varid))
       !   call check(nf90_put_att(ncid,varid,"units","mW/m2/nm/sr"))
       !   call check(nf90_put_att(ncid,varid,"name","solar-induced SIF over the OCO2 pass-time"))

       !call check(nf90_def_var(ncid,"SIF",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"units","mW/m2/nm/sr"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","solar-induced SIF"))

       !call check(nf90_def_var(ncid,"GPPpft",nf90_double,(/dimid_site,dimid_PFT,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity"))

       call check(nf90_put_att(ncid,NF90_GLOBAL,"model","Beps runs"))
       call check(nf90_put_att(ncid,NF90_GLOBAL,"institution","Nanjing University"))

       call check(nf90_def_var(ncid,"Thetam1",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil moisture in layer 1 (5 cm)"))

       call check(nf90_def_var(ncid,"LH",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Latent heat flux"))

       call check(nf90_def_var(ncid,"SH",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Sensible heat flux"))

       call check(nf90_def_var(ncid,"Trans",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Transpiration"))
       ! 2023/10/31
       call check(nf90_def_var(ncid,"ETa",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Transpiration from SPAC"))

       call check(nf90_def_var(ncid,"PWS",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","plant water storage"))

       call check(nf90_def_var(ncid,"Qupt",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Uptake soil water from SPAC"))

       call check(nf90_def_var(ncid,"Evap",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Evaporation from original Beps"))

       !call check(nf90_def_var(ncid,"Evappft",nf90_double,(/dimid_site,dimid_PFT,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"units","m/s"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","Evaporation"))

       !call check(nf90_def_var(ncid,"Thetampft",nf90_double,(/dimid_site,dimid_PFT,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       !call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","Surface soil moisture"))
       ! 2023/07/20
       call check(nf90_def_var(ncid,"fei_leaf",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","Mpa"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Leaf water potential"))

       !call check(nf90_def_var(ncid,"LAI",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_put_att(ncid,varid,"units","m2/m2"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","Leaf area index"))

       !call check(nf90_def_var(ncid,"fei_leafpft",nf90_double,(/dimid_site,dimid_PFT,dimid_time/),varid))
       !call check(nf90_put_att(ncid,varid,"units","MPa"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","Leaf water potential"))

       !
       call check(nf90_def_var(ncid,"Thetam2",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil moisture in layer 2 (15 cm) "))

       call check(nf90_def_var(ncid,"Thetam3",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil moisture in layer 3 (35 cm) "))

       call check(nf90_def_var(ncid,"Thetam4",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil moisture in layer 4 (75 cm) "))

       call check(nf90_def_var(ncid,"Thetam5",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil moisture in layer 5 (200 cm) "))

       call check(nf90_def_var(ncid,"SWP1",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","Mpa"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil water potential in layer 1 (5 cm)"))

       call check(nf90_def_var(ncid,"SWP2",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","Mpa"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil water potential in layer 2 (15 cm)"))

       call check(nf90_def_var(ncid,"SWP3",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","Mpa"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil water potential in layer 3 (35 cm)"))

       call check(nf90_def_var(ncid,"SWP4",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","Mpa"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil water potential in layer 4 (75 cm)"))

       call check(nf90_def_var(ncid,"SWP5",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","Mpa"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil water potential in layer 5 (200 cm)"))

       call check(nf90_def_var(ncid,"TS1",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","C"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil temperature in layer 1 (5 cm)"))

       call check(nf90_def_var(ncid,"TS2",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","C"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil temperature in layer 2 (15 cm)"))

       call check(nf90_def_var(ncid,"TS3",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","C"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil temperature in layer 3 (35 cm)"))

       call check(nf90_def_var(ncid,"TS4",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","C"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil temperature in layer 4 (75 cm)"))

       call check(nf90_def_var(ncid,"TS5",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","C"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil temperature in layer 5 (200 cm)"))

       !call check(nf90_def_var(ncid,"PondWater",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_put_att(ncid,varid,"units","m"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))

       !call check(nf90_def_var(ncid,"Rain_g",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_put_att(ncid,varid,"units","m/s"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))

       call check(nf90_def_var(ncid,"f_soilwater",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","dimensionless, [0-1]"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","water stress from soil moisture"))

       call check(nf90_def_var(ncid,"f_feileaf",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","dimensionless, [0-1]"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","water stress from leaf water potential"))

       !call check(nf90_def_var(ncid,"f_Tleaf",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_put_att(ncid,varid,"units","dimensionless, [0-1]"))
       !call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       !call check(nf90_put_att(ncid,varid,"name","temperature stress from leaf"))

       call check(nf90_def_var(ncid,"LHa",nf90_double,(/dimid_site,dimid_time/),varid))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Actual Latent heat flux when considering plant hydraulics"))

       call check(nf90_def_var(ncid,"GPP_o_sunlit",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity of overstory sunlit leaf"))

       call check(nf90_def_var(ncid,"GPP_o_shaded",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity of overstory shaded leaf"))

       call check(nf90_def_var(ncid,"GPP_u_sunlit",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity of understory sunlit leaf"))

       call check(nf90_def_var(ncid,"GPP_u_shaded",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity of understory shaded leaf"))

       call check(nf90_def_var(ncid,"TR_o_sunlit",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Transpiration of overstory sunlit leaf"))

       call check(nf90_def_var(ncid,"TR_o_shaded",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Transpiration of overstory shaded leaf"))

       call check(nf90_def_var(ncid,"TR_u_sunlit",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Transpiration of understory sunlit leaf"))

       call check(nf90_def_var(ncid,"TR_u_shaded",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Transpiration of understory shaded leaf"))

       call check(nf90_def_var(ncid,"Eil",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Interception of liquid water"))

       call check(nf90_def_var(ncid,"Evap_soil",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Evaporation of surface soil"))

       call check(nf90_def_var(ncid,"Evap_SW",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Evaporation of surface pond water"))

        call check(nf90_def_var(ncid,"EiS",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Evaporation of solid water (snow)"))

       call check(nf90_def_var(ncid,"Evap_SS",nf90_double,(/dimid_site,dimid_time/),varid))
       !call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Evaporation of surface snow"))

       call check(nf90_enddef(ncid))
    end if
    !! For temporary output  , Should be improved later @J.Wang
    do i = 1,nlp
       nsite(i) = i
    end do

    call check(nf90_inq_varid(ncid,"nsite",varid))
    call check(nf90_put_var(ncid,varid,nsite))

    !-- iLab::added writing of time-values (hourly output only)
    if( nhtfrq <0 ) then
       call check(nf90_inq_varid(ncid,"time",varid))
       call check(nf90_put_var(ncid,varid,(/secs_since_ref/),start=(/nt/),count=(/1/)))
    endif
    !    call check(nf90_inq_varid(ncid,"time",varid))
    !    call check(nf90_put_var(ncid,varid,nt,start=(/nt/),count=(/1/)))

    !call check(nf90_inq_varid(ncid,"NEP",varid))
    !call check(nf90_put_var(ncid,varid,NEP1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"GPP",varid))
    call check(nf90_put_var(ncid,varid,GPP1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"Net_Rad",varid))
    call check(nf90_put_var(ncid,varid,Rh1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"VOD",varid))
    !call check(nf90_put_var(ncid,varid,VOD1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"fAPAR",varid))
    !call check(nf90_put_var(ncid,varid,fAPAR1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"COS_flux",varid))
    !call check(nf90_put_var(ncid,varid,COS_flux1,start=(/1,nt/),count=(/nlp,1/)))
    !    call check(nf90_inq_varid(ncid,"SIF_sat",varid))
    !    call check(nf90_put_var(ncid,varid,SIF_sat1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"SIF",varid))
    !call check(nf90_put_var(ncid,varid,SIF1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"GPPpft",varid))
    !call check(nf90_put_var(ncid,varid,GPPpft1,start=(/1,1,nt/),count=(/nlp,PFT,1/)))
    call check(nf90_inq_varid(ncid,"Thetam1",varid))
    call check(nf90_put_var(ncid,varid,Thetam1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"LH",varid))
    call check(nf90_put_var(ncid,varid,LH1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"SH",varid))
    call check(nf90_put_var(ncid,varid,SH1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Trans",varid))
    call check(nf90_put_var(ncid,varid,Trans1,start=(/1,nt/),count=(/nlp,1/)))
    ! 2023/10/31
    call check(nf90_inq_varid(ncid,"ETa",varid))
    call check(nf90_put_var(ncid,varid,ETa1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"PWS",varid))
    call check(nf90_put_var(ncid,varid,PWS1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Qupt",varid))
    call check(nf90_put_var(ncid,varid,Qupt1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"Evap",varid))
    call check(nf90_put_var(ncid,varid,Evap1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"Thetampft",varid))
    !call check(nf90_put_var(ncid,varid,Thetampft1,start=(/1,1,nt/),count=(/nlp,PFT,1/)))
    !call check(nf90_inq_varid(ncid,"Evappft",varid))
    !call check(nf90_put_var(ncid,varid,Evappft1,start=(/1,1,nt/),count=(/nlp,PFT,1/)))
    ! 2023/07/20
    call check(nf90_inq_varid(ncid,"fei_leaf",varid))
    call check(nf90_put_var(ncid,varid,fei_leaf1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"LAI",varid))
    !call check(nf90_put_var(ncid,varid,lai1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"fei_leafpft",varid))
    !call check(nf90_put_var(ncid,varid,fei_leafpft1,start=(/1,1,nt/),count=(/nlp,PFT,1/)))

    call check(nf90_inq_varid(ncid,"Thetam2",varid))
    call check(nf90_put_var(ncid,varid,Thetam2,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Thetam3",varid))
    call check(nf90_put_var(ncid,varid,Thetam3,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Thetam4",varid))
    call check(nf90_put_var(ncid,varid,Thetam4,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Thetam5",varid))
    call check(nf90_put_var(ncid,varid,Thetam5,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"SWP1",varid))
    call check(nf90_put_var(ncid,varid,SWP1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"SWP2",varid))
    call check(nf90_put_var(ncid,varid,SWP2,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"SWP3",varid))
    call check(nf90_put_var(ncid,varid,SWP3,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"SWP4",varid))
    call check(nf90_put_var(ncid,varid,SWP4,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"SWP5",varid))
    call check(nf90_put_var(ncid,varid,SWP5,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"TS1",varid))
    call check(nf90_put_var(ncid,varid,TS1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"TS2",varid))
    call check(nf90_put_var(ncid,varid,TS2,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"TS3",varid))
    call check(nf90_put_var(ncid,varid,TS3,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"TS4",varid))
    call check(nf90_put_var(ncid,varid,TS4,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"TS5",varid))
    call check(nf90_put_var(ncid,varid,TS5,start=(/1,nt/),count=(/nlp,1/)))

    !call check(nf90_inq_varid(ncid,"PondWater",varid))
    !call check(nf90_put_var(ncid,varid,pond,start=(/1,nt/),count=(/nlp,1/)))

    !call check(nf90_inq_varid(ncid,"Rain_g",varid))
    !call check(nf90_put_var(ncid,varid,rain_g,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"f_soilwater",varid))
    call check(nf90_put_var(ncid,varid,f_soilwater1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"f_feileaf",varid))
    call check(nf90_put_var(ncid,varid,f_feileaf1,start=(/1,nt/),count=(/nlp,1/)))
    !call check(nf90_inq_varid(ncid,"f_Tleaf",varid))
    !call check(nf90_put_var(ncid,varid,f_Tleaf1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"LHa",varid))
    call check(nf90_put_var(ncid,varid,LHa1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"GPP_o_sunlit",varid))
    call check(nf90_put_var(ncid,varid,GPP_o_sunlit1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"GPP_o_shaded",varid))
    call check(nf90_put_var(ncid,varid,GPP_o_shaded1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"GPP_u_sunlit",varid))
    call check(nf90_put_var(ncid,varid,GPP_u_sunlit1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"GPP_u_shaded",varid))
    call check(nf90_put_var(ncid,varid,GPP_u_shaded1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"TR_o_sunlit",varid))
    call check(nf90_put_var(ncid,varid,TR_o_sunlit1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"TR_o_shaded",varid))
    call check(nf90_put_var(ncid,varid,TR_o_shaded1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"TR_u_sunlit",varid))
    call check(nf90_put_var(ncid,varid,TR_u_sunlit1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"TR_u_shaded",varid))
    call check(nf90_put_var(ncid,varid,TR_u_shaded1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_inq_varid(ncid,"Eil",varid))
    call check(nf90_put_var(ncid,varid,Eil1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Evap_soil",varid))
    call check(nf90_put_var(ncid,varid,Evap_soil1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Evap_SW",varid))
    call check(nf90_put_var(ncid,varid,Evap_SW1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"EiS",varid))
    call check(nf90_put_var(ncid,varid,EiS1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Evap_SS",varid))
    call check(nf90_put_var(ncid,varid,Evap_SS1,start=(/1,nt/),count=(/nlp,1/)))

    call check(nf90_close(ncid))
    !end if
    !call mpi_barrier(mpi_comm_world,ierr)

  end subroutine write_output_site


end module outputMod
