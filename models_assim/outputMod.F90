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
          !--iLab::yr,mon,day,tod now provided as arguments
          ! call get_prev_date(yr, mon, day, tod)
          !!              write(*,*) "write out data on ",yr,mon,day
          SIFpft9_sat   = SIFpft9_sat/day
          SIF9_sat      = SIF9_sat/day

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

  subroutine write_output_global(yy,mm,dd,tod)
    use netcdf
    implicit none
    !--iLab::yy,mm,dd,tod turned to arguments
    integer, intent(in) :: yy,mm,dd,tod
    real(r8),dimension(nlp)        :: NEP1,GPP1,SIF1,SIF_sat1,NPP1,temp1,Wind1,Rh1,Rain1,Snow1,Swdr1,Swdf1,lai1,LH1,SH1, &
         Trans1,Evap1,Thetam1,fAPAR1,VOD1,COS_flux1
    real(r8),dimension(nlon*nlat)  :: NEP2,GPP2,SIF2,SIF_sat2,NPP2,temp2,Wind2,Rh2,Rain2,Snow2,Swdr2,Swdf2,lai2,LH2,SH2, &
         Trans2,Evap2,Thetam2,fAPAR2,VOD2,COS_flux2
    real(r8),dimension(nlon,nlat)  :: NEP3,GPP3,SIF3,SIF_sat3,NPP3,temp3,Wind3,Rh3,Rain3,Snow3,Swdr3,Swdf3,lai3,LH3,SH3, &
         Trans3,Evap3,Thetam3,fAPAR3,VOD3,COS_flux3
    real(r8),dimension(nlp,PFT)           :: GPPpft1,LHpft1,SHpft1,Transpft1,Evappft1,Thetampft1,fAPARpft1,VODpft1,COS_fluxpft1
    real(r8),dimension(nlon*nlat,PFT)     :: GPPpft2,LHpft2,SHpft2,Transpft2,Evappft2,Thetampft2,fAPARpft2,VODpft2,COS_fluxpft2
    real(r8),dimension(nlon,nlat,PFT)     :: GPPpft3,LHpft3,SHpft3,Transpft3,Evappft3,Thetampft3,fAPARpft3,VODpft3,COS_fluxpft3

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
         Trans1,Evap1,Thetam1,fAPAR1,VOD1,COS_flux1
    real(r8),dimension(nlp,PFT)    :: GPPpft1,LHpft1,SHpft1,Transpft1,Evappft1,Thetampft1,fAPARpft1,VODpft1,COS_fluxpft1

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

    !if(myid ==0) then

    !--iLab::yy,mm,dd,tod are arguments now
    ! call get_prev_date(yy,mm,dd,tod)
    if(nhtfrq <0) then
       datestr="US_Ne1"
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
       endif
       call check(nf90_put_att(ncid,varid,"calendar","Gregorian"))

       !-- iLab::added NetCDF output initialisation with prescribed _FillValue
       !         *UPDATE*: disabled for now, since maybe not compatible on
       !                   platform which Mousong is using.
       call check(nf90_def_var(ncid,"nsite",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","-"))
       call check(nf90_put_att(ncid,varid,"long_name","site_number"))
       call check(nf90_put_att(ncid,varid,"axis","X"))

       call check(nf90_def_var(ncid,"NEP",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Net Ecosystem Productivity"))

       call check(nf90_def_var(ncid,"GPP",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity"))

       call check(nf90_def_var(ncid,"VOD",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","-"))
       call check(nf90_put_att(ncid,varid,"name","Vegetation Optical Depth"))

       call check(nf90_def_var(ncid,"fAPAR",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","-"))
       call check(nf90_put_att(ncid,varid,"name","Fraction of Absorbed Photosynthetically Active Radiation"))

       call check(nf90_def_var(ncid,"COS_flux",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"units","pmol/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","COS flux for soil and plant"))

       !   call check(nf90_def_var(ncid,"SIF_sat",nf90_double,(/dimid_site,dimid_time/),varid))
       !   call check(nf90_put_att(ncid,varid,"units","mW/m2/nm/sr"))
       !   call check(nf90_put_att(ncid,varid,"name","solar-induced SIF over the OCO2 pass-time"))

       call check(nf90_def_var(ncid,"SIF",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","mW/m2/nm/sr"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","solar-induced SIF"))

       call check(nf90_def_var(ncid,"GPPpft",nf90_double,(/dimid_site,dimid_PFT,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity"))

       call check(nf90_put_att(ncid,NF90_GLOBAL,"model","Beps runs"))
       call check(nf90_put_att(ncid,NF90_GLOBAL,"institution","Nanjing University"))

       call check(nf90_def_var(ncid,"Thetam",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Surface soil moisture"))

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

       call check(nf90_def_var(ncid,"Evap",nf90_double,(/dimid_site,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Evaporation"))

       call check(nf90_def_var(ncid,"Evappft",nf90_double,(/dimid_site,dimid_PFT,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Evaporation"))

       call check(nf90_def_var(ncid,"Thetampft",nf90_double,(/dimid_site,dimid_PFT,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Surface soil moisture"))

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

    call check(nf90_inq_varid(ncid,"NEP",varid))
    call check(nf90_put_var(ncid,varid,NEP1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"GPP",varid))
    call check(nf90_put_var(ncid,varid,GPP1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"VOD",varid))
    call check(nf90_put_var(ncid,varid,VOD1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"fAPAR",varid))
    call check(nf90_put_var(ncid,varid,fAPAR1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"COS_flux",varid))
    call check(nf90_put_var(ncid,varid,COS_flux1,start=(/1,nt/),count=(/nlp,1/)))
    !    call check(nf90_inq_varid(ncid,"SIF_sat",varid))
    !    call check(nf90_put_var(ncid,varid,SIF_sat1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"SIF",varid))
    call check(nf90_put_var(ncid,varid,SIF1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"GPPpft",varid))
    call check(nf90_put_var(ncid,varid,GPPpft1,start=(/1,1,nt/),count=(/nlp,PFT,1/)))
    call check(nf90_inq_varid(ncid,"Thetam",varid))
    call check(nf90_put_var(ncid,varid,Thetam1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"LH",varid))
    call check(nf90_put_var(ncid,varid,LH1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"SH",varid))
    call check(nf90_put_var(ncid,varid,SH1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Trans",varid))
    call check(nf90_put_var(ncid,varid,Trans1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Evap",varid))
    call check(nf90_put_var(ncid,varid,Evap1,start=(/1,nt/),count=(/nlp,1/)))
    call check(nf90_inq_varid(ncid,"Thetampft",varid))
    call check(nf90_put_var(ncid,varid,Thetampft1,start=(/1,1,nt/),count=(/nlp,PFT,1/)))
    call check(nf90_inq_varid(ncid,"Evappft",varid))
    call check(nf90_put_var(ncid,varid,Evappft1,start=(/1,1,nt/),count=(/nlp,PFT,1/)))
    call check(nf90_close(ncid))
    !end if
    !call mpi_barrier(mpi_comm_world,ierr)

  end subroutine write_output_site

end module outputMod
