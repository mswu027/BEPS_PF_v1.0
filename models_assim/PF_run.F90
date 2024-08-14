module PF_run
  use shr_kind_mod,only:r8=>shr_kind_r8
  use bepstype
  use bepstypeInit
  use controlInput_mod
  !--iLab::avoid beps_time_manager, all temporal information now passed as actual arguments
  ! use beps_time_manager
  use beps_par
  use beps_con
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

  subroutine Init_output_PF()
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
  end subroutine Init_output_PF

    !! average variables according to user's definition
  subroutine PF_weight(p,yr, mon, day, tod,kount,is_end_curr_month,is_end_week)
    implicit none
    !-- iLab::turned yr,mon,day,tod to arguments and added the further arguments
    !integer, intent(in) :: yr,mon,day,tod
    !integer, intent(in) :: kount
    !logical, intent(in) :: is_end_curr_month
    !character(len=*), intent(in) :: ref_date
    !real(r8), intent(in) :: secs_since_ref
    !integer    :: ii,iii,p
    integer, intent(in) :: kount
    logical, intent(in) :: is_end_curr_month,is_end_week
    integer, intent(in) :: yr,mon,day,tod
    integer, intent(in) :: p
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
          NEP9   = NEP9/nstpd      !! average
          GPP9   = GPP9/nstpd*3600*1000 !kg C/m2/s -----> g C/m2/h
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
             call PF_weight_update(p)

             if (is_end_week) then
                call write_PF_para(yr, mon, day, tod,p)
             end if

          else
             call PF_weight_update(p)

             if (is_end_week) then
                call write_PF_para(yr, mon, day, tod,p)
             end if

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
    else if(nhtfrq ==0) then   !!monthly output !!!!particle fliter did not consider
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
             call PF_weight_update(p)
             if (is_end_week) then
                call write_PF_para(yr, mon, day, tod,p)
             end if
          else
             call PF_weight_update(p)
             if (is_end_week) then
                call write_PF_para(yr, mon, day, tod,p)
             end if
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

  end subroutine PF_weight

  subroutine PF_weight_update(p)!假设各种平均值为都是0，误差为1
    implicit none
    integer, intent(in)   :: p
    integer   :: mean,std

    character(len=*), parameter :: sub = 'PF_weight_update'

    real(r8) :: a
    real(r8) :: likehood  !似然函数为正态分布公式
    real(r8) :: sum

    mean=0
    std=1 ! should be smaller for vod
    sum=0

    !a=-(PF_obs%obs_GPP(1)-GPP9(1)-mean)**2/(2*std**2)
    a=-(PF_obs%obs_VOD(1)-VOD9(1)-mean)**2/(2*std**2)
    !write(*,*)"a=" , a

    if (a<-200.0) then
        a=-200.0
    else
        a=a
    end if

    likehood=(1/((2*PI)**(0.5))*std)*exp(a)
    write(*,*)"likehood=" , likehood

    PF%pfweightupdate(p,1)=PF%pfweight(p,1)*likehood
    write(*,*)"PF%pfweightupdate(p,1)=" , PF%pfweightupdate(p,1)

  end subroutine PF_weight_update

  subroutine PF_weight_update_resample (p,weight)
    real(r8), intent(in) :: weight
    integer, intent(in) :: p
    integer   :: mean,std
    real(r8) :: a
    real(r8) :: likehood  !似然函数为正态分布公式
    real(r8) :: sum
    character(len=*), parameter :: sub = 'PF_weight_update_resample'

    mean=0
    std=1
    sum=0
    !a=-(PF_obs%obs_GPP(1)-GPP9(1)-mean)**2/(2*std**2)
    a=-(PF_obs%obs_VOD(1)-VOD9(1)-mean)**2/(2*std**2)
    if (a<-200.0) then
        a=-200.0
    else
        a=a
    end if
    likehood=(1/((2*PI)**(0.5))*std)*exp(a)

    PF_resample%resample_weight_update(p)=weight*likehood


  end subroutine PF_weight_update_resample


  subroutine resample(inparticles,weights)
    implicit none
    real(r8), intent(in) :: inparticles(:,:),weights(:)
    integer   ::N,a,b,i,k,kk,j
    integer, dimension (:), allocatable :: cum_copies,num_copies
    real, dimension (:,:), allocatable :: outparticles
    real, dimension (:,:), allocatable :: kkk
    real(r8) :: h,wcp

    N=size(weights,dim=1)
    a=size(inparticles,dim=1)
    b=size(inparticles,dim=2)
    write(*,*)"N" , N
    write(*,*)"a" , a
    write(*,*)"b" , b

    allocate ( cum_copies(N) )
    allocate ( num_copies(N) )
    allocate ( outparticles(a,b) )


    wcp=0
    do i=1,N
         wcp=wcp+N*weights(i)
         cum_copies(i)=NINT(wcp)

    end do

    do i =1,N
        if (i==1) then
          num_copies(i)=cum_copies(i)
        else
          num_copies(i)=cum_copies(i)-cum_copies(i-1)
        end if
    end do

    k=1
    do i=1,N
        if (num_copies(i)>0) then
           if (num_copies(i)==1) then
              outparticles(k,:)=inparticles(i,:)
              k=k+1
           else if (num_copies(i)==2) then
              outparticles(k,:)=inparticles(i,:)
              outparticles(k+1,:)=inparticles(i,:)
              k=k+2
           else if ((num_copies(i)>2) .AND. (MOD(num_copies(i),2)==0)) then
              j=int(num_copies(i)/2)
              write(*,*)"偶数" , j

              allocate(kkk(num_copies(i),b) )
              do kk =1,j-1
                  h = 10 * (exp(-((2.0 * kk - 1.0) / N*1.0) )- 0.9) * 0.01
                  kkk(kk,1:19)=inparticles(i,1:19)-h
                  kkk(kk,20)=inparticles(i,20)
              end do
              kkk(j,:)=inparticles(i,:)
              kkk(j+1,:)=inparticles(i,:)
              do kk =(j+2),num_copies(i)
                  h=10 * (exp(-((2.0 * kk - 1.0) / N*1.0)) - 0.9) * 0.01
                  kkk(kk,1:19)=inparticles(i,1:19)+h
                  kkk(kk,20)=inparticles(i,20)
              end do

              outparticles(int(k):(int(k)+int(num_copies(i))-1),:)=kkk
              k=k+num_copies(i)
              deallocate (kkk)
           else if ((num_copies(i)>2) .AND. (MOD(num_copies(i),2)/=0)) then
              j=int((num_copies(i)+1)/2)
              write(*,*)"奇数" , j
              allocate(kkk(num_copies(i),b) )
              do kk =1,j-1
                  h = 10 * (exp(-((2.0 * kk - 1.0) / N*1.0) )- 0.9) * 0.01
                  kkk(kk,1:19)=inparticles(i,1:19)-h
                  kkk(kk,20)=inparticles(i,20)
              end do
              kkk(j,:)=inparticles(i,:)
              do kk =(j+1),num_copies(i)
                  h=10 * (exp(-((2.0 * kk - 1.0) / N*1.0)) - 0.9) * 0.01
                  kkk(kk,1:19)=inparticles(i,1:19)+h
                  kkk(kk,20)=inparticles(i,20)
              end do

              outparticles(int(k):(int(k)+int(num_copies(i))-1),:)=kkk
              k=k+num_copies(i)
              deallocate (kkk)
           end if
        else
          continue
        end if
    end do

    PF_resample%outparticles=outparticles
    deallocate (cum_copies)
    deallocate (num_copies)
    deallocate (outparticles)

  end subroutine resample

  subroutine write_PF_para(yy, mm, dd, tod,p)
    implicit none
    integer, intent(in) :: yy,mm,dd,tod,p
    integer   :: fid
    integer   :: nt,status
    integer   :: ncid,dimid_site,varid

    character(len=*), parameter :: sub = 'write_PF_para'

    character(len=255)  :: fln1
    character(len=8)    :: datestr,ppp
    integer :: i
    logical :: ldebug = .False.
    ! do not use r_decay, N_leaf, agb2vod, taweff,D0,
    real(r8),dimension(nlp)    ::Vcmax,q10,VJ_slope ,VN_slope,b_h2o,sif_alpha,sif_beta,Ksat_scalar,&
                                             b_scalar,m_h2o,f_leaf,kc25,ko25,tau25,f_resp,f_decay,a,b,c,pfweight,pfweightupdate
    real(r8),dimension(nlp)        :: NEP1,GPP1,SIF1,NPP1,LH1,SH1,Trans1,Evap1,Thetam1,COS_flux1, VOD1

    Vcmax=PF%Vcmax(p,1)
    q10=PF%q10(p,1)
    VJ_slope=PF%VJ_slope(p,1)
    !N_leaf=PF%N_leaf(p,1)
    VN_slope=PF%VN_slope(p,1)
    !r_decay=PF%r_decay(p,1)
    b_h2o=PF%b_h2o(p,1)
    sif_alpha=PF%sif_alpha(p,1)
    sif_beta=PF%sif_beta(p,1)
    !taweff=PF%taweff(p,1)
    !D0=PF%D0(p,1)
    Ksat_scalar=PF%Ksat_scalar(p,1)
    b_scalar=PF%b_scalar(p,1)
    m_h2o=PF%m_h2o(p,1)
    f_leaf=PF%f_leaf(p,1)
    kc25=PF%kc25(p,1)
    ko25=PF%ko25(p,1)
    tau25=PF%tau25(p,1)
    !agb2vod=PF%agb2vod(p,1)

    f_resp=PF%f_resp(p,1)
    f_decay=PF%f_decay(p,1)
    a=PF%a(p,1)
    b=PF%b(p,1)
    c=PF%c(p,1)

    pfweight=PF%pfweight(p,1)
    pfweightupdate=PF%pfweightupdate(p,1)

    NEP1 = NEP9
    GPP1 = GPP9
    SIF1 = SIF9
    NPP1 = NPP9
    COS_flux1 = COS_flux9
    LH1 = LH9
    SH1 = SH9
    Trans1 = Trans9
    Evap1 = Evap9
    Thetam1 = Thetam9
    ! 2023/07/04
    VOD1 = VOD9

    ! 输出文件
    if(nhtfrq <0) then
       write(datestr,"(i8)") yy*10000+mm*100+dd
       nt   = (tod/3600+1)/(-nhtfrq)
       ! !-- iLab::seconds elapsed since reference time (added for time-variable output)
       ! call timemgr_diff_secs(yy_ref*10000+mm_ref*100+dd_ref, tod_ref, yy*10000+mm*100+dd, tod,&
       !      secs_since_ref(1))
    else if(nhtfrq ==0) then
       write(datestr,"(i6)") yy*100+mm
       nt   = 1
    end if

    !-- iLab::make logging output depend on flag
    if(ldebug) then
       write(*,*) "Writing out PF para file now!"
    endif
    write(ppp,"(i8)") nt
    fln1  = trim(beps_out_dir)//"beps_site_"//trim(datestr)//"_"//trim(adjustl(ppp))//".nc"
    status =  nf90_open(fln1,nf90_write,ncid)
    if(status .ne. nf90_noerr) then

       call check(nf90_create(fln1,nf90_share,ncid))
       call check(nf90_def_dim(ncid,"paraloop",parloop,dimid_site))

       !-- iLab::added NetCDF output initialisation with prescribed _FillValue
       !         *UPDATE*: disabled for now, since maybe not compatible on
       !                   platform which Mousong is using.
       !Vcmax,q10,VJ_slope ,N_leaf,r_decay,b_h2o,sif_alpha,sif_beta,taweff,D0,Ksat_scalar,&
       !b_scalar,m_h2o,f_leaf,kc25,ko25,tau25,agb2vod,pfweight,pfweightupdate
       call check(nf90_def_var(ncid,"Vcmax",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"q10",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"VJ_slope",nf90_double,(/dimid_site/),varid))
       !call check(nf90_def_var(ncid,"N_leaf",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"VN_slope",nf90_double,(/dimid_site/),varid))
       !call check(nf90_def_var(ncid,"r_decay",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"b_h2o",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"sif_alpha",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"sif_beta",nf90_double,(/dimid_site/),varid))
       !call check(nf90_def_var(ncid,"taweff",nf90_double,(/dimid_site/),varid))
       !call check(nf90_def_var(ncid,"D0",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"Ksat_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"b_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"m_h2o",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"f_leaf",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"kc25",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"ko25",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"tau25",nf90_double,(/dimid_site/),varid))
       !call check(nf90_def_var(ncid,"agb2vod",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"f_resp",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"f_decay",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"a",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"b",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"c",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"pfweight",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"pfweightupdate",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"NEP",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Net Ecosystem Productivity"))

       call check(nf90_def_var(ncid,"GPP",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","gC/m2/h"))
       call check(nf90_put_att(ncid,varid,"name","Gross Primary Productivity"))

       call check(nf90_def_var(ncid,"NPP",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","kg/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","Net Primary Productivity"))

       call check(nf90_def_var(ncid,"COS_flux",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","pmol/m2/s"))
       call check(nf90_put_att(ncid,varid,"name","COS flux for soil and plant"))

       call check(nf90_def_var(ncid,"SIF",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","mW/m2/nm/sr"))
       call check(nf90_put_att(ncid,varid,"name","solar-induced SIF"))

       call check(nf90_def_var(ncid,"Thetam",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"name","Surface soil moisture"))

       call check(nf90_def_var(ncid,"LH",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"name","Latent heat flux"))

       call check(nf90_def_var(ncid,"SH",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","W/m2"))
       call check(nf90_put_att(ncid,varid,"name","Sensible heat flux"))

       call check(nf90_def_var(ncid,"Trans",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Transpiration"))

       call check(nf90_def_var(ncid,"Evap",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","m/s"))
       call check(nf90_put_att(ncid,varid,"name","Evaporation"))

       ! 2023/07/04
       call check(nf90_def_var(ncid,"VOD",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","none"))
       call check(nf90_put_att(ncid,varid,"name","Vegetation optical depth"))

       call check(nf90_enddef(ncid))
    end if

    call check(nf90_inq_varid(ncid,"Vcmax",varid))
    call check(nf90_put_var(ncid,varid,Vcmax,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"q10",varid))
    call check(nf90_put_var(ncid,varid,q10,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"VJ_slope",varid))
    call check(nf90_put_var(ncid,varid,VJ_slope,start=(/p/),count=(/1/)))
    !call check(nf90_inq_varid(ncid,"N_leaf",varid))
    !call check(nf90_put_var(ncid,varid,N_leaf,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"VN_slope",varid))
    call check(nf90_put_var(ncid,varid,VN_slope,start=(/p/),count=(/1/)))

    !call check(nf90_inq_varid(ncid,"r_decay",varid))
    !call check(nf90_put_var(ncid,varid,r_decay,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"b_h2o",varid))
    call check(nf90_put_var(ncid,varid,b_h2o,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"sif_alpha",varid))
    call check(nf90_put_var(ncid,varid,sif_alpha,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"sif_beta",varid))
    call check(nf90_put_var(ncid,varid,sif_beta,start=(/p/),count=(/1/)))
    !call check(nf90_inq_varid(ncid,"taweff",varid))
    !call check(nf90_put_var(ncid,varid,taweff,start=(/p/),count=(/1/)))
    !call check(nf90_inq_varid(ncid,"D0",varid))
    !call check(nf90_put_var(ncid,varid,D0,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Ksat_scalar",varid))
    call check(nf90_put_var(ncid,varid,Ksat_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"b_scalar",varid))
    call check(nf90_put_var(ncid,varid,b_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"m_h2o",varid))
    call check(nf90_put_var(ncid,varid,m_h2o,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"f_leaf",varid))
    call check(nf90_put_var(ncid,varid,f_leaf,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"kc25",varid))
    call check(nf90_put_var(ncid,varid,kc25,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"ko25",varid))
    call check(nf90_put_var(ncid,varid,ko25,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"tau25",varid))
    call check(nf90_put_var(ncid,varid,tau25,start=(/p/),count=(/1/)))
    !call check(nf90_inq_varid(ncid,"agb2vod",varid))
    !call check(nf90_put_var(ncid,varid,agb2vod,start=(/p/),count=(/1/)))
    ! 2023/07/04
    call check(nf90_inq_varid(ncid,"f_resp",varid))
    call check(nf90_put_var(ncid,varid,f_resp,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"f_decay",varid))
    call check(nf90_put_var(ncid,varid,f_decay,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"a",varid))
    call check(nf90_put_var(ncid,varid,a,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"b",varid))
    call check(nf90_put_var(ncid,varid,b,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"c",varid))
    call check(nf90_put_var(ncid,varid,c,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"pfweight",varid))
    call check(nf90_put_var(ncid,varid,pfweight,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"pfweightupdate",varid))
    call check(nf90_put_var(ncid,varid,pfweightupdate,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"NEP",varid))
    call check(nf90_put_var(ncid,varid,NEP1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"GPP",varid))
    call check(nf90_put_var(ncid,varid,GPP1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"NPP",varid))
    call check(nf90_put_var(ncid,varid,NPP1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"COS_flux",varid))
    call check(nf90_put_var(ncid,varid,COS_flux1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"SIF",varid))
    call check(nf90_put_var(ncid,varid,SIF1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Thetam",varid))
    call check(nf90_put_var(ncid,varid,Thetam1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"LH",varid))
    call check(nf90_put_var(ncid,varid,LH1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"SH",varid))
    call check(nf90_put_var(ncid,varid,SH1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Trans",varid))
    call check(nf90_put_var(ncid,varid,Trans1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Evap",varid))
    call check(nf90_put_var(ncid,varid,Evap1,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"VOD",varid))
    call check(nf90_put_var(ncid,varid,VOD1,start=(/p/),count=(/1/)))

    call check(nf90_close(ncid))

  end subroutine write_PF_para

end module PF_run
