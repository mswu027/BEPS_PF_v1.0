!***************************************************************
!! This module will read control parameters,meteo input,boundary conditions,
!! yrdata,cpools etc. for beps
!! Editted by J. Wang
!! 22May 2017
!***************************************************************

module controlInput_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
  use beps_par
  use bepstype
  !--iLab::restricted use(s)
  !-- update:with required arguments in routines below, can avoid completely
  ! use beps_time_manager,only: get_curr_date, get_curr_calday
  ! use beps_time_manager,only: set_timemgr_init,timemgr_init,get_curr_date,&
  !      get_prev_date,get_curr_calday,timemgr_datediff
  use netcdf
  !--iLab::no entity of module 'esmf' used here!
  ! use esmf
  implicit none
  !include 'mpif.h'
  !--iLab::no need for include file, since we call 'use netcdf' above
  ! include 'netcdf.inc'
  save

  integer :: nlat,nlon,m
  character(len=80) :: calendar
  integer :: sim_type
  integer :: icdate
  integer :: icsec
  integer :: sim_duration
  integer :: nhtfrq,nstpd         !!nstpd determines when to output
  integer :: restart_frq
  integer :: meteo_input
  integer :: nscale,n_site
  integer :: lai_input
  character(len=255) :: meteo_path,meteo_flnm_prefix,meteo_site_flnm_prefix,prior_PF_obs_prefix
  character(len=255) :: surface_data_path
  character(len=255) :: beps_yrdata_path
  character(len=255) :: beps_site_path,site_bound_prefix,prior_para_prefix
  character(len=255) :: beps_lai_path,beps_lai_prefix,beps_lai_site_prefix
  character(len=255) :: beps_Vcmax_path,beps_Vcmax_site_path
  character(len=255) :: beps_domain
  character(len=255) :: beps_cpools
  character(len=255) :: beps_PF_obs_path
  character(len=255) :: beps_out_dir
  character(len=255) :: beps_rst_dir
  public :: rdnamelist
  private:: read_beps_domain

contains

  subroutine rdnamelist()
    implicit none
    integer ::ierr

    namelist /NLS/ nlat,nlon,nscale,calendar,icdate,icsec,sim_type,sim_duration,nhtfrq,restart_frq, &
         meteo_input,meteo_path,meteo_flnm_prefix,meteo_site_flnm_prefix,surface_data_path,&
         beps_yrdata_path,n_site,beps_site_path,site_bound_prefix,lai_input,beps_lai_path,&
         beps_lai_prefix,beps_lai_site_prefix,beps_Vcmax_path,beps_Vcmax_site_path,beps_domain,&
         prior_para_prefix,beps_cpools,beps_PF_obs_path,prior_PF_obs_prefix,beps_out_dir,beps_rst_dir

    !if(myid ==0) then
    open(5,file='beps.stdin',form='formatted')
    rewind 5
    read (5,nml=NLS,iostat=ierr)

    if (ierr >  0) then
       write(6,*)'RDNAMELIST: Namelist read returns ',ierr
       call endrun("Read Namelist Wrong!")
    end if
    if (nscale==0) then
       call read_beps_domain(beps_domain)
    end if
    !end if

    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_bcast(nlat,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(nlon,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(nscale,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(calendar,len(calendar),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(icdate,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(icsec,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(sim_type,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(sim_duration,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(nhtfrq,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(restart_frq,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(meteo_input,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(meteo_path,len(meteo_path),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(meteo_flnm_prefix,len(meteo_flnm_prefix),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(meteo_site_flnm_prefix,len(meteo_site_flnm_prefix),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(surface_data_path,len(surface_data_path),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_yrdata_path,len(beps_yrdata_path),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(n_site,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_site_path,len(beps_site_path),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(site_bound_prefix,len(site_bound_prefix),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(lai_input,1,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_lai_path,len(beps_lai_path),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_lai_prefix,len(beps_lai_prefix),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_lai_site_prefix,len(beps_lai_site_prefix),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_Vcmax_path,len(beps_Vcmax_path),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_Vcmax_site_path,len(beps_Vcmax_site_path),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_domain,len(beps_domain),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_cpools,len(beps_cpools),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_out_dir,len(beps_out_dir),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_bcast(beps_rst_dir,len(beps_rst_dir),mpi_character,0,mpi_comm_world,ierr)
    !call mpi_barrier(mpi_comm_world,ierr)

    !!  store run type
    nsrest    = sim_type

    !! setting nstpd
    if(nhtfrq <0) then    ! means hours
       nstpd  = -nhtfrq*3600/step
    end if

  end subroutine rdnamelist

  subroutine read_beps_domain(filname)
    
    implicit none
    character(len=*),intent(in):: filname
    integer:: ncid,varid(2)
    integer:: domain(nlon,nlat)
    integer:: i,j

    call check(nf90_open(trim(filname),nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"domain",varid(1)))
    call check(nf90_inq_varid(ncid,"nlp",varid(2)))
    call check(nf90_get_var(ncid,varid(1),domain))
    call check(nf90_get_var(ncid,varid(2),nlp))
    call check(nf90_close(ncid))

    allocate(stype(nlon*nlat))
    allocate(mapping(nlp))
    stype = reshape(domain,(/nlon*nlat/))

    j = 1
    do i = 1,nlon*nlat
       if(stype(i) ==1) then
          mapping(j) = i
          j = j+1
       end if
    end do
    !deallocate(stype)
    !deallocate(mapping)

  end subroutine read_beps_domain

  !***************Reading boundary data*************************

  subroutine read_boundary()
    implicit none
    integer :: i,ncid,varid(9),ierr
    integer :: lcno(nlon,nlat,PFT),stext(nlon,nlat)
    real(r8):: PCT_PFT(nlon,nlat,PFT),ci(nlon,nlat),longitude(nlon,nlat),latitude(nlon,nlat),&
         sdp(nlon,nlat),st(nlon,nlat),sw(nlon,nlat)
    integer :: lcno2(nlon*nlat,PFT),stext2(nlon*nlat)
    real(r8):: PCT_PFT2(nlon*nlat,PFT),ci2(nlon*nlat),longitude2(nlon*nlat),latitude2(nlon*nlat),&
         sdp2(nlon*nlat),st2(nlon*nlat),sw2(nlon*nlat)
    integer :: lcno3(nlp,PFT),stext3(nlp)
    real(r8):: PCT_PFT3(nlp,PFT),ci3(nlp),longitude3(nlp),latitude3(nlp),&
         sdp3(nlp),st3(nlp),sw3(nlp)
    type(surf),pointer:: p

    p=>bound

    !if(myid == 0) then
    call check(nf90_open(trim(surface_data_path),nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"PFT",varid(1)))
    call check(nf90_inq_varid(ncid,"PCT_PFT",varid(2)))
    call check(nf90_inq_varid(ncid,"ci",varid(3)))
    call check(nf90_inq_varid(ncid,"stext",varid(4)))
    call check(nf90_inq_varid(ncid,"longitude",varid(5)))
    call check(nf90_inq_varid(ncid,"latitude",varid(6)))
    call check(nf90_inq_varid(ncid,"sdp",varid(7)))
    call check(nf90_inq_varid(ncid,"st",varid(8)))
    call check(nf90_inq_varid(ncid,"sw",varid(9)))

    call check(nf90_get_var(ncid,varid(1),lcno))
    call check(nf90_get_var(ncid,varid(2),PCT_PFT))
    call check(nf90_get_var(ncid,varid(3),ci))
    call check(nf90_get_var(ncid,varid(4),stext))
    call check(nf90_get_var(ncid,varid(5),longitude))
    call check(nf90_get_var(ncid,varid(6),latitude))
    call check(nf90_get_var(ncid,varid(7),sdp))
    call check(nf90_get_var(ncid,varid(8),st))
    call check(nf90_get_var(ncid,varid(9),sw))

    call check(nf90_close(ncid))

    lcno2      =  reshape(lcno,(/nlon*nlat,PFT/))
    PCT_PFT2   =  reshape(PCT_PFT,(/nlon*nlat,PFT/))
    ci2        =  reshape(ci,(/nlon*nlat/))
    stext2     =  reshape(stext,(/nlon*nlat/))
    longitude2 =  reshape(longitude,(/nlon*nlat/))
    latitude2  =  reshape(latitude,(/nlon*nlat/))
    sdp2       =  reshape(sdp,(/nlon*nlat/))
    st2        =  reshape(st,(/nlon*nlat/))
    sw2        =  reshape(sw,(/nlon*nlat/))

    lcno3      = lcno2(mapping,:)
    PCT_PFT3   = PCT_PFT2(mapping,:)
    ci3        = ci2(mapping)
    stext3     = stext2(mapping)
    longitude3 = longitude2(mapping)
    latitude3  = latitude2(mapping)

    sdp3       = sdp2(mapping)
    st3        = st2(mapping)
    sw3        = sw2(mapping)
    !end if

    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_scatterv(ci3(1),dp,sp,MPI_real8,p%clumping(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(stext3(1),dp,sp,mpi_integer,p%stext(1),npoints,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_scatterv(longitude3(1),dp,sp,mpi_real8,p%longitude(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(latitude3(1),dp,sp,mpi_real8,p%latitude(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(sdp3(1),dp,sp,mpi_real8,p%sdp(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(st3(1),dp,sp,mpi_real8,p%st(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(sw3(1),dp,sp,mpi_real8,p%sw(1),npoints,mpi_real8,0,mpi_comm_world,ierr)


    !do i=1,PFT
    !   call mpi_scatterv(lcno3(1,i),dp,sp,mpi_integer,p%lcno(1,i),npoints,mpi_integer,0,mpi_comm_world,ierr)
    !   call mpi_scatterv(PCT_PFT3(1,i),dp,sp,mpi_real8,p%PCT_PFT(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !end do

    !call mpi_barrier(mpi_comm_world,ierr)

    p%clumping = ci3
    p%stext = stext3
    p%longitude = longitude3
    p%latitude = latitude3
    p%sdp = sdp3
    p%st = st3
    p%sw = sw3
    p%lcno = lcno3
    p%PCT_PFT = PCT_PFT3

  end subroutine read_boundary

  subroutine read_yrdata()
    implicit none
    integer  :: i,ncid,varid(2),ierr
    real(r8) :: laiyr1(nlon,nlat,PFT),nppyr1(nlon,nlat,PFT)
    real(r8) :: laiyr2(nlon*nlat,PFT),nppyr2(nlon*nlat,PFT)
    real(r8) :: laiyr3(nlp,PFT),nppyr3(nlp,PFT)
    type(surf),pointer :: p
    p=>bound

    !if(myid ==0) then
    call check(nf90_open(trim(beps_yrdata_path),nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"lai",varid(1)))
    call check(nf90_inq_varid(ncid,"npp",varid(2)))
    call check(nf90_get_var(ncid,varid(1),laiyr1))
    call check(nf90_get_var(ncid,varid(2),nppyr1))
    call check(nf90_close(ncid))

    laiyr2    = reshape(laiyr1,(/nlon*nlat,PFT/))
    nppyr2    = reshape(nppyr1,(/nlon*nlat,PFT/))

    laiyr3    = laiyr2(mapping,:)
    nppyr3    = nppyr2(mapping,:)
    !end if

    !call mpi_barrier(mpi_comm_world,ierr)

    !do i = 1,PFT
    !   call mpi_scatterv(laiyr3(1,i),dp,sp,mpi_real8,p%laiyr(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !   call mpi_scatterv(nppyr3(1,i),dp,sp,mpi_real8,p%nppyr(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !end do

    !call mpi_barrier(mpi_comm_world,ierr)
    p%laiyr = laiyr3
    p%nppyr = nppyr3

  end subroutine read_yrdata

  subroutine read_boundary_site()
    implicit none
    integer  :: i,ncid,varid(20),ierr
    integer  :: lcno(nlp,PFT),stext(nlp)
    real(r8) :: laiyr(nlp,PFT),nppyr(nlp,PFT)
    real(r8) :: PCT_PFT(nlp,PFT)
    real(r8) :: longitude(nlp),latitude(nlp),sdp(nlp),st(nlp),sw(nlp),ci(nlp)
    real(r8) :: ccd(nlp,PFT),cfmd(nlp,PFT),cfsd(nlp,PFT),cm(nlp,PFT),cp(nlp,PFT),&
         cs(nlp,PFT),csm(nlp,PFT),csmd(nlp,PFT),cssd(nlp,PFT)
    !--iLab::avoid pointer
    ! type(surf),pointer :: p
    ! p => bound
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_boundary_site'

    if( len(trim(beps_site_path))+len(trim(site_bound_prefix))+len('.nc') .gt. len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//trim(site_bound_prefix)//'.nc'
    endif
    print*, 'START read_boundary_site with *****'//trim(fname)//'*****'
    !if(myid == 0) then
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"laiyr",varid(1)))
    call check(nf90_inq_varid(ncid,"nppyr",varid(2)))
    call check(nf90_inq_varid(ncid,"lcno",varid(3)))
    call check(nf90_inq_varid(ncid,"PCT_PFT",varid(4)))
    call check(nf90_inq_varid(ncid,"ci",varid(5)))
    call check(nf90_inq_varid(ncid,"stext",varid(6)))
    call check(nf90_inq_varid(ncid,"longitude",varid(7)))
    call check(nf90_inq_varid(ncid,"latitude",varid(8)))
    call check(nf90_inq_varid(ncid,"sdp",varid(9)))
    call check(nf90_inq_varid(ncid,"st",varid(10)))
    call check(nf90_inq_varid(ncid,"sw",varid(11)))
    call check(nf90_inq_varid(ncid,"ccd",varid(12)))
    call check(nf90_inq_varid(ncid,"cfmd",varid(13)))
    call check(nf90_inq_varid(ncid,"cfsd",varid(14)))
    call check(nf90_inq_varid(ncid,"cm",varid(15)))
    call check(nf90_inq_varid(ncid,"cp",varid(16)))
    call check(nf90_inq_varid(ncid,"cs",varid(17)))
    call check(nf90_inq_varid(ncid,"csm",varid(18)))
    call check(nf90_inq_varid(ncid,"csmd",varid(19)))
    call check(nf90_inq_varid(ncid,"cssd",varid(20)))

    call check(nf90_get_var(ncid,varid(1),laiyr))
    call check(nf90_get_var(ncid,varid(2),nppyr))
    call check(nf90_get_var(ncid,varid(3),lcno))
    call check(nf90_get_var(ncid,varid(4),PCT_PFT))
    call check(nf90_get_var(ncid,varid(5),ci))
    call check(nf90_get_var(ncid,varid(6),stext))
    call check(nf90_get_var(ncid,varid(7),longitude))
    call check(nf90_get_var(ncid,varid(8),latitude))
    call check(nf90_get_var(ncid,varid(9),sdp))
    call check(nf90_get_var(ncid,varid(10),st))
    call check(nf90_get_var(ncid,varid(11),sw))
    call check(nf90_get_var(ncid,varid(12),ccd))
    call check(nf90_get_var(ncid,varid(13),cfmd))
    call check(nf90_get_var(ncid,varid(14),cfsd))
    call check(nf90_get_var(ncid,varid(15),cm))
    call check(nf90_get_var(ncid,varid(16),cp))
    call check(nf90_get_var(ncid,varid(17),cs))
    call check(nf90_get_var(ncid,varid(18),csm))
    call check(nf90_get_var(ncid,varid(19),csmd))
    call check(nf90_get_var(ncid,varid(20),cssd))
    call check(nf90_close(ncid))
    !end if
    PCT_PFT = PCT_PFT          !! convert fraction to %
    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_scatterv(ci(1),dp,sp,MPI_real8,bound%clumping(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(stext(1),dp,sp,mpi_integer,bound%stext(1),npoints,mpi_integer,0,mpi_comm_world,ierr)
    !call mpi_scatterv(longitude(1),dp,sp,mpi_real8,bound%longitude(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(latitude(1),dp,sp,mpi_real8,bound%latitude(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(sdp(1),dp,sp,mpi_real8,bound%sdp(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(st(1),dp,sp,mpi_real8,bound%st(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(sw(1),dp,sp,mpi_real8,bound%sw(1),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(name(1),dp,sp,mpi_character,bound%name(1),npoints,mpi_character,0,mpi_comm_world,ierr)
    !call mpi_barrier(mpi_comm_world,ierr)
    !do i = 1,PFT
    !    call mpi_scatterv(lcno(1,i),dp,sp,mpi_integer,bound%lcno(1,i),npoints,mpi_integer,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(PCT_PFT(1,i),dp,sp,mpi_real8,bound%PCT_PFT(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(laiyr(1,i),dp,sp,mpi_real8,bound%laiyr(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(nppyr(1,i),dp,sp,mpi_real8,bound%nppyr(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(ccd(1,i),dp,sp,MPI_real8,bound%ccd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(cfmd(1,i),dp,sp,MPI_real8,bound%cfmd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(cfsd(1,i),dp,sp,MPI_real8,bound%cfsd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(cm(1,i),dp,sp,MPI_real8,bound%cm(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(cp(1,i),dp,sp,MPI_real8,bound%cp(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(cs(1,i),dp,sp,MPI_real8,bound%cs(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(csm(1,i),dp,sp,MPI_real8,bound%csm(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(csmd(1,i),dp,sp,MPI_real8,bound%csmd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !    call mpi_scatterv(cssd(1,i),dp,sp,MPI_real8,bound%cssd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !end do
    
    !call mpi_barrier(mpi_comm_world,ierr)
    bound%clumping = ci
    bound%stext = stext
    bound%longitude = longitude
    bound%latitude = latitude
    bound%sdp = sdp
    bound%st = st
    bound%sw = sw
    bound%lcno = lcno
    bound%PCT_PFT = PCT_PFT
    bound%laiyr = laiyr
    bound%nppyr = nppyr
    bound%ccd = ccd
    bound%cfmd = cfmd
    bound%cfsd = cfsd
    bound%cm = cm
    bound%cp = cp
    bound%cs = cs
    bound%csm = csm
    bound%csmd = csmd
    bound%cssd = cssd

  end subroutine read_boundary_site

  subroutine read_meteo_daily(yr, mon, day, sec)
    implicit none
    !--iLab::made yr/mon/day/sec arguments in order to avoid 'get_curr_date'
    integer, intent(in) :: yr,mon,day,sec
    !character(len=*) :: file_path,file_flnm_prefix
    integer          :: ncid,ierr,varid(6)
    integer          :: nt
    integer          :: i
    !--iLab::avoid pointer (see also below)
    ! type(forc),pointer :: p
    real(r8)           :: T1(nlon,nlat),TMAX1(nlon,nlat),TMIN1(nlon,nlat),RH1(nlon,nlat),&
         WS1(nlon,nlat),PRCP1(nlon,nlat),SSRD1(nlon,nlat)
    real(r8)           :: T2(nlon*nlat),TMAX2(nlon*nlat),TMIN2(nlon*nlat),RH2(nlon*nlat),&
         WS2(nlon*nlat),PRCP2(nlon*nlat),SSRD2(nlon*nlat)
    real(r8)           :: T3(nlp),TMAX3(nlp),TMIN3(nlp),RH3(nlp),WS3(nlp),PRCP3(nlp),SSRD3(nlp)
    real(r8)           :: rainfall(nlp),snow(nlp)
    character(len=6)   :: monstr
    character(len=8)   :: datestr
    character(len=4)   :: yrstr
    real(r8),parameter :: density_water = 1025.
    real(r8)           :: density_new_snow
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_meteo_daily'

    !if(myid ==0) then
    !   day  = get_curr_calday()
    ! call  get_curr_date(yr,mon,day,sec)
    !   call get_prev_date(yr,mon,day,sec)   !!??
    !  write(6,*) "Reading met data on ",sec,"of", yr*10000+mon*100+day
    !  write(datestr,'(i8)')  yr*10000+mon*100+day
    write(monstr,'(i6)')   yr*100+mon
    write(yrstr,'(i4)')    yr
    !  call check(nf90_open(trim(meteo_path)//trim(yrstr)//'/'//trim(meteo_flnm_prefix)//trim(datestr)//".nc",nf90_nowrite,ncid))
    if( len(trim(meteo_path))+len(trim(meteo_flnm_prefix))+len(trim(monstr))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(meteo_path)//'/'//trim(meteo_flnm_prefix)//trim(monstr)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"TMAX",varid(1)))
    call check(nf90_inq_varid(ncid,"TMIN",varid(2)))
    call check(nf90_inq_varid(ncid,"RH",varid(3)))
    call check(nf90_inq_varid(ncid,"WS",varid(4)))
    call check(nf90_inq_varid(ncid,"PRCP",varid(5)))
    call check(nf90_inq_varid(ncid,"SSRD",varid(6)))

    !  nt    = (day-1)*24+sec/int(step)+1     !! data slide for months
    !  nt    = sec/int(step)+1
    write(*,*) "Reading meteo data on ",day,"of",monstr
    call check(nf90_get_var(ncid,varid(1),TMAX1,start=(/1,1,day/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(2),TMIN1,start=(/1,1,day/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(3),RH1,start=(/1,1,day/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(4),WS1,start=(/1,1,day/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(5),PRCP1,start=(/1,1,day/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(6),SSRD1,start=(/1,1,day/),count=(/nlon,nlat,1/)))
    call check(nf90_close(ncid))
    TMAX2       = reshape(TMAX1,(/nlon*nlat/))
    TMIN2       = reshape(TMIN1,(/nlon*nlat/))
    RH2      = reshape(RH1,(/nlon*nlat/))
    WS2      = reshape(WS1,(/nlon*nlat/))
    PRCP2    = reshape(PRCP1,(/nlon*nlat/))
    SSRD2    = reshape(SSRD1,(/nlon*nlat/))

    TMAX3       = TMAX2(mapping)    !! K
    TMIN3       = TMIN2(mapping)    !! K
    RH3      = RH2(mapping)   !!
    WS3      = WS2(mapping)   !! m/s
    PRCP3    = PRCP2(mapping)  !! m/s(consult)
    SSRD3    = SSRD2(mapping)  !! w/m2

    TMAX3       = TMAX3 -273.16     !! translate into centigrade (oC)
    TMIN3       = TMIN3 -273.16     !! translate into centigrade (oC)
    RH3      = RH3*100        !! %
    PRCP3    = max(0.,PRCP3/3600.)  !! carefull for negative prec, convert from m/hr to m/s,mousong.wu@201902

    T3       = (TMAX3 + TMIN3)/2.
    !! seperating total precip into liquid and solid according to temperature
    do i = 1,nlp
       if(T3(i)>0.) then
          rainfall(i)   = PRCP3(i)
          snow(i)       = 0.
       else
          rainfall(i)   = 0.
          density_new_snow  = 67.9+51.3*exp(T3(i)/2.6)
          snow(i)       = PRCP3(i)*density_water/density_new_snow
       end if
    end do
    !end if

    !--iLab::avoid pointer
    ! p =>clim
    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_scatterv(TMAX3(1),dp,sp,MPI_real8,clim%Tempmx(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(TMIN3(1),dp,sp,MPI_real8,clim%Tempmn(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(RH3(1),dp,sp,MPI_real8,clim%Rh(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(WS3(1),dp,sp,MPI_real8,clim%Wind(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(SSRD3(1),dp,sp,MPI_real8,clim%Srad(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(rainfall(1),dp,sp,MPI_real8,clim%Rain(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(snow(1),dp,sp,MPI_real8,clim%Snow(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_barrier(mpi_comm_world,ierr)
    clim%Tempmx = TMAX3
    clim%Tempmn = TMIN3
    clim%Rh = RH3
    clim%Wind = WS3
    clim%Srad = SSRD3
    clim%Rain = rainfall
    clim%Snow = snow
  end subroutine read_meteo_daily

  subroutine read_meteo_hourly(yr, mon, day, sec)
    implicit none
    !--iLab::made yr/mon/day/sec arguments in order to avoid 'get_curr_date'
    integer, intent(in) :: yr,mon,day,sec
    !character(len=*) :: file_path,file_flnm_prefix
    integer          :: ncid,ierr,varid(5)
    ! integer          :: yr,mon,day,sec
    integer          :: nt
    integer          :: i
    !--iLab::avoid pointer (see also below)
    ! type(forc),pointer :: p
    real(r8)           :: T1(nlon,nlat),RH1(nlon,nlat),WS1(nlon,nlat),PRCP1(nlon,nlat),SSRD1(nlon,nlat)
    real(r8)           :: T2(nlon*nlat),RH2(nlon*nlat),WS2(nlon*nlat),PRCP2(nlon*nlat),SSRD2(nlon*nlat)
    real(r8)           :: T3(nlp),RH3(nlp),WS3(nlp),PRCP3(nlp),SSRD3(nlp)
    real(r8)           :: rainfall(nlp),snow(nlp)
    !character(len=6)   :: monstr
    character(len=8)   :: datestr
    character(len=4)   :: yrstr
    real(r8),parameter :: density_water = 1025.
    real(r8)           :: density_new_snow
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_meteo_hourly'
    !if(myid ==0) then
    ! call get_curr_date(yr,mon,day,sec)

    write(datestr,'(i8)')  yr*10000+mon*100+day
    write(yrstr,'(i4)')    yr
    !  call check(nf90_open(trim(meteo_path)//trim(yrstr)//'/'//trim(meteo_flnm_prefix)//trim(datestr)//".nc",nf90_nowrite,ncid))
    if( len(trim(meteo_path))+len(trim(meteo_flnm_prefix))+len(trim(datestr))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(meteo_path)//'/'//trim(meteo_flnm_prefix)//trim(datestr)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"T",varid(1)))
    call check(nf90_inq_varid(ncid,"RH",varid(2)))
    call check(nf90_inq_varid(ncid,"WS",varid(3)))
    call check(nf90_inq_varid(ncid,"PRCP",varid(4)))
    call check(nf90_inq_varid(ncid,"SSRD",varid(5)))

    !  nt    = (day-1)*24+sec/int(step)+1     !! data slide for months
    nt    = sec/int(step)+1
    write(*,*) "Reading meteo data on ",nt,"of",datestr
    call check(nf90_get_var(ncid,varid(1),T1,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(2),RH1,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(3),WS1,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(4),PRCP1,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_get_var(ncid,varid(5),SSRD1,start=(/1,1,nt/),count=(/nlon,nlat,1/)))
    call check(nf90_close(ncid))
    T2       = reshape(T1,(/nlon*nlat/))
    RH2      = reshape(RH1,(/nlon*nlat/))
    WS2      = reshape(WS1,(/nlon*nlat/))
    PRCP2    = reshape(PRCP1,(/nlon*nlat/))
    SSRD2    = reshape(SSRD1,(/nlon*nlat/))

    T3       = T2(mapping)    !! K
    RH3      = RH2(mapping)   !!
    WS3      = WS2(mapping)   !! m/s
    PRCP3    = PRCP2(mapping)  !! m/s(consult)
    SSRD3    = SSRD2(mapping)  !! w/m2

    T3       = T3 -273.16     !! translate into centigrade (oC)
    RH3      = RH3*100        !! %
    PRCP3    = max(0.,PRCP3/3600.)  !! carefull for negative prec, convert from m/hr to m/s,mousong.wu@201902

    !! seperating total precip into liquid and solid according to temperature
    do i = 1,nlp
       if(T3(i)>0.) then
          rainfall(i)   = PRCP3(i)
          snow(i)       = 0.
       else
          rainfall(i)   = 0.
          density_new_snow  = 67.9+51.3*exp(T3(i)/2.6)
          snow(i)       = PRCP3(i)*density_water/density_new_snow
       end if

    end do
    !end if

    !--iLab::avoid pointer
    ! p =>clim

    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_scatterv(T3(1),dp,sp,MPI_real8,clim%Temp(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(RH3(1),dp,sp,MPI_real8,clim%Rh(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(WS3(1),dp,sp,MPI_real8,clim%Wind(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(SSRD3(1),dp,sp,MPI_real8,clim%Srad(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(rainfall(1),dp,sp,MPI_real8,clim%Rain(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(snow(1),dp,sp,MPI_real8,clim%Snow(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_barrier(mpi_comm_world,ierr)
    clim%Temp = T3
    clim%Rh = RH3
    clim%Wind = WS3
    clim%Srad = SSRD3
    clim%Rain = rainfall
    clim%Snow = snow

  end subroutine read_meteo_hourly
  
  subroutine read_meteo_site(nd)
    implicit none
    !character(len=*) :: file_path,file_flnm_prefix
    integer,intent(in)  :: nd
    integer          :: ncid,ierr,varid(5)
    integer          :: yr,mon,day,sec
    integer          :: nt
    integer          :: i
    !--iLab::avoid pointer (see also below)
    ! type(forc),pointer :: p
    real(r8)           :: T(nlp),RH(nlp),WS(nlp),PRCP(nlp),SSRD(nlp)
    real(r8)           :: rainfall(nlp),snow(nlp)
    !character(len=6)   :: monstr
    character(len=8)   :: datestr
    character(len=4)   :: yrstr
    real(r8),parameter :: density_water = 1025.
    real(r8)           :: density_new_snow
    logical :: ldebug
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_meteo_site'
    !--iLab::added missing initialisation
    ldebug = .false.

    !if(myid ==0) then
    !  call check(nf90_open(trim(meteo_path)//trim(yrstr)//'/'//trim(meteo_flnm_prefix)//trim(datestr)//".nc",nf90_nowrite,ncid))
    if( len(trim(beps_site_path))+len(trim(meteo_site_flnm_prefix))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//'/'//trim(meteo_site_flnm_prefix)//".nc"
    endif
    if(ldebug) then
       write(*,'(a)') 'DEBUG::'//method//':reading meteo data ***'//trim(fname)//'***'
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"T",varid(1)))
    call check(nf90_inq_varid(ncid,"RH",varid(2)))
    call check(nf90_inq_varid(ncid,"WS",varid(3)))
    call check(nf90_inq_varid(ncid,"PRCP",varid(4)))
    call check(nf90_inq_varid(ncid,"SSRD",varid(5)))
    if(ldebug) then
       write(*,'(a)') 'DEBUG::'//method//':data set identifier determined!'
    endif

    !  nt    = (day-1)*24+sec/int(step)+1     !! data slide for months
    if(ldebug) then
       write(*,'(a,2(a,i5,1x))') 'DEBUG::'//method//':reading meteo data on ',&
            "nd=",nd,"nlp=",nlp
    endif
    call check(nf90_get_var(ncid,varid(1),T,start=(/1,nd/),count=(/nlp,1/)))
    call check(nf90_get_var(ncid,varid(2),RH,start=(/1,nd/),count=(/nlp,1/)))
    call check(nf90_get_var(ncid,varid(3),WS,start=(/1,nd/),count=(/nlp,1/)))
    call check(nf90_get_var(ncid,varid(4),PRCP,start=(/1,nd/),count=(/nlp,1/)))
    call check(nf90_get_var(ncid,varid(5),SSRD,start=(/1,nd/),count=(/nlp,1/)))
    call check(nf90_close(ncid))
    if(ldebug) then
       write(*,'(a)') 'DEBUG::'//method//':NetCDF file closed again.'
    endif

    T       = T             !! translate into centigrade (oC), input is oC for site.@MOUSONG.WU
    RH      = RH*100        !! %
    PRCP    = max(0.,PRCP/3600.)  !! carefull for negative prec, convert from m/hr to m/s,mousong.wu@201902
    SSRD    = max(0.,SSRD)
    !! seperating total precip into liquid and solid according to temperature
    do i = 1,nlp
       if(T(i)>0.) then
          rainfall(i)   = PRCP(i)
          snow(i)       = 0.
       else
          rainfall(i)   = 0.
          density_new_snow  = 67.9+51.3*exp(T(i)/2.6)
          snow(i)       = PRCP(i)*density_water/density_new_snow
       end if

    end do
    !end if

    !--iLab::avoid pointer
    ! p =>clim

    !call mpi_barrier(mpi_comm_world,ierr)
    !call mpi_scatterv(T(1),dp,sp,MPI_real8,clim%Temp(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(RH(1),dp,sp,MPI_real8,clim%Rh(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(WS(1),dp,sp,MPI_real8,clim%Wind(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(SSRD(1),dp,sp,MPI_real8,clim%Srad(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(rainfall(1),dp,sp,MPI_real8,clim%Rain(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_scatterv(snow(1),dp,sp,MPI_real8,clim%Snow(1),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !call mpi_barrier(mpi_comm_world,ierr)
    clim%Temp = T
    clim%Rh = RH
    clim%Wind = WS
    clim%Srad = SSRD
    clim%Rain = rainfall
    clim%Snow = snow

  end subroutine read_meteo_site


  subroutine read_lai(yr,mn,dd,tod,day)
    implicit none
    !--iLab::made yr/mon/day/sec arguments in order to avoid 'get_curr_date'/'get_cal_day'
    integer, intent(in) :: yr,mn,dd,tod,day
    real(r8)         :: lai1(nlon,nlat,PFT),lai2(nlon*nlat,PFT),lai3(nlp,PFT)
    character(len = 4) :: datestr
    integer          :: ncid,varid,ierr,i
    !--iLab::avoid pointer (see also below)
    ! type(surf),pointer:: p
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_lai'
    !if(myid ==0) then
    ! day  = get_curr_calday()
    ! call get_curr_date(yr,mn,dd,tod)

    write(datestr,'(i4)')  yr
    write(*,*) "Reading LAI On ",yr*10000+mn*100+dd
    if( len(trim(beps_lai_path))+len(trim(beps_lai_prefix))+len(trim(datestr))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_lai_path)//trim(beps_lai_prefix)//trim(datestr)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"lai",varid))

    call check(nf90_get_var(ncid,varid,lai1,start =(/1,1,1,day/),count = (/nlon,nlat,PFT,1/)))
    call check(nf90_close(ncid))
    lai2  = reshape(lai1,(/nlon*nlat,PFT/))
    lai3  = lai2(mapping,:)
    !end if

    !--iLab::avoid pointer
    ! p => bound

    !call mpi_barrier(mpi_comm_world,ierr)
    !do i = 1,PFT
    !    call mpi_scatterv(lai3(1,i),dp,sp,mpi_real8,p%lai(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !end do
    !call mpi_barrier(mpi_comm_world,ierr)
    bound%lai = lai3

  end subroutine read_lai

  subroutine read_lai_site(day)
    implicit none
    !--iLab::no need for actual day, calendar day as input
    integer, intent(in) :: day
    ! integer          :: yr,mn,dd,tod,day
    real(r8)         :: lai1(nlp,PFT)
    integer          :: ncid,varid,ierr,i
    !--iLab::avoid pointer (see also below)
    ! type(surf),pointer:: p
    !--iLab::added to limit terminal output
    logical :: ldebug
    character(len=255) :: fname
    character(len=*), parameter :: method = "read_lai_site"
    !if(myid ==0) then
    ! day  = get_curr_calday()
    ! call get_curr_date(yr,mn,dd,tod)
    ldebug = .false.
    if(ldebug) then
       write(*,*) "INFO::"//method//":Reading LAI day=", day
    endif
    if( len(trim(beps_site_path))+len(trim(beps_lai_site_prefix))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//trim(beps_lai_site_prefix)//".nc"
    endif
    if(ldebug) then
       write(*,*) "INFO::"//method//":Reading from file ***"//trim(fname)//'***'
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"lai",varid))

    call check(nf90_get_var(ncid,varid,lai1,start =(/1,1,day/),count = (/nlp,PFT,1/)))
    call check(nf90_close(ncid))
    !end if

    !--iLab::avoid pointer (see below)
    ! p => bound

    !call mpi_barrier(mpi_comm_world,ierr)
    !do i = 1,PFT
    !    call mpi_scatterv(lai1(1,i),dp,sp,mpi_real8,p%lai(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !end do

    !call mpi_barrier(mpi_comm_world,ierr)
    bound%lai = lai1

  end subroutine read_lai_site

  !! Reading Vcmax for data assimilation use
  subroutine read_Vcmax(yr,mn,dd,tod)
    implicit none
    !--iLab::made yr/mon/day/sec arguments in order to avoid 'get_curr_date'
    integer, intent(in) :: yr,mn,dd,tod
    ! integer          :: yr,mn,dd,tod,day
    real(r8)         :: vcmax1(nlon,nlat,PFT),vcmax2(nlon*nlat,PFT),vcmax3(nlp,PFT)
    character(len = 6) :: datestr
    integer          :: ncid,varid,ierr,i
    !--iLab::avoid pointer (see also below)
    ! type(surf),pointer:: p
    character(len=255) :: fname
    character(len=*), parameter :: method = "read_Vcmax"

    !if(myid ==0) then
    !   day  = get_curr_calday()
    ! call get_curr_date(yr,mn,dd,tod)
    write(datestr,'(i6)')  yr*100+mn
    write(*,*) "Reading Vcmax on "//datestr
    if( len(trim(beps_Vcmax_path))+len('Vcmax_')+len(trim(datestr))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_Vcmax_path)//"Vcmax_"//trim(datestr)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"Vcmax",varid))

    call check(nf90_get_var(ncid,varid,vcmax1))
    call check(nf90_close(ncid))

    vcmax2  = reshape(vcmax1,(/nlon*nlat,PFT/))
    vcmax3  = vcmax2(mapping,:)
    !end if

    !--iLab::avoid pointer (see also below)
    ! p => bound

    !call mpi_barrier(mpi_comm_world,ierr)
    !do i = 1,PFT
    !    call mpi_scatterv(vcmax3(1,i),dp,sp,mpi_real8,bound%Vcmax(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !end do
    !call mpi_barrier(mpi_comm_world,ierr)
    bound%Vcmax = vcmax3

  end subroutine read_Vcmax

  subroutine read_Vcmax_site(yr,mn,dd,tod)
    implicit none
    !--iLab::made yr/mon/day/sec arguments in order to avoid 'get_curr_date'
    integer, intent(in) :: yr,mn,dd,tod
    ! integer          :: yr,mn,dd,tod,day
    real(r8)         :: vcmax(nlp,PFT)
    character(len = 6) :: datestr
    integer          :: ncid,varid,ierr,i
    !--iLab::avoid pointer (see also below)
    ! type(surf),pointer:: p
    character(len=255) :: fname
    character(len=*), parameter :: method = "read_Vcmax_site"

    !if(myid ==0) then
    !   day  = get_curr_calday()
    ! call get_curr_date(yr,mn,dd,tod)
    write(datestr,'(i6)')  yr*100+mn
    write(*,*) "Reading Vcmax on "//datestr
    if( len(trim(beps_Vcmax_site_path))+len('Site_Vcmax_')+len(trim(datestr))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_Vcmax_site_path)//"Site_Vcmax_"//trim(datestr)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"Vcmax",varid))

    call check(nf90_get_var(ncid,varid,vcmax))
    call check(nf90_close(ncid))

    !end if
    !--iLab::avoid pointer (see also below)
    ! p => bound

    !call mpi_barrier(mpi_comm_world,ierr)
    !do i = 1,PFT
    !    call mpi_scatterv(vcmax(1,i),dp,sp,mpi_real8,bound%Vcmax(1,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
    !end do
    !call mpi_barrier(mpi_comm_world,ierr)
    bound%Vcmax = vcmax

  end subroutine read_Vcmax_site


  subroutine read_cpools()
    implicit none

    integer  :: ncid,varid(9),ierr,i
    real(r8) :: ccd(nlon,nlat,PFT),cfmd(nlon,nlat,PFT),cfsd(nlon,nlat,PFT),cm(nlon,nlat,PFT),cp(nlon,nlat,PFT),&
         cs(nlon,nlat,PFT),csm(nlon,nlat,PFT),csmd(nlon,nlat,PFT),cssd(nlon,nlat,PFT)
    real(r8) :: ccd2(nlon*nlat,PFT),cfmd2(nlon*nlat,PFT),cfsd2(nlon*nlat,PFT),cm2(nlon*nlat,PFT),cp2(nlon*nlat,PFT),&
         cs2(nlon*nlat,PFT),csm2(nlon*nlat,PFT),csmd2(nlon*nlat,PFT),cssd2(nlon*nlat,PFT)
    real(r8) :: ccd3(nlp,PFT),cfmd3(nlp,PFT),cfsd3(nlp,PFT),cm3(nlp,PFT),cp3(nlp,PFT),&
         cs3(nlp,PFT),csm3(nlp,PFT),csmd3(nlp,PFT),cssd3(nlp,PFT)
    !--iLab::avoid pointer (see also below)
    ! type(surf),pointer :: p

    !--iLab::avoid pointer (see also below)
    ! p => bound

    !if( myid ==0 ) then
    call check(nf90_open(trim(beps_cpools),nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"ccd",varid(1)))
    call check(nf90_inq_varid(ncid,"cfmd",varid(2)))
    call check(nf90_inq_varid(ncid,"cfsd",varid(3)))
    call check(nf90_inq_varid(ncid,"cm",varid(4)))
    call check(nf90_inq_varid(ncid,"cp",varid(5)))
    call check(nf90_inq_varid(ncid,"cs",varid(6)))
    call check(nf90_inq_varid(ncid,"csm",varid(7)))
    call check(nf90_inq_varid(ncid,"csmd",varid(8)))
    call check(nf90_inq_varid(ncid,"cssd",varid(9)))

    call check(nf90_get_var(ncid,varid(1),ccd))
    call check(nf90_get_var(ncid,varid(2),cfmd))
    call check(nf90_get_var(ncid,varid(3),cfsd))
    call check(nf90_get_var(ncid,varid(4),cm))
    call check(nf90_get_var(ncid,varid(5),cp))
    call check(nf90_get_var(ncid,varid(6),cs))
    call check(nf90_get_var(ncid,varid(7),csm))
    call check(nf90_get_var(ncid,varid(8),csmd))
    call check(nf90_get_var(ncid,varid(9),cssd))

    call check(nf90_close(ncid))

    ccd2    = reshape(ccd ,(/nlon*nlat,PFT/))
    cfmd2   = reshape(cfmd,(/nlon*nlat,PFT/))
    cfsd2   = reshape(cfsd,(/nlon*nlat,PFT/))
    cm2     = reshape(cm  ,(/nlon*nlat,PFT/))
    cp2     = reshape(cp  ,(/nlon*nlat,PFT/))
    cs2     = reshape(cs  ,(/nlon*nlat,PFT/))
    csm2    = reshape(csm ,(/nlon*nlat,PFT/))
    csmd2   = reshape(csmd,(/nlon*nlat,PFT/))
    cssd2   = reshape(cssd,(/nlon*nlat,PFT/))

    ccd3    = ccd2(mapping,:)*1e3    !! kg->g
    cfmd3   = cfmd2(mapping,:)*1e3
    cfsd3   = cfsd2(mapping,:)*1e3
    cm3     = cm2(mapping,:)*1e3
    cp3     = cp2(mapping,:)*1e3
    cs3     = cs2(mapping,:)*1e3
    csm3    = csm2(mapping,:)*1e3
    csmd3   = csmd2(mapping,:)*1e3
    cssd3   = cssd2(mapping,:)*1e3
    !end if

    !call mpi_barrier(mpi_comm_world,ierr)
    !do i = 1,PFT
    ! call mpi_scatterv(ccd3(1,i),dp,sp,MPI_real8,bound%ccd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(cfmd3(1,i),dp,sp,MPI_real8,bound%cfmd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(cfsd3(1,i),dp,sp,MPI_real8,bound%cfsd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(cm3(1,i),dp,sp,MPI_real8,bound%cm(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(cp3(1,i),dp,sp,MPI_real8,bound%cp(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(cs3(1,i),dp,sp,MPI_real8,bound%cs(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(csm3(1,i),dp,sp,MPI_real8,bound%csm(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(csmd3(1,i),dp,sp,MPI_real8,bound%csmd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    ! call mpi_scatterv(cssd3(1,i),dp,sp,MPI_real8,bound%cssd(1,i),npoints,MPI_real8,0,mpi_comm_world,ierr)
    !end do

    !call mpi_barrier(mpi_comm_world,ierr)
    bound%ccd = ccd3
    bound%cfmd = cfmd3
    bound%cfsd = cfsd3
    bound%cm = cm3
    bound%cp = cp3
    bound%cs = cs3
    bound%csm = csm3
    bound%csmd = csmd3
    bound%cssd = cssd3

  end subroutine read_cpools

  subroutine read_prior_para()
    implicit none
    integer  :: i,ncid,varid(36),ierr
    real(r8) :: p_Vcmax(PFT,nlp),p_q10(PFT,nlp),p_VJ_slope(PFT,nlp), &
         p_N_leaf(PFT,nlp),p_r_decay(PFT,nlp),p_b_h2o(PFT,nlp),&
         p_sif_alpha(PFT,nlp),p_sif_beta(PFT,nlp),p_taweff(PFT,nlp),&
         p_D0(PFT,nlp)
    real(r8) :: p_Ksat_scalar(texture,nlp),p_b_scalar(texture,nlp)
    real(r8) :: p_m_h2o(PFT,nlp)
    real(r8) :: p_f_leaf(PFT,nlp),p_kc25(PFT,nlp),p_ko25(PFT,nlp),p_tau25(PFT,nlp),p_agb2vod(PFT,nlp)
    
    real(r8) :: u_Vcmax(PFT,nlp),u_q10(PFT,nlp),u_VJ_slope(PFT,nlp),&
         u_N_leaf(PFT,nlp),u_r_decay(PFT,nlp),u_b_h2o(PFT,nlp),&
         u_sif_alpha(PFT,nlp),u_sif_beta(PFT,nlp),u_taweff(PFT,nlp),&
         u_D0(PFT,nlp)
    real(r8) :: u_Ksat_scalar(texture,nlp),u_b_scalar(texture,nlp)
    real(r8) :: u_m_h2o(PFT,nlp)
    real(r8) :: u_f_leaf(PFT,nlp),u_kc25(PFT,nlp),u_ko25(PFT,nlp),u_tau25(PFT,nlp),u_agb2vod(PFT,nlp)

    !--iLab::avoid pointer (see also below)
    ! type(para),pointer :: p
    character(len=255) :: fname
    character(len=*), parameter :: method = "read_prior_para"

    !--iLab::avoid pointer (see also below)
    ! p => assim

    if( len(trim(beps_site_path))+len(trim(prior_para_prefix))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//trim(prior_para_prefix)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"p_Vcmax",varid(1)))
    call check(nf90_inq_varid(ncid,"p_q10",varid(2)))
    call check(nf90_inq_varid(ncid,"p_VJ_slope",varid(3)))
    call check(nf90_inq_varid(ncid,"p_N_leaf",varid(4)))
    call check(nf90_inq_varid(ncid,"p_r_decay",varid(5)))
    call check(nf90_inq_varid(ncid,"p_b_h2o",varid(6)))
    call check(nf90_inq_varid(ncid,"p_sif_alpha",varid(7)))
    call check(nf90_inq_varid(ncid,"p_sif_beta",varid(8)))
    call check(nf90_inq_varid(ncid,"p_taweff",varid(9)))
    call check(nf90_inq_varid(ncid,"p_D0",varid(10)))
    call check(nf90_inq_varid(ncid,"p_Ksat_scalar",varid(11)))
    call check(nf90_inq_varid(ncid,"p_b_scalar",varid(12)))
    call check(nf90_inq_varid(ncid,"p_m_h2o",varid(13)))
    call check(nf90_inq_varid(ncid,"p_f_leaf",varid(14)))
    call check(nf90_inq_varid(ncid,"p_kc25",varid(15)))
    call check(nf90_inq_varid(ncid,"p_ko25",varid(16)))
    call check(nf90_inq_varid(ncid,"p_tau25",varid(17)))
    !call check(nf90_inq_varid(ncid,"p_f_lr",varid(14)))
    call check(nf90_inq_varid(ncid,"p_agb2vod",varid(18)))    

    call check(nf90_inq_varid(ncid,"u_Vcmax",varid(19)))
    call check(nf90_inq_varid(ncid,"u_q10",varid(20)))
    call check(nf90_inq_varid(ncid,"u_VJ_slope",varid(21)))
    call check(nf90_inq_varid(ncid,"u_N_leaf",varid(22)))
    call check(nf90_inq_varid(ncid,"u_r_decay",varid(23)))
    call check(nf90_inq_varid(ncid,"u_b_h2o",varid(24)))
    call check(nf90_inq_varid(ncid,"u_sif_alpha",varid(25)))
    call check(nf90_inq_varid(ncid,"u_sif_beta",varid(26)))
    call check(nf90_inq_varid(ncid,"u_taweff",varid(27)))
    call check(nf90_inq_varid(ncid,"u_D0",varid(28)))
    call check(nf90_inq_varid(ncid,"u_Ksat_scalar",varid(29)))
    call check(nf90_inq_varid(ncid,"u_b_scalar",varid(30)))
    call check(nf90_inq_varid(ncid,"u_m_h2o",varid(31)))
    call check(nf90_inq_varid(ncid,"u_f_leaf",varid(32)))
    call check(nf90_inq_varid(ncid,"u_kc25",varid(33)))
    call check(nf90_inq_varid(ncid,"u_ko25",varid(34)))
    call check(nf90_inq_varid(ncid,"u_tau25",varid(35)))
    !call check(nf90_inq_varid(ncid,"u_f_lr",varid(28)))
    call check(nf90_inq_varid(ncid,"u_agb2vod",varid(36)))

    call check(nf90_get_var(ncid,varid(1),p_Vcmax))
    call check(nf90_get_var(ncid,varid(2),p_q10))
    call check(nf90_get_var(ncid,varid(3),p_VJ_slope))
    call check(nf90_get_var(ncid,varid(4),p_N_leaf))
    call check(nf90_get_var(ncid,varid(5),p_r_decay))
    call check(nf90_get_var(ncid,varid(6),p_b_h2o))
    call check(nf90_get_var(ncid,varid(7),p_sif_alpha))
    call check(nf90_get_var(ncid,varid(8),p_sif_beta))
    call check(nf90_get_var(ncid,varid(9),p_taweff)) 
    call check(nf90_get_var(ncid,varid(10),p_D0))
    call check(nf90_get_var(ncid,varid(11),p_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(12),p_b_scalar))
    call check(nf90_get_var(ncid,varid(13),p_m_h2o))
    call check(nf90_get_var(ncid,varid(14),p_f_leaf))
    call check(nf90_get_var(ncid,varid(15),p_kc25))
    call check(nf90_get_var(ncid,varid(16),p_ko25))
    call check(nf90_get_var(ncid,varid(17),p_tau25))
   ! call check(nf90_get_var(ncid,varid(14),p_f_lr))
    call check(nf90_get_var(ncid,varid(18),p_agb2vod))
    call check(nf90_get_var(ncid,varid(19),u_Vcmax))
    call check(nf90_get_var(ncid,varid(20),u_q10))
    call check(nf90_get_var(ncid,varid(21),u_VJ_slope))
    call check(nf90_get_var(ncid,varid(22),u_N_leaf))
    call check(nf90_get_var(ncid,varid(23),u_r_decay))
    call check(nf90_get_var(ncid,varid(24),u_b_h2o))
    call check(nf90_get_var(ncid,varid(25),u_sif_alpha))
    call check(nf90_get_var(ncid,varid(26),u_sif_beta))
    call check(nf90_get_var(ncid,varid(27),u_taweff))   
    call check(nf90_get_var(ncid,varid(28),u_D0))
    call check(nf90_get_var(ncid,varid(29),u_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(30),u_b_scalar))
    call check(nf90_get_var(ncid,varid(31),u_m_h2o))
    call check(nf90_get_var(ncid,varid(32),u_f_leaf))
    call check(nf90_get_var(ncid,varid(33),u_kc25))
    call check(nf90_get_var(ncid,varid(34),u_ko25))
    call check(nf90_get_var(ncid,varid(35),u_tau25))
    !call check(nf90_get_var(ncid,varid(29),u_f_lr))
    call check(nf90_get_var(ncid,varid(36),u_agb2vod))

    call check(nf90_close(ncid))

    assim%p_Vcmax = p_Vcmax
    assim%p_q10 = p_q10
    assim%p_VJ_slope = p_VJ_slope
    assim%p_N_leaf = p_N_leaf
    assim%p_r_decay = p_r_decay
    assim%p_b_h2o = p_b_h2o
    assim%p_sif_alpha = p_sif_alpha
    assim%p_sif_beta = p_sif_beta
    assim%p_taweff = p_taweff
    assim%p_D0 = p_D0
    assim%p_Ksat_scalar = p_Ksat_scalar
    assim%p_b_scalar = p_b_scalar
    assim%p_m_h2o = p_m_h2o
    assim%p_f_leaf = p_f_leaf
    assim%p_kc25 = p_kc25
    assim%p_ko25 = p_ko25
    assim%p_tau25 = p_tau25
   ! assim%p_f_lr = p_f_lr
    assim%p_agb2vod = p_agb2vod
    assim%u_Vcmax = u_Vcmax
    assim%u_q10 = u_q10
    assim%u_VJ_slope = u_VJ_slope
    assim%u_N_leaf = u_N_leaf
    assim%u_r_decay = u_r_decay
    assim%u_b_h2o = u_b_h2o
    assim%u_sif_alpha = u_sif_alpha
    assim%u_sif_beta = u_sif_beta
    assim%u_taweff = u_taweff
    assim%u_D0 = u_D0
    assim%u_Ksat_scalar = u_Ksat_scalar
    assim%u_b_scalar = u_b_scalar
    assim%u_m_h2o = u_m_h2o
    assim%u_f_leaf = u_f_leaf
    assim%u_kc25 = u_kc25
    assim%u_ko25 = u_ko25
    assim%u_tau25 = u_tau25
    !assim%u_f_lr = u_f_lr
    assim%u_agb2vod = u_agb2vod

  end subroutine read_prior_para

  real function my_random (lbound,ubound)
      real, intent(in) :: lbound,ubound
      real :: len
      !real :: my_random
      real :: t
      
      !write(*,*) 'lbound=', lbound
      !write(*,*) 'ubound=', ubound
      len=ubound-lbound  !计算范围大小
      !write(*,*) 'len=', len
      call random_number(t)  !t是0-1之间的随机数
      my_random=lbound+len*t
      !write(*,*) 'my_random=', my_random
  end function my_random

  subroutine Create_PF_para(j,ii,pl)
    implicit none
    integer  :: i,ncid,varid(36),ierr,m,n
    integer, intent(in) :: j,ii,pl

    real(r8) :: p_Vcmax(PFT,nlp),p_q10(PFT,nlp),p_VJ_slope(PFT,nlp), &
         p_N_leaf(PFT,nlp),p_r_decay(PFT,nlp),p_b_h2o(PFT,nlp),&
         p_sif_alpha(PFT,nlp),p_sif_beta(PFT,nlp),p_taweff(PFT,nlp),&
         p_D0(PFT,nlp)
    real(r8) :: p_Ksat_scalar(texture,nlp),p_b_scalar(texture,nlp)
    real(r8) :: p_m_h2o(PFT,nlp)
    real(r8) :: p_f_leaf(PFT,nlp),p_kc25(PFT,nlp),p_ko25(PFT,nlp),p_tau25(PFT,nlp),p_agb2vod(PFT,nlp)
    
    real(r8) :: u_Vcmax(PFT,nlp),u_q10(PFT,nlp),u_VJ_slope(PFT,nlp),&
         u_N_leaf(PFT,nlp),u_r_decay(PFT,nlp),u_b_h2o(PFT,nlp),&
         u_sif_alpha(PFT,nlp),u_sif_beta(PFT,nlp),u_taweff(PFT,nlp),&
         u_D0(PFT,nlp)
    real(r8) :: u_Ksat_scalar(texture,nlp),u_b_scalar(texture,nlp)
    real(r8) :: u_m_h2o(PFT,nlp)
    real(r8) :: u_f_leaf(PFT,nlp),u_kc25(PFT,nlp),u_ko25(PFT,nlp),u_tau25(PFT,nlp),u_agb2vod(PFT,nlp)

    real(r8) :: Vcmax(parloop,nlp),q10(parloop,nlp),VJ_slope(parloop,nlp),N_leaf(parloop,nlp),&
                r_decay(parloop,nlp),b_h2o(parloop,nlp),&
                sif_alpha(parloop,nlp),sif_beta(parloop,nlp),taweff(parloop,nlp),D0(parloop,nlp)
    real(r8) :: Ksat_scalar(parloop,nlp),b_scalar(parloop,nlp)
    real(r8) :: m_h2o(parloop,nlp)
    real(r8) :: f_leaf(parloop,nlp),kc25(parloop,nlp),ko25(parloop,nlp),tau25(parloop,nlp),&
                agb2vod(parloop,nlp),pfweight(parloop,nlp)
    
    !real(r8) :: Vcmax(PFT,nlp),q10(PFT,nlp),VJ_slope(PFT,nlp), &
    !     N_leaf(PFT,nlp),r_decay(PFT,nlp),b_h2o(PFT,nlp),&
    !     sif_alpha(PFT,nlp),sif_beta(PFT,nlp),taweff(PFT,nlp),&
    !     D0(PFT,nlp)
    !real(r8) :: Ksat_scalar(texture,nlp),b_scalar(texture,nlp)
    !real(r8) :: m_h2o(PFT,nlp)
    !real(r8) :: f_leaf(nlp),kc25(nlp),ko25(nlp),tau25(nlp),agb2vod(nlp)
    !--iLab::avoid pointer (see also below)
    !type(PF_para),pointer :: p
    character(len=255) :: fname
    character(len=*), parameter :: method = "Create_PF_para"

    !--iLab::avoid pointer (see also below)
    !p => PF

    if( len(trim(beps_site_path))+len(trim(prior_para_prefix))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//trim(prior_para_prefix)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
    call check(nf90_inq_varid(ncid,"p_Vcmax",varid(1)))
    call check(nf90_inq_varid(ncid,"p_q10",varid(2)))
    call check(nf90_inq_varid(ncid,"p_VJ_slope",varid(3)))
    call check(nf90_inq_varid(ncid,"p_N_leaf",varid(4)))
    call check(nf90_inq_varid(ncid,"p_r_decay",varid(5)))
    call check(nf90_inq_varid(ncid,"p_b_h2o",varid(6)))
    call check(nf90_inq_varid(ncid,"p_sif_alpha",varid(7)))
    call check(nf90_inq_varid(ncid,"p_sif_beta",varid(8)))
    call check(nf90_inq_varid(ncid,"p_taweff",varid(9)))
    call check(nf90_inq_varid(ncid,"p_D0",varid(10)))
    call check(nf90_inq_varid(ncid,"p_Ksat_scalar",varid(11)))
    call check(nf90_inq_varid(ncid,"p_b_scalar",varid(12)))
    call check(nf90_inq_varid(ncid,"p_m_h2o",varid(13)))
    call check(nf90_inq_varid(ncid,"p_f_leaf",varid(14)))
    call check(nf90_inq_varid(ncid,"p_kc25",varid(15)))
    call check(nf90_inq_varid(ncid,"p_ko25",varid(16)))
    call check(nf90_inq_varid(ncid,"p_tau25",varid(17)))
    !call check(nf90_inq_varid(ncid,"p_f_lr",varid(14)))
    call check(nf90_inq_varid(ncid,"p_agb2vod",varid(18)))    

    call check(nf90_inq_varid(ncid,"u_Vcmax",varid(19)))
    call check(nf90_inq_varid(ncid,"u_q10",varid(20)))
    call check(nf90_inq_varid(ncid,"u_VJ_slope",varid(21)))
    call check(nf90_inq_varid(ncid,"u_N_leaf",varid(22)))
    call check(nf90_inq_varid(ncid,"u_r_decay",varid(23)))
    call check(nf90_inq_varid(ncid,"u_b_h2o",varid(24)))
    call check(nf90_inq_varid(ncid,"u_sif_alpha",varid(25)))
    call check(nf90_inq_varid(ncid,"u_sif_beta",varid(26)))
    call check(nf90_inq_varid(ncid,"u_taweff",varid(27)))
    call check(nf90_inq_varid(ncid,"u_D0",varid(28)))
    call check(nf90_inq_varid(ncid,"u_Ksat_scalar",varid(29)))
    call check(nf90_inq_varid(ncid,"u_b_scalar",varid(30)))
    call check(nf90_inq_varid(ncid,"u_m_h2o",varid(31)))
    call check(nf90_inq_varid(ncid,"u_f_leaf",varid(32)))
    call check(nf90_inq_varid(ncid,"u_kc25",varid(33)))
    call check(nf90_inq_varid(ncid,"u_ko25",varid(34)))
    call check(nf90_inq_varid(ncid,"u_tau25",varid(35)))
    !call check(nf90_inq_varid(ncid,"u_f_lr",varid(28)))
    call check(nf90_inq_varid(ncid,"u_agb2vod",varid(36)))

    call check(nf90_get_var(ncid,varid(1),p_Vcmax))
    call check(nf90_get_var(ncid,varid(2),p_q10))
    call check(nf90_get_var(ncid,varid(3),p_VJ_slope))
    call check(nf90_get_var(ncid,varid(4),p_N_leaf))
    call check(nf90_get_var(ncid,varid(5),p_r_decay))
    call check(nf90_get_var(ncid,varid(6),p_b_h2o))
    call check(nf90_get_var(ncid,varid(7),p_sif_alpha))
    call check(nf90_get_var(ncid,varid(8),p_sif_beta))
    call check(nf90_get_var(ncid,varid(9),p_taweff)) 
    call check(nf90_get_var(ncid,varid(10),p_D0))
    call check(nf90_get_var(ncid,varid(11),p_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(12),p_b_scalar))
    call check(nf90_get_var(ncid,varid(13),p_m_h2o))
    call check(nf90_get_var(ncid,varid(14),p_f_leaf))
    call check(nf90_get_var(ncid,varid(15),p_kc25))
    call check(nf90_get_var(ncid,varid(16),p_ko25))
    call check(nf90_get_var(ncid,varid(17),p_tau25))
   ! call check(nf90_get_var(ncid,varid(14),p_f_lr))
    call check(nf90_get_var(ncid,varid(18),p_agb2vod))
    call check(nf90_get_var(ncid,varid(19),u_Vcmax))
    call check(nf90_get_var(ncid,varid(20),u_q10))
    call check(nf90_get_var(ncid,varid(21),u_VJ_slope))
    call check(nf90_get_var(ncid,varid(22),u_N_leaf))
    call check(nf90_get_var(ncid,varid(23),u_r_decay))
    call check(nf90_get_var(ncid,varid(24),u_b_h2o))
    call check(nf90_get_var(ncid,varid(25),u_sif_alpha))
    call check(nf90_get_var(ncid,varid(26),u_sif_beta))
    call check(nf90_get_var(ncid,varid(27),u_taweff))   
    call check(nf90_get_var(ncid,varid(28),u_D0))
    call check(nf90_get_var(ncid,varid(29),u_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(30),u_b_scalar))
    call check(nf90_get_var(ncid,varid(31),u_m_h2o))
    call check(nf90_get_var(ncid,varid(32),u_f_leaf))
    call check(nf90_get_var(ncid,varid(33),u_kc25))
    call check(nf90_get_var(ncid,varid(34),u_ko25))
    call check(nf90_get_var(ncid,varid(35),u_tau25))
    !call check(nf90_get_var(ncid,varid(29),u_f_lr))
    call check(nf90_get_var(ncid,varid(36),u_agb2vod))

    call check(nf90_close(ncid))
    
    !write(*,*) 'p_Vcmax=', p_Vcmax
    !write(*,*) 'p_Vcmax dim1=', size(p_Vcmax,dim=1)
    !write(*,*) 'p_Vcmax dim2=', size(p_Vcmax,dim=2)

    !write(*,*) 'pf%Vcmax dim1=', size(PF%Vcmax,dim=1)

    do m=1,nlp
     do n=1,parloop
      call random_seed()!库函数，生成随机数前调用
      !write(*,*) 'p_Vcmax(j,m)=', p_Vcmax(j,m)
      !write(*,*) 'u_Vcmax(j,m)=', u_Vcmax(j,m)
      !write(*,*) 'PF%Vcmax=', PF%Vcmax(n,m)
      PF%Vcmax(n,m) = my_random(p_Vcmax(j,m),u_Vcmax(j,m))
      !write(*,*) 'PF%Vcmax=', PF%Vcmax(n,m)
      PF%q10(n,m) = my_random(p_q10(j,m),+u_q10(j,m))
      PF%VJ_slope(n,m) = my_random(p_VJ_slope(j,m),u_VJ_slope(j,m))
      !write(*,*) 'p_N_leaf(j,m)=', p_N_leaf(j,m)
      !write(*,*) 'u_N_leaf(j,m)=', u_N_leaf(j,m)
      PF%N_leaf(n,m) = my_random(p_N_leaf(j,m),u_N_leaf(j,m))
      PF%r_decay(n,m) =  my_random(p_r_decay(j,m),u_r_decay(j,m))
      PF%b_h2o(n,m) = my_random(p_b_h2o(j,m),u_b_h2o(j,m))
      PF%sif_alpha(n,m) = my_random(p_sif_alpha(j,m),u_sif_alpha(j,m))
      PF%sif_beta(n,m) =  my_random(p_sif_beta(j,m),u_sif_beta(j,m))
      PF%taweff(n,m) =  my_random(p_taweff(j,m),u_taweff(j,m))
      PF%D0(n,m) = my_random(p_D0(j,m),u_D0(j,m))
      PF%Ksat_scalar(n,m) = my_random(p_Ksat_scalar(ii,m),u_Ksat_scalar(ii,m))
      PF%b_scalar(n,m) = my_random(p_b_scalar(ii,m),u_b_scalar(ii,m))
      PF%m_h2o(n,m) =  my_random(p_m_h2o(j,m),u_m_h2o(j,m))
      PF%f_leaf(n,m) =  my_random(p_f_leaf(j,m),u_f_leaf(j,m))
      PF%kc25(n,m) =  my_random(p_kc25(j,m),u_kc25(j,m))
      PF%ko25(n,m) = my_random(p_ko25(j,m),u_ko25(j,m))
      PF%tau25(n,m) =  my_random(p_tau25(j,m),u_tau25(j,m))
      ! assim%p_f_lr = p_f_lr
      PF%agb2vod(n,m) =  my_random(p_agb2vod(j,m),u_agb2vod(j,m))
      PF%pfweight(n,m) =  1./parloop
     end do
    end do

    !PF%Vcmax = my_random(p_Vcmax(j)-u_Vcmax(j),p_Vcmax(j)+u_Vcmax(j))
    !PF%q10 = my_random(p_q10(j)-u_q10(j),p_q10(j)+u_q10(j))
    !PF%VJ_slope = my_random(p_VJ_slope(j)-u_VJ_slope(j),p_VJ_slope(j)+u_VJ_slope(j))
    !PF%N_leaf = my_random(p_N_leaf(j)-u_N_leaf(j),p_N_leaf(j)+u_N_leaf(j))
    !PF%r_decay = my_random(p_r_decay(j)-u_r_decay(j),p_r_decay(j)+u_r_decay(j))
    !PF%b_h2o = my_random(p_b_h2o(j)-u_b_h2o(j),p_b_h2o(j)+u_b_h2o(j))
    !PF%sif_alpha = my_random(p_sif_alpha(j)-u_sif_alpha(j),p_sif_alpha(j)+u_sif_alpha(j))
    !PF%sif_beta = my_random(p_sif_beta(j)-u_sif_beta(j),p_sif_beta(j)+u_sif_beta(j))
    !PF%taweff = my_random(p_taweff(j)-u_taweff(j),p_taweff(j)+u_taweff(j))
    !PF%D0 = my_random(p_D0(j)-u_D0(j),p_D0(j)+u_D0(j))
    !PF%Ksat_scalar = my_random(p_Ksat_scalar(ii)-u_Ksat_scalar(ii),p_Ksat_scalar(ii)+u_Ksat_scalar(ii))
    !PF%b_scalar = my_random(p_b_scalar(ii)-u_b_scalar(ii),p_b_scalar(ii)+u_b_scalar(ii))
    !PF%m_h2o = my_random(p_m_h2o(j)-u_m_h2o(j),p_m_h2o(j)+u_m_h2o(j))
    !PF%f_leaf = my_random(p_f_leaf(j)-u_f_leaf(j),p_f_leaf(j)+u_f_leaf(j))
    !PF%kc25 = my_random(p_kc25(j)-u_kc25(j),p_kc25(j)+u_kc25(j))
    !PF%ko25 = my_random(p_ko25(j)-u_ko25(j),p_ko25(j)+u_ko25(j))
    !PF%tau25 = my_random(p_tau25(j)-u_tau25(j),p_tau25(j)+u_tau25(j))
   ! assim%p_f_lr = p_f_lr
    !PF%agb2vod = my_random(p_agb2vod-u_agb2vod,p_agb2vod+u_agb2vod)

  end subroutine Create_PF_para

  subroutine read_PF_obs(nd)
    implicit none

    integer,intent(in)  :: nd
    integer          :: ncid,ierr,varid
    integer          :: yr,mon,day,sec
    integer          :: nt
    integer          :: i
    !--iLab::avoid pointer (see also below)
    ! type(forc),pointer :: p
    real(r8)           :: obs_GPP(nlp)
    !character(len=6)   :: monstr
    character(len=8)   :: datestr
    character(len=4)   :: yrstr
    logical :: ldebug
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_PF_obs'
    !--iLab::added missing initialisation
    ldebug = .false.

    !--iLab::avoid pointer (see also below)
    ! p => assim

    if( len(trim(beps_PF_obs_path))+len(trim(prior_PF_obs_prefix))+len('.nc') .gt.&
        len(fname) ) then
        write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
        stop
    else
        fname = trim(beps_PF_obs_path)//trim(prior_PF_obs_prefix)//".nc"
    endif
    !write(*,*) "nf90_open"
    call check(nf90_open(fname,nf90_nowrite,ncid))
    !write(*,*) "nf90_open"
    call check(nf90_inq_varid(ncid,"obs_GPP",varid))
    !write(*,*) "nf90_inq_varid"
    call check(nf90_get_var(ncid,varid,obs_GPP,start=(/1,nd/),count=(/nlp,1/)))
    !write(*,*) "nf90_get_var"
    call check(nf90_close(ncid))

    PF_obs%obs_GPP = obs_GPP !g C/m2/h

  end subroutine read_PF_obs


  subroutine check(status)
    implicit none
    integer, intent (in) :: status
    !write(*,*) 'get rid of check status'
    if(status /= nf90_noerr) then
       call endrun("Reading nc file is wrong!!")
    end if
  end subroutine check


end module controlInput_mod
