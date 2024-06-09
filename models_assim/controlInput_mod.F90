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
  character(len=255) :: meteo_path,meteo_flnm_prefix,meteo_site_flnm_prefix,prior_PF_obs_prefix,prior_MC_obs_prefix
  character(len=255) :: surface_data_path
  character(len=255) :: beps_yrdata_path
  character(len=255) :: beps_site_path,site_bound_prefix,prior_para_prefix,PF_prior_para_prefix,MC_prior_para_prefix
  character(len=255) :: beps_lai_path,beps_lai_prefix,beps_lai_site_prefix
  character(len=255) :: beps_Vcmax_path,beps_Vcmax_site_path
  character(len=255) :: beps_domain
  character(len=255) :: beps_cpools
  character(len=255) :: beps_PF_obs_path, beps_MC_obs_path
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
         prior_para_prefix,PF_prior_para_prefix,MC_prior_para_prefix,beps_cpools,beps_PF_obs_path,prior_PF_obs_prefix, &
         beps_MC_obs_path,prior_MC_obs_prefix,beps_out_dir,beps_rst_dir

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
    integer  :: i,ncid,varid(26),ierr
    integer  :: lcno(nlp,PFT),stext(nlp)
    real(r8) :: laiyr(nlp,PFT),nppyr(nlp,PFT)
    real(r8) :: PCT_PFT(nlp,PFT)
    real(r8) :: longitude(nlp),latitude(nlp),sdp(nlp),st(nlp),sw(nlp),ci(nlp)
    real(r8) :: ccd(nlp,PFT),cfmd(nlp,PFT),cfsd(nlp,PFT),cm(nlp,PFT),cp(nlp,PFT),&
         cs(nlp,PFT),csm(nlp,PFT),csmd(nlp,PFT),cssd(nlp,PFT)

    real(r8) :: tt_veg(nlp,PFT),tt_rep(nlp,PFT),phot_type(nlp,PFT),emer_doy(nlp,PFT),har_doy(nlp,PFT)
    real(r8) :: Hp(nlp)
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
    !print *, 'check(nf90_inq_varid(ncid,"cssd",varid(20)))!'
    !*****************************************************these are C4 Crop parameters  Xiuli*********************************
    call check(nf90_inq_varid(ncid,"tt_veg",varid(21)))   !crop module 	 thermal requirement of stage 1 of crop development (degree days).
    !print *, 'check(nf90_inq_varid(ncid,"tt_veg",varid(21)))!'
    call check(nf90_inq_varid(ncid,"tt_rep",varid(22)))   !crop module 	 thermal requirement of stage 2 of crop development (degree days).
    call check(nf90_inq_varid(ncid,"phot_type",varid(23)))  !crop module 	 photoperiod yes=1 not=0
    call check(nf90_inq_varid(ncid,"emer_doy",varid(24)))  !crop module emergence doy
    call check(nf90_inq_varid(ncid,"har_doy",varid(25)))  !crop module harvest doy
    ! 2024/03/29 canopy height
    call check(nf90_inq_varid(ncid,"Hp",varid(26)))
    !*****************************************************these are C4 Crop parameters  Xiuli*********************************

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
    !print *, 'check(nf90_get_var(ncid,varid(20),cssd))!'
    !*****************************************************these are C4 Crop parameters  Xiuli*********************************
    call check(nf90_get_var(ncid,varid(21),tt_veg))
    !print *, 'check(nf90_get_var(ncid,varid(21),tt_veg))!'
    call check(nf90_get_var(ncid,varid(22),tt_rep))
    call check(nf90_get_var(ncid,varid(23),phot_type))
    call check(nf90_get_var(ncid,varid(24),emer_doy))
    call check(nf90_get_var(ncid,varid(25),har_doy))
    call check(nf90_get_var(ncid,varid(26),Hp))
    !*****************************************************these are C4 Crop parameters  Xiuli*********************************
    call check(nf90_close(ncid))
    !end if
    PCT_PFT = PCT_PFT * 1.0         !! convert fraction to %
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
    ! 2024/03/29
    bound%HeightC = Hp

    !*****************************************************these are C4 Crop parameters  Xiuli*********************************
    pgdd%tt_veg = tt_veg
    pgdd%tt_rep= tt_rep
    pgdd%phot_type= phot_type
    pgdd%emer_doy= emer_doy
    pgdd%har_doy= har_doy

    !*****************************************************these are C4 Crop parameters  Xiuli*********************************

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

!..... copy from VOD-version ...................
 subroutine read_meteo_site_reftime()
    implicit none
    integer          :: ncid,timevar_id
    logical :: ldebug
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_meteo_site_reftime'
    character(len=*), parameter :: time_unit_expected = "hours since YYYY-MM-DD"
    character(len=128) :: time_unit
    ldebug = .false.

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
    call check(nf90_inq_varid(ncid,"time",timevar_id))
    if(ldebug) then
       write(*,'(a)') 'DEBUG::'//method//':data set identifier determined!'
    endif

    !-- iLab::get reference time in meteorological forcing
    call check(nf90_get_att(ncid, timevar_id, "units", time_unit))
    if(len(trim(time_unit)).ne.len(time_unit_expected)) then
       write(*, '(a)') ' FATAL::'//method//': unexpected unit of time in met forcing '//&
            '***'//trim(time_unit)//'***'
       stop
    else
       !-- expected: "hours since yyyy-mm-dd"
       clim%meteo_ref_yyyymmdd = time_unit(13:22)
    endif

    call check(nf90_close(ncid))
    if(ldebug) then
       write(*,'(a)') 'DEBUG::'//method//':NetCDF file closed again.'
    endif
  end subroutine read_meteo_site_reftime
 !............................................................

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
    logical :: exist
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

	! copy from VOD-version 203/06/30
	!-- iLab::added: stop in case meteo file is not present
    inquire(FILE=fname, exist=exist)
    if (.not.exist) then
       write(*, '(a)') ' FATAL::'//method//': file ***'//trim(fname)//'*** does NOT exist'
       stop
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
    logical :: exist
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

    ! copy from vod-version 2023/06/30
    !-- iLab::added: stop in case LAI file is not present
    inquire(FILE=fname, exist=exist)
    if (.not.exist) then
       write(*, '(a)') ' FATAL::'//method//': file ***'//trim(fname)//'*** does NOT exist'
       stop
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
    integer  :: i,ncid,varid(84),ierr  !  2024/03/17 42 parameters
    ! p means prior values; u means uncertainty
    real(r8) :: p_Vcmax(PFT,nlp),p_VJ_slope(PFT,nlp),p_VN_slope(PFT,nlp),p_b_h2o(PFT,nlp),&
         p_m_h2o(PFT,nlp),p_f_leaf(PFT,nlp),p_kc25(PFT,nlp),p_ko25(PFT,nlp),p_tau25(PFT,nlp)

    real(r8) :: p_sif_alpha(PFT,nlp),p_sif_beta(PFT,nlp)

    real(r8) :: p_q10(PFT,nlp),p_f_resp(PFT,nlp)

    real(r8) :: p_f_decay(PFT,nlp),p_Ksat_scalar(texture,nlp),p_b_scalar(texture,nlp),&
                p_porosity_scalar(texture,nlp),p_vfc_scalar(texture,nlp),p_vwp_scalar(texture,nlp),&
                p_psisat_scalar(texture,nlp),p_drainage_scalar(PFT,nlp)

    real(r8) :: p_vod_a(PFT,nlp),p_vod_b(PFT,nlp),p_vod_c(PFT,nlp)
    real(r8) :: p_theta_Amin(PFT,nlp),p_pox(PFT,nlp),p_fei_c(PFT,nlp),p_spac_p1(PFT,nlp),p_spac_p2(PFT,nlp),&
       p_tWA(PFT,nlp),p_tWB(PFT,nlp),p_Ttrig(PFT,nlp),p_r_xylem(PFT,nlp),p_r_r(PFT,nlp),p_Lr(PFT,nlp),&
       p_deltal_min(PFT,nlp),p_deltal_max(PFT,nlp),p_p_delta(PFT,nlp),p_ppslh(PFT,nlp),p_fei_min(PFT,nlp),&
       p_fei_th(PFT,nlp),p_p_excess(PFT,nlp)
    !real(r8) :: p_Tleaf_H(PFT,nlp),p_Tleaf_L(PFT,nlp),p_Tleaf_O(PFT,nlp)

    real(r8) :: u_Vcmax(PFT,nlp),u_VJ_slope(PFT,nlp),u_VN_slope(PFT,nlp),u_b_h2o(PFT,nlp),&
         u_m_h2o(PFT,nlp),u_f_leaf(PFT,nlp),u_kc25(PFT,nlp),u_ko25(PFT,nlp),u_tau25(PFT,nlp)

    real(r8) :: u_sif_alpha(PFT,nlp),u_sif_beta(PFT,nlp)

    real(r8) :: u_q10(PFT,nlp),u_f_resp(PFT,nlp)

    real(r8) :: u_f_decay(PFT,nlp),u_Ksat_scalar(texture,nlp),u_b_scalar(texture,nlp),&
                u_porosity_scalar(texture,nlp),u_vfc_scalar(texture,nlp),u_vwp_scalar(texture,nlp),&
                u_psisat_scalar(texture,nlp),u_drainage_scalar(PFT,nlp)

    real(r8) :: u_vod_a(PFT,nlp),u_vod_b(PFT,nlp),u_vod_c(PFT,nlp)
    real(r8) :: u_theta_Amin(PFT,nlp),u_pox(PFT,nlp),u_fei_c(PFT,nlp),u_spac_p1(PFT,nlp),u_spac_p2(PFT,nlp),&
       u_tWA(PFT,nlp),u_tWB(PFT,nlp),u_Ttrig(PFT,nlp),u_r_xylem(PFT,nlp),u_r_r(PFT,nlp),u_Lr(PFT,nlp),&
       u_deltal_min(PFT,nlp),u_deltal_max(PFT,nlp),u_p_delta(PFT,nlp),u_ppslh(PFT,nlp),u_fei_min(PFT,nlp),&
       u_fei_th(PFT,nlp),u_p_excess(PFT,nlp)
    !real(r8) :: u_Tleaf_H(PFT,nlp),u_Tleaf_L(PFT,nlp),u_Tleaf_O(PFT,nlp)

    !--iLab::avoid pointer (see also below)
    ! type(para),pointer :: p
    character(len=255) :: fname
    character(len=*), parameter :: method = "read_prior_para"
    logical :: exist  !copy from vod-version

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
 ! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
    call check(nf90_inq_varid(ncid,"p_Vcmax",varid(1)))
    call check(nf90_inq_varid(ncid,"p_VJ_slope",varid(2)))
    call check(nf90_inq_varid(ncid,"p_VN_slope",varid(3)))
    call check(nf90_inq_varid(ncid,"p_b_h2o",varid(4)))
    call check(nf90_inq_varid(ncid,"p_m_h2o",varid(5)))
    call check(nf90_inq_varid(ncid,"p_f_leaf",varid(6)))
    call check(nf90_inq_varid(ncid,"p_kc25",varid(7)))
    call check(nf90_inq_varid(ncid,"p_ko25",varid(8)))
    call check(nf90_inq_varid(ncid,"p_tau25",varid(9)))
 ! parameters for SIF modelling
    call check(nf90_inq_varid(ncid,"p_sif_alpha",varid(10)))
    call check(nf90_inq_varid(ncid,"p_sif_beta",varid(11)))
! parameters related to plant respiration: plant_resp
    call check(nf90_inq_varid(ncid,"p_q10",varid(12)))
! parameters related to soil respiration: soil_resp
    call check(nf90_inq_varid(ncid,"p_f_resp",varid(13)))
! soil parameters
    call check(nf90_inq_varid(ncid,"p_f_decay",varid(14)))
    call check(nf90_inq_varid(ncid,"p_Ksat_scalar",varid(15)))
    call check(nf90_inq_varid(ncid,"p_b_scalar",varid(16)))
    call check(nf90_inq_varid(ncid,"p_porosity_scalar",varid(17)))
    call check(nf90_inq_varid(ncid,"p_vfc_scalar",varid(18)))
    call check(nf90_inq_varid(ncid,"p_vwp_scalar",varid(19)))
    call check(nf90_inq_varid(ncid,"p_psisat_scalar",varid(20)))
    call check(nf90_inq_varid(ncid,"p_drainage_scalar",varid(21)))
 ! parameters for vod modelling
    call check(nf90_inq_varid(ncid,"p_vod_a",varid(22)))
    call check(nf90_inq_varid(ncid,"p_vod_b",varid(23)))
    call check(nf90_inq_varid(ncid,"p_vod_c",varid(24)))
! parameters for SPAC in soilwateruptake
    call check(nf90_inq_varid(ncid,"p_theta_Amin",varid(25)))
    call check(nf90_inq_varid(ncid,"p_pox",varid(26)))
    call check(nf90_inq_varid(ncid,"p_fei_c",varid(27)))
    call check(nf90_inq_varid(ncid,"p_spac_p1",varid(28)))
    call check(nf90_inq_varid(ncid,"p_spac_p2",varid(29)))
    call check(nf90_inq_varid(ncid,"p_tWA",varid(30)))
    call check(nf90_inq_varid(ncid,"p_tWB",varid(31)))
    call check(nf90_inq_varid(ncid,"p_Ttrig",varid(32)))
    call check(nf90_inq_varid(ncid,"p_r_xylem",varid(33)))
    call check(nf90_inq_varid(ncid,"p_r_r",varid(34)))
    call check(nf90_inq_varid(ncid,"p_Lr",varid(35)))
    call check(nf90_inq_varid(ncid,"p_deltal_min",varid(36)))
    call check(nf90_inq_varid(ncid,"p_deltal_max",varid(37)))
    call check(nf90_inq_varid(ncid,"p_p_delta",varid(38)))
    call check(nf90_inq_varid(ncid,"p_ppslh",varid(39)))
    call check(nf90_inq_varid(ncid,"p_fei_min",varid(40)))
    call check(nf90_inq_varid(ncid,"p_fei_th",varid(41)))
    call check(nf90_inq_varid(ncid,"p_p_excess",varid(42)))
    !call check(nf90_inq_varid(ncid,"p_Tleaf_H",varid(43)))
    !call check(nf90_inq_varid(ncid,"p_Tleaf_L",varid(44)))
    !call check(nf90_inq_varid(ncid,"p_Tleaf_O",varid(45)))

    call check(nf90_inq_varid(ncid,"u_Vcmax",varid(43)))
    call check(nf90_inq_varid(ncid,"u_VJ_slope",varid(44)))
    call check(nf90_inq_varid(ncid,"u_VN_slope",varid(45)))
    call check(nf90_inq_varid(ncid,"u_b_h2o",varid(46)))
    call check(nf90_inq_varid(ncid,"u_m_h2o",varid(47)))
    call check(nf90_inq_varid(ncid,"u_f_leaf",varid(48)))
    call check(nf90_inq_varid(ncid,"u_kc25",varid(49)))
    call check(nf90_inq_varid(ncid,"u_ko25",varid(50)))
    call check(nf90_inq_varid(ncid,"u_tau25",varid(51)))

    call check(nf90_inq_varid(ncid,"u_sif_alpha",varid(52)))
    call check(nf90_inq_varid(ncid,"u_sif_beta",varid(53)))

    call check(nf90_inq_varid(ncid,"u_q10",varid(54)))
    call check(nf90_inq_varid(ncid,"u_f_resp",varid(55)))

    call check(nf90_inq_varid(ncid,"u_f_decay",varid(56)))
    call check(nf90_inq_varid(ncid,"u_Ksat_scalar",varid(57)))
    call check(nf90_inq_varid(ncid,"u_b_scalar",varid(58)))
    call check(nf90_inq_varid(ncid,"u_porosity_scalar",varid(59)))
    call check(nf90_inq_varid(ncid,"u_vfc_scalar",varid(60)))
    call check(nf90_inq_varid(ncid,"u_vwp_scalar",varid(61)))
    call check(nf90_inq_varid(ncid,"u_psisat_scalar",varid(62)))
    call check(nf90_inq_varid(ncid,"u_drainage_scalar",varid(63)))

    call check(nf90_inq_varid(ncid,"u_vod_a",varid(64)))
    call check(nf90_inq_varid(ncid,"u_vod_b",varid(65)))
    call check(nf90_inq_varid(ncid,"u_vod_c",varid(66)))

    call check(nf90_inq_varid(ncid,"u_theta_Amin",varid(67)))
    call check(nf90_inq_varid(ncid,"u_pox",varid(68)))
    call check(nf90_inq_varid(ncid,"u_fei_c",varid(69)))
    call check(nf90_inq_varid(ncid,"u_spac_p1",varid(70)))
    call check(nf90_inq_varid(ncid,"u_spac_p2",varid(71)))
    call check(nf90_inq_varid(ncid,"u_tWA",varid(72)))
    call check(nf90_inq_varid(ncid,"u_tWB",varid(73)))
    call check(nf90_inq_varid(ncid,"u_Ttrig",varid(74)))
    call check(nf90_inq_varid(ncid,"u_r_xylem",varid(75)))
    call check(nf90_inq_varid(ncid,"u_r_r",varid(76)))
    call check(nf90_inq_varid(ncid,"u_Lr",varid(77)))
    call check(nf90_inq_varid(ncid,"u_deltal_min",varid(78)))
    call check(nf90_inq_varid(ncid,"u_deltal_max",varid(79)))
    call check(nf90_inq_varid(ncid,"u_p_delta",varid(80)))
    call check(nf90_inq_varid(ncid,"u_ppslh",varid(81)))
    call check(nf90_inq_varid(ncid,"u_fei_min",varid(82)))
    call check(nf90_inq_varid(ncid,"u_fei_th",varid(83)))
    call check(nf90_inq_varid(ncid,"u_p_excess",varid(84)))
 !   call check(nf90_inq_varid(ncid,"u_Tleaf_H",varid(88)))
 !   call check(nf90_inq_varid(ncid,"u_Tleaf_L",varid(89)))
  !  call check(nf90_inq_varid(ncid,"u_Tleaf_O",varid(90)))


    call check(nf90_get_var(ncid,varid(1),p_Vcmax))
    call check(nf90_get_var(ncid,varid(2),p_VJ_slope))
    call check(nf90_get_var(ncid,varid(3),p_VN_slope))
    call check(nf90_get_var(ncid,varid(4),p_b_h2o))
    call check(nf90_get_var(ncid,varid(5),p_m_h2o))
    call check(nf90_get_var(ncid,varid(6),p_f_leaf))
    call check(nf90_get_var(ncid,varid(7),p_kc25))
    call check(nf90_get_var(ncid,varid(8),p_ko25))
    call check(nf90_get_var(ncid,varid(9),p_tau25))

    call check(nf90_get_var(ncid,varid(10),p_sif_alpha))
    call check(nf90_get_var(ncid,varid(11),p_sif_beta))

    call check(nf90_get_var(ncid,varid(12),p_q10))
    call check(nf90_get_var(ncid,varid(13),p_f_resp))

    call check(nf90_get_var(ncid,varid(14),p_f_decay))
    call check(nf90_get_var(ncid,varid(15),p_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(16),p_b_scalar))
    call check(nf90_get_var(ncid,varid(17),p_porosity_scalar))
    call check(nf90_get_var(ncid,varid(18),p_vfc_scalar))
    call check(nf90_get_var(ncid,varid(19),p_vwp_scalar))
    call check(nf90_get_var(ncid,varid(20),p_psisat_scalar))
    call check(nf90_get_var(ncid,varid(21),p_drainage_scalar))

    call check(nf90_get_var(ncid,varid(22),p_vod_a))
    call check(nf90_get_var(ncid,varid(23),p_vod_b))
    call check(nf90_get_var(ncid,varid(24),p_vod_c))

    call check(nf90_get_var(ncid,varid(25),p_theta_Amin))
    call check(nf90_get_var(ncid,varid(26),p_pox))
    call check(nf90_get_var(ncid,varid(27),p_fei_c))
    call check(nf90_get_var(ncid,varid(28),p_spac_p1))
    call check(nf90_get_var(ncid,varid(29),p_spac_p2))
    call check(nf90_get_var(ncid,varid(30),p_tWA))
    call check(nf90_get_var(ncid,varid(31),p_tWB))
    call check(nf90_get_var(ncid,varid(32),p_Ttrig))
    call check(nf90_get_var(ncid,varid(33),p_r_xylem))
    call check(nf90_get_var(ncid,varid(34),p_r_r))
    call check(nf90_get_var(ncid,varid(35),p_Lr))
    call check(nf90_get_var(ncid,varid(36),p_deltal_min))
    call check(nf90_get_var(ncid,varid(37),p_deltal_max))
    call check(nf90_get_var(ncid,varid(38),p_p_delta))
    call check(nf90_get_var(ncid,varid(39),p_ppslh))
    call check(nf90_get_var(ncid,varid(40),p_fei_min))
    call check(nf90_get_var(ncid,varid(41),p_fei_th))
    call check(nf90_get_var(ncid,varid(42),p_p_excess))
    !call check(nf90_get_var(ncid,varid(43),p_Tleaf_H))
    !call check(nf90_get_var(ncid,varid(44),p_Tleaf_L))
    !call check(nf90_get_var(ncid,varid(45),p_Tleaf_O))

    call check(nf90_get_var(ncid,varid(43),u_Vcmax))
    call check(nf90_get_var(ncid,varid(44),u_VJ_slope))
    call check(nf90_get_var(ncid,varid(45),u_VN_slope))
    call check(nf90_get_var(ncid,varid(46),u_b_h2o))
    call check(nf90_get_var(ncid,varid(47),u_m_h2o))
    call check(nf90_get_var(ncid,varid(48),u_f_leaf))
    call check(nf90_get_var(ncid,varid(49),u_kc25))
    call check(nf90_get_var(ncid,varid(50),u_ko25))
    call check(nf90_get_var(ncid,varid(51),u_tau25))

    call check(nf90_get_var(ncid,varid(52),u_sif_alpha))
    call check(nf90_get_var(ncid,varid(53),u_sif_beta))

    call check(nf90_get_var(ncid,varid(54),u_q10))
    call check(nf90_get_var(ncid,varid(55),u_f_resp))

    call check(nf90_get_var(ncid,varid(56),u_f_decay))
    call check(nf90_get_var(ncid,varid(57),u_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(58),u_b_scalar))
    call check(nf90_get_var(ncid,varid(59),u_porosity_scalar))
    call check(nf90_get_var(ncid,varid(60),u_vfc_scalar))
    call check(nf90_get_var(ncid,varid(61),u_vwp_scalar))
    call check(nf90_get_var(ncid,varid(62),u_psisat_scalar))
    call check(nf90_get_var(ncid,varid(63),u_drainage_scalar))

    call check(nf90_get_var(ncid,varid(64),u_vod_a))
    call check(nf90_get_var(ncid,varid(65),u_vod_b))
    call check(nf90_get_var(ncid,varid(66),u_vod_c))

    call check(nf90_get_var(ncid,varid(67),u_theta_Amin))
    call check(nf90_get_var(ncid,varid(68),u_pox))
    call check(nf90_get_var(ncid,varid(69),u_fei_c))
    call check(nf90_get_var(ncid,varid(70),u_spac_p1))
    call check(nf90_get_var(ncid,varid(71),u_spac_p2))
    call check(nf90_get_var(ncid,varid(72),u_tWA))
    call check(nf90_get_var(ncid,varid(73),u_tWB))
    call check(nf90_get_var(ncid,varid(74),u_Ttrig))
    call check(nf90_get_var(ncid,varid(75),u_r_xylem))
    call check(nf90_get_var(ncid,varid(76),u_r_r))
    call check(nf90_get_var(ncid,varid(77),u_Lr))
    call check(nf90_get_var(ncid,varid(78),u_deltal_min))
    call check(nf90_get_var(ncid,varid(79),u_deltal_max))
    call check(nf90_get_var(ncid,varid(80),u_p_delta))
    call check(nf90_get_var(ncid,varid(81),u_ppslh))
    call check(nf90_get_var(ncid,varid(82),u_fei_min))
    call check(nf90_get_var(ncid,varid(83),u_fei_th))
    call check(nf90_get_var(ncid,varid(84),u_p_excess))
    !call check(nf90_get_var(ncid,varid(88),u_Tleaf_H))
    !call check(nf90_get_var(ncid,varid(89),u_Tleaf_L))
    !call check(nf90_get_var(ncid,varid(90),u_Tleaf_O))

    call check(nf90_close(ncid))

    assim%p_Vcmax   = p_Vcmax
    assim%p_VJ_slope   = p_VJ_slope
    assim%p_VN_slope   = p_VN_slope
    assim%p_b_h2o = p_b_h2o
    assim%p_m_h2o = p_m_h2o
    assim%p_f_leaf   = p_f_leaf
    assim%p_kc25   = p_kc25
    assim%p_ko25   = p_ko25
    assim%p_tau25   = p_tau25

    assim%p_sif_alpha   = p_sif_alpha
    assim%p_sif_beta   = p_sif_beta

    assim%p_q10     = p_q10
    assim%p_f_resp   = p_f_resp

    assim%p_f_decay   = p_f_decay
    assim%p_Ksat_scalar   = p_Ksat_scalar
    assim%p_b_scalar  = p_b_scalar
    assim%p_porosity_scalar  = p_porosity_scalar
    assim%p_vfc_scalar  = p_vfc_scalar
    assim%p_vwp_scalar  = p_vwp_scalar
    assim%p_psisat_scalar  = p_psisat_scalar
    assim%p_drainage_scalar  = p_drainage_scalar

    assim%p_vod_a   = p_vod_a
    assim%p_vod_b   = p_vod_b
    assim%p_vod_c   = p_vod_c

    assim%p_theta_Amin = p_theta_Amin
    assim%p_pox = p_pox
    assim%p_fei_c = p_fei_c
    assim%p_spac_p1 = p_spac_p1
    assim%p_spac_p2 = p_spac_p2
    assim%p_tWA = p_tWA
    assim%p_tWB = p_tWB
    assim%p_Ttrig = p_Ttrig
    assim%p_r_xylem = p_r_xylem
    assim%p_r_r = p_r_r
    assim%p_Lr = p_Lr
    assim%p_deltal_min = p_deltal_min
    assim%p_deltal_max = p_deltal_max
    assim%p_p_delta = p_p_delta
    assim%p_ppslh = p_ppslh
    assim%p_fei_min = p_fei_min
    assim%p_fei_th = p_fei_th
    assim%p_p_excess = p_p_excess
    !assim%p_Tleaf_H = p_Tleaf_H
    !assim%p_Tleaf_L = p_Tleaf_L
    !assim%p_Tleaf_O = p_Tleaf_O

    assim%u_Vcmax   = u_Vcmax
    assim%u_VJ_slope   = u_VJ_slope
    assim%u_VN_slope   = u_VN_slope
    assim%u_b_h2o = u_b_h2o
    assim%u_m_h2o = u_m_h2o
    assim%u_f_leaf   = u_f_leaf
    assim%u_kc25   = u_kc25
    assim%u_ko25   = u_ko25
    assim%u_tau25   = u_tau25

    assim%u_sif_alpha   = u_sif_alpha
    assim%u_sif_beta   = u_sif_beta

    assim%u_q10     = u_q10
    assim%u_f_resp   = u_f_resp

    assim%u_f_decay   = u_f_decay
    assim%u_Ksat_scalar   = u_Ksat_scalar
    assim%u_b_scalar  = u_b_scalar
    assim%u_porosity_scalar  = u_porosity_scalar
    assim%u_vfc_scalar  = u_vfc_scalar
    assim%u_vwp_scalar  = u_vwp_scalar
    assim%u_psisat_scalar  = u_psisat_scalar
    assim%u_drainage_scalar  = u_drainage_scalar

    assim%u_vod_a   = u_vod_a
    assim%u_vod_b   = u_vod_b
    assim%u_vod_c   = u_vod_c

    assim%u_theta_Amin = u_theta_Amin
    assim%u_pox = u_pox
    assim%u_fei_c = u_fei_c
    assim%u_spac_p1 = u_spac_p1
    assim%u_spac_p2 = u_spac_p2
    assim%u_tWA = u_tWA
    assim%u_tWB = u_tWB
    assim%u_Ttrig = u_Ttrig
    assim%u_r_xylem = u_r_xylem
    assim%u_r_r = u_r_r
    assim%u_Lr = u_Lr
    assim%u_deltal_min = u_deltal_min
    assim%u_deltal_max = u_deltal_max
    assim%u_p_delta = u_p_delta
    assim%u_ppslh = u_ppslh
    assim%u_fei_min = u_fei_min
    assim%u_fei_th = u_fei_th
    assim%u_p_excess = u_p_excess
    !assim%u_Tleaf_H = u_Tleaf_H
    !assim%u_Tleaf_L = u_Tleaf_L
    !assim%u_Tleaf_O = u_Tleaf_O

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
    integer  :: i,ncid,varid(90),ierr,m,n  ! 45*2 parameters
    integer, intent(in) :: j,ii,pl ! j, land cover index; ii, soil texture index
  ! ub means upper boundary; lb means lower boundary
  ! nlp means the number of sites
    real(r8) :: ub_Vcmax(PFT,nlp),ub_VJ_slope(PFT,nlp),ub_VN_slope(PFT,nlp),ub_b_h2o(PFT,nlp),&
         ub_m_h2o(PFT,nlp),ub_f_leaf(PFT,nlp),ub_kc25(PFT,nlp),ub_ko25(PFT,nlp),ub_tau25(PFT,nlp)

    real(r8) :: ub_sif_alpha(PFT,nlp),ub_sif_beta(PFT,nlp)

    real(r8) :: ub_q10(PFT,nlp),ub_f_resp(PFT,nlp)

    real(r8) :: ub_f_decay(PFT,nlp),ub_Ksat_scalar(texture,nlp),ub_b_scalar(texture,nlp),&
                ub_porosity_scalar(texture,nlp),ub_vfc_scalar(texture,nlp),ub_vwp_scalar(texture,nlp),&
                ub_psisat_scalar(texture,nlp),ub_drainage_scalar(PFT,nlp)

    real(r8) :: ub_vod_a(PFT,nlp),ub_vod_b(PFT,nlp),ub_vod_c(PFT,nlp)
    real(r8) :: ub_theta_Amin(PFT,nlp),ub_pox(PFT,nlp),ub_fei_c(PFT,nlp),ub_spac_p1(PFT,nlp),ub_spac_p2(PFT,nlp),&
       ub_tWA(PFT,nlp),ub_tWB(PFT,nlp),ub_Ttrig(PFT,nlp),ub_r_xylem(PFT,nlp),ub_r_r(PFT,nlp),ub_Lr(PFT,nlp),&
       ub_deltal_min(PFT,nlp),ub_deltal_max(PFT,nlp),ub_p_delta(PFT,nlp),ub_ppslh(PFT,nlp),ub_fei_min(PFT,nlp),&
       ub_fei_th(PFT,nlp),ub_p_excess(PFT,nlp),ub_Tleaf_H(PFT,nlp),ub_Tleaf_L(PFT,nlp),ub_Tleaf_O(PFT,nlp)

    real(r8) :: lb_Vcmax(PFT,nlp),lb_VJ_slope(PFT,nlp),lb_VN_slope(PFT,nlp),lb_b_h2o(PFT,nlp),&
         lb_m_h2o(PFT,nlp),lb_f_leaf(PFT,nlp),lb_kc25(PFT,nlp),lb_ko25(PFT,nlp),lb_tau25(PFT,nlp)

    real(r8) :: lb_sif_alpha(PFT,nlp),lb_sif_beta(PFT,nlp)

    real(r8) :: lb_q10(PFT,nlp),lb_f_resp(PFT,nlp)

    real(r8) :: lb_f_decay(PFT,nlp),lb_Ksat_scalar(texture,nlp),lb_b_scalar(texture,nlp),&
                lb_porosity_scalar(texture,nlp),lb_vfc_scalar(texture,nlp),lb_vwp_scalar(texture,nlp),&
                lb_psisat_scalar(texture,nlp),lb_drainage_scalar(PFT,nlp)

    real(r8) :: lb_vod_a(PFT,nlp),lb_vod_b(PFT,nlp),lb_vod_c(PFT,nlp)
    real(r8) :: lb_theta_Amin(PFT,nlp),lb_pox(PFT,nlp),lb_fei_c(PFT,nlp),lb_spac_p1(PFT,nlp),lb_spac_p2(PFT,nlp),&
       lb_tWA(PFT,nlp),lb_tWB(PFT,nlp),lb_Ttrig(PFT,nlp),lb_r_xylem(PFT,nlp),lb_r_r(PFT,nlp),lb_Lr(PFT,nlp),&
       lb_deltal_min(PFT,nlp),lb_deltal_max(PFT,nlp),lb_p_delta(PFT,nlp),lb_ppslh(PFT,nlp),lb_fei_min(PFT,nlp),&
       lb_fei_th(PFT,nlp),lb_p_excess(PFT,nlp),lb_Tleaf_H(PFT,nlp),lb_Tleaf_L(PFT,nlp),lb_Tleaf_O(PFT,nlp)

    !--iLab::avoid pointer (see also below)
    !type(PF_para),pointer :: p
    character(len=255) :: fname
    character(len=*), parameter :: method = "Create_PF_para"

    !--iLab::avoid pointer (see also below)
    !p => PF

    if( len(trim(beps_site_path))+len(trim(PF_prior_para_prefix))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//trim(PF_prior_para_prefix)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
    call check(nf90_inq_varid(ncid,"ub_Vcmax",varid(1)))
    call check(nf90_inq_varid(ncid,"ub_VJ_slope",varid(2)))
    call check(nf90_inq_varid(ncid,"ub_VN_slope",varid(3)))
    call check(nf90_inq_varid(ncid,"ub_b_h2o",varid(4)))
    call check(nf90_inq_varid(ncid,"ub_m_h2o",varid(5)))
    call check(nf90_inq_varid(ncid,"ub_f_leaf",varid(6)))
    call check(nf90_inq_varid(ncid,"ub_kc25",varid(7)))
    call check(nf90_inq_varid(ncid,"ub_ko25",varid(8)))
    call check(nf90_inq_varid(ncid,"ub_tau25",varid(9)))
 ! parameters for SIF modelling
    call check(nf90_inq_varid(ncid,"ub_sif_alpha",varid(10)))
    call check(nf90_inq_varid(ncid,"ub_sif_beta",varid(11)))
! parameters related to plant respiration: plant_resp
    call check(nf90_inq_varid(ncid,"ub_q10",varid(12)))
! parameters related to soil respiration: soil_resp
    call check(nf90_inq_varid(ncid,"ub_f_resp",varid(13)))
! soil parameters
    call check(nf90_inq_varid(ncid,"ub_f_decay",varid(14)))
    call check(nf90_inq_varid(ncid,"ub_Ksat_scalar",varid(15)))
    call check(nf90_inq_varid(ncid,"ub_b_scalar",varid(16)))
    call check(nf90_inq_varid(ncid,"ub_porosity_scalar",varid(17)))
    call check(nf90_inq_varid(ncid,"ub_vfc_scalar",varid(18)))
    call check(nf90_inq_varid(ncid,"ub_vwp_scalar",varid(19)))
    call check(nf90_inq_varid(ncid,"ub_psisat_scalar",varid(20)))
    call check(nf90_inq_varid(ncid,"ub_drainage_scalar",varid(21)))
 ! parameters for vod modelling
    call check(nf90_inq_varid(ncid,"ub_vod_a",varid(22)))
    call check(nf90_inq_varid(ncid,"ub_vod_b",varid(23)))
    call check(nf90_inq_varid(ncid,"ub_vod_c",varid(24)))
! parameters for SPAC in soilwateruptake
    call check(nf90_inq_varid(ncid,"ub_theta_Amin",varid(25)))
    call check(nf90_inq_varid(ncid,"ub_pox",varid(26)))
    call check(nf90_inq_varid(ncid,"ub_fei_c",varid(27)))
    call check(nf90_inq_varid(ncid,"ub_spac_p1",varid(28)))
    call check(nf90_inq_varid(ncid,"ub_spac_p2",varid(29)))
    call check(nf90_inq_varid(ncid,"ub_tWA",varid(30)))
    call check(nf90_inq_varid(ncid,"ub_tWB",varid(31)))
    call check(nf90_inq_varid(ncid,"ub_Ttrig",varid(32)))
    call check(nf90_inq_varid(ncid,"ub_r_xylem",varid(33)))
    call check(nf90_inq_varid(ncid,"ub_r_r",varid(34)))
    call check(nf90_inq_varid(ncid,"ub_Lr",varid(35)))
    call check(nf90_inq_varid(ncid,"ub_deltal_min",varid(36)))
    call check(nf90_inq_varid(ncid,"ub_deltal_max",varid(37)))
    call check(nf90_inq_varid(ncid,"ub_p_delta",varid(38)))
    call check(nf90_inq_varid(ncid,"ub_ppslh",varid(39)))
    call check(nf90_inq_varid(ncid,"ub_fei_min",varid(40)))
    call check(nf90_inq_varid(ncid,"ub_fei_th",varid(41)))
    call check(nf90_inq_varid(ncid,"ub_p_excess",varid(42)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_H",varid(43)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_L",varid(44)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_O",varid(45)))

    call check(nf90_inq_varid(ncid,"lb_Vcmax",varid(46)))
    call check(nf90_inq_varid(ncid,"lb_VJ_slope",varid(47)))
    call check(nf90_inq_varid(ncid,"lb_VN_slope",varid(48)))
    call check(nf90_inq_varid(ncid,"lb_b_h2o",varid(49)))
    call check(nf90_inq_varid(ncid,"lb_m_h2o",varid(50)))
    call check(nf90_inq_varid(ncid,"lb_f_leaf",varid(51)))
    call check(nf90_inq_varid(ncid,"lb_kc25",varid(52)))
    call check(nf90_inq_varid(ncid,"lb_ko25",varid(53)))
    call check(nf90_inq_varid(ncid,"lb_tau25",varid(54)))

    call check(nf90_inq_varid(ncid,"lb_sif_alpha",varid(55)))
    call check(nf90_inq_varid(ncid,"lb_sif_beta",varid(56)))

    call check(nf90_inq_varid(ncid,"lb_q10",varid(57)))
    call check(nf90_inq_varid(ncid,"lb_f_resp",varid(58)))

    call check(nf90_inq_varid(ncid,"lb_f_decay",varid(59)))
    call check(nf90_inq_varid(ncid,"lb_Ksat_scalar",varid(60)))
    call check(nf90_inq_varid(ncid,"lb_b_scalar",varid(61)))
    call check(nf90_inq_varid(ncid,"lb_porosity_scalar",varid(62)))
    call check(nf90_inq_varid(ncid,"lb_vfc_scalar",varid(63)))
    call check(nf90_inq_varid(ncid,"lb_vwp_scalar",varid(64)))
    call check(nf90_inq_varid(ncid,"lb_psisat_scalar",varid(65)))
    call check(nf90_inq_varid(ncid,"lb_drainage_scalar",varid(66)))

    call check(nf90_inq_varid(ncid,"lb_vod_a",varid(67)))
    call check(nf90_inq_varid(ncid,"lb_vod_b",varid(68)))
    call check(nf90_inq_varid(ncid,"lb_vod_c",varid(69)))

    call check(nf90_inq_varid(ncid,"lb_theta_Amin",varid(70)))
    call check(nf90_inq_varid(ncid,"lb_pox",varid(71)))
    call check(nf90_inq_varid(ncid,"lb_fei_c",varid(72)))
    call check(nf90_inq_varid(ncid,"lb_spac_p1",varid(73)))
    call check(nf90_inq_varid(ncid,"lb_spac_p2",varid(74)))
    call check(nf90_inq_varid(ncid,"lb_tWA",varid(75)))
    call check(nf90_inq_varid(ncid,"lb_tWB",varid(76)))
    call check(nf90_inq_varid(ncid,"lb_Ttrig",varid(77)))
    call check(nf90_inq_varid(ncid,"lb_r_xylem",varid(78)))
    call check(nf90_inq_varid(ncid,"lb_r_r",varid(79)))
    call check(nf90_inq_varid(ncid,"lb_Lr",varid(80)))
    call check(nf90_inq_varid(ncid,"lb_deltal_min",varid(81)))
    call check(nf90_inq_varid(ncid,"lb_deltal_max",varid(82)))
    call check(nf90_inq_varid(ncid,"lb_p_delta",varid(83)))
    call check(nf90_inq_varid(ncid,"lb_ppslh",varid(84)))
    call check(nf90_inq_varid(ncid,"lb_fei_min",varid(85)))
    call check(nf90_inq_varid(ncid,"lb_fei_th",varid(86)))
    call check(nf90_inq_varid(ncid,"lb_p_excess",varid(87)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_H",varid(88)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_L",varid(89)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_O",varid(90)))

    call check(nf90_get_var(ncid,varid(1),ub_Vcmax))
    call check(nf90_get_var(ncid,varid(2),ub_VJ_slope))
    call check(nf90_get_var(ncid,varid(3),ub_VN_slope))
    call check(nf90_get_var(ncid,varid(4),ub_b_h2o))
    call check(nf90_get_var(ncid,varid(5),ub_m_h2o))
    call check(nf90_get_var(ncid,varid(6),ub_f_leaf))
    call check(nf90_get_var(ncid,varid(7),ub_kc25))
    call check(nf90_get_var(ncid,varid(8),ub_ko25))
    call check(nf90_get_var(ncid,varid(9),ub_tau25))

    call check(nf90_get_var(ncid,varid(10),ub_sif_alpha))
    call check(nf90_get_var(ncid,varid(11),ub_sif_beta))

    call check(nf90_get_var(ncid,varid(12),ub_q10))
    call check(nf90_get_var(ncid,varid(13),ub_f_resp))

    call check(nf90_get_var(ncid,varid(14),ub_f_decay))
    call check(nf90_get_var(ncid,varid(15),ub_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(16),ub_b_scalar))
    call check(nf90_get_var(ncid,varid(17),ub_porosity_scalar))
    call check(nf90_get_var(ncid,varid(18),ub_vfc_scalar))
    call check(nf90_get_var(ncid,varid(19),ub_vwp_scalar))
    call check(nf90_get_var(ncid,varid(20),ub_psisat_scalar))
    call check(nf90_get_var(ncid,varid(21),ub_drainage_scalar))

    call check(nf90_get_var(ncid,varid(22),ub_vod_a))
    call check(nf90_get_var(ncid,varid(23),ub_vod_b))
    call check(nf90_get_var(ncid,varid(24),ub_vod_c))

    call check(nf90_get_var(ncid,varid(25),ub_theta_Amin))
    call check(nf90_get_var(ncid,varid(26),ub_pox))
    call check(nf90_get_var(ncid,varid(27),ub_fei_c))
    call check(nf90_get_var(ncid,varid(28),ub_spac_p1))
    call check(nf90_get_var(ncid,varid(29),ub_spac_p2))
    call check(nf90_get_var(ncid,varid(30),ub_tWA))
    call check(nf90_get_var(ncid,varid(31),ub_tWB))
    call check(nf90_get_var(ncid,varid(32),ub_Ttrig))
    call check(nf90_get_var(ncid,varid(33),ub_r_xylem))
    call check(nf90_get_var(ncid,varid(34),ub_r_r))
    call check(nf90_get_var(ncid,varid(35),ub_Lr))
    call check(nf90_get_var(ncid,varid(36),ub_deltal_min))
    call check(nf90_get_var(ncid,varid(37),ub_deltal_max))
    call check(nf90_get_var(ncid,varid(38),ub_p_delta))
    call check(nf90_get_var(ncid,varid(39),ub_ppslh))
    call check(nf90_get_var(ncid,varid(40),ub_fei_min))
    call check(nf90_get_var(ncid,varid(41),ub_fei_th))
    call check(nf90_get_var(ncid,varid(42),ub_p_excess))
    call check(nf90_get_var(ncid,varid(43),ub_Tleaf_H))
    call check(nf90_get_var(ncid,varid(44),ub_Tleaf_L))
    call check(nf90_get_var(ncid,varid(45),ub_Tleaf_O))

    call check(nf90_get_var(ncid,varid(46),lb_Vcmax))
    call check(nf90_get_var(ncid,varid(47),lb_VJ_slope))
    call check(nf90_get_var(ncid,varid(48),lb_VN_slope))
    call check(nf90_get_var(ncid,varid(49),lb_b_h2o))
    call check(nf90_get_var(ncid,varid(50),lb_m_h2o))
    call check(nf90_get_var(ncid,varid(51),lb_f_leaf))
    call check(nf90_get_var(ncid,varid(52),lb_kc25))
    call check(nf90_get_var(ncid,varid(53),lb_ko25))
    call check(nf90_get_var(ncid,varid(54),lb_tau25))

    call check(nf90_get_var(ncid,varid(55),lb_sif_alpha))
    call check(nf90_get_var(ncid,varid(56),lb_sif_beta))

    call check(nf90_get_var(ncid,varid(57),lb_q10))
    call check(nf90_get_var(ncid,varid(58),lb_f_resp))

    call check(nf90_get_var(ncid,varid(59),lb_f_decay))
    call check(nf90_get_var(ncid,varid(60),lb_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(61),lb_b_scalar))
    call check(nf90_get_var(ncid,varid(62),lb_porosity_scalar))
    call check(nf90_get_var(ncid,varid(63),lb_vfc_scalar))
    call check(nf90_get_var(ncid,varid(64),lb_vwp_scalar))
    call check(nf90_get_var(ncid,varid(65),lb_psisat_scalar))
    call check(nf90_get_var(ncid,varid(66),lb_drainage_scalar))

    call check(nf90_get_var(ncid,varid(67),lb_vod_a))
    call check(nf90_get_var(ncid,varid(68),lb_vod_b))
    call check(nf90_get_var(ncid,varid(69),lb_vod_c))

    call check(nf90_get_var(ncid,varid(70),lb_theta_Amin))
    call check(nf90_get_var(ncid,varid(71),lb_pox))
    call check(nf90_get_var(ncid,varid(72),lb_fei_c))
    call check(nf90_get_var(ncid,varid(73),lb_spac_p1))
    call check(nf90_get_var(ncid,varid(74),lb_spac_p2))
    call check(nf90_get_var(ncid,varid(75),lb_tWA))
    call check(nf90_get_var(ncid,varid(76),lb_tWB))
    call check(nf90_get_var(ncid,varid(77),lb_Ttrig))
    call check(nf90_get_var(ncid,varid(78),lb_r_xylem))
    call check(nf90_get_var(ncid,varid(79),lb_r_r))
    call check(nf90_get_var(ncid,varid(80),lb_Lr))
    call check(nf90_get_var(ncid,varid(81),lb_deltal_min))
    call check(nf90_get_var(ncid,varid(82),lb_deltal_max))
    call check(nf90_get_var(ncid,varid(83),lb_p_delta))
    call check(nf90_get_var(ncid,varid(84),lb_ppslh))
    call check(nf90_get_var(ncid,varid(85),lb_fei_min))
    call check(nf90_get_var(ncid,varid(86),lb_fei_th))
    call check(nf90_get_var(ncid,varid(87),lb_p_excess))
    call check(nf90_get_var(ncid,varid(88),lb_Tleaf_H))
    call check(nf90_get_var(ncid,varid(89),lb_Tleaf_L))
    call check(nf90_get_var(ncid,varid(90),lb_Tleaf_O))

    call check(nf90_close(ncid))

    !write(*,*) 'pf%Vcmax dim1=', size(PF%Vcmax,dim=1)

    do m=1,nlp
     do n=1,parloop
      call random_seed()!Library function，generating random number
      ! currently do not considering parameters related to photosynthesis and respiration
      PF%Vcmax(n,m) = ub_Vcmax(j,m)
      PF%VJ_slope(n,m) = ub_VJ_slope(j,m)
      PF%VN_slope(n,m) = ub_VN_slope(j,m)
      PF%b_h2o(n,m) = ub_b_h2o(j,m)
      PF%m_h2o(n,m) =  ub_m_h2o(j,m)
      PF%f_leaf(n,m) =  ub_f_leaf(j,m)
      PF%kc25(n,m) =  ub_kc25(j,m)
      PF%ko25(n,m) = ub_ko25(j,m)
      PF%tau25(n,m) =  ub_tau25(j,m)

      PF%sif_alpha(n,m) = ub_sif_alpha(j,m)
      PF%sif_beta(n,m) =  ub_sif_beta(j,m)

      PF%q10(n,m) = ub_q10(j,m)
      PF%f_resp(n,m) = ub_f_resp(j,m)
      !...............................................................................

      PF%f_decay(n,m) = ub_f_decay(j,m)
      ! for sm modeling
      PF%Ksat_scalar(n,m) = 10.**my_random(lb_Ksat_scalar(ii,m),ub_Ksat_scalar(ii,m))
      PF%b_scalar(n,m) = my_random(lb_b_scalar(ii,m),ub_b_scalar(ii,m))
      PF%porosity_scalar(n,m)=my_random(lb_porosity_scalar(ii,m),ub_porosity_scalar(ii,m))
      PF%vfc_scalar(n,m)=my_random(lb_vfc_scalar(ii,m),ub_vfc_scalar(ii,m))
      PF%vwp_scalar(n,m)=my_random(lb_vwp_scalar(ii,m),ub_vwp_scalar(ii,m))
      PF%psisat_scalar(n,m)=my_random(lb_psisat_scalar(ii,m),ub_psisat_scalar(ii,m))
      PF%drainage_scalar(n,m)=my_random(lb_drainage_scalar(j,m),ub_drainage_scalar(j,m))

      PF%vod_a(n,m) = ub_vod_a(j,m)
      PF%vod_b(n,m) = ub_vod_b(j,m)
      PF%vod_c(n,m) = ub_vod_c(j,m)

      !PF%theta_Amin(n,m) = my_random(lb_theta_Amin(j,m),ub_theta_Amin(j,m))
      !PF%pox(n,m) = my_random(lb_pox(j,m),ub_pox(j,m))
      !PF%fei_c(n,m) = my_random(lb_fei_c(j,m),ub_fei_c(j,m))
      !PF%spac_p1(n,m) = my_random(lb_spac_p1(j,m),ub_spac_p1(j,m))
      !PF%spac_p2(n,m) = my_random(lb_spac_p2(j,m),ub_spac_p2(j,m))
      !PF%tWA(n,m) = my_random(lb_tWA(j,m),ub_tWA(j,m))
      !PF%tWB(n,m) = my_random(lb_tWB(j,m),ub_tWB(j,m))
      !PF%Ttrig(n,m) = my_random(lb_Ttrig(j,m),ub_Ttrig(j,m))
      !PF%r_xylem(n,m) = my_random(lb_r_xylem(j,m),ub_r_xylem(j,m))
      !PF%r_r(n,m) = 1000.*10.**my_random(lb_r_r(j,m),ub_r_r(j,m)) ! 1000.*10^b, b [-1,1]]
      !PF%Lr(n,m) = 0.1*10.**my_random(lb_Lr(j,m),ub_Lr(j,m))
      !PF%deltal_min(n,m) = my_random(lb_deltal_min(j,m),ub_deltal_min(j,m))
      !PF%deltal_max(n,m) = my_random(lb_deltal_max(j,m),ub_deltal_max(j,m))
      !PF%p_delta(n,m) = my_random(lb_p_delta(j,m),ub_p_delta(j,m))
      !PF%ppslh(n,m) = my_random(lb_ppslh(j,m),ub_ppslh(j,m))
      !PF%fei_min(n,m) = my_random(lb_fei_min(j,m),ub_fei_min(j,m))
      !PF%fei_th(n,m) = my_random(lb_fei_th(j,m),ub_fei_th(j,m))
      !PF%p_excess(n,m) = my_random(lb_p_excess(j,m),ub_p_excess(j,m))
      !PF%Tleaf_H(n,m) = my_random(lb_Tleaf_H(j,m),ub_Tleaf_H(j,m))
      !PF%Tleaf_L(n,m) = my_random(lb_Tleaf_L(j,m),ub_Tleaf_L(j,m))
      !PF%Tleaf_O(n,m) = my_random(lb_Tleaf_O(j,m),ub_Tleaf_O(j,m))

      PF%theta_Amin(n,m) = ub_theta_Amin(j,m)
      PF%pox(n,m) = ub_pox(j,m)
      PF%fei_c(n,m) = ub_fei_c(j,m)
      PF%spac_p1(n,m) = ub_spac_p1(j,m)
      PF%spac_p2(n,m) = ub_spac_p2(j,m)
      PF%tWA(n,m) = ub_tWA(j,m)
      PF%tWB(n,m) = ub_tWB(j,m)
      PF%Ttrig(n,m) = ub_Ttrig(j,m)
      PF%r_xylem(n,m) = ub_r_xylem(j,m)
      PF%r_r(n,m) = ub_r_r(j,m) ! 1000.*10^b, b [-1,1]]
      PF%Lr(n,m) = ub_Lr(j,m)
      PF%deltal_min(n,m) = ub_deltal_min(j,m)
      PF%deltal_max(n,m) = ub_deltal_max(j,m)
      PF%p_delta(n,m) = ub_p_delta(j,m)
      PF%ppslh(n,m) = ub_ppslh(j,m)
      PF%fei_min(n,m) = ub_fei_min(j,m)
      PF%fei_th(n,m) = ub_fei_th(j,m)
      PF%p_excess(n,m) = ub_p_excess(j,m)
      PF%Tleaf_H(n,m) = ub_Tleaf_H(j,m)
      PF%Tleaf_L(n,m) = ub_Tleaf_L(j,m)
      PF%Tleaf_O(n,m) = ub_Tleaf_O(j,m)

      PF%pfweight(n,m) =  1./parloop

     end do

    end do

  end subroutine Create_PF_para


  subroutine Create_MC_para(j,ii,MCnum)

    implicit none
    integer  :: i,ncid,varid(90),ierr,m,n  ! for LWP analysis, 45*2 parameters
    integer, intent(in) :: j,ii,MCnum ! j-> land cover index; ii->soil texture index
  ! ub means upper boundary; lb means lower boundary
    real(r8) :: ub_Vcmax(PFT,nlp),ub_VJ_slope(PFT,nlp),ub_VN_slope(PFT,nlp),ub_b_h2o(PFT,nlp),&
         ub_m_h2o(PFT,nlp),ub_f_leaf(PFT,nlp),ub_kc25(PFT,nlp),ub_ko25(PFT,nlp),ub_tau25(PFT,nlp)

    real(r8) :: ub_sif_alpha(PFT,nlp),ub_sif_beta(PFT,nlp)

    real(r8) :: ub_q10(PFT,nlp),ub_f_resp(PFT,nlp)

    real(r8) :: ub_f_decay(PFT,nlp),ub_Ksat_scalar(texture,nlp),ub_b_scalar(texture,nlp),&
                ub_porosity_scalar(texture,nlp),ub_vfc_scalar(texture,nlp),ub_vwp_scalar(texture,nlp),&
                ub_psisat_scalar(texture,nlp),ub_drainage_scalar(PFT,nlp)

    real(r8) :: ub_vod_a(PFT,nlp),ub_vod_b(PFT,nlp),ub_vod_c(PFT,nlp)
    real(r8) :: ub_theta_Amin(PFT,nlp),ub_pox(PFT,nlp),ub_fei_c(PFT,nlp),ub_spac_p1(PFT,nlp),ub_spac_p2(PFT,nlp),&
       ub_tWA(PFT,nlp),ub_tWB(PFT,nlp),ub_Ttrig(PFT,nlp),ub_r_xylem(PFT,nlp),ub_r_r(PFT,nlp),ub_Lr(PFT,nlp),&
       ub_deltal_min(PFT,nlp),ub_deltal_max(PFT,nlp),ub_p_delta(PFT,nlp),ub_ppslh(PFT,nlp),ub_fei_min(PFT,nlp),&
       ub_fei_th(PFT,nlp),ub_p_excess(PFT,nlp),ub_Tleaf_H(PFT,nlp),ub_Tleaf_L(PFT,nlp),ub_Tleaf_O(PFT,nlp)

    real(r8) :: lb_Vcmax(PFT,nlp),lb_VJ_slope(PFT,nlp),lb_VN_slope(PFT,nlp),lb_b_h2o(PFT,nlp),&
         lb_m_h2o(PFT,nlp),lb_f_leaf(PFT,nlp),lb_kc25(PFT,nlp),lb_ko25(PFT,nlp),lb_tau25(PFT,nlp)

    real(r8) :: lb_sif_alpha(PFT,nlp),lb_sif_beta(PFT,nlp)

    real(r8) :: lb_q10(PFT,nlp),lb_f_resp(PFT,nlp)

    real(r8) :: lb_f_decay(PFT,nlp),lb_Ksat_scalar(texture,nlp),lb_b_scalar(texture,nlp),&
                lb_porosity_scalar(texture,nlp),lb_vfc_scalar(texture,nlp),lb_vwp_scalar(texture,nlp),&
                lb_psisat_scalar(texture,nlp),lb_drainage_scalar(PFT,nlp)

    real(r8) :: lb_vod_a(PFT,nlp),lb_vod_b(PFT,nlp),lb_vod_c(PFT,nlp)
    real(r8) :: lb_theta_Amin(PFT,nlp),lb_pox(PFT,nlp),lb_fei_c(PFT,nlp),lb_spac_p1(PFT,nlp),lb_spac_p2(PFT,nlp),&
       lb_tWA(PFT,nlp),lb_tWB(PFT,nlp),lb_Ttrig(PFT,nlp),lb_r_xylem(PFT,nlp),lb_r_r(PFT,nlp),lb_Lr(PFT,nlp),&
       lb_deltal_min(PFT,nlp),lb_deltal_max(PFT,nlp),lb_p_delta(PFT,nlp),lb_ppslh(PFT,nlp),lb_fei_min(PFT,nlp),&
       lb_fei_th(PFT,nlp),lb_p_excess(PFT,nlp),lb_Tleaf_H(PFT,nlp),lb_Tleaf_L(PFT,nlp),lb_Tleaf_O(PFT,nlp)

    character(len=255) :: fname
    character(len=*), parameter :: method = "Create_MC_para"

    !--iLab::avoid pointer (see also below)
    !p => PF

    if( len(trim(beps_site_path))+len(trim(MC_prior_para_prefix))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//trim(MC_prior_para_prefix)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
 ! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
    call check(nf90_inq_varid(ncid,"ub_Vcmax",varid(1)))
    call check(nf90_inq_varid(ncid,"ub_VJ_slope",varid(2)))
    call check(nf90_inq_varid(ncid,"ub_VN_slope",varid(3)))
    call check(nf90_inq_varid(ncid,"ub_b_h2o",varid(4)))
    call check(nf90_inq_varid(ncid,"ub_m_h2o",varid(5)))
    call check(nf90_inq_varid(ncid,"ub_f_leaf",varid(6)))
    call check(nf90_inq_varid(ncid,"ub_kc25",varid(7)))
    call check(nf90_inq_varid(ncid,"ub_ko25",varid(8)))
    call check(nf90_inq_varid(ncid,"ub_tau25",varid(9)))
 ! parameters for SIF modelling
    call check(nf90_inq_varid(ncid,"ub_sif_alpha",varid(10)))
    call check(nf90_inq_varid(ncid,"ub_sif_beta",varid(11)))
! parameters related to plant respiration: plant_resp
    call check(nf90_inq_varid(ncid,"ub_q10",varid(12)))
! parameters related to soil respiration: soil_resp
    call check(nf90_inq_varid(ncid,"ub_f_resp",varid(13)))
! soil parameters
    call check(nf90_inq_varid(ncid,"ub_f_decay",varid(14)))
    call check(nf90_inq_varid(ncid,"ub_Ksat_scalar",varid(15)))
    call check(nf90_inq_varid(ncid,"ub_b_scalar",varid(16)))
    call check(nf90_inq_varid(ncid,"ub_porosity_scalar",varid(17)))
    call check(nf90_inq_varid(ncid,"ub_vfc_scalar",varid(18)))
    call check(nf90_inq_varid(ncid,"ub_vwp_scalar",varid(19)))
    call check(nf90_inq_varid(ncid,"ub_psisat_scalar",varid(20)))
    call check(nf90_inq_varid(ncid,"ub_drainage_scalar",varid(21)))
 ! parameters for vod modelling
    call check(nf90_inq_varid(ncid,"ub_vod_a",varid(22)))
    call check(nf90_inq_varid(ncid,"ub_vod_b",varid(23)))
    call check(nf90_inq_varid(ncid,"ub_vod_c",varid(24)))
! parameters for SPAC in soilwateruptake
    call check(nf90_inq_varid(ncid,"ub_theta_Amin",varid(25)))
    call check(nf90_inq_varid(ncid,"ub_pox",varid(26)))
    call check(nf90_inq_varid(ncid,"ub_fei_c",varid(27)))
    call check(nf90_inq_varid(ncid,"ub_spac_p1",varid(28)))
    call check(nf90_inq_varid(ncid,"ub_spac_p2",varid(29)))
    call check(nf90_inq_varid(ncid,"ub_tWA",varid(30)))
    call check(nf90_inq_varid(ncid,"ub_tWB",varid(31)))
    call check(nf90_inq_varid(ncid,"ub_Ttrig",varid(32)))
    call check(nf90_inq_varid(ncid,"ub_r_xylem",varid(33)))
    call check(nf90_inq_varid(ncid,"ub_r_r",varid(34)))
    call check(nf90_inq_varid(ncid,"ub_Lr",varid(35)))
    call check(nf90_inq_varid(ncid,"ub_deltal_min",varid(36)))
    call check(nf90_inq_varid(ncid,"ub_deltal_max",varid(37)))
    call check(nf90_inq_varid(ncid,"ub_p_delta",varid(38)))
    call check(nf90_inq_varid(ncid,"ub_ppslh",varid(39)))
    call check(nf90_inq_varid(ncid,"ub_fei_min",varid(40)))
    call check(nf90_inq_varid(ncid,"ub_fei_th",varid(41)))
    call check(nf90_inq_varid(ncid,"ub_p_excess",varid(42)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_H",varid(43)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_L",varid(44)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_O",varid(45)))

    call check(nf90_inq_varid(ncid,"lb_Vcmax",varid(46)))
    call check(nf90_inq_varid(ncid,"lb_VJ_slope",varid(47)))
    call check(nf90_inq_varid(ncid,"lb_VN_slope",varid(48)))
    call check(nf90_inq_varid(ncid,"lb_b_h2o",varid(49)))
    call check(nf90_inq_varid(ncid,"lb_m_h2o",varid(50)))
    call check(nf90_inq_varid(ncid,"lb_f_leaf",varid(51)))
    call check(nf90_inq_varid(ncid,"lb_kc25",varid(52)))
    call check(nf90_inq_varid(ncid,"lb_ko25",varid(53)))
    call check(nf90_inq_varid(ncid,"lb_tau25",varid(54)))

    call check(nf90_inq_varid(ncid,"lb_sif_alpha",varid(55)))
    call check(nf90_inq_varid(ncid,"lb_sif_beta",varid(56)))

    call check(nf90_inq_varid(ncid,"lb_q10",varid(57)))
    call check(nf90_inq_varid(ncid,"lb_f_resp",varid(58)))

    call check(nf90_inq_varid(ncid,"lb_f_decay",varid(59)))
    call check(nf90_inq_varid(ncid,"lb_Ksat_scalar",varid(60)))
    call check(nf90_inq_varid(ncid,"lb_b_scalar",varid(61)))
    call check(nf90_inq_varid(ncid,"lb_porosity_scalar",varid(62)))
    call check(nf90_inq_varid(ncid,"lb_vfc_scalar",varid(63)))
    call check(nf90_inq_varid(ncid,"lb_vwp_scalar",varid(64)))
    call check(nf90_inq_varid(ncid,"lb_psisat_scalar",varid(65)))
    call check(nf90_inq_varid(ncid,"lb_drainage_scalar",varid(66)))

    call check(nf90_inq_varid(ncid,"lb_vod_a",varid(67)))
    call check(nf90_inq_varid(ncid,"lb_vod_b",varid(68)))
    call check(nf90_inq_varid(ncid,"lb_vod_c",varid(69)))

    call check(nf90_inq_varid(ncid,"lb_theta_Amin",varid(70)))
    call check(nf90_inq_varid(ncid,"lb_pox",varid(71)))
    call check(nf90_inq_varid(ncid,"lb_fei_c",varid(72)))
    call check(nf90_inq_varid(ncid,"lb_spac_p1",varid(73)))
    call check(nf90_inq_varid(ncid,"lb_spac_p2",varid(74)))
    call check(nf90_inq_varid(ncid,"lb_tWA",varid(75)))
    call check(nf90_inq_varid(ncid,"lb_tWB",varid(76)))
    call check(nf90_inq_varid(ncid,"lb_Ttrig",varid(77)))
    call check(nf90_inq_varid(ncid,"lb_r_xylem",varid(78)))
    call check(nf90_inq_varid(ncid,"lb_r_r",varid(79)))
    call check(nf90_inq_varid(ncid,"lb_Lr",varid(80)))
    call check(nf90_inq_varid(ncid,"lb_deltal_min",varid(81)))
    call check(nf90_inq_varid(ncid,"lb_deltal_max",varid(82)))
    call check(nf90_inq_varid(ncid,"lb_p_delta",varid(83)))
    call check(nf90_inq_varid(ncid,"lb_ppslh",varid(84)))
    call check(nf90_inq_varid(ncid,"lb_fei_min",varid(85)))
    call check(nf90_inq_varid(ncid,"lb_fei_th",varid(86)))
    call check(nf90_inq_varid(ncid,"lb_p_excess",varid(87)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_H",varid(88)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_L",varid(89)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_O",varid(90)))


    call check(nf90_get_var(ncid,varid(1),ub_Vcmax))
    call check(nf90_get_var(ncid,varid(2),ub_VJ_slope))
    call check(nf90_get_var(ncid,varid(3),ub_VN_slope))
    call check(nf90_get_var(ncid,varid(4),ub_b_h2o))
    call check(nf90_get_var(ncid,varid(5),ub_m_h2o))
    call check(nf90_get_var(ncid,varid(6),ub_f_leaf))
    call check(nf90_get_var(ncid,varid(7),ub_kc25))
    call check(nf90_get_var(ncid,varid(8),ub_ko25))
    call check(nf90_get_var(ncid,varid(9),ub_tau25))

    call check(nf90_get_var(ncid,varid(10),ub_sif_alpha))
    call check(nf90_get_var(ncid,varid(11),ub_sif_beta))

    call check(nf90_get_var(ncid,varid(12),ub_q10))
    call check(nf90_get_var(ncid,varid(13),ub_f_resp))

    call check(nf90_get_var(ncid,varid(14),ub_f_decay))
    call check(nf90_get_var(ncid,varid(15),ub_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(16),ub_b_scalar))
    call check(nf90_get_var(ncid,varid(17),ub_porosity_scalar))
    call check(nf90_get_var(ncid,varid(18),ub_vfc_scalar))
    call check(nf90_get_var(ncid,varid(19),ub_vwp_scalar))
    call check(nf90_get_var(ncid,varid(20),ub_psisat_scalar))
    call check(nf90_get_var(ncid,varid(21),ub_drainage_scalar))

    call check(nf90_get_var(ncid,varid(22),ub_vod_a))
    call check(nf90_get_var(ncid,varid(23),ub_vod_b))
    call check(nf90_get_var(ncid,varid(24),ub_vod_c))

    call check(nf90_get_var(ncid,varid(25),ub_theta_Amin))
    call check(nf90_get_var(ncid,varid(26),ub_pox))
    call check(nf90_get_var(ncid,varid(27),ub_fei_c))
    call check(nf90_get_var(ncid,varid(28),ub_spac_p1))
    call check(nf90_get_var(ncid,varid(29),ub_spac_p2))
    call check(nf90_get_var(ncid,varid(30),ub_tWA))
    call check(nf90_get_var(ncid,varid(31),ub_tWB))
    call check(nf90_get_var(ncid,varid(32),ub_Ttrig))
    call check(nf90_get_var(ncid,varid(33),ub_r_xylem))
    call check(nf90_get_var(ncid,varid(34),ub_r_r))
    call check(nf90_get_var(ncid,varid(35),ub_Lr))
    call check(nf90_get_var(ncid,varid(36),ub_deltal_min))
    call check(nf90_get_var(ncid,varid(37),ub_deltal_max))
    call check(nf90_get_var(ncid,varid(38),ub_p_delta))
    call check(nf90_get_var(ncid,varid(39),ub_ppslh))
    call check(nf90_get_var(ncid,varid(40),ub_fei_min))
    call check(nf90_get_var(ncid,varid(41),ub_fei_th))
    call check(nf90_get_var(ncid,varid(42),ub_p_excess))
    call check(nf90_get_var(ncid,varid(43),ub_Tleaf_H))
    call check(nf90_get_var(ncid,varid(44),ub_Tleaf_L))
    call check(nf90_get_var(ncid,varid(45),ub_Tleaf_O))

    call check(nf90_get_var(ncid,varid(46),lb_Vcmax))
    call check(nf90_get_var(ncid,varid(47),lb_VJ_slope))
    call check(nf90_get_var(ncid,varid(48),lb_VN_slope))
    call check(nf90_get_var(ncid,varid(49),lb_b_h2o))
    call check(nf90_get_var(ncid,varid(50),lb_m_h2o))
    call check(nf90_get_var(ncid,varid(51),lb_f_leaf))
    call check(nf90_get_var(ncid,varid(52),lb_kc25))
    call check(nf90_get_var(ncid,varid(53),lb_ko25))
    call check(nf90_get_var(ncid,varid(54),lb_tau25))

    call check(nf90_get_var(ncid,varid(55),lb_sif_alpha))
    call check(nf90_get_var(ncid,varid(56),lb_sif_beta))

    call check(nf90_get_var(ncid,varid(57),lb_q10))
    call check(nf90_get_var(ncid,varid(58),lb_f_resp))

    call check(nf90_get_var(ncid,varid(59),lb_f_decay))
    call check(nf90_get_var(ncid,varid(60),lb_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(61),lb_b_scalar))
    call check(nf90_get_var(ncid,varid(62),lb_porosity_scalar))
    call check(nf90_get_var(ncid,varid(63),lb_vfc_scalar))
    call check(nf90_get_var(ncid,varid(64),lb_vwp_scalar))
    call check(nf90_get_var(ncid,varid(65),lb_psisat_scalar))
    call check(nf90_get_var(ncid,varid(66),lb_drainage_scalar))

    call check(nf90_get_var(ncid,varid(67),lb_vod_a))
    call check(nf90_get_var(ncid,varid(68),lb_vod_b))
    call check(nf90_get_var(ncid,varid(69),lb_vod_c))

    call check(nf90_get_var(ncid,varid(70),lb_theta_Amin))
    call check(nf90_get_var(ncid,varid(71),lb_pox))
    call check(nf90_get_var(ncid,varid(72),lb_fei_c))
    call check(nf90_get_var(ncid,varid(73),lb_spac_p1))
    call check(nf90_get_var(ncid,varid(74),lb_spac_p2))
    call check(nf90_get_var(ncid,varid(75),lb_tWA))
    call check(nf90_get_var(ncid,varid(76),lb_tWB))
    call check(nf90_get_var(ncid,varid(77),lb_Ttrig))
    call check(nf90_get_var(ncid,varid(78),lb_r_xylem))
    call check(nf90_get_var(ncid,varid(79),lb_r_r))
    call check(nf90_get_var(ncid,varid(80),lb_Lr))
    call check(nf90_get_var(ncid,varid(81),lb_deltal_min))
    call check(nf90_get_var(ncid,varid(82),lb_deltal_max))
    call check(nf90_get_var(ncid,varid(83),lb_p_delta))
    call check(nf90_get_var(ncid,varid(84),lb_ppslh))
    call check(nf90_get_var(ncid,varid(85),lb_fei_min))
    call check(nf90_get_var(ncid,varid(86),lb_fei_th))
    call check(nf90_get_var(ncid,varid(87),lb_p_excess))
    call check(nf90_get_var(ncid,varid(88),lb_Tleaf_H))
    call check(nf90_get_var(ncid,varid(89),lb_Tleaf_L))
    call check(nf90_get_var(ncid,varid(90),lb_Tleaf_O))

    call check(nf90_close(ncid))

 !only soil hydraulic parameters were created by Monte Carlo

    do m=1,nlp
     do n=1,MCnum
      call random_seed()!Library function, generating random numbers

      MC%Vcmax(n,m) = ub_Vcmax(j,m)
      MC%VJ_slope(n,m) = ub_VJ_slope(j,m)
      MC%VN_slope(n,m) = ub_VN_slope(j,m)
      MC%b_h2o(n,m) = ub_b_h2o(j,m)
      MC%m_h2o(n,m) =  ub_m_h2o(j,m)
      MC%f_leaf(n,m) =  ub_f_leaf(j,m)
      MC%kc25(n,m) =  ub_kc25(j,m)
      MC%ko25(n,m) = ub_ko25(j,m)
      MC%tau25(n,m) =  ub_tau25(j,m)

      MC%sif_alpha(n,m) = ub_sif_alpha(j,m)
      MC%sif_beta(n,m) =  ub_sif_beta(j,m)

      MC%q10(n,m) = ub_q10(j,m)
      MC%f_resp(n,m) = ub_f_resp(j,m)

      MC%f_decay(n,m) = my_random(lb_f_decay(j,m),ub_f_decay(j,m))
      MC%Ksat_scalar(n,m) = 10.**my_random(lb_Ksat_scalar(ii,m),ub_Ksat_scalar(ii,m))
      MC%b_scalar(n,m) = my_random(lb_b_scalar(ii,m),ub_b_scalar(ii,m))
      MC%porosity_scalar(n,m)=my_random(lb_porosity_scalar(ii,m),ub_porosity_scalar(ii,m))
      MC%vfc_scalar(n,m)=my_random(lb_vfc_scalar(ii,m),ub_vfc_scalar(ii,m))
      MC%vwp_scalar(n,m)=my_random(lb_vwp_scalar(ii,m),ub_vwp_scalar(ii,m))
      MC%psisat_scalar(n,m)=my_random(lb_psisat_scalar(ii,m),ub_psisat_scalar(ii,m))
      MC%drainage_scalar(n,m)=my_random(lb_drainage_scalar(j,m),ub_drainage_scalar(j,m))

      MC%vod_a(n,m) = ub_vod_a(j,m)
      MC%vod_b(n,m) = ub_vod_b(j,m)
      MC%vod_c(n,m) = ub_vod_c(j,m)

      MC%theta_Amin(n,m) = my_random(lb_theta_Amin(j,m),ub_theta_Amin(j,m))
      MC%pox(n,m) = my_random(lb_pox(j,m),ub_pox(j,m))
      MC%fei_c(n,m) = my_random(lb_fei_c(j,m),ub_fei_c(j,m))
      MC%spac_p1(n,m) = my_random(lb_spac_p1(j,m),ub_spac_p1(j,m))
      MC%spac_p2(n,m) = my_random(lb_spac_p2(j,m),ub_spac_p2(j,m))
      MC%tWA(n,m) = my_random(lb_tWA(j,m),ub_tWA(j,m))
      MC%tWB(n,m) = my_random(lb_tWB(j,m),ub_tWB(j,m))
      MC%Ttrig(n,m) = my_random(lb_Ttrig(j,m),ub_Ttrig(j,m))
      MC%r_xylem(n,m) = my_random(lb_r_xylem(j,m),ub_r_xylem(j,m))
      MC%r_r(n,m) = 1000.*10.**my_random(lb_r_r(j,m),ub_r_r(j,m)) ! 1000.*10^b, b [-1,1]]
      MC%Lr(n,m) = 0.1*10.**my_random(lb_Lr(j,m),ub_Lr(j,m))
      MC%deltal_min(n,m) = my_random(lb_deltal_min(j,m),ub_deltal_min(j,m))
      MC%deltal_max(n,m) = my_random(lb_deltal_max(j,m),ub_deltal_max(j,m))
      MC%p_delta(n,m) = my_random(lb_p_delta(j,m),ub_p_delta(j,m))
      MC%ppslh(n,m) = my_random(lb_ppslh(j,m),ub_ppslh(j,m))
      MC%fei_min(n,m) = my_random(lb_fei_min(j,m),ub_fei_min(j,m))
      MC%fei_th(n,m) = my_random(lb_fei_th(j,m),ub_fei_th(j,m))
      MC%p_excess(n,m) = my_random(lb_p_excess(j,m),ub_p_excess(j,m))
      MC%Tleaf_H(n,m) = my_random(lb_Tleaf_H(j,m),ub_Tleaf_H(j,m))
      MC%Tleaf_L(n,m) = my_random(lb_Tleaf_L(j,m),ub_Tleaf_L(j,m))
      MC%Tleaf_O(n,m) = my_random(lb_Tleaf_O(j,m),ub_Tleaf_O(j,m))

     end do

    end do

  end subroutine Create_MC_para

  subroutine Create_PF_para_normal(j,ii,pl)
    implicit none
    integer  :: i,ncid,varid(90),ierr,m,n
    integer, intent(in) :: j,ii,pl
    real(r8) :: pi,temp_num, mean=0.5, sd=0.25

  ! ub means upper boundary; lb means lower boundary
  ! nlp means the number of sites
    real(r8) :: ub_Vcmax(PFT,nlp),ub_VJ_slope(PFT,nlp),ub_VN_slope(PFT,nlp),ub_b_h2o(PFT,nlp),&
         ub_m_h2o(PFT,nlp),ub_f_leaf(PFT,nlp),ub_kc25(PFT,nlp),ub_ko25(PFT,nlp),ub_tau25(PFT,nlp)

    real(r8) :: ub_sif_alpha(PFT,nlp),ub_sif_beta(PFT,nlp)

    real(r8) :: ub_q10(PFT,nlp),ub_f_resp(PFT,nlp)

    real(r8) :: ub_f_decay(PFT,nlp),ub_Ksat_scalar(texture,nlp),ub_b_scalar(texture,nlp),&
                ub_porosity_scalar(texture,nlp),ub_vfc_scalar(texture,nlp),ub_vwp_scalar(texture,nlp),&
                ub_psisat_scalar(texture,nlp),ub_drainage_scalar(PFT,nlp)

    real(r8) :: ub_vod_a(PFT,nlp),ub_vod_b(PFT,nlp),ub_vod_c(PFT,nlp)
    real(r8) :: ub_theta_Amin(PFT,nlp),ub_pox(PFT,nlp),ub_fei_c(PFT,nlp),ub_spac_p1(PFT,nlp),ub_spac_p2(PFT,nlp),&
       ub_tWA(PFT,nlp),ub_tWB(PFT,nlp),ub_Ttrig(PFT,nlp),ub_r_xylem(PFT,nlp),ub_r_r(PFT,nlp),ub_Lr(PFT,nlp),&
       ub_deltal_min(PFT,nlp),ub_deltal_max(PFT,nlp),ub_p_delta(PFT,nlp),ub_ppslh(PFT,nlp),ub_fei_min(PFT,nlp),&
       ub_fei_th(PFT,nlp),ub_p_excess(PFT,nlp),ub_Tleaf_H(PFT,nlp),ub_Tleaf_L(PFT,nlp),ub_Tleaf_O(PFT,nlp)

    real(r8) :: lb_Vcmax(PFT,nlp),lb_VJ_slope(PFT,nlp),lb_VN_slope(PFT,nlp),lb_b_h2o(PFT,nlp),&
         lb_m_h2o(PFT,nlp),lb_f_leaf(PFT,nlp),lb_kc25(PFT,nlp),lb_ko25(PFT,nlp),lb_tau25(PFT,nlp)

    real(r8) :: lb_sif_alpha(PFT,nlp),lb_sif_beta(PFT,nlp)

    real(r8) :: lb_q10(PFT,nlp),lb_f_resp(PFT,nlp)

    real(r8) :: lb_f_decay(PFT,nlp),lb_Ksat_scalar(texture,nlp),lb_b_scalar(texture,nlp),&
                lb_porosity_scalar(texture,nlp),lb_vfc_scalar(texture,nlp),lb_vwp_scalar(texture,nlp),&
                lb_psisat_scalar(texture,nlp),lb_drainage_scalar(PFT,nlp)

    real(r8) :: lb_vod_a(PFT,nlp),lb_vod_b(PFT,nlp),lb_vod_c(PFT,nlp)
    real(r8) :: lb_theta_Amin(PFT,nlp),lb_pox(PFT,nlp),lb_fei_c(PFT,nlp),lb_spac_p1(PFT,nlp),lb_spac_p2(PFT,nlp),&
       lb_tWA(PFT,nlp),lb_tWB(PFT,nlp),lb_Ttrig(PFT,nlp),lb_r_xylem(PFT,nlp),lb_r_r(PFT,nlp),lb_Lr(PFT,nlp),&
       lb_deltal_min(PFT,nlp),lb_deltal_max(PFT,nlp),lb_p_delta(PFT,nlp),lb_ppslh(PFT,nlp),lb_fei_min(PFT,nlp),&
       lb_fei_th(PFT,nlp),lb_p_excess(PFT,nlp),lb_Tleaf_H(PFT,nlp),lb_Tleaf_L(PFT,nlp),lb_Tleaf_O(PFT,nlp)

    real(r8) :: array(parloop)

    character(len=255) :: fname
    character(len=*), parameter :: method = "Create_PF_para"

    !--iLab::avoid pointer (see also below)
    !p => PF

    if( len(trim(beps_site_path))+len(trim(PF_prior_para_prefix))+len('.nc') .gt.&
         len(fname) ) then
       write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
       stop
    else
       fname = trim(beps_site_path)//trim(PF_prior_para_prefix)//".nc"
    endif
    call check(nf90_open(fname,nf90_nowrite,ncid))
 ! parameters related to plant photosynthesis: Vcmax_Jmax, photosynthesis
    call check(nf90_inq_varid(ncid,"ub_Vcmax",varid(1)))
    call check(nf90_inq_varid(ncid,"ub_VJ_slope",varid(2)))
    call check(nf90_inq_varid(ncid,"ub_VN_slope",varid(3)))
    call check(nf90_inq_varid(ncid,"ub_b_h2o",varid(4)))
    call check(nf90_inq_varid(ncid,"ub_m_h2o",varid(5)))
    call check(nf90_inq_varid(ncid,"ub_f_leaf",varid(6)))
    call check(nf90_inq_varid(ncid,"ub_kc25",varid(7)))
    call check(nf90_inq_varid(ncid,"ub_ko25",varid(8)))
    call check(nf90_inq_varid(ncid,"ub_tau25",varid(9)))
 ! parameters for SIF modelling
    call check(nf90_inq_varid(ncid,"ub_sif_alpha",varid(10)))
    call check(nf90_inq_varid(ncid,"ub_sif_beta",varid(11)))
! parameters related to plant respiration: plant_resp
    call check(nf90_inq_varid(ncid,"ub_q10",varid(12)))
! parameters related to soil respiration: soil_resp
    call check(nf90_inq_varid(ncid,"ub_f_resp",varid(13)))
! soil parameters
    call check(nf90_inq_varid(ncid,"ub_f_decay",varid(14)))
    call check(nf90_inq_varid(ncid,"ub_Ksat_scalar",varid(15)))
    call check(nf90_inq_varid(ncid,"ub_b_scalar",varid(16)))
    call check(nf90_inq_varid(ncid,"ub_porosity_scalar",varid(17)))
    call check(nf90_inq_varid(ncid,"ub_vfc_scalar",varid(18)))
    call check(nf90_inq_varid(ncid,"ub_vwp_scalar",varid(19)))
    call check(nf90_inq_varid(ncid,"ub_psisat_scalar",varid(20)))
    call check(nf90_inq_varid(ncid,"ub_drainage_scalar",varid(21)))
 ! parameters for vod modelling
    call check(nf90_inq_varid(ncid,"ub_vod_a",varid(22)))
    call check(nf90_inq_varid(ncid,"ub_vod_b",varid(23)))
    call check(nf90_inq_varid(ncid,"ub_vod_c",varid(24)))
! parameters for SPAC in soilwateruptake
    call check(nf90_inq_varid(ncid,"ub_theta_Amin",varid(25)))
    call check(nf90_inq_varid(ncid,"ub_pox",varid(26)))
    call check(nf90_inq_varid(ncid,"ub_fei_c",varid(27)))
    call check(nf90_inq_varid(ncid,"ub_spac_p1",varid(28)))
    call check(nf90_inq_varid(ncid,"ub_spac_p2",varid(29)))
    call check(nf90_inq_varid(ncid,"ub_tWA",varid(30)))
    call check(nf90_inq_varid(ncid,"ub_tWB",varid(31)))
    call check(nf90_inq_varid(ncid,"ub_Ttrig",varid(32)))
    call check(nf90_inq_varid(ncid,"ub_r_xylem",varid(33)))
    call check(nf90_inq_varid(ncid,"ub_r_r",varid(34)))
    call check(nf90_inq_varid(ncid,"ub_Lr",varid(35)))
    call check(nf90_inq_varid(ncid,"ub_deltal_min",varid(36)))
    call check(nf90_inq_varid(ncid,"ub_deltal_max",varid(37)))
    call check(nf90_inq_varid(ncid,"ub_p_delta",varid(38)))
    call check(nf90_inq_varid(ncid,"ub_ppslh",varid(39)))
    call check(nf90_inq_varid(ncid,"ub_fei_min",varid(40)))
    call check(nf90_inq_varid(ncid,"ub_fei_th",varid(41)))
    call check(nf90_inq_varid(ncid,"ub_p_excess",varid(42)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_H",varid(43)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_L",varid(44)))
    call check(nf90_inq_varid(ncid,"ub_Tleaf_O",varid(45)))

    call check(nf90_inq_varid(ncid,"lb_Vcmax",varid(46)))
    call check(nf90_inq_varid(ncid,"lb_VJ_slope",varid(47)))
    call check(nf90_inq_varid(ncid,"lb_VN_slope",varid(48)))
    call check(nf90_inq_varid(ncid,"lb_b_h2o",varid(49)))
    call check(nf90_inq_varid(ncid,"lb_m_h2o",varid(50)))
    call check(nf90_inq_varid(ncid,"lb_f_leaf",varid(51)))
    call check(nf90_inq_varid(ncid,"lb_kc25",varid(52)))
    call check(nf90_inq_varid(ncid,"lb_ko25",varid(53)))
    call check(nf90_inq_varid(ncid,"lb_tau25",varid(54)))

    call check(nf90_inq_varid(ncid,"lb_sif_alpha",varid(55)))
    call check(nf90_inq_varid(ncid,"lb_sif_beta",varid(56)))

    call check(nf90_inq_varid(ncid,"lb_q10",varid(57)))
    call check(nf90_inq_varid(ncid,"lb_f_resp",varid(58)))

    call check(nf90_inq_varid(ncid,"lb_f_decay",varid(59)))
    call check(nf90_inq_varid(ncid,"lb_Ksat_scalar",varid(60)))
    call check(nf90_inq_varid(ncid,"lb_b_scalar",varid(61)))
    call check(nf90_inq_varid(ncid,"lb_porosity_scalar",varid(62)))
    call check(nf90_inq_varid(ncid,"lb_vfc_scalar",varid(63)))
    call check(nf90_inq_varid(ncid,"lb_vwp_scalar",varid(64)))
    call check(nf90_inq_varid(ncid,"lb_psisat_scalar",varid(65)))
    call check(nf90_inq_varid(ncid,"lb_drainage_scalar",varid(66)))

    call check(nf90_inq_varid(ncid,"lb_vod_a",varid(67)))
    call check(nf90_inq_varid(ncid,"lb_vod_b",varid(68)))
    call check(nf90_inq_varid(ncid,"lb_vod_c",varid(69)))

    call check(nf90_inq_varid(ncid,"lb_theta_Amin",varid(70)))
    call check(nf90_inq_varid(ncid,"lb_pox",varid(71)))
    call check(nf90_inq_varid(ncid,"lb_fei_c",varid(72)))
    call check(nf90_inq_varid(ncid,"lb_spac_p1",varid(73)))
    call check(nf90_inq_varid(ncid,"lb_spac_p2",varid(74)))
    call check(nf90_inq_varid(ncid,"lb_tWA",varid(75)))
    call check(nf90_inq_varid(ncid,"lb_tWB",varid(76)))
    call check(nf90_inq_varid(ncid,"lb_Ttrig",varid(77)))
    call check(nf90_inq_varid(ncid,"lb_r_xylem",varid(78)))
    call check(nf90_inq_varid(ncid,"lb_r_r",varid(79)))
    call check(nf90_inq_varid(ncid,"lb_Lr",varid(80)))
    call check(nf90_inq_varid(ncid,"lb_deltal_min",varid(81)))
    call check(nf90_inq_varid(ncid,"lb_deltal_max",varid(82)))
    call check(nf90_inq_varid(ncid,"lb_p_delta",varid(83)))
    call check(nf90_inq_varid(ncid,"lb_ppslh",varid(84)))
    call check(nf90_inq_varid(ncid,"lb_fei_min",varid(85)))
    call check(nf90_inq_varid(ncid,"lb_fei_th",varid(86)))
    call check(nf90_inq_varid(ncid,"lb_p_excess",varid(87)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_H",varid(88)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_L",varid(89)))
    call check(nf90_inq_varid(ncid,"lb_Tleaf_O",varid(90)))

    call check(nf90_get_var(ncid,varid(1),ub_Vcmax))
    call check(nf90_get_var(ncid,varid(2),ub_VJ_slope))
    call check(nf90_get_var(ncid,varid(3),ub_VN_slope))
    call check(nf90_get_var(ncid,varid(4),ub_b_h2o))
    call check(nf90_get_var(ncid,varid(5),ub_m_h2o))
    call check(nf90_get_var(ncid,varid(6),ub_f_leaf))
    call check(nf90_get_var(ncid,varid(7),ub_kc25))
    call check(nf90_get_var(ncid,varid(8),ub_ko25))
    call check(nf90_get_var(ncid,varid(9),ub_tau25))

    call check(nf90_get_var(ncid,varid(10),ub_sif_alpha))
    call check(nf90_get_var(ncid,varid(11),ub_sif_beta))

    call check(nf90_get_var(ncid,varid(12),ub_q10))
    call check(nf90_get_var(ncid,varid(13),ub_f_resp))

    call check(nf90_get_var(ncid,varid(14),ub_f_decay))
    call check(nf90_get_var(ncid,varid(15),ub_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(16),ub_b_scalar))
    call check(nf90_get_var(ncid,varid(17),ub_porosity_scalar))
    call check(nf90_get_var(ncid,varid(18),ub_vfc_scalar))
    call check(nf90_get_var(ncid,varid(19),ub_vwp_scalar))
    call check(nf90_get_var(ncid,varid(20),ub_psisat_scalar))
    call check(nf90_get_var(ncid,varid(21),ub_drainage_scalar))

    call check(nf90_get_var(ncid,varid(22),ub_vod_a))
    call check(nf90_get_var(ncid,varid(23),ub_vod_b))
    call check(nf90_get_var(ncid,varid(24),ub_vod_c))

    call check(nf90_get_var(ncid,varid(25),ub_theta_Amin))
    call check(nf90_get_var(ncid,varid(26),ub_pox))
    call check(nf90_get_var(ncid,varid(27),ub_fei_c))
    call check(nf90_get_var(ncid,varid(28),ub_spac_p1))
    call check(nf90_get_var(ncid,varid(29),ub_spac_p2))
    call check(nf90_get_var(ncid,varid(30),ub_tWA))
    call check(nf90_get_var(ncid,varid(31),ub_tWB))
    call check(nf90_get_var(ncid,varid(32),ub_Ttrig))
    call check(nf90_get_var(ncid,varid(33),ub_r_xylem))
    call check(nf90_get_var(ncid,varid(34),ub_r_r))
    call check(nf90_get_var(ncid,varid(35),ub_Lr))
    call check(nf90_get_var(ncid,varid(36),ub_deltal_min))
    call check(nf90_get_var(ncid,varid(37),ub_deltal_max))
    call check(nf90_get_var(ncid,varid(38),ub_p_delta))
    call check(nf90_get_var(ncid,varid(39),ub_ppslh))
    call check(nf90_get_var(ncid,varid(40),ub_fei_min))
    call check(nf90_get_var(ncid,varid(41),ub_fei_th))
    call check(nf90_get_var(ncid,varid(42),ub_p_excess))
    call check(nf90_get_var(ncid,varid(43),ub_Tleaf_H))
    call check(nf90_get_var(ncid,varid(44),ub_Tleaf_L))
    call check(nf90_get_var(ncid,varid(45),ub_Tleaf_O))

    call check(nf90_get_var(ncid,varid(46),lb_Vcmax))
    call check(nf90_get_var(ncid,varid(47),lb_VJ_slope))
    call check(nf90_get_var(ncid,varid(48),lb_VN_slope))
    call check(nf90_get_var(ncid,varid(49),lb_b_h2o))
    call check(nf90_get_var(ncid,varid(50),lb_m_h2o))
    call check(nf90_get_var(ncid,varid(51),lb_f_leaf))
    call check(nf90_get_var(ncid,varid(52),lb_kc25))
    call check(nf90_get_var(ncid,varid(53),lb_ko25))
    call check(nf90_get_var(ncid,varid(54),lb_tau25))

    call check(nf90_get_var(ncid,varid(55),lb_sif_alpha))
    call check(nf90_get_var(ncid,varid(56),lb_sif_beta))

    call check(nf90_get_var(ncid,varid(57),lb_q10))
    call check(nf90_get_var(ncid,varid(58),lb_f_resp))

    call check(nf90_get_var(ncid,varid(59),lb_f_decay))
    call check(nf90_get_var(ncid,varid(60),lb_Ksat_scalar))
    call check(nf90_get_var(ncid,varid(61),lb_b_scalar))
    call check(nf90_get_var(ncid,varid(62),lb_porosity_scalar))
    call check(nf90_get_var(ncid,varid(63),lb_vfc_scalar))
    call check(nf90_get_var(ncid,varid(64),lb_vwp_scalar))
    call check(nf90_get_var(ncid,varid(65),lb_psisat_scalar))
    call check(nf90_get_var(ncid,varid(66),lb_drainage_scalar))

    call check(nf90_get_var(ncid,varid(67),lb_vod_a))
    call check(nf90_get_var(ncid,varid(68),lb_vod_b))
    call check(nf90_get_var(ncid,varid(69),lb_vod_c))

    call check(nf90_get_var(ncid,varid(70),lb_theta_Amin))
    call check(nf90_get_var(ncid,varid(71),lb_pox))
    call check(nf90_get_var(ncid,varid(72),lb_fei_c))
    call check(nf90_get_var(ncid,varid(73),lb_spac_p1))
    call check(nf90_get_var(ncid,varid(74),lb_spac_p2))
    call check(nf90_get_var(ncid,varid(75),lb_tWA))
    call check(nf90_get_var(ncid,varid(76),lb_tWB))
    call check(nf90_get_var(ncid,varid(77),lb_Ttrig))
    call check(nf90_get_var(ncid,varid(78),lb_r_xylem))
    call check(nf90_get_var(ncid,varid(79),lb_r_r))
    call check(nf90_get_var(ncid,varid(80),lb_Lr))
    call check(nf90_get_var(ncid,varid(81),lb_deltal_min))
    call check(nf90_get_var(ncid,varid(82),lb_deltal_max))
    call check(nf90_get_var(ncid,varid(83),lb_p_delta))
    call check(nf90_get_var(ncid,varid(84),lb_ppslh))
    call check(nf90_get_var(ncid,varid(85),lb_fei_min))
    call check(nf90_get_var(ncid,varid(86),lb_fei_th))
    call check(nf90_get_var(ncid,varid(87),lb_p_excess))
    call check(nf90_get_var(ncid,varid(88),lb_Tleaf_H))
    call check(nf90_get_var(ncid,varid(89),lb_Tleaf_L))
    call check(nf90_get_var(ncid,varid(90),lb_Tleaf_O))

    call check(nf90_close(ncid))

    !write(*,*) 'p_Vcmax=', p_Vcmax
    !write(*,*) 'p_Vcmax dim1=', size(p_Vcmax,dim=1)
    !write(*,*) 'p_Vcmax dim2=', size(p_Vcmax,dim=2)

    !write(*,*) 'pf%Vcmax dim1=', size(PF%Vcmax,dim=1)

    call random_seed() ! Library function,for generating random number
    call random_number(array) ! Uniform distribution

    ! Now convert to normal distribution
    DO i = 1, n-1, 2
       temp_num = sd * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
       array(i+1) = sd * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
       array(i) = temp_num
    END DO

    do m=1,nlp
     do n=1,parloop

      PF%Vcmax(n,m) = lb_Vcmax(j,m)+(ub_Vcmax(j,m)-lb_Vcmax(j,m))*array(n)
      PF%VJ_slope(n,m) = lb_VJ_slope(j,m)+(ub_VJ_slope(j,m)-lb_VJ_slope(j,m))*array(n)
      PF%VN_slope(n,m) = lb_VN_slope(j,m)+(ub_VN_slope(j,m)-lb_VN_slope(j,m))*array(n)
      PF%b_h2o(n,m) = lb_b_h2o(j,m)+(ub_b_h2o(j,m)-lb_b_h2o(j,m))*array(n)
      PF%m_h2o(n,m) =  lb_m_h2o(j,m)+(ub_m_h2o(j,m)-lb_m_h2o(j,m))*array(n)
      PF%f_leaf(n,m) =  lb_f_leaf(j,m)+(ub_f_leaf(j,m)-lb_f_leaf(j,m))*array(n)
      PF%kc25(n,m) =  lb_kc25(j,m)+(ub_kc25(j,m)-lb_kc25(j,m))*array(n)
      PF%ko25(n,m) = lb_ko25(j,m)+(ub_ko25(j,m)-lb_ko25(j,m))*array(n)
      PF%tau25(n,m) =  lb_tau25(j,m)+(ub_tau25(j,m)-lb_tau25(j,m))*array(n)

      PF%sif_alpha(n,m) = lb_sif_alpha(j,m)+(ub_sif_alpha(j,m)-lb_sif_alpha(j,m))*array(n)
      PF%sif_beta(n,m) =  lb_sif_beta(j,m)+(ub_sif_beta(j,m)-lb_sif_beta(j,m))*array(n)

      PF%q10(n,m) = lb_q10(j,m)+(ub_q10(j,m)-lb_q10(j,m))*array(n)
      PF%f_resp(n,m) = lb_f_resp(j,m)+(ub_f_resp(j,m)-lb_f_resp(j,m))*array(n)

      PF%f_decay(n,m) = lb_f_decay(j,m)+(ub_f_decay(j,m)-lb_f_decay(j,m))*array(n)
      PF%Ksat_scalar(n,m) = 10.**(lb_Ksat_scalar(ii,m)+(ub_Ksat_scalar(ii,m)-lb_Ksat_scalar(ii,m))*array(n))
      PF%b_scalar(n,m) = lb_b_scalar(ii,m)+(ub_b_scalar(ii,m)-lb_b_scalar(ii,m))*array(n)
      PF%porosity_scalar(n,m)=lb_porosity_scalar(ii,m)+(ub_porosity_scalar(ii,m)-lb_porosity_scalar(ii,m))*array(n)
      PF%vfc_scalar(n,m) = lb_vfc_scalar(ii,m)+(ub_vfc_scalar(ii,m)-lb_vfc_scalar(ii,m))*array(n)
      PF%vwp_scalar(n,m) = lb_vwp_scalar(ii,m)+(ub_vwp_scalar(ii,m)-lb_vwp_scalar(ii,m))*array(n)
      PF%psisat_scalar(n,m) = lb_psisat_scalar(ii,m)+(ub_psisat_scalar(ii,m)-lb_psisat_scalar(ii,m))*array(n)
      PF%drainage_scalar(n,m) = lb_drainage_scalar(ii,m)+(ub_drainage_scalar(ii,m)-lb_drainage_scalar(ii,m))*array(n)

      PF%vod_a(n,m) = lb_vod_a(j,m)+(ub_vod_a(j,m)-lb_vod_a(j,m))*array(n)
      PF%vod_b(n,m) = lb_vod_b(j,m)+(ub_vod_b(j,m)-lb_vod_b(j,m))*array(n)
      PF%vod_c(n,m) = lb_vod_c(j,m)+(ub_vod_c(j,m)-lb_vod_c(j,m))*array(n)

      PF%theta_Amin(n,m) = lb_theta_Amin(j,m)+(ub_theta_Amin(j,m)-lb_theta_Amin(j,m))*array(n)
      PF%pox(n,m) = lb_pox(j,m)+(ub_pox(j,m)-lb_pox(j,m))*array(n)
      PF%fei_c(n,m) = lb_fei_c(j,m)+(ub_fei_c(j,m)-lb_fei_c(j,m))*array(n)
      PF%spac_p1(n,m) = lb_spac_p1(j,m)+(ub_spac_p1(j,m)-lb_spac_p1(j,m))*array(n)
      PF%spac_p2(n,m) = lb_spac_p2(j,m)+(ub_spac_p2(j,m)-lb_spac_p2(j,m))*array(n)
      PF%tWA(n,m) = lb_tWA(j,m)+(ub_tWA(j,m)-lb_tWA(j,m))*array(n)
      PF%tWB(n,m) = lb_tWB(j,m)+(ub_tWB(j,m)-lb_tWB(j,m))*array(n)
      PF%Ttrig(n,m) = lb_Ttrig(j,m)+(ub_Ttrig(j,m)-lb_Ttrig(j,m))*array(n)
      PF%r_xylem(n,m) = lb_r_xylem(j,m)+(ub_r_xylem(j,m)-lb_r_xylem(j,m))*array(n)
      PF%r_r(n,m) = 1000.*10.**(lb_r_r(j,m)+(ub_r_r(j,m)-lb_r_r(j,m))*array(n)) ! 1000.*10^b, b [-1,1]]
      PF%Lr(n,m) = 0.1*10.**(lb_Lr(j,m)+(ub_Lr(j,m)-lb_Lr(j,m))*array(n))
      PF%deltal_min(n,m) = lb_deltal_min(j,m)+(ub_deltal_min(j,m)-lb_deltal_min(j,m))*array(n)
      PF%deltal_max(n,m) = lb_deltal_max(j,m)+(ub_deltal_max(j,m)-lb_deltal_max(j,m))*array(n)
      PF%p_delta(n,m) = lb_p_delta(j,m)+(ub_p_delta(j,m)-lb_p_delta(j,m))*array(n)
      PF%ppslh(n,m) = lb_ppslh(j,m)+(ub_ppslh(j,m)-lb_ppslh(j,m))*array(n)
      PF%fei_min(n,m) = lb_fei_min(j,m)+(ub_fei_min(j,m)-lb_fei_min(j,m))*array(n)
      PF%fei_th(n,m) = lb_fei_th(j,m)+(ub_fei_th(j,m)-lb_fei_th(j,m))*array(n)
      PF%p_excess(n,m) = lb_p_excess(j,m)+(ub_p_excess(j,m)-lb_p_excess(j,m))*array(n)
      PF%Tleaf_H(n,m) = lb_Tleaf_H(j,m)+(ub_Tleaf_H(j,m)-lb_Tleaf_H(j,m))*array(n)
      PF%Tleaf_L(n,m) = lb_Tleaf_L(j,m)+(ub_Tleaf_L(j,m)-lb_Tleaf_L(j,m))*array(n)
      PF%Tleaf_O(n,m) = lb_Tleaf_O(j,m)+(ub_Tleaf_O(j,m)-lb_Tleaf_O(j,m))*array(n)

      PF%pfweight(n,m) =  1./parloop
     end do
    end do
  end subroutine Create_PF_para_normal

  subroutine read_PF_obs(nd)
    implicit none

    integer,intent(in)  :: nd  ! nd is the time index
    integer          :: ncid,ierr,varid
    integer          :: yr,mon,day,sec
    integer          :: nt
    integer          :: i
    !--iLab::avoid pointer (see also below)
    ! type(forc),pointer :: p
    real(r8)           :: obs_var(nlp) ! changed to obs_VOD from obs_GPP
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
    call check(nf90_inq_varid(ncid,"obs_SM",varid)) ! var name is defined by user
    !write(*,*) "nf90_inq_varid"
    call check(nf90_get_var(ncid,varid,obs_var,start=(/1,nd/),count=(/nlp,1/)))
    !write(*,*) "nf90_get_var"
    call check(nf90_close(ncid))

    ! PF_obs%obs_GPP = obs_GPP
    PF_obs%obs_var = obs_var

  end subroutine read_PF_obs

  subroutine read_MC_obs(nd)
    ! nd is the time index, this function read one value for each time
    implicit none

    integer,intent(in)  :: nd
    integer          :: ncid,ierr,varid(2)
    integer          :: yr,mon,day,sec
    integer          :: nt
    integer          :: i
    !--iLab::avoid pointer (see also below)
    ! type(forc),pointer :: p
    real(r8)           :: obs_SM(nlp) ! for sm analysis
    real(r8)           :: obs_LWP(nlp) ! for leaf water potential analysis
    !character(len=6)   :: monstr
    character(len=8)   :: datestr
    character(len=4)   :: yrstr
    logical :: ldebug
    character(len=255) :: fname
    character(len=*), parameter :: method = 'read_MC_obs'
    !--iLab::added missing initialisation
    ldebug = .false.

    !--iLab::avoid pointer (see also below)
    ! p => assim

    if( len(trim(beps_MC_obs_path))+len(trim(prior_MC_obs_prefix))+len('.nc') .gt.&
        len(fname) ) then
        write(*, '(a)') 'FATAL::INTERNAL-ERROR::'//method//':file name too long!'
        stop
    else
        fname = trim(beps_MC_obs_path)//trim(prior_MC_obs_prefix)//".nc"
    endif
    !write(*,*) "nf90_open"
    call check(nf90_open(fname,nf90_nowrite,ncid))
    !write(*,*) "nf90_open"
    call check(nf90_inq_varid(ncid,"obs_SM",varid(1)))
    call check(nf90_inq_varid(ncid,"obs_LWP",varid(2)))
    !write(*,*) "nf90_inq_varid"
    call check(nf90_get_var(ncid,varid(1),obs_SM,start=(/1,nd/),count=(/nlp,1/)))
    call check(nf90_get_var(ncid,varid(2),obs_LWP,start=(/1,nd/),count=(/nlp,1/)))
    !write(*,*) "nf90_get_var"
    call check(nf90_close(ncid))

    MC_obs%obs_SM = obs_SM
    MC_obs%obs_LWP = obs_LWP

  end subroutine read_MC_obs

  subroutine check(status)
    implicit none
    integer, intent (in) :: status
    !write(*,*) 'get rid of check status'
    if(status /= nf90_noerr) then
       call endrun("Reading nc file is wrong!!")
    end if
  end subroutine check


end module controlInput_mod
