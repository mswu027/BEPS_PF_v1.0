

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file bepsfunc_setup.F90
!> \brief provides interfaces for setup and configuration of running
!>        the functional implementation of BEPS as
!>        suitable for use within sensitivity and optimisation framework.
!>
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  January 2020
!> \last  April 2021
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     initf
!
!> @brief determines sizes of control vector (n) and overall size of
!>        simulation vector (m)
!>        In addition necessary initialisations to actually run selected observational operator(s)
!>        need to be performed.
!
!> @details 
!
!> @param[out]  n  overall length of control vector
!> @param[out]  m  overall length of simulation vector
!
!> \authors Michael Vossbbeck, The Inversion Lab
!> \date    January 2020
subroutine initf(n, m)
  use shr_kind_mod, only: r8 =>shr_kind_r8
  use controlInput_mod
  use beps_time_manager, only:set_timemgr_init, timemgr_init, get_curr_date, &
       advance_timestep, is_last_step, get_doys, get_curr_calday, get_nstep, &
       get_start_date, timemgr_diff_secs
  use bepstype
  use bepstypeInit
  use restart
  use esmf, only:ESMF_INITIALIZE
  use mo_bepsfunc_ctl
  implicit none
  !- arguments
  integer, intent(out) :: n, m
  !-local variables
  character(len=*), parameter :: sub = 'initf::'
  logical :: ldebug
  integer :: min
  character(len=32) :: fmt

  !-- init
  ldebug = .false.
  n = -1
  m = -1

  !! setting up MPI enviroments with the namelists
  !   call Initmpi()

  call rdnamelist()
  if(nscale == 1) then
     nlp = n_site
     npoints = nlp
     write(*, '(a,2(a,i3,1x))') ' INFO::'//sub,'nlp=',nlp,'npoints=',npoints
     write(*,*) 'site points check', npoints
  end if

  !! Initialize the beps types
  call Initbepstype()

  !! Initialize output
  call Init_output

  !--------------------------------------------------------------------------
  ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
  ! because it is needed for the time manager, even if the ESMF_INTERFACE
  ! is not used.
  !--------------------------------------------------------------------------
  call ESMF_Initialize()

  !! setting time manager
  if(nsrest == nsrStartup) then
!!! calling time_manager set init
     call set_timemgr_init(calendar_in  = calendar,&
          start_ymd_in = icdate  ,&
          start_tod_in = icsec   ,&
          nelapse_in   = sim_duration,&
          dtime_in     = step)

  else if(nsrest == nsrContinue) then
     call restart_io("read")

     call set_timemgr_init(calendar_in  = calendar,   &
          start_ymd_in = rst_icdate, &
          start_tod_in = rst_icsec , &
          nelapse_in   = sim_duration,&
          dtime_in     = step)
  end if

  call timemgr_init()

  if (nscale == 0) then     ! nscale = 0 for global simulation, 1 for site simulation
     write(*, '(a)') ' INFO::'//sub//'BEPS run at global scale!'
     !! Reading boundary fields, yearly data,and soil Cpools for BEPS @J.Wang (note: if yearly and C pools data fields will change year by year, these datasets should be read in the time looping
     call read_boundary()
     call read_yrdata()
     call read_cpools()
  else
     write(*, '(a)') ' INFO::'//sub//'BEPS run at site scale!'
     call read_boundary_site()   ! read site data, including yrdata, boundary data, and carbon pools
  end if

  !-- get number of time steps
  call beps_time_setup(ntp)
  !-- iLab::we do have 'nsimvar' simulated variables,
  !         'ntp' timesteps and 'nlp' land points (i.e. sites)
  !         toplevel function 'bepsf_timesum' (see bepsfunc.F90) build temporal sum over
  !         variables:
  m = ntp*nlp*nsimvar
  !
  fmt = '(a,'//trim(ifmt(ntp))//')'
  write(*,fmt) ' INFO::'//sub//' ntp=',ntp
  fmt = '(a,'//trim(ifmt(nlp))//')'
  write(*,fmt) ' INFO::'//sub//' nlp=',nlp
  fmt = '(a,'//trim(ifmt(nsimvar))//')'
  write(*,fmt) ' INFO::'//sub//' nsimvar=',nsimvar
  fmt = '(a,1(a,'//trim(ifmt(m))//'),a)'
  write(*, fmt) ' INFO::'//sub//' ...dependents DONE. (','m=',m,')'

  !MOUSONG WU,2020-09-21
  ! 1. the parameters to be optimized have been modified, thes parameters include p_Ksat,p_b,these two are soil texture differentiated, in total 2*11,
  ! p_Vcmax,p_VJ_slope,sif_alpha,sif_beta,p_q10,p_D0,p_taueff, these are PFT differentiated, in total 7*9
  ! p_f_leaf,p_kc25,p_ko25,p_tau25,p_agbvod, in total 5*1
  ! Totally, 90 parameters are optimized in this new version. In practice, when we only focus on carbon fluxes, there might not be so many parameters to be optimized, for example,
  ! the parameters related to SIF, VOD will not be optimized. this means that we have a maximum of 90 parameters to be optimized in this model.
  ! these parameters have been assigned with prior values and uncertainties.
  !-- iLab::**7** PFT differentiated parameter
  !         (Vcmax,VJ_slope,Q10,sif_alpha,sif_beta,D0,taueff)
  !         **2** texture differentiated parameter (Ksat, b)
  !         **5** global parameter (f_leaf,kc25,ko25,tau25,agb2vod)
  write(*, '(a)') ' INFO::'//sub//' start determining number of parameter...'
  n = 7*PFT + 2*texture + 5
  !
  fmt = '(a,1(a,'//trim(ifmt(n))//'),a)'
  write(*, fmt) ' INFO::'//sub//' ...parameter DONE. (','n=',n,')'
  
  call readobs(m, min,'obs.b')
  
contains
  character(len=4) function ifmt(n) result(fmt)
    implicit none
    integer, intent(in) :: n
    if( n.lt.10 ) then
       fmt = 'i1'
    else if( n.lt.100 ) then
       fmt = 'i2'
    else if( n.lt.1000 ) then
       fmt = 'i3'
    else if( n.lt.10000 ) then
       fmt = 'i4'
    else if( n.lt.100000 ) then
       fmt = 'i5'
    else if( n.lt.1000000 ) then
       fmt = 'i6'
    else if( n.lt.10000000 ) then
       fmt = 'i7'
    else
       fmt = 'i15'
    endif
  end function ifmt

  subroutine beps_time_setup(nt)
    implicit none
    ! local declarations
    integer, intent(out) :: nt
    integer :: yr,mn,dy,tod
    integer :: yr_ref,mn_ref,dy_ref,tod_ref
    integer, allocatable :: ymds(:,:)
    integer :: it
    integer :: caldy, doys, kount
    character(len=*), parameter :: sub = 'beps_time_setup'

    !-- set timer to startup (very likely not necessary here)
    call timemgr_init()
    if( ldebug ) then
       open(1, file='beps_times_1.asc', form='formatted', action='write')
    endif

    !-- take the first
    
    !-- determine number of time-steps
    nt = 0
    timeloop: do
       nt = nt + 1
       call get_curr_date(yr,mn,dy,tod)
       !-- take the first date as reference time/date
       if( nt.eq.1 ) then
          if( tod.ne.0 ) then
             write(*, '(1x,a)') 'FATAL::'//sub//'tod was expected to be zero for reference date!'
             stop 1
          else
             yr_ref  = yr
             mn_ref  = mn
             dy_ref  = dy
             tod_ref = tod
             write(ref_date(1:4),   '(i4.4)') yr_ref
             write(ref_date(5:5),   '(a)')    '-'
             write(ref_date(6:7),   '(i2.2)') mn_ref
             write(ref_date(8:8),   '(a)')    '-'
             write(ref_date(9:10),  '(i2.2)') dy_ref
             write(ref_date(11:19), '(a)')    'T00:00:00'
          endif
       endif
       if( ldebug ) then
          caldy = get_curr_calday()
          doys = get_doys(yr)
          kount = get_nstep()
          write(1, '(a,i3,1x,a,i4,1x,2(a,i2,1x),a,i5,1x,2(a,i3,1x),a,i3)') 'nt=',nt,&
               'yr=',yr,'mn=',mn,'dy=',dy,'tod=',tod,&
               'caldy=',caldy,'doys=',doys,'kount=',kount
       endif
       !! advance time
       call advance_timestep()
       if(is_last_step()) exit timeloop
    end do timeloop !! end time loop
    if( ldebug ) then
       close(1)
    endif

    !-- allocate time-point array
    !   NOTE:we store one time-point more for restart purpose
    allocate( time_points(6,nt+1), seconds_since_ref(nt) )

    !-- reset timer to first time-step
    call timemgr_init()
    if( ldebug ) then
       open(1, file='beps_times_2.asc', form='formatted', action='write')
    endif
    it = 0
    ttloop: do
       it = it + 1
       call get_curr_date(yr,mn,dy,tod)
       if( ldebug ) then
          caldy = get_curr_calday()
          doys = get_doys(yr)
          kount = get_nstep()
          write(1, '(a,i3,1x,a,i4,1x,2(a,i2,1x),a,i5,1x,2(a,i3,1x),a,i3)') 'nt=',it,&
               'yr=',yr,'mn=',mn,'dy=',dy,'tod=',tod,&
               'caldy=',caldy,'doys=',doys,'kount=',kount
       endif
       time_points(1,it) = yr
       time_points(2,it) = mn
       time_points(3,it) = dy
       time_points(4,it) = tod
       time_points(5,it) = get_curr_calday()
       time_points(6,it) = get_doys(yr)
       !-- seconds elapsed since reference date/time
       call timemgr_diff_secs(yr_ref*10000+mn_ref*100+dy_ref, tod_ref, yr*10000+mn*100+dy, tod,&
            seconds_since_ref(it))
       !! advance time
       call advance_timestep()
       if(is_last_step()) then !-- save first time-point *after* simulation period (restart)
          call get_curr_date(yr,mn,dy,tod)
          time_points(1,it+1) = yr
          time_points(2,it+1) = mn
          time_points(3,it+1) = dy
          time_points(4,it+1) = tod
          time_points(5,it+1) = get_curr_calday()
          time_points(6,it+1) = get_doys(yr)
          exit ttloop
       end if
    enddo ttloop

    if( ldebug ) then
       close(1)
    endif

    !-- reset timer to first time-step
    call timemgr_init()

  end subroutine beps_time_setup

end subroutine initf


!***********************************************************
!     initx
!
!> @brief method to access the prior control vector in normalised units
!>        and the prior uncertainty in physical units
!>
!> @param[in]   n   overall length of control vector
!> @param[out]  x   control vector (in normalised units)
!> @param[out]  sx  uncertainty of control vector elements (physical units)
!
!> \authors Michael Vossbbeck, The Inversion Lab
!> \date    January 2020
subroutine initx(n, x, sx)
  use mo_prior
  use bepstype, only:assim
  use beps_par, only:PFT,texture
  use controlInput_mod, only:read_prior_para
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n), sx(n)
  ! local variables
  real(kind=8) :: xphys(n)
  integer :: i1,i2,ii
  character(len=2) :: str2

  !-- allocate prior module
  allocate( x_pr(n), x_sigma(n), x_prname(n) )

  call read_prior_para()       ! put the parameters to be optimized in a NETCDF
  !file and read them as well their
  !uncertainties,@MOUSONG.WU,2019-11

  !-- assign/mapping
  !-- Vcmax
  i1 = 1
  i2 = i1+PFT-1
  xphys(i1:i2)   = assim%p_Vcmax
  x_sigma(i1:i2) = assim%u_Vcmax
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'Vcmax_PFT'//str2
  enddo
  !-- VJ_slope
  i1 = i2+1
  i2 = i1+PFT-1
  xphys(i1:i2)   = assim%p_VJ_slope
  x_sigma(i1:i2) = assim%u_VJ_slope
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'VJ_slope_PFT'//str2
  enddo
  !-- Q10
  !   (parameter to adjust 'q10' in calculation plant respiration)
  i1 = i2+1
  i2 = i1+PFT-1
  xphys(i1:i2)   = assim%p_q10
  x_sigma(i1:i2) = assim%u_q10
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'Q10_PFT'//str2
  enddo
  !-- sif_alpha
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  xphys(i1:i2)   = assim%p_sif_alpha
  x_sigma(i1:i2) = assim%u_sif_alpha
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'sif_alpha_PFT'//str2
  enddo
  !-- sif_beta
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  xphys(i1:i2)   = assim%p_sif_beta
  x_sigma(i1:i2) = assim%u_sif_beta
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'sif_beta_PFT'//str2
  enddo
  !-- D0
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  xphys(i1:i2)   = assim%p_D0
  x_sigma(i1:i2) = assim%u_D0
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'D0_PFT'//str2
  enddo
  !-- TAWEFF
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  xphys(i1:i2)   = assim%p_taweff
  x_sigma(i1:i2) = assim%u_taweff
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'taweff_PFT'//str2
  enddo
  !-- Ksat
  !   (parameter to adjust 'Ksat', the saturated hydraulic conductivity, in beps_soilMod.F90)
  i1 = i2+1
  i2 = i1+texture-1
  xphys(i1:i2)   = assim%p_Ksat_scalar
  x_sigma(i1:i2) = assim%u_Ksat_scalar
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'Ksat_Txt'//str2
  enddo
  !-- u
  !   (parameter to adjust 'b', the parameters for calculating soil water characteristic, in beps_soilMode.F90)
  i1 = i2+1
  i2 = i1+texture-1
  xphys(i1:i2)   = assim%p_b_scalar
  x_sigma(i1:i2) = assim%u_b_scalar
  do ii=i1,i2
     write(str2, '(i2.2)') ii-i1+1
     x_prname(ii) = 'b_Txt'//str2
  enddo
  !-- f_leaf
  i1 = i2+1
  xphys(i1)    = assim%p_f_leaf
  x_sigma(i1)  = assim%u_f_leaf
  x_prname(i1) = 'f_leaf'
  !-- kc25
  i1 = i1+1
  xphys(i1)    = assim%p_kc25
  x_sigma(i1)  = assim%u_kc25
  x_prname(i1) = 'kc25'
  !-- ko25
  i1 = i1+1
  xphys(i1)    = assim%p_ko25
  x_sigma(i1)  = assim%u_ko25
  x_prname(i1) = 'ko25'
  !-- tau25
  i1 = i1+1
  xphys(i1)    = assim%p_tau25
  x_sigma(i1)  = assim%u_tau25
  x_prname(i1) = 'tau25'
  !-- agb2vod
  i1 = i1+1
  xphys(i1)    = assim%p_agb2vod
  x_sigma(i1)  = assim%u_agb2vod
  x_prname(i1) = 'agb2vod'

  !-- save normalised prior module
  x_pr = xphys/x_sigma

  !-- assign to output arguments
  x = x_pr
  sx = x_sigma

  !-- write list of parameter names (to ease consistency of post-processing)
  open(18, file='ctlvec_scaled.txt', form='formatted')
  write(str2, '(i2.2)') PFT
  write(18, '(a)') '#npft='//str2
  write(str2, '(i2.2)') texture
  write(18, '(a)') '#ntexture='//str2
  write(18, '(a,a31,a20,a20)') '#','name','scaled prior','prior uncertainty'
  do ii=1,n
     write(18,'(a32,e20.6,e20.6)') trim(x_prname(ii)), x_pr(ii), x_sigma(ii)
  enddo
  close(18)
end subroutine initx


character(len=32) function pname(i)
  use mo_prior, only:x_prname 
  implicit none
  integer, intent(in) :: i

  pname = x_prname(i)
end function pname


!***********************************************************
!     finishf
!
!> @brief Postprocessing method after having run the model (evalf),
!>        finalisation and potentially perform diagnostic or I/O,
!         e.g. could write controlvector and/or simulation to file...
!
!> \authors Michael Vossbbeck, The Inversion Lab
!> \date    January 2020
subroutine finishf(n, x, m, y)
  use esmf, only:ESMF_Finalize
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in) :: x(n), y(m)
  ! local declarations
  character(len=*), parameter :: sub = 'finishf'

  !! For clean up
  call ESMF_Finalize()

  write(*, '(a)') ' INFO::this is routine ***'//sub//'***'//&
       ' based on user demands add diagnostic/postprocessing after model evaluation here'
end subroutine finishf


!***********************************************************
!     finishc
!
!> @brief Postprocessing method after having run the cost function/optimisation,
!>        finalisation and potentially perform diagnostic or I/O,
!         e.g. could write controlvector and/or simulation to file...
!
!> \authors Michael Vossbbeck, Thomas Kaminski, The Inversion Lab
!> \date    April 2021
subroutine finishc(n,x,m,f,ok)
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in) :: x(n), f
  logical, intent(in) :: ok
  ! local declarations
  character(len=*), parameter :: sub = 'finishc'

  write(*, '(a)') ' INFO::this is routine ***'//sub//'***'//&
       ' based on user demands add diagnostic/postprocessing after '//&
       'cost function or optimisation here'

end subroutine finishc


!***********************************************************
!     get_nmeteo
!
!> @brief procedure to determine for a given point in time
!         (yr,mn,dy,tod) the corresponding index in the
!         meteorological forcing file.
!         The implementation to determine the index applies the
!         approach made by MSWU in driver.f90 (bepspkg_2020-09-21_essi)
!
!         ATTENTION::in the actual implementation the reference time
!                    is hard-coded to '2010-01-01' which equals the
!                    reference time of the meteorological forcing data
!                    (inputdata/beps_site/Site_meteo_2010_2015_hourly.nc)
!
!> \authors Michael Vossbbeck, The Inversion Lab
!> \date    April 2021
subroutine get_nmeteo(yr, mn, dy, tod, n_meteo)
  use shr_kind_mod, only: r8 =>shr_kind_r8
  use beps_time_manager, only:timemgr_diff_secs
  implicit none
  !-- arguments
  integer, intent(in) :: yr, mn, dy, tod
  integer, intent(out) :: n_meteo
  !-- local
  character(len=*), parameter :: meteo_refdate = '2013-01-01' !-- to be updated
  real(r8) :: secs_meteo
  integer :: metyr, metmn, metdy
  integer :: iostat

  !-- get year/month/day from reference date
  call str2int(meteo_refdate(1:4),  metyr, iostat)
  call str2int(meteo_refdate(6:7),  metmn, iostat)
  call str2int(meteo_refdate(9:10), metdy, iostat)

  call timemgr_diff_secs(metyr*10000+metmn*100+metdy,0,yr*10000+mn*100+dy,tod,&
       secs_meteo)
  n_meteo = int(secs_meteo/3600 + 1)
contains 

  elemental subroutine str2int(str,int,iostat)
    implicit none
    !-- arguments
    character(len=*), intent(in) :: str
    integer, intent(out)         :: int
    integer, intent(out)         :: iostat

    read(str,*,iostat=iostat)  int
  end subroutine str2int

end subroutine get_nmeteo



!***********************************************************
!     wrt_vec1d_bin
!
!> @brief Simple I/O wrapper for writing a 1D vector to binary file
!
!> \authors Michael Vossbbeck, The Inversion Lab
!> \date    April 2020
subroutine wrt_vec1d_bin(fname, n, v)
  implicit none
  !-- arguments
  character(len=*), intent(in) :: fname
  integer, intent(in) :: n
  real(kind=8), intent(in) :: v(n)
  !-- local declarations
  character(len=*), parameter :: sub = 'wrt_vec1d_bin'
  integer :: iostat

  !-- writing 1-D vector to binary file
  write(*, '(a)') ' INFO::'//sub//':start writing file ***'//trim(fname)//'***...'
  open(unit=1, file=fname, action='write', form='unformatted', iostat=iostat)
  call wrt_iocheck(iostat, 'opening ***'//trim(fname)//'*** did not succeed')
  write(1, iostat=iostat) n
  call wrt_iocheck(iostat, 'writing vector size did not succeed')
  write(1, iostat=iostat) v
  call wrt_iocheck(iostat, 'writing vector did not succeed')
  close(1, iostat=iostat)
  call wrt_iocheck(iostat, 'error message on closing ***'//trim(fname)//'***')
  write(*, '(a)') ' INFO::'//sub//'...writing DONE.'

contains
  subroutine wrt_iocheck(iostat, msg)
    implicit none
    integer, intent(in) :: iostat
    character(len=*), intent(in) :: msg
    if( iostat.ne.0 ) then
       write(*, '(a)') ' ERROR::'//sub//':'//trim(msg)
       stop
    endif
  end subroutine wrt_iocheck
end subroutine wrt_vec1d_bin


!***********************************************************
!     rd_vec1d_bin
!
!> @brief Simple I/O wrapper for reading a 1D vector from binary file
!
!> \authors Michael Vossbbeck, The Inversion Lab
!> \date    April 2020
subroutine rd_vec1d_bin(fname, n, v)
  implicit none
  !-- arguments
  character(len=*), intent(in) :: fname
  integer, intent(in) :: n
  real(kind=8), intent(out) :: v(n)
  !-- local declarations
  character(len=*), parameter :: sub = 'rd_vec1d_bin'
  integer :: iostat
  integer :: nn

  !-- writing 1-D vector to binary file
  write(*, '(a)') 'INFO::'//sub//':start reading file ***'//trim(fname)//'***...'
  open(unit=1, file=fname, action='read', form='unformatted', iostat=iostat)
  call rd_iocheck(iostat, 'opening ***'//trim(fname)//'*** did not succeed')
  read(1, iostat=iostat) nn
  call rd_iocheck(iostat, 'reading simulation vector size did not succeed')
  if( nn.ne.n ) then
     write(*, '(a,1x,2(a,i5,1x))') ' ERROR::detected inconsistent vector sizes:',&
          'expected=',n,'got=',nn
     stop
  endif
  read(1, iostat=iostat) v
  call rd_iocheck(iostat, 'reading simulation vector did not succeed')
  close(1, iostat=iostat)
  call rd_iocheck(iostat, 'error message on closing ***'//trim(fname)//'***')
  write(*, '(a)') 'INFO::'//sub//'...writing DONE.'

contains
  subroutine rd_iocheck(iostat, msg)
    implicit none
    integer, intent(in) :: iostat
    character(len=*), intent(in) :: msg
    if( iostat.ne.0 ) then
       write(*, '(a)') ' ERROR::'//sub//':'//trim(msg)
       stop
    endif
  end subroutine rd_iocheck
end subroutine rd_vec1d_bin
