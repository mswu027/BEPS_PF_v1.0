module obs_netcdf
  use netcdf
  
contains
  subroutine nc_check(status)
    implicit none
    integer, intent (in) :: status
    character(len=*), parameter :: sub = 'nc_check'

    if(status .ne. nf90_noerr) then
       call endrun(' ERROR::'//sub//':: !!!while writing NetCDF output!!!'//'('//trim(nf90_strerror(status))//')')
    endif
  end subroutine nc_check

  subroutine update_missing_value(ncid, varid, data_buffer)
    use obs, only: fill_value
    implicit none
    !- arguments
    integer, intent(in) :: ncid, varid
    real(kind=8), intent(inout) :: data_buffer(:,:)
    !- local declarations
    real(kind=8) :: read_miss_value
    integer :: rc
    rc = nf90_get_att(ncid, varid, 'missing_value', read_miss_value)
    if ( rc.eq.nf90_noerr ) then
       !nothing todo, no missing value specified for current NetCDF dataset
    else
       !-- map missing value from NetCDF file to fill value used in assimilation system
       where(data_buffer.eq.read_miss_value)
          data_buffer = fill_value
       end where
    endif
  end subroutine update_missing_value

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

end module obs_netcdf


!***********************************************************
!     ncreadobs
!
!> @brief reading observed SIF, surface soil moisture, and COSflux
!         and their uncertainties from NetCDF file for assimilation purposes
!
!> @comment data variables read are passed to 1D observation buffers
!           yobs and syobs in module 'obs'
!
!> \authors Michael Vossbeck, The Inversion Lab
!> \date    June 2022
subroutine ncreadobs(m, filename)
  use obs_netcdf
  use obs, only: fill_value, yobs, syobs, mmax
  use mo_bepsfunc_ctl
  use beps_par, only: nlp
  implicit none
  !-- arguments
  integer, intent(in) :: m !-- expected length of simulation vector
  character(len=*), intent(in) :: filename
  !-- local declarations
  character(len=*), parameter :: sub = 'ncreadobs'
  integer, parameter :: ntc = 4  !-- year,month,day,hour
  integer, parameter :: zlev = 4 !-- gzip-like compression level
  integer, parameter :: chunklim = 2*1000000 !-- 2MB [bytes]
  real(kind=8), allocatable :: sif(:,:), thetam(:,:), cosflux(:,:)
  real(kind=8), allocatable :: sif_unc(:,:), thetam_unc(:,:), cosflux_unc(:,:)
  integer(kind=4), allocatable :: yyyymmddss(:,:)
  !-- NetCDF stuff
  integer :: ncid, rc
  integer :: dimid_nlp, dimid_time, dimid_ntc
  integer :: nlp_read, ntp_read, ntc_read
  integer :: varid
  integer :: read_start(2), read_count(2)
  !--
  logical :: exist
  integer :: cnt(1)
  integer :: itp, it_s, it_e
  character(len=32) :: fmt
  !-- activate in order to trace problems
  logical :: ldebug = .false.

  if( allocated(yobs) .or. allocated(syobs) ) then
     write(*, '(a)') ' FATAL::'//sub//': yobs/syobs are not expected to be already allocated here!'
     stop
  else
     !-- 'm' is expected size of simulation/obs vector
     allocate( yobs(m), syobs(m) )
  endif

  inquire(FILE=filename, EXIST=exist)
  if( .not.exist ) then
     yobs  = 0._8
     syobs = 1._8
     write(*, '(a)') ' WARNING::'//sub//': observation file ***'//trim(filename)//&
          '*** NOT found, running with adhoc values yobs=0 syobs=1'
     return
  endif
 
  !-- open file
  rc = nf90_open(filename, nf90_nowrite, ncid)
  call nc_check(rc)
  if(ldebug) then
     write(*, '(a)') ' DEBUG::'//sub//': opened for writing ***'//trim(filename)//'***'
  endif

  !-- read dimension sizes
  !>> #land points
  call nc_check(nf90_inq_dimid(ncid, "nlp", dimid_nlp))
  call nc_check(nf90_inquire_dimension(ncid, dimid_nlp, len=nlp_read))
  if( nlp.ne.nlp_read ) then
     write(*,'(a,2(a,i6,1x))') ' FATAL::'//sub//': unexpected number of land points, ', &
          'nlp,expected=',nlp,'nlp, in file=',nlp_read
     stop
  endif
  !>> #time points
  call nc_check(nf90_inq_dimid(ncid, "time", dimid_time))
  call nc_check(nf90_inquire_dimension(ncid, dimid_time, len=ntp_read))
  if( ntp.gt.ntp_read ) then
     write(*,'(a,2(a,i6,1x))') ' FATAL::'//sub//': too few time points in observational file, ', &
          'ntp,expected=',ntp,'ntp, in file=',ntp_read
     stop
  endif
  !>> # calender componnts
  call nc_check(nf90_inq_dimid(ncid, "ntc", dimid_ntc))
  call nc_check(nf90_inquire_dimension(ncid, dimid_ntc, len=ntc_read))
  if( ntc.ne.ntc_read ) then
     write(*,'(a,2(a,i6,1x))') ' FATAL::'//sub//': unexpected number of calendar components, ', &
          'ntc,expected=',ntc,'ntc, in file=',ntc_read
     stop
  endif

  !-- ensure suitable size of simulation vector
  if( m.ne.nsimvar*ntp*nlp ) then
     write(*,'(a,2(a,i6,1x))') ' FATAL::'//sub//': unexpected size of observation vector, ', &
          'expected=',nsimvar*ntp*nlp, 'got=',m
     stop
  else if( m.gt.mmax ) then
     write(*,'(a,2(a,i6,1x))') ' FATAL::'//sub//&
          ': expected size of observation vector exceeds buffer, ', &
          'expected=',m, 'mmax=',mmax
     stop
  endif

  !-- initialise observation buffer
  yobs  = fill_value
  syobs = fill_value

  !-- allocate data variables
  allocate(sif(nlp,ntp), thetam(nlp,ntp), cosflux(nlp,ntp))
  allocate(sif_unc(nlp,ntp), thetam_unc(nlp,ntp), cosflux_unc(nlp,ntp))
  allocate(yyyymmddss(ntc,ntp_read))

  !-- read data variables
  !>> calendar
  call nc_check(nf90_inq_varid(ncid, "yyyymmddss", varid))
  call nc_check(nf90_get_var(ncid, varid, yyyymmddss))
  !-- check temporal coverage of observations fully covers simulation period,
  !   stop processing in case of mismatch
  it_s = -1
  it_e = -1
  obsloop:do itp=1,ntp_read
     !-- compare against first simulation time
     if( maxval(abs(yyyymmddss(1:ntc,itp)-time_points(1:ntc,1))).eq.0 ) then
        it_s = itp
     endif
     !-- comapre against last simulation time
     if( maxval(abs(yyyymmddss(1:ntc,itp)-time_points(1:ntc,ntp))).eq.0 ) then
        it_e = itp
     endif
     if( it_s.ne.-1 .and. it_e.ne.-1 ) then
        exit obsloop
     end if
  enddo obsloop
  if( it_s.eq.-1 .or. it_e.eq.-1 ) then
     write(*,'(a,2(a,i6,1x))') ' FATAL::'//sub//&
          ': observations do not cover range of simulation! ',&
          'it_start=',it_s, 'it_end=',it_e
     stop
  else if( it_e-it_s+1.ne.ntp ) then
     write(*,'(a,3(a,i6,1x))') ' FATAL::'//sub//&
          ': observations do not cover range of simulation! ',&
          'it_start=',it_s, 'it_end=',it_e, 'ntp=',ntp
     stop
  endif
  cnt = count(reshape(yyyymmddss(:,it_s:it_e),(/ntp*ntc/)).ne.reshape(time_points(1:ntc,1:ntp),(/ntp*ntc/)))
  if( cnt(1).ne.0 ) then
     write(*,'(a,1(a,i6,1x))') ' FATAL::'//sub//&
          ': expected same time-points for simulation setup and observations. ', &
          'detected #differing time points',cnt(1)
     stop
  endif
  !-- logging
  write(*, '(a,a,i4,1x,i2,1x,i2,1x,i5,1x,a,i6)') ' INFO::'//sub//': ',&
       'first simulation time point ', time_points(1:ntc,1), &
       'discovered in observational data at time-index it_start=', it_s
  write(*, '(a,a,i4,1x,i2,1x,i2,1x,i5,1x,a,i6)') ' INFO::'//sub//': ',&
       'last simulation time point ', time_points(1:ntc,ntp), &
       'discovered in observational data at time-index it_end=', it_e

  !-- range for reading data from file
  !   *full* spatial range 1:nlp
  !   temporal range of simulation it_s:it_e = it_s:it_s+ntp-1
  read_start = (/1, it_s/)
  read_count = (/nlp,ntp/)

  !>> SIF
  call nc_check(nf90_inq_varid(ncid, "sif", varid))
  call nc_check(nf90_get_var(ncid, varid, sif, start=read_start, count=read_count))
  call update_missing_value(ncid, varid, sif)
  !>> SIF uncertainty
  call nc_check(nf90_inq_varid(ncid, "sif_unc", varid))
  call nc_check(nf90_get_var(ncid, varid, sif_unc, start=read_start, count=read_count))
  call update_missing_value(ncid, varid, sif_unc)
  !>> Thetam
  call nc_check(nf90_inq_varid(ncid, "sm", varid))
  call nc_check(nf90_get_var(ncid, varid, thetam, start=read_start, count=read_count))
  call update_missing_value(ncid, varid, thetam)
  !>> Thetam uncertainty
  call nc_check(nf90_inq_varid(ncid, "sm_unc", varid))
  call nc_check(nf90_get_var(ncid, varid, thetam_unc, start=read_start, count=read_count))
  call update_missing_value(ncid, varid, thetam_unc)
  !>> COSflux
  call nc_check(nf90_inq_varid(ncid, "COS_flux", varid))
  call nc_check(nf90_get_var(ncid, varid, cosflux, start=read_start, count=read_count))
  call update_missing_value(ncid, varid, cosflux)
  !>> COSflux uncertainty
  call nc_check(nf90_inq_varid(ncid, "COS_flux_unc", varid))
  call nc_check(nf90_get_var(ncid, varid, cosflux_unc, start=read_start, count=read_count))
  call update_missing_value(ncid, varid, cosflux_unc)

  !-- map to 1D observation vector(s)
  !-- 1D ordering slowest-to-fastest varying: time->location->variable,
  !   with sif     -> first position
  !        thetam  -> second position
  !        cosflux -> third position
  yobs(1:m:nsimvar)  = reshape(sif,         (/nlp*ntp/))
  yobs(2:m:nsimvar)  = reshape(thetam,      (/nlp*ntp/))
  yobs(3:m:nsimvar)  = reshape(cosflux,     (/nlp*ntp/))
  syobs(1:m:nsimvar) = reshape(sif_unc,     (/nlp*ntp/))
  syobs(2:m:nsimvar) = reshape(thetam_unc,  (/nlp*ntp/))
  syobs(3:m:nsimvar) = reshape(cosflux_unc, (/nlp*ntp/))

  !-- properly close file handle
  call nc_check(nf90_close(ncid))

  cnt = count(yobs(1:m).eq.fill_value)
  !--
  fmt = '(a,a,'//trim(ifmt(m))//',1x,a,'//trim(ifmt(cnt(1)))//')'
  write(*, fmt) ' INFO::'//sub//': 1D observation vector with ','size=',m,' nmissing=',cnt(1)
  !-- dispose memory
  deallocate(sif, thetam, cosflux)
  deallocate(sif_unc, thetam_unc, cosflux_unc)
  deallocate(yyyymmddss)

end subroutine ncreadobs


!***********************************************************
!     ncwriteobs
!
!> @brief writing pseudo-observations of SIF, surface soil moisture, and COSflux
!         to NetCDF file suitable as input for identical twin experiment(s)
!
!> @comment the input vectors yobs and syobs are distributed to specific variable,
!           land point and time according to the ordering
!           time->landpoints->variable (slowest-to-fastest)
!
!> \authors Michael Vossbeck, The Inversion Lab
!> \date    June 2022
subroutine ncwriteobs(m, yobs, syobs)
  use obs_netcdf
  use mo_bepsfunc_ctl
  use bepstype, only: bound
  use beps_par, only: nlp
  use obs, only: fill_value
  implicit none
  !-- arguments
  integer, intent(in) :: m
  real(kind=8), intent(in) :: yobs(m), syobs(m)
  !-- local declarations
  character(len=*), parameter :: sub = 'ncwriteobs'
  character(len=*), parameter :: outname = 'obs.nc' !-- file to be generated
  integer, parameter :: ntc = 4  !-- year,month,day,hour
  integer, parameter :: zlev = 4 !-- gzip-like compression level
  integer, parameter :: chunklim = 2*1000000 !-- 2MB [bytes]
  real(kind=8), allocatable :: data3D(:,:,:), dataunc3D(:,:,:)
  !-- NetCDF stuff
  integer :: ncid, rc
  integer :: dimid_nlp, dimid_time, dimid_ntc
  integer :: varid
  integer :: var_chunksizes(2) !-- for data variables only dimensioned with (nlp,ntp)
  integer :: chunk_ntp
  !--
  logical :: ldebug = .false.

  
  !-- basic consistency
  if( m.ne.nsimvar*ntp*nlp ) then
     write(*,'(a,2(a,i6,1x))') ' FATAL::'//sub//': unexpected size of observation vector, ', &
          'expected=',nsimvar*ntp*nlp, 'got=',m
  endif

  allocate(data3D(nsimvar,nlp,ntp), dataunc3D(nsimvar,nlp,ntp))

  !-- disassemble 1D obs vector to separate data variables depending on space and time
  !-- 1D ordering, slowest-to-fastest varying: time->location->variable
  data3D    = reshape(yobs(1:m), (/nsimvar,nlp,ntp/))
  dataunc3D = reshape(syobs(1:m), (/nsimvar,nlp,ntp/))

  !-- open file
  rc = nf90_create(outname, nf90_netcdf4, ncid)
  call nc_check(rc)
  if(ldebug) then
     write(*, '(a)') ' DEBUG::'//sub//': opened for writing ***'//trim(outname)//'***'
  endif
  !-- define dimensions
  call nc_check(nf90_def_dim(ncid, "nlp",  nlp, dimid_nlp))
  call nc_check(nf90_def_dim(ncid, "time", ntp, dimid_time))
  call nc_check(nf90_def_dim(ncid, "ntc",  ntc, dimid_ntc))   !yyyymmdd

  !-- dimension variables
  !>> calendar
  call nc_check(nf90_def_var(ncid, "yyyymmddss", nf90_int,(/dimid_ntc,dimid_time/), varid))
  call nc_check(nf90_put_att(ncid, varid, "comment", "year/month/day/time-of-day[s]"))
  !>> longitude
  call nc_check(nf90_def_var(ncid, 'lon', nf90_double, (/dimid_nlp/), varid))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "longitude"))
  call nc_check(nf90_put_att(ncid, varid, "units", "degrees_east"))
  !>> latitude
  call nc_check(nf90_def_var(ncid, 'lat', nf90_double, (/dimid_nlp/), varid))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "latitude"))
  call nc_check(nf90_put_att(ncid, varid, "units", "degrees_north"))

  !-- data variables
  if( ntp*nlp*8 .gt. chunklim ) then
     chunk_ntp = chunklim/(nlp*8) !-- (nlp*8) because all variable data expected as kind=8.
  else
     chunk_ntp = ntp
  endif
  if(ldebug) then
     write(*, '(a,2(a,i6,1x))') ' DEBUG::'//sub//': chunksizes explicitly set to ', &
          'chunk_nlp=',nlp, 'chunk_ntp=',chunk_ntp
  endif
  var_chunksizes = (/nlp,chunk_ntp/) !-- always *all* land points
  if(ldebug) then
     write(*, '(a,2(i6,1x))') ' DEBUG::'//sub//': var_chunksizes = ',&
          var_chunksizes(1), var_chunksizes(2)
  endif
  !>> SIF
  call nc_check(nf90_def_var(ncid, "sif", nf90_double, (/dimid_nlp,dimid_time/), varid, &
       chunksizes=var_chunksizes, shuffle=.true., deflate_level=zlev))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "Solar induced fluorescence"))
  call nc_check(nf90_put_att(ncid, varid, "units", "mW/m2/nm/sr"))
  call nc_check(nf90_put_att(ncid, varid, "missing_value", fill_value))
  !>> SIF uncertainty
  call nc_check(nf90_def_var(ncid, "sif_unc", nf90_double, (/dimid_nlp,dimid_time/), varid, &
       chunksizes=var_chunksizes, shuffle=.true., deflate_level=zlev))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "Uncertainty of Solar induced fluorescence"))
  call nc_check(nf90_put_att(ncid, varid, "units", "mW/m2/nm/sr"))
  call nc_check(nf90_put_att(ncid, varid, "missing_value", fill_value))
  !>> Thetam
  call nc_check(nf90_def_var(ncid, "sm", nf90_double, (/dimid_nlp,dimid_time/), varid, &
       chunksizes=var_chunksizes, shuffle=.true., deflate_level=zlev))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "fractional surface layer soil moisture"))
  call nc_check(nf90_put_att(ncid, varid, "units", ""))
  call nc_check(nf90_put_att(ncid, varid, "missing_value", fill_value))
  !>> Thetam uncertainty
  call nc_check(nf90_def_var(ncid, "sm_unc", nf90_double, (/dimid_nlp,dimid_time/), varid, &
       chunksizes=var_chunksizes, shuffle=.true., deflate_level=zlev))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "Uncertainty of fractional surface layer soil moisture"))
  call nc_check(nf90_put_att(ncid, varid, "units", ""))
  call nc_check(nf90_put_att(ncid, varid, "missing_value", fill_value))
  !>> cosflux
  call nc_check(nf90_def_var(ncid, "COS_flux", nf90_double, (/dimid_nlp,dimid_time/), varid, &
       chunksizes=var_chunksizes, shuffle=.true., deflate_level=zlev))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "COS flux for soil and plant"))
  call nc_check(nf90_put_att(ncid, varid, "units", "pmol/m2/s"))
  call nc_check(nf90_put_att(ncid, varid, "missing_value", fill_value))
  !>> cosflux uncertainty
  call nc_check(nf90_def_var(ncid, "COS_flux_unc", nf90_double, (/dimid_nlp,dimid_time/), varid, &
       chunksizes=var_chunksizes, shuffle=.true., deflate_level=zlev))
  call nc_check(nf90_put_att(ncid, varid, "long_name", "Uncertainty of COS flux for soil and plant"))
  call nc_check(nf90_put_att(ncid, varid, "units", "pmol/m2/s"))
  call nc_check(nf90_put_att(ncid, varid, "missing_value", fill_value))

  !-- terminate definition mode
  call nc_check(nf90_enddef(ncid))

  !-- start writing
  !>> calendar
  call nc_check(nf90_inq_varid(ncid, "yyyymmddss", varid))
  call nc_check(nf90_put_var(ncid, varid, time_points(1:ntc,1:ntp)))
  !>> longitude
  call nc_check(nf90_inq_varid(ncid, "lon", varid))
  call nc_check(nf90_put_var(ncid, varid, bound%longitude(:)))
  !>> latitude
  call nc_check(nf90_inq_varid(ncid, "lat", varid))
  call nc_check(nf90_put_var(ncid, varid, bound%latitude(:)))
  !>> SIF
  call nc_check(nf90_inq_varid(ncid, "sif", varid))
  call nc_check(nf90_put_var(ncid, varid, data3D(1,:,:)))
  !>> SIF uncertainty
  call nc_check(nf90_inq_varid(ncid, "sif_unc", varid))
  call nc_check(nf90_put_var(ncid, varid, dataunc3D(1,:,:)))
  !>> Thetam
  call nc_check(nf90_inq_varid(ncid, "sm", varid))
  call nc_check(nf90_put_var(ncid, varid, data3D(2,:,:)))
  !>> Thetam uncertainty
  call nc_check(nf90_inq_varid(ncid, "sm_unc", varid))
  call nc_check(nf90_put_var(ncid, varid, dataunc3D(2,:,:)))
  !>> COSflux
  call nc_check(nf90_inq_varid(ncid, "COS_flux", varid))
  call nc_check(nf90_put_var(ncid, varid, data3D(3,:,:)))
  !>> COSflux uncertainty
  call nc_check(nf90_inq_varid(ncid, "COS_flux_unc", varid))
  call nc_check(nf90_put_var(ncid, varid, dataunc3D(3,:,:)))

  !-- properly close file handle
  call nc_check(nf90_close(ncid))

  !-- dispose memory
  deallocate(data3D, dataunc3D)

end subroutine ncwriteobs
