!       module to handle observations in assimilation system
!       ilab march 2021
!       last: 2022/06
module obs
  use mo_bepsfunc_ctl, only: nsimvar
  real(kind=8),    parameter :: fill_value = -99999._8
  integer(kind=4), parameter :: nday_max = 730  !-- around 2 years simulation period
  integer(kind=4), parameter :: nlp_max  = 500  !-- maximal number of land points (i.e. site locs)
  integer(kind=4), parameter :: mmax = nday_max*24*nlp_max*nsimvar
  real(kind=8), allocatable :: yobs(:), syobs(:)
  !-- activate to trace potential problems
  logical, parameter :: debug = .false.
end module obs
 

subroutine getobs(m,yobs_out,syobs_out,obs_missing_value)
  use obs
  implicit none
  ! arguments
  integer(kind=4), intent(in) :: m
  real(kind=8), intent(out)   :: yobs_out(m), syobs_out(m)
  real(kind=8), intent(out)   :: obs_missing_value

  yobs_out = yobs(1:m)
  syobs_out = syobs(1:m)
  obs_missing_value = fill_value

  if (debug) then
     print*, 'getobs m = ', m
     print*, 'yobs = ', yobs(1:m)
     print*, 'syobs = ', syobs(1:m)
     print*, 'missing_value=', obs_missing_value
  endif
end subroutine getobs


!/////////////////////////////////////////////////////////////////////
!
! NOTE:
! - routines 'readobs' and 'writeobs' below implemented a Fortran binary
!   interface for I/O of (pseudo-) observations.
!   This I/O has been converted to NetCDF format in the current
!   implementation, the routines below are now ! D E P R E C E A T E D !
!
subroutine readobs (m,mread,filename)
  use obs
  implicit none
  ! arguments
  integer ( kind = 4 ), intent(in) :: m
  character(*),         intent(in)  :: filename
  integer ( kind = 4 ), intent(out) :: mread
  ! local
  character(len=*), parameter :: sub = 'readobs'
  logical :: exist

  if( allocated(yobs) .or. allocated(syobs) ) then
     write(*, '(a)') ' FATAL::'//sub//': yobs/syobs are not expected to be already allocated here!'
     stop
  else
     !-- 'm' is expected size of simulation/obs vector
     allocate( yobs(m), syobs(m) )
  endif
  inquire(FILE=filename,EXIST=exist) 
  if (exist) then
     open (unit=1, file=filename, form='unformatted')
     read (1) mread
     read (1) yobs(1:mread)
     read (1) syobs(1:mread)
     close (1)
     print*, 'read ', mread,' pseudo observations from file obs.b '
  else
     mread = m
     yobs  = 0._8
     syobs = 1._8
  endif
  if (debug) then
     print*, 'read mread = ', mread
     print*, 'yobs = ', yobs(1:mread)
     print*, 'syobs = ', syobs(1:mread)
     print*, 'done readobs'
  endif
end subroutine readobs

subroutine writeobs (m, yobs, syobs)
  implicit none
  ! arguments
  integer ( kind = 4 ) :: m
  real ( kind = 8 ) :: yobs(m), syobs(m)
  ! local
  syobs = 1. ! to be refined
  open (unit=1, file='obs.b', form='unformatted')
  write (1) m
  write (1) yobs
  write (1) syobs
  close (1)
  print*, 'wrote ', m,' pseudo observations to file obs.b '
end subroutine writeobs


