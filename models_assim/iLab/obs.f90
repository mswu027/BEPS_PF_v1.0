!       module to handle observations in assimilation system
!       ilab march 2021
module obs
  !-- maximal length of simulation vector
  !   (must be large enough to capture 3 output variables for #sites * #hours)
  integer(kind=4), parameter :: nday_max = 730  !-- around 2 years simulation period
  integer(kind=4), parameter :: nlp_max  = 2    !-- maximal number of land points (i.e. site locs)
  integer(kind=4), parameter :: nsimvar  = 3    !-- SIF,Thetam,COS_flux
  integer ( kind = 4 ), parameter :: mmax = nday_max*24*nlp_max*nsimvar
  real ( kind = 8 ) :: yobs(mmax), syobs(mmax)
  logical, parameter :: debug = .false.
end module obs
 
subroutine putobs(m,yobs_in,syobs_in)
  use obs
  implicit none
  ! arguments
  integer ( kind = 4 ) :: m
  real ( kind = 8 ) :: yobs_in(m), syobs_in(m)
  if(m.gt.mmax) then 
     print*, 'putobs: m.gt.mmax = ', mmax
     stop
  endif
  yobs(1:m) = yobs_in
  syobs(1:m) = syobs_in
end subroutine putobs

subroutine getobs(m,yobs_out,syobs_out)
  use obs
  implicit none
  ! arguments
  integer ( kind = 4 ) :: m
  real ( kind = 8 ) :: yobs_out(m), syobs_out(m)
  yobs_out = yobs(1:m)
  syobs_out = syobs(1:m)

  if (debug) then
     print*, 'getobs m = ', m
     print*, 'yobs = ', yobs(1:m)
     print*, 'syobs = ', syobs(1:m)
  endif
end subroutine getobs

subroutine readobs (m,min,filename)
  use obs
  implicit none
  ! arguments
  integer ( kind = 4 ), intent(in) :: m
  character(*),         intent(in)  :: filename
  integer ( kind = 4 ), intent(out) :: min
  ! local
  logical :: exist
  inquire(FILE=filename,EXIST=exist) 
  if (exist) then
     open (unit=1, file=filename, form='unformatted')
     read (1) min
     read (1) yobs(1:min)
     read (1) syobs(1:min)
     close (1)
     print*, 'read ', min,' pseudo observations from file obs.b '
  else
     min = m
     yobs = 0.
     syobs = 1.
  endif
  if (debug) then
     print*, 'read m = ', min
     print*, 'yobs = ', yobs(1:min)
     print*, 'syobs = ', syobs(1:min)
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

