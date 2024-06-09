!       simple function
!       iLab

subroutine devprior (n, x, priordiff)
  implicit none
  ! arguments
  integer ( kind = 4 ) :: n
  real ( kind = 8 ) :: x(n), priordiff(n)
  ! local
  real ( kind = 8 ) :: sx(n), x0(n) 
  logical :: mask(n)
  ! get prior
  call getprior(n, x0, sx, mask) 
  where(mask)
     priordiff = (x-x0) ! x and x0 already normalised by sx
  elsewhere
     priordiff = 0.
  endwhere
end subroutine devprior
