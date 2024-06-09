subroutine getprior (n, x, sx, mask)
  use mo_prior
  implicit none
  ! arguments
  integer ( kind = 4 ), intent(in)  :: n
  real ( kind = 8 ), intent(out)    :: x(n), sx(n)
  logical, intent(out) :: mask(n)
  ! local
  ! initialise parameters
  x = x_pr(1:n)
  sx = x_sigma(1:n)
!  mask = .true. ! activate prior in setup with real data
  mask = .false. ! only for demonstration of package in identical twin setup with pseudodata
end subroutine getprior
