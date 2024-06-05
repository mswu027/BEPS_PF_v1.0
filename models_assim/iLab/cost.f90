!       module to compute cost function
!       ilab march 2021
subroutine cost ( n, x, m, f)
  implicit none
  ! arguments
  integer(kind=4), intent(in) :: n, m
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: f
  ! local
  real(kind=8) :: obsdiff(m), priordiff(n)
  
  call misfit(n,x,m,obsdiff)
  call devprior(n,x,priordiff)
  f = 0.5 * (sum(obsdiff**2) + sum(priordiff**2))
end subroutine cost
