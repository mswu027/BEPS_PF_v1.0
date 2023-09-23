!       module to compute model-data misfit
!       ilab march 2021
subroutine misfit(n,x,m,obsdiff)
  implicit none
  ! arguments
  integer ( kind = 4 ) n, m, j
  real ( kind = 8 ) :: x(n), obsdiff(m)
  ! local
  real ( kind = 8 ) :: y(m), yobs(m), syobs(m)
  ! read obs
  call getobs(m,yobs,syobs)
  ! simulate obs
  call evalf(n,x,m,y)
  ! difference
  obsdiff = (y-yobs)/syobs
end subroutine misfit
