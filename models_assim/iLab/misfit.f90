!       module to compute model-data misfit
!       ilab march 2021
!       last: 06/2022
subroutine misfit(n,x,m,obsdiff)
  implicit none
  ! arguments
  integer(kind=4), intent(in) :: n, m
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: obsdiff(m)
  ! local
  real(kind=8) :: y(m), yobs(m), syobs(m), obs_missing_value

  ! read obs
  call getobs(m,yobs,syobs,obs_missing_value)
  ! simulate obs
  call evalf(n,x,m,y)
  ! difference
  where(yobs.ne.obs_missing_value)
     obsdiff = (y-yobs)/syobs
  elsewhere
     obsdiff = 0.
  endwhere
end subroutine misfit
