!       simple driver to simulate observations
!       ilab march 2021

program simobs
  implicit none
  integer (kind = 4) :: n, m
  integer :: i
  real ( kind = 8 ), allocatable :: x(:), sx(:), y(:), sy(:)
  character(len=32), external :: pname
  ! set dimensions and allocate
  call initf(n,m)
  allocate (x(n),sx(n))
  allocate (y(m),sy(m))
  ! init unknowns
  call initx(n,x,sx)
  ! simulate obs
  call evalf(n,x,m,y)
  ! output result
  call writeobs(m,y,sy)
  call finishf(n,x,m,y)
end program simobs
