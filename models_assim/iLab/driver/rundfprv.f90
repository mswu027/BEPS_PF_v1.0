!       driver calling nr optimisation routine
!       with reverse mode derivative
!       ilab 2017
module sizes
  integer nn, mm
end module sizes
program runopti
  use sizes
  implicit none
  integer ( kind = 4 )           :: iter, n, m, itmax
  logical :: ok
  real ( kind = 8 )              :: f, gtol, pert
  real ( kind = 8 ), allocatable :: x(:), sx(:)
  external func, dfunc
  ! init configuration
  call initf_bw(n,m)
  nn = n
  mm = m
  ! init unknowns
  allocate (x(n),sx(n))
  call initx_cd(n,x,sx)
  ! perturb initial parameter values
  pert = 0.5
  x = x * (1+pert)
  ! intialise optimisation parameters
  iter = 0
  gtol = 1.e-8
  f = 0.
  itmax=200
  ! output table
  print '(3a20)', 'f', 'norm g', 'x'
  call dfpmin(x,n,gtol,iter,f,func,dfunc)
  print*, 'dfpmin, iterations: ', iter
  ok = (iter.lt.itmax)
  call finishc(n,x,f,ok)
end program runopti

function func(x)
  use sizes
  implicit none
  real ( kind = 8 ) :: x(nn), func, f
  call cost_cd(nn,x,mm, f)
  func=f
  print '(e20.6,a20,6e20.6)', f, '-', x(1:min(3,nn))
end function func

subroutine dfunc(x,fd)
  use sizes
  implicit none
  real ( kind = 8 ) :: x(nn), fd(nn), f, fb 
  fd = 0.
  fb = 1.
  call COST_bw(nn, x, fd, mm, f, fb)
  print '(a20,6e20.6)', '-', sqrt(sum(fd**2)), x(1:min(3,nn))  ! reverse model provides no value of f
  write(3,*)  sqrt(sum(fd**2)), x(1:min(100,nn))
end subroutine dfunc
