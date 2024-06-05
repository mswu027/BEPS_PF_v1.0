!       driver calling nr optimisation routine
!       with reverse mode derivative
!       ilab 2023
module sizes
  real(kind=8), parameter :: missing_value = -99999._8
  integer nn, mm
  integer(kind=4)      :: ifunc=0, igrad=0
  logical              :: write_iter = .false. !-- whether to write control/gradient per iteration
  character(7+8)       :: fileout = 'output/x-0000.b'
  character(7+8)       :: gradout = 'output/g-0000.b'
  !                                  123456789
  integer, parameter :: ichar = 10 !-- position of first digit in 4 digit counter
  real ( kind = 8 ), allocatable :: x0(:), sx(:)
end module sizes
program runopti
  use sizes
  implicit none
  character(len=*), parameter :: prog = 'runopti'
  character(len=*), parameter :: trace_fname = 'rundfprv-trace.asc'
  integer(kind=4)             :: iter, n, m, itmax, iostat
  real (kind=8)               :: f, gtol, pert
  real (kind=8), allocatable  :: x(:) !-- actual control vector
  logical :: prior_enabled = .true.  !-- whether prior term is active in cost function
  logical :: ok
  logical :: exist
  logical :: debug
  external func, dfunc

  ! handling cmdline
  integer::narg,cptArg !#of arg & counter of arg
  integer :: nconsumed
  character(len=32) :: argname !Arg name
  character(len=32) :: argval

  !-- default settings
  !-- perturbation of prior control vector for initial guess
  pert = 0._8

  !Check if any arguments are found
  narg=command_argument_count()
  !Loop over the arguments
  if(narg>0)then
     !loop across options
     nconsumed=0
     arg: do cptArg=1,narg
        if( nconsumed.gt.0 ) then
           nconsumed = nconsumed - 1
           cycle arg
        endif
        call get_command_argument(cptArg,argname)
        select case(adjustl(argname))
        case("--help","-h")
           write(*,'(/)')
           write(*,*)" This is program '"//prog//"'"
           call dump_options()
           stop 0
           exit arg
        case('--pert')
           call get_command_argument(cptArg+1, argval)
           nconsumed = 1
           read(argval, *) pert
        case('--disable_prior')
           prior_enabled = .false.
        case('--write_iter')
           write_iter = .true.
        case default
           write(*,*)"Option ",adjustl(argname),"unknown"
           write(*,'(a)') "Available options:"
           call dump_options()
           stop 0
        end select
     end do arg
  end if

  ! opti iteration files should go to directory 'output',
  ! ensure directory exists (command suitable for Linux platform, only)
  if( write_iter ) then
     call system('mkdir -p '//'output')
  endif

  ! disable NetCDF output and generation of restart files
  call disable_netcdf_output()
  call disable_write_restart()
  
  ! init configuration
  call initf_fwd(n,m)
  nn = n
  mm = m
  ! init unknowns
  allocate (x0(n), x(n), sx(n))

  !-- read prior
  call initx_cd(n,x0,sx)
  !-- initial control vector
  inquire(file='x.b', exist=exist)
  if( exist ) then
     write(*, '(a)') ' INFO::'//prog//': initial control vector is read from file ***'//&
          'x.b'//'***'
     open(unit=1, file='x.b', form='unformatted')
     read(1) x
     close(1)
  else
     write(*, '(a,e10.4,a)') ' INFO::'//prog//': apply perturbation [frac of x0] of ', pert, &
          ' to prior control vector'
     x = x0 * (1+pert)
  endif

  if( prior_enabled ) then
     write(*, '(a)') ' INFO::'//prog//': prior-term active in cost function.'
     call enable_prior()
  else
     write(*, '(a)') ' INFO::'//prog//': prior-term de-activated in cost function.'
     call disable_prior()
  endif
  
  ! intialise optimisation parameters
  iter = 0
  gtol = 1.e-8
  f = 0.
  itmax=200
  open(unit=3, file=trace_fname, form='formatted', action='write', iostat=iostat)
  if( iostat.ne.0 ) then
     write(*, '(a)') ' FATAL::'//prog//': trace file ***'//trace_fname//'*** could not be '//&
          'opened for writing!'
     stop
  endif
  ! output table
  print '(2a10,5a20)', 'g_calls', 'f_calls', 'f', 'norm_g', 'x'
  call dfpmin(x,n,gtol,iter,f,func,dfunc)
  print*, 'dfpmin, iterations: ', iter
  call woptimum(n,x0,sx,x)
  !-- close trace file
  close(3)
  write(*, '(a)') ' INFO::'//prog//': generated trace file ***'//trace_fname//'***'
  ok = (iter.lt.itmax)
  call finishc(n,x,f,ok)

  !-- dispose resources
  deallocate(x0, x, sx)
contains

  subroutine dump_options()
    implicit none
    character(len=15) :: option
    write(*,'(a)') '=============================='
    write(*,*) " Options:"
    write(*,'(/)')
    option = "--pert"
    write(*,'(2x, a15,2x,a,e10.4,a)') option,&
         &"perturbation of prior to set initial control vector "&
         &//"(default:",pert,")"
    option = "--disable_prior"
    write(*,'(2x, a15, 2x, a)') option,&
         &"whether to deactivate prior-term in cost function."
    option = "--write_iter"
    write(*,'(2x, a15, 2x, a)') option,&
         &"whether to write control vector and gradient to binary file per evaluation."
    write(*,'(2/)')
  end subroutine dump_options
end program runopti

function func(x)
  use sizes
  implicit none
  real ( kind = 8 ) :: x(nn), func, f
  if( write_iter ) then
     ! save control vector used for evaluation
     call writeiter(1,ifunc+1,ichar,nn,x,fileout)
  endif
  call cost_cd(nn,x,mm, f)
  func=f
  ifunc = ifunc + 1
  print '(2i10,e20.6,a20,6e20.6)', igrad, ifunc, f, '-', x(1:min(3,nn))
  write(3,*)  f, missing_value, x(1:nn)
end function func

subroutine dfunc(x,fd)
  use sizes
  implicit none
  real ( kind = 8 ) :: x(nn), fd(nn), f, fb 
  fd = 0.
  fb = 1.
  call COST_bw(nn, x, fd, mm, f, fb)
  !-- increment counter, save gradient
  igrad=igrad+1
  print '(2i10,a20,6e20.6)', igrad, ifunc, '-', sqrt(sum(fd**2)), x(1:min(3,nn))  ! reverse model pr  write(3,*)  missing_value, sqrt(sum(fd**2)), x(1:min(100,nn))
  if( write_iter ) then
     call writeiter(1,igrad,ichar,nn,fd,gradout)
  endif
end subroutine dfunc
