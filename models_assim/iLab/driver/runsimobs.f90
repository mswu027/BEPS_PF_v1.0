!       simple driver to simulate observations
!       ilab march 2021

program simobs
  implicit none
  character(len=*), parameter :: prog = 'simobs'
  integer (kind = 4) :: n, m
  integer :: i
  real ( kind = 8 ), allocatable :: x0(:), x(:), sx(:), y(:), sy(:)
  character(len=32), external :: pname
  !-- settings
  real(kind=8), parameter :: syobs_default = 1._8
  real(kind=8), parameter :: unc_floor = 1.e-4
  real(kind=8) :: relunc = -1._8 !-- by default uniform absolute uncertainty is applied
  !-- local
  logical :: exist

  ! handling cmdline
  integer::narg,cptArg !#of arg & counter of arg
  integer :: nconsumed
  character(len=32) :: argname !Arg name
  character(len=32) :: argval

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
        case('--runc')
           call get_command_argument(cptArg+1, argval)
           nconsumed = 1
           read(argval, *) relunc
        case default
           write(*,*)"Option ",adjustl(argname),"unknown"
           write(*,'(a)') "Available options:"
           call dump_options()
           stop 0
        end select
     end do arg
  end if


  ! set dimensions and allocate
  call initf(n,m)
  allocate (x0(n),x(n),sx(n))
  allocate (y(m),sy(m))
  !-- read prior
  call initx(n,x0,sx)
  !-- initial control vector
  inquire(file='x.b', exist=exist)
  if( exist ) then
     write(*, '(a)') ' INFO::'//prog//': initial control vector is read from file ***'//&
          'x.b'//'***'
     open(unit=1, file='x.b', form='unformatted')
     read(1) x
     close(1)
  else
     x = x0 !-- simulation at prior
  endif
  ! simulated pseudo observations
  call evalf(n,x,m,y)
  ! output result
  if( relunc.lt.0._8 ) then
     write(*, '(a)') ' INFO::'//prog//': pseudo observation uncertainty uniformly set to 1'
     sy = 1._8 !TODO: to be refined
  else
     write(*, '(a,e10.4)') ' INFO::'//prog//': pseudo observation generated with relative uncertainty of ', relunc
     sy = relunc*abs(y)
  endif
  write(*, '(a,e10.4)') ' INFO::'//prog//': applied uncertainty floor value of ',unc_floor
  sy = max(sy, unc_floor)

  call ncwriteobs(m,y,sy)
  call finishf(n,x,m,y)


contains

  subroutine dump_options()
    implicit none
    character(len=15) :: option
    write(*,'(a)') '=============================='
    write(*,*) " Options:"
    write(*,'(/)')
    option = "--runc"
    write(*,'(2x, a15,2x,a,e10.4,a)') option,&
         &"apply relative uncertainty with positive value (instead of uniform absolute uncertainty of 1) "&
         &//"(default:",relunc,")"
    write(*,'(2/)')
  end subroutine dump_options
end program simobs
