!       simple driver to run cost function
!       ilab september 2023

program runcost
  implicit none
  character(len=*), parameter :: prog = 'runcost'
  integer(kind=4) :: n, m
  integer :: i
  real(kind=8), allocatable :: x0(:), x(:), sx(:)
  real(kind=8) :: fcost
  real(kind=8) :: pert = 0._8
  character(len=32), external :: pname
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
        case("--pert")
           print*, 'NOT YET IMPLEMENTED'
           stop
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
     x = x0 !-- evaluation at prior
  endif
  ! simulated pseudo observations
  call cost(n,x,m,fcost)

  write(*, '(a,e25.16)') ' INFO::'//prog//': fcost = ', fcost
contains

  subroutine dump_options()
    implicit none
    character(len=15) :: option
    write(*,'(a)') '=============================='
    write(*,*) " Options:"
    write(*,'(/)')
    option = "--pert"
    write(*,'(2x, a15,2x,a,e10.4,a)') option,&
         &"apply relative perturbation of control vector "&
         &//"(default:",pert,")"
    write(*,'(2/)')
  end subroutine dump_options
end program runcost
