


subroutine wrf_error_fatal(msg)
!  use shr_sys_mod, only: shr_sys_abort
  implicit none
  character(len=*), intent(in) :: msg
  !iLab-added:function with optional argument requires an explicit interface!
  interface
     subroutine endrun(msg)
       implicit none
       character(len=*), intent(in), optional :: msg
     end subroutine endrun
  end interface
  write(6,*) 'wrf_error_fatal: ',trim(msg)
!  call shr_sys_abort( msg )
  call endrun(msg)
end subroutine wrf_error_fatal

