subroutine endrun(msg)
!Abort the model for abnormal termination
implicit none

character(len=*), intent(in), optional :: msg    ! string to be printed
!include 'mpif.h'

if (present (msg)) then
   write(6,*)'ENDRUN:', msg
else
   write(6,*)'ENDRUN IS CALLED' 
end if

!call mpi_abort (MPI_COMM_WORLD, 1)
return

end subroutine 


