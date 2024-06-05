module restart
use shr_kind_mod,only:r8=>shr_kind_r8
use beps_par
use outputMod
use bepstype
!--iLab::can avoid beps_time_manager, date information now passed as actual arguments
! use beps_time_manager, only:get_curr_date
use beps_soilMod
implicit none

integer    :: rst_icdate
integer    :: rst_icsec

contains

subroutine restart_io(flag,yy,mm,dd,tod)
implicit none
character(len=*) :: flag
!-- iLab::yy,mm,dd,tod turned to input arguments
!         (can be omitted when in  flag=='read')
integer, intent(in), optional :: yy,mm,dd,tod
!--iLab::added flag, original log message below only on activated flag
logical :: log_flag
real(r8)         :: v2last1(nlp,0:40,PFT)

integer          :: i,j,ierr
character(len=255) :: fln1,fln2
! integer          :: yy,mm,dd,tod
character(len=8) :: datestr

log_flag = .false.
if(trim(flag) == "write") then
   if( .not.present(yy) ) then
      write(*, '(a)') ' FATAL::restart_io::yy must be given in write mode!'
      stop
   else if( .not.present(mm) ) then
      write(*, '(a)') ' FATAL::restart_io::mm must be given in write mode!'
      stop
   else if( .not.present(dd) ) then
      write(*, '(a)') ' FATAL::restart_io::dd must be given in write mode!'
      stop
   else if( .not.present(tod))  then
      write(*, '(a)') ' FATAL::restart_io::tod must be given in write mode!'
      stop
   endif
!   call mpi_barrier(mpi_comm_world,ierr)
!   do i = 1,PFT
!      do j = 0,40
!         call mpi_gatherv(v2last(1,j,i),npoints,mpi_real8,v2last1(1,j,i),dp,sp,mpi_real8,0,mpi_comm_world,ierr)
!      end do
!   end do
!   call mpi_barrier(mpi_comm_world,ierr)

 v2last1 = v2last
  
   ! call get_curr_date(yy,mm,dd,tod)
   write(datestr,'(i8)') yy*10000+mm*100+dd
   !   if(myid ==0)
   !--iLab::log message only on activated flag
   if( log_flag ) write(*,*) "restart file on ",tod
   
!   if(myid ==0) then
      fln1 = trim(beps_rst_dir)//"beps.restart."//trim(datestr)

       rst_icdate  = yy*10000+mm*100+dd
       rst_icsec   = tod
       open(80,file=fln1,form = "unformatted",status='unknown')
       write(80) rst_icdate,rst_icsec,v2last1
       close(80)

       fln2 = trim(beps_rst_dir)//"rpointer"
       open(81,file=fln2,form = "formatted",status = 'REPLACE')
       write(81,'(A)') fln1
       close(81)
!    end if
!    call mpi_barrier(mpi_comm_world,ierr)

else if(trim(flag) == "read") then
!   if(myid ==0) then
      fln2 = trim(beps_rst_dir)//"rpointer"
      open(81,file= fln2,form = "formatted")
      read(81,'(A)') fln1
      close(81)

      open(80,file=trim(fln1),form="unformatted")
      read(80) rst_icdate,rst_icsec,v2last1
      close(80)
!   end if
  
!   call mpi_barrier(mpi_comm_world,ierr)
!   do i = 1,PFT
!     do j = 0,40
!        call mpi_scatterv(v2last1(1,j,i),dp,sp,mpi_real8,v2last(1,j,i),npoints,mpi_real8,0,mpi_comm_world,ierr)
!     end do
!   end do

!   call mpi_bcast(rst_icdate,1,mpi_integer,0,mpi_comm_world,ierr)
!   call mpi_bcast(rst_icsec,1,mpi_integer,0,mpi_comm_world,ierr)
!   call mpi_barrier(mpi_comm_world,ierr)
v2last = v2last1
end if

end subroutine
end module
