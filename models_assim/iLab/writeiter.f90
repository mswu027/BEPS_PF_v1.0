!***********************************************************
!     writeiter
!
!> @brief writes control vector or cost function gradient associated to
!>        a particular iteration to binary output file.
!>        need to be performed.
!
!> @details 
!
!> @param[in]  iunit    unit number to be associated to binary file
!> @param[in]  icount   iteration counter
!> @param[in]  ichar    position of first digit in 4 digit counter within filename
!> @param[in]  n        length of (1D) array
!> @param[in]  x(n)     1D array (control vector or cost function gradient)
!> @param[in]  filename 
!
!> \authors Michael Vossbeck, The Inversion Lab
!> \date    September 2023
subroutine writeiter(iunit,icount,ichar,n, x,filename)
  implicit none
  !-- arguments
  integer, intent(in)          :: iunit, n
  integer, intent(in)          :: icount    ! counter
  integer, intent(in)          :: ichar     ! first character in filename to contain counter
  real(kind=8), intent(in)     :: x(n)
  character(*), intent (inout) :: filename 
  !-- local declarations
  character(len=*), parameter :: sub = 'writeiter'
  integer :: iostat
  if( icount.gt.9999 ) then
     write(*, '(a,i8)') ' FATAL::'//sub//': icount exceeds limit=9999, icount=', icount
     stop
  endif
  ! save control vector
  if(icount.gt.    0) write(filename(ichar+3:ichar+3),'(i1)') icount
  if(icount.gt.    9) write(filename(ichar+2:ichar+3),'(i2)') icount
  if(icount.gt.   99) write(filename(ichar+1:ichar+3),'(i3)') icount
  if(icount.gt.  999) write(filename(ichar+0:ichar+3),'(i4)') icount
  open (unit=iunit, file=filename, form='unformatted', action='write', iostat=iostat)
  if( iostat.ne.0 ) then
     write(*, '(a)') ' FATAL::'//sub//': file ***'//trim(filename)//'***'// &
          ' was not opened properly!'
     stop
  else
     write (iunit) x
     close (iunit)
  endif
end subroutine writeiter
