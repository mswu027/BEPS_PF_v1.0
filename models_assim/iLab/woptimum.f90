!***********************************************************
!     woptimum
!
!> @brief writes control vector (in normalised units) to binary file,
!         writes parameter uncertainties (physical units) to binary file
!
!> @details 
!
!> @param[in]  n  overall length of control vector
!> @param[in]  x0(n) priori control vector
!> @param[in]  sx(n) prior control vector uncertainties
!> @param[in]  x(n)  posterior control vector
!
!> \authors Michael Vossbeck, The Inversion Lab
!> \date    January 2020
subroutine woptimum(n,x0,sx,x)
  implicit none
  !-- arguments
  integer(kind=4), intent(in) ::  n
  real(kind=8),    intent(in) :: x(n), x0(n), sx(n)
  !-- local
  character(len=*), parameter :: sub = 'woptimum'
  character(len=*), parameter :: ctl_outname = 'xopt.b'
  character(len=*), parameter :: ctlunc_outname = 'sx.b'
  integer(kind=4) :: i
  real(kind=8) :: eps

  !-- initialise
  eps = epsilon(eps)*10

  write ( *, '(a)' ) ' '
  write ( *, '(a9x,4a20)' ) 'i', 'Prior', 'Posterior', 'Change [% of x0]', 'Chg [% of sx0]'
  write ( *, '(a)' ) ' '
  do i = 1, n
     if(x0(i).lt.eps) then
        print '(i9x,2e20.6,2f20.5)', &
             i, sx(i)*x0(i), sx(i)*x(i), 100*(x(i)-x0(i)), 100*(x(i)-x0(i))
     else
        print '(i9x,2e20.6,2f20.5)', &
             i, sx(i)*x0(i), sx(i)*x(i), 100*(x(i)-x0(i))/x0(i), 100*(x(i)-x0(i))
     endif
  enddo

  open (unit=1, file=ctl_outname, form='unformatted')
  write (1) x
  close (1)
  write(*, '(a)') ' INFO::'//sub//': writing final control vector to ***'//trim(ctl_outname)//'***'

  open (unit=1, file=ctlunc_outname, form='unformatted')
  write (1) sx
  close (1)
  write(*, '(a)') ' INFO::'//sub//': writing prior control vector uncertainty to ***'//&
       trim(ctlunc_outname)//'***'
end subroutine woptimum
