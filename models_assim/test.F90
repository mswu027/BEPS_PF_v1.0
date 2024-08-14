 program test

 use,intrinsic :: iso_fortran_env, only: wp => real64
 use pyplot_module

 implicit none

 real(wp),dimension(100) :: x,sx
 type(pyplot) :: plt
 integer :: i

 !generate some data:
 x = [(real(i,wp), i=0,size(x)-1)]/5.0_wp
 sx = sin(x)

 !plot it:
 call plt%initialize(grid=.true.,xlabel='angle (rad)',&
                     title='Plot of $\sin(x)$',legend=.true.)
 call plt%add_plot(x,sx,label='$\sin(x)$',linestyle='b-o',markersize=5,linewidth=2)
 call plt%savefig('sinx.png', pyfile='sinx.py')

 end program test