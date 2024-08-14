!******************************************
!Function: parameters in BEPS
!Created : Jun Wang
!Date    : 2016/12/2
!******************************************

module beps_par
use shr_kind_mod,only: r8=>shr_kind_r8
implicit none

! For plant
integer,parameter :: PFT  = 10 !previous 9, with c4 grass and crop added
 !(/1,2,6,9,10,13,14,15,40,41/) =>(/conifer evergreen,conifer decidous,broadleaf decidous,broadleaf evergreen,mix,shrub,grass,crop,c4grass,c4crop/)
integer,parameter :: texture = 11 ! 11 soil texture classes
! For Soil
integer :: nparameters=20000
integer :: parloop=200
integer :: layer = 5
integer :: FW_VERSION  = 1    ! 1 =>soil water uptake using R*fpsisr  other val=>soil water uptake using R
integer,parameter :: MAX_LAYERS  = 10   ! max layers
integer :: DEPTH_F     = 6    ! ??

! For timestep
#ifdef COUP_CSM
   integer :: step       =  1200     ! time step  (s)
   integer :: kstep      =  300        ! loop in timestep (units:s)
   integer :: kloop      =  4        ! loop counts
#else
   integer :: step       = 3600
   integer :: kstep      = 360
   integer :: kloop      = 10
#endif

! CO2 concentration
real(r8):: CO2_air    = 380   !ppm
real(r8):: COS_air    = 450   !ppt
! FOR MPI
integer   :: nlp              ! num of land points in total
integer,allocatable :: stype(:) ! land mask 1=>land
integer,allocatable :: mapping(:)  !map land points
integer,allocatable :: dp(:)  ! ! number points per processor
integer,allocatable :: sp(:)  ! ! start point
integer,allocatable :: ep(:)  ! ! end point

! For each processor
integer   :: myid,nproc
integer   :: npoints


! For initial/restart/branch
integer, public :: nsrest             = 0
integer, public, parameter :: nsrStartup  = 0
integer, public, parameter :: nsrContinue = 1
integer, public, parameter :: nsrBranch   = 2

end module

