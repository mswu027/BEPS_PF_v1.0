subroutine beps_phenology(lc,daylen,dt,theta,trans,lai)

!***********************************************************
!* WOK, 2008-07-30
!* IMPLIFIED PHENOLOGY MODEL
!* simplified phenology model
!***********************************************************
  use shr_kind_mod, only: r8 =>shr_kind_r8
  USE beps_helper, ONLY : errf, mins, maxs, minx, maxx, fominef_ss, fomaxef_ss

  IMPLICIT NONE

!  REAL, ALLOCATABLE, DIMENSION (:,:,:) :: mlai ! monthly LAI fields from external data
  integer,intent(in)  :: lc
  real,intent(in)     :: daylen
  real,intent(in)     :: dt
  real,intent(in)     :: theta
  real,intent(in)     :: trans
  real,intent(inout)  :: lai
  real                :: tmpm                  ! air-temperature memory [deg C]
  real                :: laim                  ! water limited LAI memory
  real                :: laihi                 ! highest recorded LAI (with a decay rate, for setting 'zfc')
! WOK-ADD-070723 litter production to be calculated directly in phenology (not in cbalance indirectly)
  real                :: tmpmmult, laimmult    ! auxiliary fields
  real                :: leafshed              ! output field
  real                :: laihimult
  real, PARAMETER     :: taulaihi = 5.0        ! memory time for updating fractional cover
!  REAL, PARAMETER :: laimin = 1e-6            ! minimum LAI for pot. transpiration per LAI estimates
  real, PARAMETER     :: eta = 0.99999         ! curvature parameter for mins/maxs
! WOK-ADD-070723 the list of controlling parameters
! FREE PARAMETERS
  real                :: plaimax(1:10)               ! maximum LAI
  real                :: rootdepth(1:10)             ! rootdepth
  real                :: ptphen(1:10)                ! leaf onset temperature [deg C]
  real                :: ptphenr(1:10)               ! range of leaf onset temperature [1/deg C]
  real                :: pdphen                ! leaf shedding daylength [hours]
  real                :: pdphenr               ! range of leaf shedding daylength [hours]
!  REAL, ALLOCATABLE, DIMENSION (:) :: PTSHD   ! leaf shedding temperature [deg C]
!  REAL, ALLOCATABLE, DIMENSION (:) :: PTSHDS  ! spread of leaf shedding temperature [1/deg C]
  real                :: plgr                  ! leaf growth factor [1/days]
  real                :: pkl(1:10)                   ! inverse leaf longevity from start of senescense [1/days]
  real                :: ptauw(1:10)                 ! target survival time at current soil moisture [days]
! PARAMETERS LEFT FIXED
  real                :: pks                   ! inverse memory time for soil moisture-limited LAI [1/days]
  real                :: pkm                   ! inverse memory time for air temperature [1/days]
  real                :: pasm
  real                :: zfc, zlai
  real                :: fcmax0, lailim0, cdrm
  real                :: xdtmp, lait, laiw, fx, t0, ts, ft, fd, fg
  real                :: lailast
  real                :: laimaxw, laimax, r, lailim, wai
  real                :: dptrp
  integer             :: plt
! WOK-090309 'ph' is now only used in nscale
! xph     1: warm-evergreen; 2: cold-evergreen; 3: summergreen; 4: raingreen; 5: grass; 6: annual crop;
 ! INTEGER, DIMENSION (0:13), PARAMETER :: xph= &
!PFT:0  1  2  3  4  5  6  7  8  9 10 11 12 13
 !  (/5, 1, 4, 1, 3, 2, 3, 1, 4, 5, 5, 2, 5, 6/)
  real                :: sla(1:10)

select case (lc)
   case (1)    !conifer evergreen
    plt = 1
   case(2)      !conifer decidous
    plt = 2
   case(6)      !broadleaf decidous
    plt = 3
   case(9)      !broadleaf evergreen
    plt = 4
   case(10)     !mix
    plt = 5
   case(13)     !shrub
    plt = 6
   case(14)     ! grass
    plt = 7
   case(15)     ! crop
    plt = 8
   case(40)     ! C4 grass
    plt = 9
   case(41)     ! C4 crop
    plt = 10
   case(1001)     ! C4 crop
    plt = 10
   case(1002)     
    plt = 8
   case(1003)     
    plt = 8
   case(1004)     
    plt = 8
end select

  sla  = (/4.1, 11.3, 12.8, 7.8, 9.0, 9.2, 16.9, 25.3,16.9,16.9/)
!sla(0:13)=(/0., 9.9, 14.1, 5.7, 11.5, 4.1, 11.3, 6.9, &
!		& 11.5, 16.9, 16.9, 6.9, 16.9, 25.3/)
  ptphen = (/10.0, 10.0, 5.0, 0., 5., 4.0, 2.0, 15.0, 2.0,2.0/)
!  ptphen(0:13)=(/0.,0.,0., 0., 10.0, 10.0, 10.0, 0., 8.0, &
!       & 2.0, 2.0, 2.0, 2.0, 15.0/)

  ptphenr = (/2.0, 2.0, 2., 2.,2.,2.0, 2., 2., 2.,2./)
!  ptphenr(0:13)=(/0.,0., 0., 0., 2.0, 2.0, 2.0, 0., 2.0, &
!        & 2.0, 2.0, 2.0, 2.0, 2.0/)
  pdphen = 10.5
  pdphenr = 0.5
  plgr = 0.5

  pkl  = (/0.1, 0., 5.e-3, 0.1,0.1, 0.1, 5.e-3, 0.1, 5.e-3,5.e-3/)
 ! pkl(0:13) = (/0.,0.1, 0., 0.1, 5.e-3, 0.1, 0., 0.1, &
 !       & 0.1, 0.1, 5.e-3, 0.1, 0.1, 0.1/)

  ptauw  = (/30., 30., 30., 30., 30., 30., 30., 30., 30.,30./)
!  ptauw(1:13) = (/0., 30., 30., 30., 30., 30., 30., &
!        & 30., 30., 30., 30., 30., 30./)
!  LIST OF PFTs in BETHY:
!  1:  tropical broadleaf evergreen tree
!  2:  tropical broadleaf deciduous tree
!  3:  temperate broadleaf evergreen tree
!  4:  temperate broadleaf deciduous tree
!  5:  evergreen coniferous tree
!  6:  deciduous coniferous tree
!  7:  evergreen shrub
!  8:  deciduous shrub
!  9:  C3 grass
! 10:  C4 grass
! 11:  tundra
! 12:  swamp
! 13:  arable crop
  rootdepth = (/0.6, 0.6, 0.8, 0.8, 0.7, 0.5, 0.3, 0.3,0.3,0.3/)
  plaimax = (/4.5, 4.5, 4.5, 4.5, 4.5, 3.3, 3.0, 4.5, 3.0,3.0/)
  fcmax0    = 1.0
  lailim0   = 3.0
  cdrm      = 0.45
  pkm = 1. / 30.
  pks = 1. / 30.

!   multiplier for advancing temperature memory by one day
    tmpmmult = exp (-pkm)
!   multiplier for advancing soil-water limited LAI memory by one day
    laimmult = exp (-pks)
!   decay multiplier for evergreen LAI
!    laimult = exp (-pkl)
!   the air-temperature memory
    tmpm = 0.
!   the water stress index memory
    laim = 0.
    lai  = plaimax(plt)
!   decay multiplier for maximum LAI used to set fractional cover
    laihimult = exp (-1./(taulaihi*365.))
!   control for fractional cover 'zfc'
    laihi = 0.
    zfc = fcmax0
!   output field
    leafshed = 0.
!   spin-up of temperature memory
    tmpm = dt * (1. - tmpmmult) + tmpm * tmpmmult
!   calculate plant available soil moisture and daily potential transpiration, @MOUSONG.WU, 201905
    pasm = theta*rootdepth(plt)*1000. ! convert from m3/m3 to mm
    dptrp = trans*86400.*1000.       ! convert from m/s to mm

    lailast = lai
!------------------------------------------------------------------
! advances LAI and fractional cover by one day
! from its current state to the state at day 'iday'
!------------------------------------------------------------------
DO plt = 1,9
!      IF (ph(k)==1.or.ph(k)==4) THEN ! warm-evergreen and warm-deciduous phenology
      IF (plt==3 .or. plt==4 .or. plt==5 .or. plt==6) THEN ! warm-evergreen and warm-deciduous phenology
        ! effective maximum LAI, taking into account structural limiations
!        laimax = plaimax(k) * (1. - exp(-laimaxw/plaimax(k)))
        !   initialize LAI
        laimax = pasm * lailast / ptauw(plt) / maxx (dptrp, 1.e-3, 2.e-2)
!        laimax = mins (laimax, plaimax(k), 0.9)
        laimax = fominef_ss (laimax, plaimax(plt), 2.e-1)+1.e-1 !snb, test
        ! update water limited LAI memory
        laim = laimax * (1. - laimmult) + laim * laimmult
        ! rate of change of LAI towards limit
        r = plgr
        ! limit LAI
        lailim = laim
        ! update LAI
!        lai(k) = lailim - (lailim - lai(k)) * exp (-r)

!      ELSE IF (ph(k)==2.or.ph(k)==3) THEN ! cold-evergreen and cold-deciduous phenology
      ELSE IF (plt==1 .or. plt==2) THEN ! cold-evergreen and cold-deciduous phenology
        ! update memory of daily mean temperature
        tmpm = dt * (1. - tmpmmult) + tmpm * tmpmmult
        ! fraction of vegetation above temperature threshold
        ft = errf((tmpm-ptphen(plt))/ptphenr(plt))
        ! fraction of vegetation above daylength threshold
        fd = errf((daylen-pdphen)/pdphenr)
        r = ft * fd * plgr + (1. - ft * fd) * pkl(plt) + 1.e-9
        lailim = maxx (ft * fd * plgr * plaimax(plt) / r, 1.e-9, 5.e-3)
!        lai(k) = plaimax(k) - (plaimax(k) - lai(k)) * exp (-r)
!        lai(k) = lailim - (lailim - lai(k)) * exp (-r)

      ELSE ! grass and annual crop phenology
        tmpm = dt * (1. - tmpmmult) + tmpm * tmpmmult
        ft = errf((tmpm-ptphen(plt))/ptphenr(plt))
        laimax = pasm * lailast / ptauw(plt) / maxx (dptrp, 1.e-3, 2.e-2)
!        laimax = mins (laimax, plaimax(k), 0.9)
!        if(pft(k)==9) print '(a,8g30.14)','PFT9-a: iday,tmpm(k),ft,laimax,pasm(k),zlai(k),dptrp(k),lai(k)',iday,tmpm(k),ft,laimax,pasm(k),zlai(k),dptrp(k),lai(k)
        laimax = fominef_ss (laimax, plaimax(plt), 2.e-1)+1.e-1 !snb, test
        laim = laimax * (1. - laimmult) + laim * laimmult
        r = ft * plgr + (1. - ft) * pkl(plt) + 1.e-9
        lailim = maxx (ft * plgr * laim / r, 1.e-9, 5.e-3)
!        lai(k) = lailim - (lailim - lai(k)) * exp (-r)
      ENDIF
!      leafshed(k) = maxx (lailast - lai(k), 0., 1e-3) / sla(k) * 1000. * cdrm
      leafshed = maxx ((lailim-lai)*(1.-exp(-r)), 0., 1.e-3) / sla(plt) * 1000. * cdrm
      lai = lailim - (lailim - lai) * exp (-r)
!      if(pft(k)==9) print '(a,8g30.14)','PFT9-b: lai(k),r,lailim',lai(k),r,lailim
END DO
      ! set fractional cover
!      laihi(k) = maxs (lai(k), laihi(k), eta) * laihimult
      laihi = fomaxef_ss (lai, laihi, 2.e-6 ) * laihimult
!      zfc(k) = maxs (laihi(k) / lailim0, lai(k) / lailim0, eta)
      zfc = fomaxef_ss (laihi / lailim0, lai / lailim0, 2.e-6)
!      zfc(k) = mins ( zfc(k), 1., eta) * fcmax0
      zfc = minx( zfc, 1., 2.e-6) * fcmax0 ! snb, test

end subroutine

