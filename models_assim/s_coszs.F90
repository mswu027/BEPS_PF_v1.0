!************************************************
! original code was written by W. Ju on July 2004
! in C program.
! This Fortran version is updated by Jun Wang on
!  5/12/2016
!************************************************

subroutine s_coszs(yr,mn,dd,tod,day,doy,lat,lon,CosZs,hr,Hsolar1)
use shr_kind_mod,only: r8=>shr_kind_r8
use beps_con,only: PI
!--iLab::restrict use of beps_time_manager to required entities
!--iLab-update::extended argument list allows to remove beps_time_manager
! use beps_time_manager, only:get_calendar, NO_LEAP_C, GREGORIAN_C
implicit none
!--iLab::current date/day now input in order to avoid invoking beps_time_manager
integer, intent(in) :: yr,mn,dd,tod
integer, intent(in) :: day !-- calendar day
integer, intent(in) :: doy !-- #days in year

!!integer,intent(in)  :: day,hour    ! the hour in day and the day in year (hour begins from 0, ends at 23)
real(r8),intent(in) :: lat,lon
real(r8),intent(out):: CosZs       ! Solar Zenith
real(r8),intent(out):: hr
real(r8),intent(out):: Hsolar1

real(r8)            :: lon1        ! change 0.5~359.5 into -179.5~179.5
integer             :: hour
! integer             :: yr,mn,dd,tod
!iLab::length should (potentially) fit ESMF_MAXSTR,
!      but we don't want to have 'use EMSF' statement here and apply the value
character(len = 128)  :: calendar
! character(len = 30)  :: calendar
real(r8)  :: Delta,Lat_arc         !  Hsolar1

! call get_curr_date(yr,mn,dd,tod)
hour   = tod/3600


Delta = 0.006918-0.399912*cos(day*2.0*PI/doy)+0.070257*sin(day*2.0*PI/doy) &
        -0.006758*cos(day*4.0*PI/doy)+0.000907*sin(day*4.0*PI/doy)
!delta is the declination angle of sun.

!! longitude 0.5~359.5  => -179.5~179.5

!if(lon > 180 .and. lon < 360) then
!     lon1 = lon - 360.
!else
!     lon1 = lon
!end if
lon1=lon
hr=hour
!hr    = hour + lon1/15.0   !convert local time to UTC, but don"t need for xiuli
if(hr > 24) hr = 24 - hr
if(hr < 0 ) hr = 24 + hr

 Lat_arc   = PI*lat/180.0
 Hsolar1   = (hr - 12.0)*2.0*PI/24.0  !local hour angle in arc.
 CosZs     = cos(Delta)*cos(Lat_arc)*cos(Hsolar1)+sin(Delta)*sin(Lat_arc)

end subroutine

