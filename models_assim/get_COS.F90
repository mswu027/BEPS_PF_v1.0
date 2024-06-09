subroutine get_COS_concentration(yr,COSA)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none

integer,intent(in)     :: yr
real(r8),intent(inout) :: COSA

select case (yr)
    case(2000)
        COSA = 490.56
    case (2001)
        COSA = 497.21
    case (2002)
        COSA = 488.86
    case (2003)
        COSA = 495.70
    case (2004)
        COSA = 494.16
    case (2005)
        COSA = 496.59
    case (2006)
        COSA = 497.51
    case (2007)
       COSA = 503.40                  ! ppt
    case (2008)
       COSA = 500.40
    case (2009)
       COSA = 492.99
    case (2010)
       COSA = 498.47
    case (2011)
       COSA = 498.28
    case (2012)
       COSA = 504.20
    case (2013)
       COSA = 504.85
    case (2014)
       COSA = 502.03
    case (2015)
       COSA = 500.89
    case (2016)
       COSA = 508.41
    case (2017)
       COSA = 499.40
    case (2018)
       COSA = 491.96
    case (2019)
       COSA = 493.07
end select

end
