

!! The CO2 concentration is referred to www.esrl.noaa.gov/gmd/ccgg/trends/

subroutine get_CO2_concentration(yr,CO2)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none

integer,intent(in)     :: yr
real(r8),intent(inout) :: CO2   ! from NOAA ESRL global annual mean data

select case (yr)
    case (1980)
        CO2 = 338.80
    case (1981)
        CO2 = 340.00
    case (1982)
        CO2 = 340.76
    case (1983)
        CO2 = 342.44
    case (1984)
        CO2 = 343.99
    case (1985)
        CO2 = 345.47
    case (1986)
        CO2 = 346.87
    case (1987)
        CO2 = 348.62
    case (1988)
        CO2 = 351.15
    case (1989)
        CO2 = 352.79
    case (1990)
        CO2 = 353.98
    case (1991)
        CO2 = 355.29
    case (1992)
        CO2 = 355.99
    case (1993)
        CO2 = 356.71
    case (1994)
        CO2 = 358.21
    case (1995)
        CO2 = 360.05
    case (1996)
        CO2 = 361.79
    case (1997)
        CO2 = 362.90
    case (1998)
        CO2 = 365.54
    case (1999)
        CO2 = 367.64
    case (2000)
        CO2 = 368.84
    case (2001)
        CO2 = 370.41
    case (2002)
        CO2 = 372.42
    case (2003)
        CO2 = 374.96
    case (2004)
        CO2 = 376.80
    case (2005)
        CO2 = 378.81
    case (2006)
        CO2 = 380.94
    case (2007)
       CO2 = 382.68
    case (2008)
       CO2 = 384.78
    case (2009)
       CO2 = 386.29
    case (2010)
       CO2 = 388.57
    case (2011)
       CO2 = 390.45
    case (2012)
       CO2 = 392.46
    case (2013)
       CO2 = 395.19
    case (2014)
       CO2 = 397.12
    case (2015)
       CO2 = 399.41
    case (2016)
       CO2 = 402.86
    case (2017)
       CO2 = 405.00
    case (2018)
       CO2 = 407.37
end select

end
