!! The CO2 concentration is referred to www.esrl.noaa.gov/gmd/ccgg/trends/
!! update @ Lu Hu: 1979-2023, 2024/03/30, https://gml.noaa.gov/ccgg/trends/global.html
subroutine get_CO2_concentration(yr,CO2)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none

integer,intent(in)     :: yr
real(r8),intent(inout) :: CO2   ! from NOAA ESRL global annual mean data

select case (yr)
    case (1979)
        CO2 = 336.85
    case (1980)
        CO2 = 338.91
    case (1981)
        CO2 = 340.11
    case (1982)
        CO2 = 340.86
    case (1983)
        CO2 = 342.53
    case (1984)
        CO2 = 344.07
    case (1985)
        CO2 = 345.54
    case (1986)
        CO2 = 346.97
    case (1987)
        CO2 = 348.68
    case (1988)
        CO2 = 351.16
    case (1989)
        CO2 = 352.79
    case (1990)
        CO2 = 354.06
    case (1991)
        CO2 = 355.39
    case (1992)
        CO2 = 356.09
    case (1993)
        CO2 = 356.83
    case (1994)
        CO2 = 358.33
    case (1995)
        CO2 = 360.17
    case (1996)
        CO2 = 361.93
    case (1997)
        CO2 = 363.05
    case (1998)
        CO2 = 365.70
    case (1999)
        CO2 = 367.79
    case (2000)
        CO2 = 368.96
    case (2001)
        CO2 = 370.57
    case (2002)
        CO2 = 372.58
    case (2003)
        CO2 = 375.14
    case (2004)
        CO2 = 376.95
    case (2005)
        CO2 = 378.98
    case (2006)
        CO2 = 381.15
    case (2007)
        CO2 = 382.90
    case (2008)
        CO2 = 385.02
    case (2009)
        CO2 = 386.50
    case (2010)
        CO2 = 388.76
    case (2011)
        CO2 = 390.63
    case (2012)
        CO2 = 392.65
    case (2013)
        CO2 = 395.40
    case (2014)
        CO2 = 397.34
    case (2015)
        CO2 = 399.65
    case (2016)
        CO2 = 403.07
    case (2017)
        CO2 = 405.22
    case (2018)
        CO2 = 407.61
    case (2019)
        CO2 = 410.07
    case (2020)
        CO2 = 412.44
    case (2021)
        CO2 = 414.70
    case (2022)
        CO2 = 417.07
    case (2023)
        CO2 = 419.25
end select

end
