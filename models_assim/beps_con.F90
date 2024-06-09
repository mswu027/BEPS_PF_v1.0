!*********************************************
! Function: constants used in BEPS model
! Created : Jun Wang
! Date    : 2016/12/2
!*********************************************

module beps_con
use shr_kind_mod, only: r8=>shr_kind_r8
implicit none

real(r8),parameter:: PI    = 3.1415926
real(r8),parameter:: rho_a = 1.292              ! density of air at 0C
real(r8),parameter:: rho_w = 1025.0             ! density of water kg/m3
real(r8),parameter:: cp_ice= 2228.261           ! J/kg/K
real(r8),parameter:: latent_fusion=3.34*1000000.
real(r8),parameter:: PI2   = 6.283185307        ! 2 times PI
real(r8),parameter:: zero  = 0.0000000001

real(r8),parameter:: rugc  = 8.314              ! J mole-1 K-1
real(r8),parameter:: hkin  = 200000.0           ! enthalpy term, J mol-1
real(r8),parameter:: ejm   = 55000.0            !  // activation energy for electron transport, J mol-1
real(r8),parameter:: evc   = 55000.0            ! activation energy for carboxylation, J mol-1
real(r8),parameter:: kc25  = 274.6              ! kinetic coef for CO2 at 25 C, microbar
real(r8),parameter:: ko25  = 419.8              ! kinetic coef for O2 at 25C,  millibars
real(r8),parameter:: o2    = 210.0              ! oxygen concentration  mmol mol-1
real(r8),parameter:: tau25 = 2904.12            ! tau coefficient

!Arrhenius constants
!Eact for Michaelis-Menten const. for KC, KO and dark respiration
!These values are from Harley
real(r8),parameter:: ekc   = 80500.0            ! Activation energy for K of CO2; J mol-1
real(r8),parameter:: eko   = 14500.0            ! Activation energy for K of O2, J mol-1
real(r8),parameter:: erd   = 38000.0            ! activation energy for dark respiration, eg Q10=2
real(r8),parameter:: ektau = -29000.0           ! J mol-1 (Jordan and Ogren, 1984)
real(r8),parameter:: tk_25 = 298.16             ! absolute temperature at 25 C
real(r8),parameter:: toptvc= 301.0              ! optimum temperature for maximum carboxylation
real(r8),parameter:: toptjm= 301.0              ! optimum temperature for maximum electron transport

real(r8),parameter:: sigma = 5.67e-08           ! Stefan-Boltzmann constant W M-2 K-4
real(r8),parameter:: mass_air  = 29.0           ! air g mmol-1
real(r8),parameter:: mass_CO2  = 44.0           ! CO2 g mol-1

contains

!*******************************************************
!*   SUBROUTINE abort                                   *
!*   Forces the program to exit with a non-zero status *
!*******************************************************
  subroutine abortf
    implicit none
    print *,'============================================================'
    print *,'Forcing an error condition by opening a non-existing file...'
    print *,'============================================================'
    open(1,file='/../inexistent',status='old') ! This forces the exit
    stop ! This is to make flow analysis easier
  end subroutine abortf

end module
