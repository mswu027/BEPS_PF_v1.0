!--iLab::convf no longer required
! subroutine cos_plant(lc,cosa,convf,g_sw,g_b,vcmax,ffpa,f_soilwater,cos_assim)
subroutine cos_plant(lc,cosa,g_sw,g_b,vcmax,ffpa,f_soilwater,cos_assim)
    use shr_kind_mod, only: r8 =>shr_kind_r8
    implicit none

    !Input Variables
    ! real(r8) :: vcmax, g_sw, g_b,convf,f_soilwater,cos_assim
    real(r8) :: vcmax, g_sw, g_b,f_soilwater,cos_assim
    integer  :: lc
    !Local Variables
    real(r8) :: c4flag,ffpa
    real(r8) :: cosa   ! CAS COS concentration (pmol COS/mol air)
    real(r8) :: gcosm
    real(r8) :: gtcos

    c4flag = 0.

    if (lc==40 .or. lc==41) c4flag = 1.0
    g_sw = max(1.e-6,g_sw)
    g_b = max(1.e-6,g_b)
    gcosm = 1.40e3 * vcmax * 1.0e-6 * (1.0 + 5.33*c4flag) * ffpa * f_soilwater     ! mol/m2/s
    gcosm = max(gcosm,1.e-6)
    !gtcos = 1.0/(1.94/(g_sw*convf) + 1.0/(g_b*convf) + 1.0/(gcosm*convf))  ! m/s
    gtcos = 1.0/(1.94/(g_sw) + 1.0/(g_b) + 1.0/(gcosm))  ! mol/m2/s
    cos_assim = gtcos * cosa                                               ! pmol/m2/s

end subroutine
