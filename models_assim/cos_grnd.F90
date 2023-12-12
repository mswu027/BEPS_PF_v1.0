subroutine cos_grnd(soilp,cos_soil)
    use shr_kind_mod, only: r8 =>shr_kind_r8
    use beps_soilMod
    implicit none

    !Input Variables
    type(soil), intent(in)           ::soilp

    !Local Variables
    real(r8) :: cos_soil    ! local ground COS flux (pmol/m2/sec)
    real(r8) :: F_opt, S_opt, F_g, S_g, a
    real(r8), parameter :: k_cos_soil = 1.2E-4

    !Misc Variables
    integer :: j
    real(r8):: soil_T, soil_S, dsoil, soil_ice
    real(r8):: cos_soil_abiotic, cos_soil_biotic
    intrinsic log,dble,exp

    soil_T    = 0.
    soil_S    = 0.
    soil_ice  = 0.
    dsoil     = 0.

    !...ground uptake of COS, calculated from Whelan et al., 2016, ACP. calculate the abiotic and biotic part of ground uptake separately.
    do j = 1,3
        soil_T = soil_T + soilp%temp_soil_c(j-1)*soilp%d_soil(j-1)
        soil_S = soil_S + soilp%thetam(j-1)*soilp%d_soil(j-1)
        soil_ice = soil_ice + soilp%ice_ratio(j-1)*soilp%d_soil(j-1)
        dsoil = dsoil + soilp%d_soil(j-1)
    end do

    soil_T = soil_T/dsoil
    soil_ice = soil_ice/dsoil

    cos_soil_abiotic = 0.437 *exp(0.0984 * soil_T)

    F_opt = -0.00986 * soil_T * soil_T + 0.197 * soil_T - 9.32
    S_opt = 0.28 * soil_T + 14.5
    F_g = -0.0119 * soil_T * soil_T + 0.110 * soil_T -1.18
    S_g  = 35.0
    a = log(F_opt/F_g) * (log(S_opt/S_g) + (S_g/S_opt - 1.))**(-1)

    cos_soil_biotic = F_opt * (soil_S/S_opt)**a * exp(-a * (soil_S/S_opt -  1.))

    cos_soil = - (cos_soil_abiotic + cos_soil_biotic)    ! uptake is positive, so emission is negative

end subroutine
