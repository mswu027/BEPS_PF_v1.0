

! this module will calculate sensible heat from overstorey,understorey, and ground editted by XZ Luo, May23,2015
! inputs: 
! temperature of sunlit and shaded leaves from other storey (leaf temperature module)
! temperature of air,relative humidity
! temperature of ground (soil heat flux module)
! aerodynamic heat conductance of sunlit shaded leaves from over/understorey
! aerodynamic heat conductance of ground
! lAI sunlit and shaded, over/understorey (LAI module)
! Outputs:
! sensible heat from over/understorey and ground

subroutine SensibleHeat(tempL_o_sunlit,tempL_o_shaded,tempL_u_sunlit,tempL_u_shaded,temp_g,temp_air,rh_air, &
                       Gheat_o_sunlit,Gheat_o_shaded,Gheat_u_sunlit,Gheat_u_shaded,Gheat_g,lai_o_sunlit, &
                       lai_o_shaded,lai_u_sunlit,lai_u_shaded,SH_o,SH_u,SH_g)
use meteoMod
use shr_kind_mod, only: r8=>shr_kind_r8
implicit none
real(r8),intent(in)  :: tempL_o_sunlit,tempL_o_shaded,tempL_u_sunlit,tempL_u_shaded,temp_g,temp_air,rh_air
real(r8),intent(in)  :: Gheat_o_sunlit,Gheat_o_shaded,Gheat_u_sunlit,Gheat_u_shaded,Gheat_g
real(r8),intent(in)  :: lai_o_sunlit,lai_o_shaded,lai_u_sunlit,lai_u_shaded
real(r8),intent(out) :: SH_o,SH_u,SH_g

real(r8)  :: SH_o_sunlit,SH_o_shaded,SH_u_sunlit,SH_u_shaded

call meteo_pack(temp_air,rh_air)

SH_o_sunlit  = (tempL_o_sunlit - temp_air)*density_air*cp_air*Gheat_o_sunlit
SH_o_shaded  = (tempL_o_shaded - temp_air)*density_air*cp_air*Gheat_o_shaded

SH_u_sunlit  = (tempL_u_sunlit - temp_air)*density_air*cp_air*Gheat_u_sunlit
SH_u_shaded  = (tempL_u_shaded - temp_air)*density_air*cp_air*Gheat_u_shaded

SH_o         = SH_o_sunlit*lai_o_sunlit + SH_o_shaded*lai_o_shaded
SH_u         = SH_u_sunlit*lai_u_sunlit + SH_u_shaded*lai_u_shaded

SH_o         = max(-200.,SH_o)
SH_u         = max(-200.,SH_u)
SH_g         = (temp_g -temp_air)*density_air*cp_air*Gheat_g

return
end subroutine
