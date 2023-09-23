

module mid_results
use shr_kind_mod,only: r8=>shr_kind_r8
implicit none

type :: results
   real(r8)  :: gpp_o_sunlit
   real(r8)  :: gpp_u_sunlit
   real(r8)  :: gpp_o_shaded
   real(r8)  :: gpp_u_shaded
   real(r8)  :: plant_resp
   real(r8)  :: npp_o
   real(r8)  :: npp_u
   real(r8)  :: GPP
   real(r8)  :: SIF
   real(r8)  :: NPP
   real(r8)  :: NEP
   real(r8)  :: soil_resp
   real(r8)  :: Net_Rad
   real(r8)  :: SH
   real(r8)  :: LH
   real(r8)  :: Trans
   real(r8)  :: Evap
   real(r8)  :: thetam_surf
   real(r8)  :: COS_flux
   real(r8)  :: lai
   real(r8)  :: lai_old
   real(r8)  :: lai_new
   real(r8)  :: COS_plant
   real(r8)  :: COS_grnd
!   real(r8)  :: COS_grnd2
   real(r8)  :: fAPAR
   real(r8)  :: PWS
   real(r8)  :: ETa
   real(r8)  :: VOD
   real(r8)  :: fei_leaf
end type results

type :: climatedata
   real(r8) :: temp      ! temperatre  (oC)
   real(r8)  :: Srad     !solar radiation
   real(r8)  :: LR       !downward longwave radiation
   real(r8)  :: rainfall !liquid water rainfall
   real(r8)  :: snow     !snow
   real(r8)  :: S_dff    !diffuse solar radiation
   real(r8)  :: S_dir    !direct solar radiation
   real(r8)  :: rh
   real(r8)  :: wind
   real(r8)  :: tempmx
   real(r8)  :: tempmn
end type climatedata


!! These types are adapted from DB.h, but many variables are not used @J.Wang
!! Importantly, these three types are only used by photosynthesis module.
type :: meteorology
   real(r8)  :: ustar                 !friction velocity, m s-1
   real(r8)  :: ustarnew              !updated friction velocity with new H, m s-1
   real(r8)  :: rhova_g               !absolute humidity, g m-3
   real(r8)  :: rhova_kg              !absolute humidity, kg m-3
   real(r8)  :: sensible_heat_flux    !sensible heat flux, W M-2
   real(r8)  :: H_old                 !old sensible heat flux, W m-2
   real(r8)  :: air_density           !air density, kg m-3
   real(r8)  :: T_Kelvin              !absolute air temperature, K
   real(r8)  :: press_kpa             !station pressure, kPa
   real(r8)  :: press_bars            !station pressure, bars
   real(r8)  :: press_Pa              !pressure, Pa
   real(r8)  :: pstat273              !gas constant computations
   real(r8)  :: air_density_mole      !air density, mole m-3
   real(r8)  :: relative_humidity     !relative humidity, ea/es(T)
   real(r8)  :: vpd                   !vapor pressure deficit
   real(r8)  :: ir_in                 !infrared flux density
end type meteorology

type  :: factors
   real(r8) :: latent                 !latent heat of vaporization, J kg-1
   real(r8) :: latent18               !latent heat of vaporization times molecular mass of vapor, 18 g mol-1
   real(r8) :: heatcoef               !factor for sensible heat flux density
   real(r8) :: a_filt                 !filter coefficients
   real(r8) :: b_filt                 !filter coefficients
   real(r8) :: co2                    !CO2 factor, ma/mc * rhoa (mole m-3)
end type factors

type :: boundary_layer_resistances
   real(r8) :: vapor                  !resistance for water vapor, s/m
   real(r8) :: heat                   !resistance for heat, s/m
   real(r8) :: co2                    !resistance for CO2, s/m
end type boundary_layer_resistances

contains
  !-- iLab::added for debugging purpose
  subroutine midres_dump(midres, fmt)
    implicit none
    type(results), intent(in) :: midres
    character(len=*), intent(in), optional :: fmt
    !-- local
    character(len=16) :: r8_fmt
    if( present(fmt) ) then
       r8_fmt = '(a15,'//trim(fmt)//')'
    else
       r8_fmt = '(a15,e25.16)'
    endif
    write(*, r8_fmt) 'gpp_o_sunlit=',midres%gpp_o_sunlit
    write(*, r8_fmt) 'gpp_u_sunlit=',midres%gpp_u_sunlit
    write(*, r8_fmt) 'gpp_o_shaded=',midres%gpp_o_shaded
    write(*, r8_fmt) 'gpp_u_shaded=',midres%gpp_u_shaded
    write(*, r8_fmt) 'plant_resp=',midres%plant_resp
    write(*, r8_fmt) 'npp_o=',midres%npp_o
    write(*, r8_fmt) 'npp_u=',midres%npp_u
    write(*, r8_fmt) 'GPP=',midres%GPP
    write(*, r8_fmt) 'SIF=',midres%SIF
    write(*, r8_fmt) 'NPP=',midres%NPP
    write(*, r8_fmt) 'NEP=',midres%NEP
    write(*, r8_fmt) 'soil_resp=',midres%soil_resp
    write(*, r8_fmt) 'Net_Rad=',midres%Net_Rad
    write(*, r8_fmt) 'SH=',midres%SH
    write(*, r8_fmt) 'LH=',midres%LH
    write(*, r8_fmt) 'Trans=',midres%Trans
    write(*, r8_fmt) 'Evap=',midres%Evap
    write(*, r8_fmt) 'thetam_surf=',midres%thetam_surf
    write(*, r8_fmt) 'COS_flux=',midres%COS_flux
    write(*, r8_fmt) 'lai=',midres%lai
    write(*, r8_fmt) 'lai_old=',midres%lai_old
    write(*, r8_fmt) 'lai_new=',midres%lai_new
    write(*, r8_fmt) 'COS_plant=',midres%COS_plant
    write(*, r8_fmt) 'COS_grnd=',midres%COS_grnd
    !write(*, r8_fmt) 'COS_grnd2=',midres%COS_grnd2
    write(*, r8_fmt) 'fapar=',midres%fapar
  end subroutine midres_dump
end module

