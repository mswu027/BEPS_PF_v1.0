subroutine flux_at_S(soil_S_in, F_opt_in, F_g_in, S_g_in, S_opt_in, result)
    real, intent(in) :: soil_S_in, F_opt_in, F_g_in, S_g_in, S_opt_in
    real, intent(out) :: result
    real :: a
	
    a = log(F_opt_in/F_g_in) * (log(S_opt_in/S_g_in) + (S_g_in/S_opt_in - 1.))**(-1)
    result = F_opt_in * (soil_S_in / S_opt_in) ** a * exp(-a * ((soil_S_in / S_opt_in) - 1))
	
end subroutine


subroutine cos_grnd(soilp,cos_soil,LC)
    use shr_kind_mod, only: r8 =>shr_kind_r8
    use beps_soilMod
    implicit none

    !Input Variables
    type(soil), intent(in)   ::soilp
	integer,intent(in)   :: LC                                          ! landcover type

    !Local Variables
    real(r8) :: cos_soil                                                ! local ground COS flux (pmol/m2/sec)
    real(r8) :: F_opt, S_opt, F_g, S_g, alpha , beta, depth_ratio       ! add alpha, beta, depth_ratio   2023/12/16  Huajie Zhu          delete a  2024/01/18 Huajie Zhu  
	real(r8) :: F_opt_T, F_g_T                                          ! add F_opt_T and F_g_T for temperature correction        2024/01/18   Huajie Zhu
    real(r8), parameter :: k_cos_soil = 1.2E-4

    !Misc Variables
    integer :: j
    real(r8):: soil_T, soil_S, dsoil, soil_ice
    real(r8):: cos_soil_abiotic, cos_soil_biotic
    intrinsic log,dble,exp

    !...ground uptake of COS, calculated from Whelan et al., 2016, ACP. calculate the abiotic and biotic part of ground uptake separately.
    !   the paramerization scheme of the abiotic part was taken from Abadie et al. 2022, BG (see Table A2, and Table 2 in Whelan et al, 2016, ACP for details).   
    !   the paramerization scheme of the biotic part was taken from Whelan et al. 2022, JGRB (see Sect. 2.2 and Table 3 for details).   
    !   Referring to Abadie et al. 2022 and Whelan et al. 2022, here we only calculated the COS soil flux from the top 9 cm ( Top layer (5cm) + 0.4 * second layer (10cm) in BEPS)

    soil_T    = 0.
    soil_S    = 0.
    ! soil_ice  = 0.
    ! dsoil     = 0.

    ! do j = 1,3
        ! soil_T = soil_T + soilp%temp_soil_c(j-1)*soilp%d_soil(j-1)
        ! soil_S = soil_S + soilp%thetam(j-1)*soilp%d_soil(j-1)
        ! soil_ice = soil_ice + soilp%ice_ratio(j-1)*soilp%d_soil(j-1)
        ! dsoil = dsoil + soilp%d_soil(j-1)
    ! end do

    ! soil_T = soil_T/dsoil
    ! soil_ice = soil_ice/dsoil
	
	depth_ratio = 0.4                                      ! top 9 cm ( Top layer (5cm) + 0.4 * second layer (10cm) in BEPS)
    soil_T = (soilp%temp_soil_c(0)*soilp%d_soil(0)+soilp%temp_soil_c(1)*soilp%d_soil(1)*depth_ratio)
    soil_T = soil_T/(soilp%d_soil(0)+soilp%d_soil(1)*depth_ratio)               ! unit: ℃
    soil_S = (soilp%thetam(0)*soilp%d_soil(0)+soilp%thetam(1)*soilp%d_soil(1)*depth_ratio)
    soil_S = soil_S/(soilp%d_soil(0)+soilp%d_soil(1)*depth_ratio) * 100         ! unit: %
	
	select case (LC)
		case (1)       ! conifer evergreen
			alpha = -7.77
			beta = 0.119
			
		case (2)       ! conifer decidous
			alpha = -7.77
			beta = 0.119		
			
		case (6)       ! broadleaf decidous
			alpha = -7.77
			beta = 0.119

		case (9)       !  broadleaf evergreen 
			alpha = -7.77
			beta = 0.119			

		case (10)       ! mix forest
			alpha = -7.77
			beta = 0.119

		case (13)       ! shrub                    
			alpha = -7.77
			beta = 0.119	

		case (14)       ! grass
			alpha = -9.54
			beta = 0.108	

		case (15)       ! crop
			alpha = -6.12
			beta = 0.108				

		case (40)       ! C4 grass
			alpha = -9.54
			beta = 0.108		

		case (41)       ! C4 crop
			alpha = -6.12
			beta = 0.108				
	end select

    cos_soil_abiotic = exp(alpha + beta * soil_T)

    ! F_opt = -0.00986 * soil_T * soil_T + 0.197 * soil_T - 9.32
    ! S_opt = 0.28 * soil_T + 14.5
    ! F_g = -0.0119 * soil_T * soil_T + 0.110 * soil_T -1.18
    ! S_g  = 35.0
	
	if(soil_T>0)then
		select case (LC)
			case (1)       ! conifer evergreen
				S_opt = 12.5                                   ! accroding to Broadleaf or needleleaf 
				F_opt = -18                                    ! In whelan's COS soil model, negetive indicate COS uptake.
				S_g = 19.3
				F_g = -5.9
				call flux_at_S(soil_T, F_opt, -12.0, 35.0, 28.0, F_opt_T)           ! apply tempreture correction with a same model of moisture correction, here soil_T must be positive  
				call flux_at_S(soil_T, F_g, -3.8, 35.0, 28.0, F_g_T)              ! apply tempreture correction with a same model of moisture correction
				call flux_at_S(soil_S, F_opt_T, F_g_T, S_g, S_opt, cos_soil_biotic)         ! apply moisture correction based on F_opt_T and F_g_T
				
			case (2)       ! conifer decidous
				S_opt = 12.5                                  
				F_opt = -18                               
				S_g = 19.3
				F_g = -5.9
				call flux_at_S(soil_T, F_opt, -12.0, 35.0, 28.0, F_opt_T)    
				call flux_at_S(soil_T, F_g, -3.8, 35.0, 28.0, F_g_T)      
				call flux_at_S(soil_S, F_opt_T, F_g_T, S_g, S_opt, cos_soil_biotic)		
				
			case (6)       ! broadleaf decidous
				S_opt = 24.6       
				F_opt = -12.6
				S_g = 51
				F_g = -0.18 * soil_T + 0.48
				call flux_at_S(soil_S, F_opt, F_g, S_g, S_opt, cos_soil_biotic)             ! note, here F_opt and F_g was used.

			case (9)       !  broadleaf evergreen 
				S_opt = 24.6       
				F_opt = -12.6
				S_g = 51
				F_g = -0.18 * soil_T + 0.48			
				call flux_at_S(soil_S, F_opt, F_g, S_g, S_opt, cos_soil_biotic)
				
			case (14)       ! grass
				S_opt = 12.5
				F_opt = -4.5
				S_g = 26.9
				F_g = -2.3
				call flux_at_S(soil_T, F_opt, -1.5, 25.0, 10.9, F_opt_T)          
				call flux_at_S(soil_T, F_g, -1.3, 25.0, 14.8, F_g_T)              
				call flux_at_S(soil_S, F_opt_T, F_g_T, S_g, S_opt, cos_soil_biotic)
				
			case (15)       ! crop
				S_opt = 17.7
				F_opt = -9.7
				S_g = 22
				F_g = -5.36	
				call flux_at_S(soil_S, F_opt, F_g, S_g, S_opt, cos_soil_biotic)

            case (40)       ! C4 grass
				S_opt = 12.5
				F_opt = -4.5
				S_g = 26.9
				F_g = -2.3	
				call flux_at_S(soil_T, F_opt, -1.5, 25.0, 10.9, F_opt_T)           
				call flux_at_S(soil_T, F_g, -1.3, 25.0, 14.8, F_g_T)          
				call flux_at_S(soil_S, F_opt_T, F_g_T, S_g, S_opt, cos_soil_biotic)
				
            case (41)       ! C4 crop
				S_opt = 17.7
				F_opt = -9.7
				S_g = 22
				F_g = -5.36
	            call flux_at_S(soil_S, F_opt, F_g, S_g, S_opt, cos_soil_biotic)
        end select
    else
        cos_soil_biotic = 0
    end if

    cos_soil = -(cos_soil_abiotic + cos_soil_biotic)         ! In our model, positive indicate COS uptake                  

end subroutine
