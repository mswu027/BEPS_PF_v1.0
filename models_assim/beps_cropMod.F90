!*******************************************************
! Function: crop module( a  phenology  module  and  a  dry  matter  allocation  module .)
! C code: Xiaorong Wang
! Created : Xiuli Xing
! Date    : 2023/3/7
!*******************************************************
module beps_cropMod
use shr_kind_mod, only: r8=>shr_kind_r8
use beps_par
use bepstype
use bepstypeInit
use meteoMod
use beps_con
use mid_results
implicit none

!*************************************************parameters/
type,public::param_vars

    real(r8) :: t_bse  ! Base temp for development (K)
    real(r8) :: t_opt  ! Optimum temp for development (K).
    real(r8) :: t_max  ! Maximum temp for development (K).
    real(r8) :: tt_emr  ! Thermal time for emergence (degree days).
    real(r8) :: crit_pp  ! Critical daylength for photoperiod sensitivity (hours).
    real(r8) :: pp_sens  ! Sensitivity to daylength (hours-1).
    real(r8) :: rt_dir  
    real(r8) :: alpha1  ! Alpha for root growth direction
    real(r8) :: alpha2
    real(r8) :: alpha3
    real(r8) :: beta1
    real(r8) :: beta2
    real(r8) :: beta3  ! Coefficients to calculate partition coefficients
    real(r8) :: r_gamma
    real(r8) :: delta  ! Coefficients for sla calculation (m2 kg-1).
    real(r8) :: remob  ! Remobilisation factor
    real(r8) :: cfrac_s  ! Carbon fraction of stems
    real(r8) :: cfrac_r  ! Carbon fraction of roots
    real(r8) :: cfrac_l  ! Carbon fraction of leaves
    real(r8) :: allo1  ! Allometric coefficients for stemc <% canht
    real(r8) :: allo2  ! Allometric coefficients for stemc <% canht
    real(r8) :: mu  ! Coefficient for senescence calculation
    real(r8) :: nu  ! Coefficient for senescence calculation
    real(r8) :: yield_frac  ! Fraction of the harv carbon pool converted to yield carbon
    real(r8) :: initial_carbon  ! Carbon in crops at DVI=initial_c_dvi (kgC/m2).
    real(r8) :: initial_c_dvi  ! DVI at which total crop carbon is set to initial_carbon. Should be emergence (0.0) or shortly after.
    real(r8) :: sen_dvi  ! DVI at which leaf senescence begins. 
    real(r8) :: t_mort  ! Soil temperature (second level) at which to kill crop if dvi>1 (K).
    real(r8) :: croplai_min  ! minimum value for crop LAI
    real(r8) :: cropcanht_min  ! minimum value for crop canopy height
    real(r8) :: croprootc_min  ! minimum value for crop root carbon
    real(r8) :: cropharvc_min  ! minimum value for carbon in harvested crop parts
    real(r8) :: cropreservec_min  ! minimum value for crop stem reserve pool
    real(r8) :: dvi_at_sla_min  ! minimum value for crop SLA
    real(r8) :: cropstemc_min  ! minimum value for crop stem carbon
    real(r8) :: cropleafc_min  ! minimum value for crop leaf carbon 

end type param_vars

!*************************************************phenp/
type,public::Phen

    real(r8) :: npp_ft_acc  ! Accumulated NPP (kg m-2).
    real(r8) :: gpp_ft_acc  ! Accumulated GPP (kg m-2).
    real(r8) :: dvi  ! crop development index
    real(r8) :: rootc  ! root biomass (kg m-2)
    real(r8) :: harvc  ! crop yield (kg m-2)
    real(r8) :: reservec  ! stem reserve pool (kg m-2)
    real(r8) :: nonyield_diag  ! carbon that is leaving the crop model, which is not yield (kg m-2)
    real(r8) :: leafc  ! leaf carbon pool (kg m-2)
    real(r8) :: stemc  ! stem carbon pool (kg m-2)
    real(r8) :: f_root  ! fraction to roots
    real(r8) :: f_stem  ! fraction to stem
    real(r8) :: f_leaf  ! fraction to leaf
    real(r8) :: f_harv  ! fraction to harvested parts
    real(r8) :: extrac  ! carbon used to keep pools above their minima (kg m-2)
    real(r8) :: sen_fac  ! fraction of leafc to move to harvc for senescence
    real(r8) :: lai  ! crop leaf area index
    real(r8) :: sla  ! specific leaf area
    real(r8) :: canht  ! crop canopy height
    real(r8) :: teff  ! effective temperature (K)
    real(r8) :: rpe  ! relative photoperiod effect
    real(r8) :: phot ! photoperiod (hours)
    real(r8) :: temp_daily  ! daily mean temperature (K)

end type Phen   ! crop phenp

public :: Init_crop_parameters,  &       ! initial
          Init_Phen_parameters,  &
          beps_crop_Development, &
          beps_crop_pool_alloc,  &
          no_carbon_pools_below_minimum
          
!private :: carbon_fraction_from_dvi, &
!           calc_cropleafc_min, &
!           calc_cropstemc_min,  &
!           leafc_from_prognostics, &
!           stemc_from_prognostics

contains
! Initialization process
 subroutine Init_crop_parameters(LC,param)     
 implicit none
 integer,intent(in)  :: LC
 type(param_vars)           :: param

 select case (LC)
 case(1001)   ! corn
    param%t_bse= 8+273.15                                                                    
			! Base temp for development (K).
		param%t_opt= 30+273.15                                                                    
			! Optimum temp for development (K).
		param%t_max= 42+273.15                                                                    
			! Maximum temp for development (K).
		param%tt_emr=80+273.15                                                                   
			! Thermal time for emergence (degree days).
		param%crit_pp= 0.                                                                 
			! Critical daylength for photoperiod sensitivity (hours).
		param%pp_sens= 0.                                                                 
		  ! Sensitivity to daylength (hours-1).
		param%rt_dir= 0.                                                                  
		  ! Alpha for root growth direction
		param%alpha1= 13.5                                                                   
		param%alpha2= 12.5                                                                   
		param%alpha3= 13.0                                                                   
		param%beta1= -15.5                                                                    
		param%beta2= -12.5                                                                    
		param%beta3= -14.0                                                                    
		  ! Coefficients to calculate partition coefficients
		param%r_gamma= 22.5                                                                  
		param%delta= -0.2587                                                                    
		  ! Coefficients for sla calculation (m2 kg-1).
		param%remob= 0.35                                                                    
		  ! Remobilisation factor
		param%cfrac_s= 0.5                                                                  
		  ! Carbon fraction of stems
		param%cfrac_r= 0.5                                                                  
		  ! Carbon fraction of roots
		param%cfrac_l= 0.5                                                                  
		  ! Carbon fraction of leaves
		param%allo1= 3.5                                                                    
		param%allo2= 0.4                                                                    
		  ! Allometric coefficients for stemc <% canht
		param%croplai_min= 0.001                                                        
		  ! minimum value for crop LAI
		param%cropcanht_min= 0.001                                                      
		  ! minimum value for crop canopy height
		param%croprootc_min= 0.001                                                      
		  ! minimum value for crop root carbon
		param%cropharvc_min= 0.                                                        
		  ! minimum value for carbon in harvested crop parts
		param%cropreservec_min= 0.
		  ! minimum value for crop stem reserve pool
		param%dvi_at_sla_min= 2 
		  ! minimum value for crop SLA
		
		
		  !Note: parameters below are not used in C language code,here allocate zero to them
	    !      if you want to use them ,please refer T.Osborne paper 2016 and JULES official technical guide online
			
		param%mu= 0.                                                                      
			! Coefficient for senescence calculation
		param%nu= 0.                                                                      
		  ! Coefficient for senescence calculation
		param%yield_frac= 0.                                                              
		  ! Fraction of the harv carbon pool converted to yield carbon
		param%initial_carbon= 0.                                                          
		  ! Carbon in crops at DVI=initial_c_dvi (kgC/m2).
		param%initial_c_dvi= 0.                                                           
		  ! DVI at which total crop carbon is set to initial_carbon. Should be 
		  ! emergence (0.0) or shortly after.
		param%sen_dvi= 0.                                                                 
		  ! DVI at which leaf senescence begins. 
		param%t_mort= 0.
		  ! Soil temperature (second level) at which to kill crop if dvi>1 (K).
 case(1002)    !soybean
    	param%t_bse= 5+273.15                                                                    
			  ! Base temp for development (K).
			param%t_opt= 27+273.15                                                                    
			  ! Optimum temp for development (K).
			param%t_max= 40+273.15                                                                    
			  ! Maximum temp for development (K).
			param%tt_emr= 35                                                                   
			  ! Thermal time for emergence (degree days).
			param%crit_pp= 0.                                                                 
			  ! Critical daylength for photoperiod sensitivity (hours).
			param%pp_sens= 0.                                                                 
			  ! Sensitivity to daylength (hours-1).
			param%rt_dir= 0.                                                                  
			  ! Alpha for root growth direction
			param%alpha1= 20                                                                   
			param%alpha2= 18.5                                                                   
			param%alpha3= 19.5                                                                   
			param%beta1= -16.5                                                                    
			param%beta2= -14.5                                                                    
			param%beta3= -15.0                                                                    
			  ! Coefficients to calculate partition coefficients
			param%r_gamma= 25.9                                                                  
			param%delta= -0.1451                                                                    
			  ! Coefficients for sla calculation (m2 kg-1).
			param%remob= 0.18                                                                    
			  ! Remobilisation factor
			param%cfrac_s= 0.5                                                                  
			  ! Carbon fraction of stems
			param%cfrac_r= 0.5                                                                  
			  ! Carbon fraction of roots
			param%cfrac_l= 0.5                                                                  
			  ! Carbon fraction of leaves
			param%allo1= 1.6                                                                    
			param%allo2= 0.4                                                                    
			  ! Allometric coefficients for stemc <% canht
			param%croplai_min= 0.001                                                        
			  ! minimum value for crop LAI
			param%cropcanht_min= 0.001                                                      
			  ! minimum value for crop canopy height
			param%croprootc_min= 0.001                                                      
			  ! minimum value for crop root carbon
			param%cropharvc_min= 0.                                                        
			  ! minimum value for carbon in harvested crop parts
			param%cropreservec_min= 0.
			  ! minimum value for crop stem reserve pool
			param%dvi_at_sla_min= 2 
			  ! minimum value for crop SLA
		
		    !Note: parameters below are not used in C language code,here allocate zero to them
			!      if you want to use them ,please refer T.Osborne paper 2016 and JULES official technical guide online
			
			param%mu= 0.                                                                      
			  ! Coefficient for senescence calculation
			param%nu= 0.                                                                      
			  ! Coefficient for senescence calculation
			param%yield_frac= 0.                                                              
			  ! Fraction of the harv carbon pool converted to yield carbon
			param%initial_carbon= 0.                                                          
			  ! Carbon in crops at DVI=initial_c_dvi (kgC/m2).
			param%initial_c_dvi= 0.                                                           
			  ! DVI at which total crop carbon is set to initial_carbon. Should be 
			  ! emergence (0.0) or shortly after.
			param%sen_dvi= 0.                                                                 
			  ! DVI at which leaf senescence begins. 
			param%t_mort= 0.
			  ! Soil temperature (second level) at which to kill crop if dvi>1 (K).
 case(1003)   ! spring wheat
    	param%t_bse= 0+273.15                                                                    
			  ! Base temp for development (K).
			param%t_opt= 20+273.15                                                                    
			  ! Optimum temp for development (K).
			param%t_max= 30+273.15                                                                    
			  ! Maximum temp for development (K).
			param%tt_emr= 35+273.15                                                                   
			  ! Thermal time for emergence (degree days).
			param%crit_pp= 0.                                                                 
			  ! Critical daylength for photoperiod sensitivity (hours).
			param%pp_sens= 0.                                                                 
			  ! Sensitivity to daylength (hours-1).
			param%rt_dir= 0.                                                                  
			  ! Alpha for root growth direction
			param%alpha1= 18.5                                                                   
			param%alpha2= 16.0                                                                   
			param%alpha3= 18.0                                                                   
			param%beta1= -20.0                                                                    
			param%beta2= -15.0                                                                    
			param%beta3= -18.5                                                                    
			  ! Coefficients to calculate partition coefficients
			param%r_gamma= 27.3                                                                  
			param%delta= -0.0507                                                                    
			  ! Coefficients for sla calculation (m2 kg-1).
			param%remob= 0.40                                                                    
			  ! Remobilisation factor
			param%cfrac_s= 0.5                                                                  
			  ! Carbon fraction of stems
			param%cfrac_r= 0.5                                                                  
			  ! Carbon fraction of roots
			param%cfrac_l= 0.5                                                                  
			  ! Carbon fraction of leaves
			param%allo1= 1.4                                                                    
			param%allo2= 0.4                                                                    
			  ! Allometric coefficients for stemc <% canht
			param%croplai_min= 0.001                                                        
			  ! minimum value for crop LAI
			param%cropcanht_min= 0.001                                                      
			  ! minimum value for crop canopy height
			param%croprootc_min= 0.001                                                      
			  ! minimum value for crop root carbon
			param%cropharvc_min= 0.                                                        
			  ! minimum value for carbon in harvested crop parts
			param%cropreservec_min= 0.
			  ! minimum value for crop stem reserve pool
			param%dvi_at_sla_min= 2 
			  ! minimum value for crop SLA
		
		    !Note: parameters below are not used in C language code,here allocate zero to them
			!      if you want to use them ,please refer T.Osborne paper 2016 and JULES official technical guide online
			
			param%mu= 0.                                                                      
			  ! Coefficient for senescence calculation
			param%nu= 0.                                                                      
			  ! Coefficient for senescence calculation
			param%yield_frac= 0.                                                              
			  ! Fraction of the harv carbon pool converted to yield carbon
			param%initial_carbon= 0.                                                          
			  ! Carbon in crops at DVI=initial_c_dvi (kgC/m2).
			param%initial_c_dvi= 0.                                                           
			  ! DVI at which total crop carbon is set to initial_carbon. Should be 
			  ! emergence (0.0) or shortly after.
			param%sen_dvi= 0.                                                                 
			  ! DVI at which leaf senescence begins. 
			param%t_mort= 0.
			  ! Soil temperature (second level) at which to kill crop if dvi>1 (K).
 case(1004)    !rice
      param%t_bse= 8+273.15                                                                    
			  ! Base temp for development (K).
			param%t_opt= 30+273.15                                                                    
			  ! Optimum temp for development (K).
			param%t_max= 42+273.15                                                                    
			  ! Maximum temp for development (K).
			param%tt_emr= 60+273.15                                                                   
			  ! Thermal time for emergence (degree days).
			param%crit_pp= 0.                                                                 
			  ! Critical daylength for photoperiod sensitivity (hours).
			param%pp_sens= 0.                                                                 
			  ! Sensitivity to daylength (hours-1).
			param%rt_dir= 0.                                                                  
			  ! Alpha for root growth direction
			param%alpha1= 18.5                                                                   
			param%alpha2= 19.0                                                                   
			param%alpha3= 19.5                                                                   
			param%beta1= -19.0                                                                    
			param%beta2= -17.0                                                                    
			param%beta3= -18.5                                                                    
			  ! Coefficients to calculate partition coefficients
			param%r_gamma= 20.9                                                                  
			param%delta= -0.2724                                                                    
			  ! Coefficients for sla calculation (m2 kg-1).
			param%remob= 0.25                                                                    
			  ! Remobilisation factor
			param%cfrac_s= 0.5                                                                  
			  ! Carbon fraction of stems
			param%cfrac_r= 0.5                                                                  
			  ! Carbon fraction of roots
			param%cfrac_l= 0.5                                                                  
			  ! Carbon fraction of leaves
			param%allo1= 1.4                                                                    
			param%allo2= 0.4                                                                    
			  ! Allometric coefficients for stemc <% canht
			param%croplai_min= 0.001                                                        
			  ! minimum value for crop LAI
			param%cropcanht_min= 0.001                                                      
			  ! minimum value for crop canopy height
			param%croprootc_min= 0.001                                                      
			  ! minimum value for crop root carbon
			param%cropharvc_min= 0.                                                        
			  ! minimum value for carbon in harvested crop parts
			param%cropreservec_min= 0.
			  ! minimum value for crop stem reserve pool
			param%dvi_at_sla_min= 2 
			  ! minimum value for crop SLA
		
		
		    !Note: parameters below are not used in C language code,here allocate zero to them
			!      if you want to use them ,please refer T.Osborne paper 2016 and JULES official technical guide online
			
			param%mu= 0.                                                                      
			  ! Coefficient for senescence calculation
			param%nu= 0.                                                                      
			  ! Coefficient for senescence calculation
			param%yield_frac= 0.                                                              
			  ! Fraction of the harv carbon pool converted to yield carbon
			param%initial_carbon= 0.                                                          
			  ! Carbon in crops at DVI=initial_c_dvi (kgC/m2).
			param%initial_c_dvi= 0.                                                           
			  ! DVI at which total crop carbon is set to initial_carbon. Should be 
			  ! emergence (0.0) or shortly after.
			param%sen_dvi= 0.                                                                 
			  !! DVI at which leaf senescence begins. 
			param%t_mort= 0.
			  !! Soil temperature (second level) at which to kill crop if dvi>1 (K).
 case default
      param%t_bse= 0.                                                                   
	    !! Base temp for development (K).
	    param%t_opt= 0.                                                                   
	    !! Optimum temp for development (K).
	    param%t_max= 0.                                                                   
	    !! Maximum temp for development (K).
	    param%tt_emr= 0.                                                                  
	    !! Thermal time for emergence (degree days).
	    param%crit_pp= 0.                                                                 
	    !! Critical daylength for photoperiod sensitivity (hours).
	    param%pp_sens= 0.                                                                 
	    !! Sensitivity to daylength (hours-1).
	    param%rt_dir= 0.                                                                  
	    !! Alpha for root growth direction
	    param%alpha1= 0.                                                                  
	    param%alpha2= 0.                                                                  
	    param%alpha3= 0.                                                                  
	    param%beta1= 0.                                                                   
	    param%beta2= 0.                                                                   
	    param%beta3= 0.                                                                   
	    !! Coefficients to calculate partition coefficients
	    param% r_gamma= 0.                                                                 
	    param%delta= 0.                                                                   
	    !! Coefficients for sla calculation (m2 kg-1).
	    param%remob= 0.                                                                   
	    !! Remobilisation factor
	    param%cfrac_s= 0.                                                                 
	    !! Carbon fraction of stems
	    param%cfrac_r= 0.                                                                 
	    !! Carbon fraction of roots
	    param%cfrac_l= 0.                                                                 
	    !! Carbon fraction of leaves
	    param%allo1= 0.                                                                   
	    param%allo2= 0.                                                                   
	    !! Allometric coefficients for stemc <% canht
	    param%mu= 0.                                                                      
	    !! Coefficient for senescence calculation
	    param%nu= 0.                                                                      
	    !! Coefficient for senescence calculation
	    param%yield_frac= 0.                                                              
	    !! Fraction of the harv carbon pool converted to yield carbon
	    param%initial_carbon= 0.                                                          
	    !! Carbon in crops at DVI=initial_c_dvi (kgC/m2).
	    param%initial_c_dvi= 0.                                                           
	    !! DVI at which total crop carbon is set to initial_carbon. Should be 
	    !! emergence (0.0) or shortly after.
	    param%sen_dvi= 0.                                                                 
	    !! DVI at which leaf senescence begins. 
	    param%t_mort= 0.
	    ! Soil temperature (second level) at which to kill crop if dvi>1 (K).
	    param%croplai_min= 0.                                                       
	    ! minimum value for crop LAI
	    param%cropcanht_min= 0.                                                     
	    ! minimum value for crop canopy height
	    param%croprootc_min= 0.                                                     
	    ! minimum value for crop root carbon
	    param%cropharvc_min= 0.                                                        
	    ! minimum value for carbon in harvested crop parts
	    param%cropreservec_min= 0.
	    ! minimum value for crop stem reserve pool
	    param%dvi_at_sla_min= 0.
	    ! minimum value for crop SLA
	  
	    param%cropstemc_min= 0.  ! minimum value for crop stem carbon
	    param%cropleafc_min= 0.  ! minimum value for crop leaf carbon  
 end select
 return
 end subroutine Init_crop_parameters


 subroutine Init_Phen_parameters(phenp)     
 		implicit none
 		type(Phen)         :: phenp

    phenp%npp_ft_acc = 0.   	 ! Accumulated NPP (kg m-2).
		phenp%gpp_ft_acc = 0.   	 ! Accumulated GPP (kg m-2).
		phenp% dvi = 0.            ! crop development index
		phenp% rootc = 0.          ! root biomass (kg m-2)
		phenp% harvc = 0.          ! crop yield (kg m-2)
		phenp% reservec = 0.       ! stem reserve pool (kg m-2)
		phenp% nonyield_diag = 0.  ! carbon that is leaving the crop model, which is not yield (kg m-2)
	  

		phenp% leafc = 0.        ! leaf carbon pool (kg m-2)
		phenp% stemc = 0.        ! stem carbon pool (kg m-2)
	
		phenp% f_root = 0.         ! fraction to roots
		phenp% f_stem = 0.         ! fraction to stem
		phenp% f_leaf = 0.         ! fraction to leaf
		phenp% f_harv = 0.         ! fraction to harvested parts
		phenp% extrac = 0.         ! carbon used to keep pools above their minima (kg m-2)
		phenp% sen_fac = 0.        ! fraction of leafc to move to harvc for senescence
	
		phenp%lai = 0.   ! crop leaf area index
    phenp%sla = 0.   ! specific leaf area
    phenp%canht = 0. ! crop canopy height
	
		phenp%teff = 0.  ! effective temperature (K)
		phenp%rpe  = 0.  ! relative photoperiod effect
		phenp%phot = 0.  ! photoperiod (hours)
	
		phenp%temp_daily=0.  ! daily mean temperature (K)
 return
 end subroutine Init_Phen_parameters

 subroutine beps_crop_Development(i,j,param,phenp)
	
	!-----------------------------------------------------------------------------
	! Description:
	!   This subroutine increments the plant development index using the thermal
	!   time accumulation.
	!   Increases dvi by a rate determined by thermal time accumulation. dvi
	!   increases from 0 to 2.
	!
	! Method:
	!   Crop should have already emerged when this subroutine is called.
	!   An effective temperature (Teff) is calculated from the tile temperature
	!   and the increase in crop development index (DVI) is then found.
	!   During the first stage of development, the rate of increase in DVI
	!   is slowed (for some plants) by multiplying by RPE (the relative
	!   photoperiod effect).
	!
	!   See JULES-crop technical documentation (Tom Osborne, Josh Hooker).
	!
	! Code Owner: Please refer to ModuleLeaders.txt
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
    !
	! Code update:xiaorongwang rewrite to C language, 2021/06
	! 0:emerge  1:flowering  2:maturity
	!-----------------------------------------------------------------------------	

  implicit none
	integer, INTENT(IN)    :: i,j
  type(param_vars)  :: param
  type(Phen)  :: phenp
  

  if ( phenp%temp_daily <= param%t_opt .AND. phenp%temp_daily >= param%t_bse ) then

		phenp%teff = phenp%temp_daily - param%t_bse

	else if ( phenp%temp_daily > param%t_opt .AND. phenp%temp_daily < param%t_max ) then

		phenp%teff = param%t_opt - param%t_bse - ( ( phenp%temp_daily - param%t_opt ) &
		             / ( param%t_max - param%t_opt ) ) * ( param%t_opt - param%t_bse )                    
	 
	else if ( phenp%temp_daily < param%t_bse .OR. phenp%temp_daily >= param%t_max ) then

		phenp%teff = 0.0
  end if
	
	if ( pgdd%phot_type(i,j)==1.0 .AND. phenp%dvi < 1.0 ) then  !! Photo-period effect on development only during first phase 
		
			phenp%rpe = 1.0 - ( phenp%phot - param%crit_pp ) * param%pp_sens 
			phenp%rpe = min( 1.0, phenp%rpe ) 
			phenp%rpe = max( 0.0, phenp%rpe ) 
			phenp%dvi = phenp%dvi + phenp%rpe * phenp%teff / pgdd%tt_veg(i,j) 
	
	else  !! No photoperiod effect
		
			phenp%dvi = phenp%dvi + phenp%teff / pgdd%tt_rep(i,j) 
		
  end if

 return
 end subroutine beps_crop_Development

! !------------------------------------------------------------------------------------------------
!!                                 Function calculate photoperiod
!!------------------------------------------------------------------------------------------------  
!void photoperiod()
!{
!			!pgdd%phot_type  yes=1 no=0  read from input file
!		    !Note: photoperiod are not considered in C language code
!			!      if you want to use it ,please refer T.Osborne paper 2016 and JULES official technical guide online
	 
!}


!------------------------------------------------------------------------------------------------
!                                 Do the daily Partition coefficients and carbon allocation 
!------------------------------------------------------------------------------------------------

 subroutine beps_crop_pool_alloc(param,phenp)
	!-----------------------------------------------------------------------------
	! Description:
	!   Partitions crop tile npp to carbon stores
	!
	! Method:
	!   Updates crop carbon pools using crop DVI and NPP.
	!   Checks crops is established.
	!   Includes leaf senescence.
	!   Doesn't allow carbon pools to drop below zero.
	!
	!   See JULES-crop technical documentation (Tom Osborne, Josh Hooker).
	!
	! Code Owner: Please refer to ModuleLeaders.txt
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
	!
	! Code update:xiaorongwang rewrite to C language,2021/06
	!-----------------------------------------------------------------------------	
  implicit none
  type(param_vars)  :: param
  type(Phen)  :: phenp

  ! Calculate Partition coefficients
	call carbon_fraction_from_dvi(param,phenp) 


	
	!-----------------------------------------------------------------------------
	! Update carbon pools according to DVI
	! Partition according to crop type and development stage
	! units are kgC/m^2
	! NB. Reserves taken from STEM growth rate
	!-----------------------------------------------------------------------------
   
	phenp%rootc    = phenp%rootc + phenp%npp_ft_acc * phenp%f_root 

	phenp%leafc    = phenp%leafc + phenp%npp_ft_acc * phenp%f_leaf 

	phenp%stemc    = phenp%stemc + (phenp%npp_ft_acc * phenp%f_stem * ( 1.0 - param%remob ) ) 

	phenp%harvc    = phenp%harvc + phenp%npp_ft_acc * phenp%f_harv 

	phenp%reservec = phenp%reservec + ( phenp%npp_ft_acc * phenp%f_stem * param%remob ) 

	
	
	
	!-----------------------------------------------------------------------------
	! Re-allocation e.g. from reserves to grain once
	! stems have stopped growing
	! NB. 0.1 (0.9) is hard-wired
	!-----------------------------------------------------------------------------
	
	
	if ( phenp%f_stem < 0.01 ) then
		
			phenp%harvc    = phenp%harvc + ( 0.1 * phenp%reservec ) 
			phenp%reservec = phenp%reservec * 0.9 

	end if

	!-----------------------------------------------------------------------------
	! Leaf senescence
	!-----------------------------------------------------------------------------
	
	if ( phenp%dvi > param%sen_dvi ) then
		
			phenp%sen_fac = param%mu * ( phenp%dvi - param%sen_dvi )** param%nu 
			phenp%sen_fac = min(phenp%sen_fac, 1.0) 

			phenp%harvc   = phenp%harvc + (phenp%sen_fac * phenp%leafc) 
			phenp%leafc   = phenp%leafc * (1.0 - phenp%sen_fac) 
	end if
	
	if ( phenp%dvi > 1.5 ) then !here use the value of sen_dvi in T.Osborne 2016 directly
		
			phenp%harvc   = phenp%harvc + (0.05 * phenp%leafc) 
			phenp%leafc   = phenp%leafc * (1.0 - 0.05) 
		
	end if

	
	!-----------------------------------------------------------------------------
	! don't allow carbon pools to drop below their minimum values
	! ----------------------------------------------------------------------------
	
	
	!no_carbon_pools_below_minimum(param,phenp) 
	!phenp%nonyield_diag = phenp%nonyield_diag - phenp%extrac 

 return
 end subroutine beps_crop_pool_alloc


 subroutine carbon_fraction_from_dvi(param,phenp)
	!-----------------------------------------------------------------------------
	! Description:
	!   Calculates carbon fraction to root, stem, leaf and harvested part of crop
	!   using the crop development index
	!
	! Method:
	!   See JULES-crop technical documentation (Tom Osborne, Josh Hooker).
	!
	! Code Owner: Please refer to ModuleLeaders.txt
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
	!
	! Code update:xiaorongwang rewrite to C language, 2021/06
	!-----------------------------------------------------------------------------
  implicit none
  type(param_vars)  :: param
  type(Phen)  :: phenp


  !! Local variables
	REAL(r8) :: s1
  REAL(r8) :: s2
  REAL(r8) :: s3
  REAL(r8) :: denom
	!!-----------------------------------------------------------------------------
	
	
	s1 = exp( param%alpha1 + ( param%beta1 * phenp%dvi ) )  !root
	s2 = exp( param%alpha2 + ( param%beta2 * phenp%dvi ) )  !stem
	s3 = exp( param%alpha3 + ( param%beta3 * phenp%dvi ) )  !leaf

	denom = s1 + s2 + s3 + 1.0 

	phenp%f_root = s1 / denom 
	phenp%f_stem = s2 / denom 
	phenp%f_leaf = s3 / denom 
	phenp%f_harv =  1.0 / denom 

 return
 end subroutine carbon_fraction_from_dvi

 subroutine no_carbon_pools_below_minimum(param,phenp)

	!-----------------------------------------------------------------------------
	! Description:
	!   Makes sure the carbon pools do not drop below their lower limit
	!
	! Code Owner: Please refer to ModuleLeaders.txt
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
	! Code update:xiaorongwang rewrite to C language, 2021/06
	!-----------------------------------------------------------------------------
	
  implicit none
  type(param_vars)  :: param
  type(Phen)  :: phenp

  phenp%extrac = 0.0 

	if ( phenp%rootc < param%croprootc_min ) then
		
			phenp%extrac = phenp%extrac + param%croprootc_min - phenp%rootc 
			phenp%rootc  = param%croprootc_min 
	end if

	if ( phenp%harvc < param%cropharvc_min ) then
		
			phenp%extrac = phenp%extrac + param%cropharvc_min - phenp%harvc 
			phenp%harvc  = param%cropharvc_min 
	end if

	if ( phenp%reservec < param%cropreservec_min ) then
		
			phenp%extrac   = phenp%extrac + param%cropreservec_min - phenp%reservec 
			phenp%reservec = param%cropreservec_min 
	end if

	call calc_cropleafc_min(param,phenp)

	if ( phenp%leafc < param%cropleafc_min ) then
		
			phenp%extrac = phenp%extrac + param%cropleafc_min - phenp%leafc 
			phenp%leafc  = param%cropleafc_min 
	end if
	
	call calc_cropstemc_min(param,phenp) 

	if ( phenp%stemc < param%cropstemc_min ) then
		
			phenp%extrac = phenp%extrac + param%cropstemc_min - phenp%stemc 
			phenp%stemc  = param%cropstemc_min 
	end if 

 return
 end subroutine no_carbon_pools_below_minimum

 subroutine calc_cropleafc_min(param,phenp)

	!-----------------------------------------------------------------------------
	! Description:
	!   Calculates lower limit for leaf carbon using minimum lai i.e. ensure that
	!   lai never goes below croplai_min, whatever the value of DVI is.
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
	!
	! Code update:xiaorongwang rewrite to C language, 2021/06
	!-----------------------------------------------------------------------------

  implicit none
  type(param_vars)  :: param
  type(Phen)  :: phenp

  !param%dvi_at_sla_min = 2.0   
	call leafc_from_prognostics(param,phenp)
	param%cropleafc_min =  phenp%leafc

 return
 end subroutine calc_cropleafc_min

 subroutine calc_cropstemc_min(param,phenp)
 
	!-----------------------------------------------------------------------------
	! Description:
	!   Calculates lower limit for stem carbon using minimum canopy height.
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
	!-----------------------------------------------------------------------------
	
  implicit none
  type(param_vars)  :: param
  type(Phen)  :: phenp

  call stemc_from_prognostics(param,phenp) 
	param%cropstemc_min = phenp%stemc

 return
 end subroutine calc_cropstemc_min

 subroutine leafc_from_prognostics(param,phenp)
 
	!-----------------------------------------------------------------------------
	! Description:
	!   Calculates leaf carbon from the DVI (crop development index)
	!   and LAI (leaf area index), which are prognostics, and some crop-specific
	!   allometric variables.
	!
	! Method:
	!   See JULES-crop technical documentation (Tom Osborne, Josh Hooker).
	!
	! Code Owner: Please refer to ModuleLeaders.txt
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
	!
	! Code update:xiaorongwang rewrite to C language, 2021/06
	!-----------------------------------------------------------------------------
	
  implicit none
  type(param_vars)  :: param
  type(Phen)  :: phenp

  if ( phenp%dvi >= 0.0 ) then
		
			phenp%sla = param%r_gamma * ( phenp%dvi + 0.06 )**param%delta 
			phenp%leafc = ( phenp%lai / phenp%sla ) * param%cfrac_l 
	end if
	!else 
			!phenp%leafc = calc_cropleafc_min(n) 
		!}

 return
 end subroutine leafc_from_prognostics


 subroutine stemc_from_prognostics(param,phenp)

	!-----------------------------------------------------------------------------
	! Description:
	!   Calculates stem carbon from the canopy height, which is a prognostic,
	!   and some crop-specific allometric variables.
	!
	! Method:
	!   See JULES-crop technical documentation (Tom Osborne, Josh Hooker).
	!
	! Code Owner: Please refer to ModuleLeaders.txt
	!
	! Code Description:
	!   Language: Fortran 90.
	!   This code is written to JULES coding standards v1.
	!
	! Code update:xiaorongwang rewrite to C language, 2021/06
	!-----------------------------------------------------------------------------
  implicit none
  type(param_vars)  :: param
  type(Phen)  :: phenp

  phenp%stemc = param%cfrac_s * ( ( phenp%canht / param%allo1 )**(1.0 / param%allo2)) 

 return
 end subroutine stemc_from_prognostics

end module beps_cropMod
