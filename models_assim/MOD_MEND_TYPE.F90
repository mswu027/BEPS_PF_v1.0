MODULE MOD_MEND_TYPE    
! File:   STRUCT_MEND.F90
! Author: GANGSHENG WANG @ ORNL
! Updated: July 15, 2015
! Created on February 26, 2013, 9:30 AM
    use shr_kind_mod, only: r8=>shr_kind_r8
     ! ``STRUCTURE /name/ ... END STRUCTURE'' becomes
     ! ``TYPE name ... END TYPE''
     ! ``RECORD /name/ variable'' becomes ``TYPE(name) variable''
    INTEGER, PARAMETER:: const_nVARc = 13                           !# of state variables
    INTEGER, PARAMETER:: const_nFLXc = 27                           !# of flux variables
    
    INTEGER, PARAMETER:: const_nPOC = 2                             !# of particular organic carbon pools
    INTEGER, PARAMETER:: const_nISO = 2                             !# of isotopes, e.g., C12, C14, C13
    
    real(r8), PARAMETER:: const_R = 8.314                          ![J/mol/K],universal gas constant
    real(r8), PARAMETER:: const_tmp_C2K = 273.15                  !conversion of degree C to degree K
!    real(8), PARAMETER:: const_Rstd(const_nISO - 1) = (/1d-12/)     !C14/C12
    real(r8), PARAMETER:: const_Rstd(const_nISO - 1) = (/0.0112372/)   !C13/C12, Pee Dee Belemnite (PDB) Standard
    real(r8), PARAMETER:: const_cm2MPa = 98e-6                       !1cm water column 
    
    real(r8), PARAMETER:: const_Tref = 20.                         ![degree c],reference temperature
    real(r8), PARAMETER:: const_SWPmin = -3.0e4                        ![MPa], lowest SWP
    !real(r8), PARAMETER:: const_FillValue = -999                   !filled value
    
    !character, PARAMETER :: cBackspace = char(13)
      
!-----------------------------------------------------------------------------!   
!Adsorption and Desorption of DOC
    TYPE sSORP_PAR
        REAL(r8) Qmax                        ![mg C/g soil], adsorption capacity
        REAL(r8) Kads                        ![mg C/mg C/h], specific adsorption rate
        REAL(r8) Kdes                        ![mg C/g soil/h], desorption rate
    END TYPE sSORP_PAR
    
    TYPE sSORP_INP
        REAL(r8) adsorbate                   !e.g., DOC
        REAL(r8) adsorbent                   !e.g., QOC (MOC)
    END TYPE sSORP_INP
    
    TYPE sSORP_OUT
        REAL(r8) ads                         ![mg C/g soil/h], adsorption flux
        REAL(r8) des                         ![mg C/g soil/h], desorption flux
        REAL(r8) ads_net                     ![mg C/g soil/h], net adsorption flux = ads - des
    END TYPE sSORP_OUT
!-----------------------------------------------------------------------------!   
!Michaelis-Menten Kinetics
     TYPE sMM_PAR
         REAL(r8) vm                         ![mg C/mg C/h], specific enzyme activity
         REAL(r8) km                         ![mg C/g soil], half-saturation constant
     END TYPE sMM_PAR
     
     TYPE sMM_INP
         REAL(r8) substrate                  ![mg C/g soil], substrate concentration
         REAL(r8) enzyme                     ![mg C/g soil], enzyme concentration
     END TYPE sMM_INP

!=============================================================================!
!Microbial-Enzyme-mediated Nitrification-Denitrification-Decomposition (MEND)
!-----------------------------------------------------------------------------!              
     TYPE sMEND_CPOOL
         REAL(r8) SOC                       ![mg C/g soil], total SOC = POC + MOC + QOC
         REAL(r8) POC(const_nPOC)           ![mg C/g soil],Particulate Organic Carbon, size of {} determined by nPOC
         REAL(r8) MOC                       ![mg C/g soil],Mineral Associate Organic Carbon
         REAL(r8) QOC                       ![mg C/g soil],adsorbed phase of DOC
         REAL(r8) DOC                       ![mg C/g soil],Dissolved Organic Carbon
         REAL(r8) MBC                       ![mg C/g soil],Microbial Biomass Carbon
         REAL(r8) MBCA                      ![mg C/g soil],Active Microbial Biomass Carbon
         REAL(r8) MBCD                      ![mg C/g soil],Dormant Microbial Biomass Carbon
         REAL(r8) ENZP(const_nPOC)          ![mg C/g soil],ENZyme for POC
         REAL(r8) ENZM                      ![mg C/g soil],Enzyme for MAOC     
         REAL(r8) CO2                       ![mg C/g soil],cumulative CO2 in the soil
     END TYPE sMEND_CPOOL
!-----------------------------------------------------------------------------!     
     TYPE sMEND_CFLUX
         REAL(r8) POCdec(const_nPOC)        ![mg C/g soil/h],decomposition of POC
         REAL(r8) POCdec_to_DOC(const_nPOC) ![mg C/g soil/h],decomposition of POC allocated to DOC
         REAL(r8) POCdec_to_MOC(const_nPOC) ![mg C/g soil/h],decomposition of POC allocated to MOC
         REAL(r8) MOC_to_DOC                ![mg C/g soil/h],decomposition of MOC
         REAL(r8) QOC_to_DOC                ![mg C/g soil/h], desorption flux
         REAL(r8) DOC_to_QOC                ![mg C/g soil/h], adsorption flux
         REAL(r8) DOC_to_QOC_net            ![mg C/g soil/h], net adsorption flux = DOC_to_QOC - QOC_to_DOC
         REAL(r8) DOC_to_MBC                ![mg C/g soil/h], uptake of DOC by MBC
         REAL(r8) CO2_growth                ![mg C/g soil/h], growth respiration rate
         REAL(r8) CO2_maintn                ![mg C/g soil/h], maintenance respiration rate
         REAL(r8) CO2_maintn_dorm           ![mg C/g soil/h], maintenance respiration rate of dormant microbes
         REAL(r8) CO2_gm                    ![mg C/g soil/h], growth + maintenance respiration rate
         REAL(r8) MBC_mortality             ![mg C/g soil/h], microbial mortality
         REAL(r8) MBC_to_ENZP(const_nPOC)   ![mg C/g soil/h], enzyme-POC production
         REAL(r8) MBC_to_ENZM               ![mg C/g soil/h], enzyme-MOC production
         REAL(r8) MBC_PM                    ![mg C/g soil/h], microbial mortality + enzyme production
         REAL(r8) ENZP_to_DOC(const_nPOC)   ![mg C/g soil/h], enzyme-POC turnover
         REAL(r8) ENZM_to_DOC               ![mg C/g soil/h], enzyme-MOC turnover         
         REAL(r8) MBC_to_DOC                ![mg C/g soil/h], turnover of MBC allocated to DOC
         REAL(r8) MBC_to_POC                ![mg C/g soil/h], turnover of MBC allocated to POC (lignin)
         REAL(r8) MBCA_to_MBCD              ![mg C/g soil/h], dormancy
         REAL(r8) MBCD_to_MBCA              ![mg C/g soil/h], reactivation
     END TYPE sMEND_CFLUX
 !-----------------------------------------------------------------------------!      
     TYPE sMEND_CADD
         REAL(r8) POCadd(const_nPOC)        ![mg POC/g soil/h],inputs to POC
         REAL(r8) DOCadd                    ![mg DOC/g soil/h], inputs to DOC
     END TYPE sMEND_CADD
 !-----------------------------------------------------------------------------!                       
     TYPE sMEND_PAR
!         INTEGER nPar
!         INTEGER iKinetics                 !decomposition kinetics: 0-Michaelis-Menten, 1-First Order, 2-Second Order
         REAL(r8) VdPOC(const_nPOC)         ![mg POC/mg ENZP/h], maximum reaction rate for conversion of POC by ENZP
         REAL(r8) VdMOC                     ![mg MOC/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC
         REAL(r8) KsPOC(const_nPOC)         ![mg POC/cm3], half-saturation constant for conversion of POC by ENZP
         REAL(r8) KsMOC                     ![mg MOC/cm3], half-saturation constant for conversion of MAOC by ENZMAOC
         REAL(r8) Qmax                      ![mg C/g soil], adsorption capacity
         REAL(r8) Kba                       ![mg C/g soil/h], binding affinity
         REAL(r8) Kdes                      ![mg DOC/h],desorption rate constant
         REAL(r8) Kads                      ![mg DOC/mg DOC/h], adsorption rate constant = Kba*Kdes
         REAL(r8) rENZP(const_nPOC)         ![1/h],turnover rate of ENZP
         REAL(r8) rENZM                     ![1/h], turnover rate of ENZMAOC
         REAL(r8) pENZP                     ![mg ENZP/mg MBC/h], production rate of ENZP
         REAL(r8) pENZM                     ![mg ENZM/mg MBC/h], production rate of ENZMAOC
         REAL(r8) frPOC2DOC                 ![-], fraction of decomposed POC allocated to DOC
         REAL(r8) frMBC2DOC                 ![-], fraction of dead MBC allocated to DOC
         REAL(r8) Vg                        ![mg DOC/mg MBC/h], maximum uptake rate of DOC by MB
         REAL(r8) alpha                     ![-], = Vm/(Vg+Vm)
         REAL(r8) Vm                        ![1/h], specific microbial maintenance rate
         REAL(r8) KsDOC                     ![mg DOC/cm3], half-saturation constant for uptake of DOC by MB
         REAL(r8) Yg                        ![-], carbon use efficiency in uptake of DOC by MB
         REAL(r8) wdie                      ![-], exponential in SWP scalar (fSWP0) for microbial mortality
         REAL(r8) gamma                     ![-], rMORT = gama*Vm
         REAL(r8) rMORT                     ![1/h], microbial mortality
         REAL(r8) beta                      ![-], VmD = beta*Vm = dormant maintenance rate
         REAL(r8) VmD                       ![1/h], VmD = beta*Vm, specific microbial maintenance rate for dormant microbes
         REAL(r8) VmA2D                     ![1/h], dormancy rate
         REAL(r8) VmD2A                     ![1/h], resuscitation rate
         REAL(r8) SWP_A2D                   ![MPa], negative value, threshold SWP for microbial dormancy (e.g., -0.4 MPa)
         REAL(r8) tau                       ![-], SWP_D2A = tau*SWP_A2D
         REAL(r8) SWP_D2A                   ![MPa], negative value, threshold SWP for microbial resuscitation (e.g., 1/4*SWP_A2D)
         REAL(r8) wdorm                     ![-], exponential in SWP scalar for A2D and D2A         
      END TYPE sMEND_PAR
!-----------------------------------------------------------------------------!       
     TYPE sMEND_INP
         REAL(r8) SIN                                    ![mgC/cm3/h],SOC input, e.g., litter
         REAL(r8) tmp                                    ![C],temperature
         REAL(r8) pH                                     ![-],pH value
         REAL(r8) SWP                                    ![MPa],soil water potential [-], water-filled pore space
         REAL(r8) dt                                     ![h], time-step, dt = 0.5 h for CLM 
!         INTEGER iKinetics                              !decomposition kinetics: 0-Michaelis-Menten, 1-First Order, 2-Second Order
         TYPE(sMEND_CADD) CADD                          !external carbon input
         TYPE(sMEND_CADD) CADDI(const_nISO)             !external carbon input (isotopes)
         TYPE(sMEND_CPOOL) CPOOL                        !carbon pools
         TYPE(sMEND_CPOOL) CPOOLI(const_nISO)           !isotopic carbon pools
         TYPE(sMEND_CPOOL) CPOOLIFR(const_nISO)         !isotopic fraction in carbon pools
         TYPE(sMEND_CPOOL) CPOOLI_SIG(const_nISO - 1)   ![‰],isotopic signature in carbon pools
     END TYPE sMEND_INP
!-----------------------------------------------------------------------------!             
      TYPE sMEND_OUT
         TYPE(sMEND_CPOOL) CPOOL                        !carbon pools
         TYPE(sMEND_CPOOL) CPOOLI(const_nISO)           !isotopic carbon pools
         TYPE(sMEND_CPOOL) CPOOLIFR(const_nISO)         !isotopic fraction in carbon pools
         TYPE(sMEND_CPOOL) CPOOLI_SIG(const_nISO - 1)   ![‰],isotopic signature in carbon pools
         TYPE(sMEND_CFLUX) CFLUX                        !CARBON FLUXES
         REAL(8) CO2_gm_iso                             !C14 or C13 flux
         REAL(8) TOCinp                                 !total external C input
         REAL(8) TOCout                                 !total output: respiration
         REAL(8) TOCbeg                                 !total mass at the beginning
         REAL(8) TOCend                                 !total mass at the end
         REAL(8) RE                                     ![-],default = 0%, mass balance check = (OC2 - OC1)/OC1*100% 
      END TYPE sMEND_OUT
      
!-----------------------------------------------------------------------------! 
       TYPE sMEND_INI
         INTEGER nCase                                  !# of cases, e.g., # of treatments combined together for calibration
         !INTEGER iFout_SIM_obs_cases                    !file unit for response variables matching observations for all cases
         !INTEGER nVARopt_cases                          !# of observed response variables for optimization
         !REAL(r8) VARobjw_cases(10)                      !VAR objective function weighting factor (any number, will be normalized)
         !CHARACTER(len=4) VARobj_cases(10)              !VAR objective function, e.g.,  "NSEC", "MARE"
         !REAL(r8) rOBJ_cases(10)                         !Nash-Sutcliffe Efficency Coefficient, rNSE(1): mean NSEC, rNSE(2:nObs_var+1): individual NSEC for nObs_var 
         !REAL(r8) rOBJw_cases(10)                        !obj weighting factor
         !CHARACTER(len=20),ALLOCATABLE:: CASEname(:)    !case names  
!         INTEGER, ALLOCATABLE:: VARopt_int_cases(:, :)  !nrow = nVARopt, ncol = 2 (index of the VAR for opt, no. of data, tstep); 
                 
         TYPE(sMEND_INP) sINP                           !INPUT
         CHARACTER(len=10) SITE                         !!Site Name, used for prefix of output files
         CHARACTER(len=3) BIOME                         !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
         CHARACTER(len=3) SOM                           !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
         !CHARACTER(len=100) dirinp                      !!input dir
         !CHARACTER(len=100) dirout                      !!output dir
         !CHARACTER(len=100) dirinp_case                 !!input dir for a case
         !character(len=8) sDate_beg_all                 !!"yyyymmdd": available input data: begin date for optimizaton
         !character(len=8) sDate_end_all                 !!"yyyymmdd": available input data: end date for optimizaton
         !character(len=8) sDate_beg_sim                 !!"yyyymmdd": simulation: begin date 
         !character(len=8) sDate_end_sim                 !!"yyyymmdd": simulation: end date 
         !character(len=8) sDate_beg_inp2                !!"yyyymmdd": begin date for available constant input data 
         !character(len=8) sDate_end_inp2                !!"yyyymmdd": end date for available constant input data 
         !CHARACTER(len=20) SOIL_INI_file               !soil initialization file
         !CHARACTER(len=20),ALLOCATABLE:: VARfile(:)     !VAR data file
         !CHARACTER(len=4), ALLOCATABLE:: VARobj(:)      !VAR objective function type (NSEC or MARE)  
         INTEGER iModel                                 !0: run model; 1: run model optimization; 2: sensitivity        
         INTEGER nPar                                   !# of parameters to be optimized
         INTEGER iKinetics                              !decomposition kinetics: 0-Michaelis-Menten, 1-First Order, 2-Second Order
         !INTEGER nVARopt                                !# of observed response variables for optimization
         !INTEGER nOBS_tot                               !total # of observations for those output variables, = sum(VARopt_int(:, 2) )
         !INTEGER nYear                                  !simulation years (no-use) 
         !INTEGER nHour                                  !available input data period length
         !INTEGER nHour_sim                              !simulation period length
         !INTEGER nOutStep                               !output time interval (step)        
         !INTEGER iFout_SIM_obs                          !file unit for response variables matching observations
         !INTEGER iFout_SIM_day                          !file unit for response variables, mean daily output
         !INTEGER iFout_SIM_mon                          !file unit for response variables, mean monthly output
         !INTEGER iFout_VAR_hour                          !file unit for all response variables, user-defined output time-step
         !INTEGER iFout_FLX_hour                          !file unit for all flux variables, user-defined output time-step
         !INTEGER iFout_RATE_hour                        !file unit for derived rates, e.g., 1st-order decomposition rate, active fraction, phi
         !INTEGER iFout_ITW_hour                         !file unit for Input (e.g., litter), Temperature, Water content & potential 
         !INTEGER iFout_PAR_hour                         !file unit for hourly parameters modified by Temperature, Water potential, and other factors 
         
         !INTEGER iFout_UQvar                            !file unit for response variables in COFI uncertainty
         REAL(r8) r0                                     !initial active-microbe fraction
         REAL(r8) LCI0                                   !initial lignocellulose index = Lignin/(Lignin+Cellulose)!file unit for output
         REAL(r8) soilDepth                              ![cm],soil depth
         REAL(r8) SIN_C12_C14                            !ratio of C12 to C14 in SOC input
         !REAL(r8), ALLOCATABLE:: rOBJ(:)                 !Nash-Sutcliffe Efficency Coefficient, rNSE(1): mean NSEC, rNSE(2:nObs_var+1): individual NSEC for nObs_var 
         !REAL(r8), ALLOCATABLE:: rOBJw(:)                !obj weighting factor
         !REAL(r8), ALLOCATABLE:: dOBS_opt(:, :)          !date,obs,iVARopt
         !REAL(r8), ALLOCATABLE:: dSIM_opt(:, :)          !date,sim,sim_sd
         REAL(r8) :: STP                  ![C], soil temperature
         REAL(r8) :: SWC                  ![fraction],soil water content
         REAL(r8) :: SWP                  ![MPa], soil water potential 
         REAL(r8) :: SpH                  !soil pH
         REAL(r8) :: SIN                  ![mgC/cm3/h],SOC input, e.g., litter
         REAL(r8) SIN_other(2,3)                         ![mgC/cm3/h | mgC/cm3], other constant inputs to 3 pools, e.g., coarse wood, roots, SIN_other(1,:): between sDate_beg_all & sDate_end_all; SIN_other(2,:): between sDate_beg_inp2 & sDate_end_inp2; 
         REAL(r8) SIN_frac(3)                            !fraction of SOC input to 3 pools (POC1, POC2, & DOC)
         REAL(r8) SIN_Multiplier                         !multiplier for litter input during post-data period simulation: SIN*SIN_multiplier
         !REAL(r8), ALLOCATABLE:: VARobjw(:)              !VAR objective function weighting factor (any number, will be normalized) 
         !INTEGER, ALLOCATABLE:: VARopt(:)               !if use observed VAR to calibrate model, see MEND.ini
         !INTEGER, ALLOCATABLE:: VARstep(:)              !time-step, 0(hourly),1(daily),2(monthly),3(yearly)
         !INTEGER, ALLOCATABLE:: VARcol(:)               !VAR data in which column of input file
         !INTEGER, ALLOCATABLE:: VARopt_int(:, :)        !nrow = nVARopt, ncol = 3 (index of the VAR for opt, no. of data, tstep); tstep = 0(hourly),1(daily),2(monthly),3(yearly)
      END TYPE sMEND_INI

!-----------------------------------------------------------------------------! 
END MODULE MOD_MEND_TYPE 

