MODULE MOD_MEND  
! File:   TMEND.F90
! Author: GANGSHENG WANG @ ORNL
! Updated: April 28, 2015
! Created on February 26, 2013, 11:02 AM
! Updated on 7/15/2015
    
    USE MOD_MEND_TYPE
    USE MOD_USRFS
    use shr_kind_mod, only: r8=>shr_kind_r8
    IMPLICIT NONE

    !PRIVATE:: subMEND_INI
    !PRIVATE:: subMEND_output
    !PRIVATE:: subMEND_output_rate
    !PRIVATE:: subMEND_RUN
    !PRIVATE:: subMEND
    !PRIVATE:: subMEND_PAR
    !PRIVATE:: sOUT_OPT_h
    
    !PUBLIC :: sINP_Read
    !PUBLIC :: sOUT_Day2Mon
    !PUBLIC :: sOUT_OPT
    !PUBLIC :: sOUT_ALL_tscale
    
    !PUBLIC:: fMEND_OBJ
    PUBLIC:: fV_MM     !!Michaelis-Menten Kinetics
    PUBLIC:: fAds      !!adsorption-desorption
    
    PUBLIC:: fTArrhenius
    PUBLIC:: fTArh0
    PUBLIC:: fTArh
    PUBLIC:: fT_Linear
    PUBLIC:: fT_CUE
    PUBLIC:: fTAG
    
    PUBLIC:: fSWP2SWC
    PUBLIC:: fSWC2SWP
    PUBLIC:: fSWP0
    PUBLIC:: fSWP
    PUBLIC:: fSWP_Death
    PUBLIC:: fSWP_OPT
    PUBLIC:: fSWP_CLM
    PUBLIC:: fSWPsat
    PUBLIC:: fSWP_A2D
    PUBLIC:: fSWP_D2A
    PUBLIC:: fpH
    PUBLIC:: fpH0
    PUBLIC:: fPermil
    
    
    CONTAINS
!!---------------------------------------------------------------------------- 
!! [1] MEND MODEL: INITILIZATION, MODEL, RUN CONTROL, COMPUTATION OF OBF
!! [2] MODEL KINETICS
!! [3] SOIL TEMPERATURE SCALAR
!! [4] SOIL WATER SCALAR
!! [5] SOIL pH SCALAR
!! [6] ISOTOPES
!! [7] OUTPUT PROCESSING: VARIABLE EXTRACTION, TIME SCALE CONVERSION
!!----------------------------------------------------------------------------

!=============================================================================
!MEND: Microbial-Enzyme-mediated Nitrification Denitrification & Decomposition
!-----------------------------------------------------------------------------
SUBROUTINE subMEND_INI(sINI)
    !model initialization

    !!ARGUMENTS:
    TYPE(sMEND_INI), intent(inout) :: sINI
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_INP):: sINP
    INTEGER, PARAMETER :: nProp = 13 !# of soil properties
    REAL(r8) dINI(nProp)
    !INTEGER nObs_time, nObs_var, nObs_col != nSoil*nSubstrate*nObs_var  !# of columns in observation data file

    INTEGER i, j, k!, lp
    !INTEGER ifini!, ifobs  
    REAL(r8) frISO(const_nISO), frISOadd(const_nISO)
    REAL(r8) rLig_Cel(2) !proportion of Lignin and Cellulose in POC
    REAL(r8) rQOC !proportion of QOC in total MOC

    !CHARACTER(len = 100) sRead!!propName(nProp) !name of 10 properties 
    !INTEGER iRead  !!iTemp
    !REAL(r8) rRead  !!rTemp
    
    rLig_Cel = (/sINI%LCI0, 1.0 - sINI%LCI0/) !proportion of Lignin and Cellulose in POC
    rQOC = 0.01 !proportion of QOC in total MOC
    dINI = (/7.98, 0.88, 3.32, 4.65, 0.17, 0.10, 0.0003, 0.0003, 0.0003, 510., 350., 140., 10./)
    ! SOC, TON, POC, MOC, MBC, DOC, EP1, EP2, EM, Sand, Silt, Clay, Depth
    ! units: clay, sand, silt, 1/1000; SOC, TON, POC, MOC, MBC, mgC/g soil; Depth, cm
    frISO(2) = const_Rstd(1)/(1.0 + const_Rstd(1)) !standard ratio of C14/C12 - 1e-12
    frISO(1) = 1. - frISO(2) !e.g., C12 
    sINI % r0 = 0.1

!! sINP%CPOOL        !carbon pools
    sINI % soilDepth = dINI(13)  !![cm], soil depth
    sINP % CPOOL % POC = dINI(3) * rLig_Cel ![mg C/g soil],Particulate Organic Carbon, size of {} determined by nPOC, !(/DBLE(0.958), DBLE(2.290)/)  !reshape((/DBLE(0.958), DBLE(2.290)/),shape(sINP%CPOOL%POC))
    sINP % CPOOL % MOC = dINI(4)*(1.0 - rQOC) ![mg C/g soil],Mineral Associate Organic Carbon  !DBLE(27.868) 
    sINP % CPOOL % QOC = dINI(4) * rQOC ![mg C/g soil],adsorbed phase of DOC  !DBLE(0.100)
    sINP % CPOOL % MBC = dINI(5) ![mg C/g soil],Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % DOC = dINI(6) ![mg C/g soil],Dissolved Organic Carbon     !DBLE(0.210)  
    sINP % CPOOL % MBCA = sINP % CPOOL % MBC * sINI % r0 ![mg C/g soil],Active Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % MBCD = sINP % CPOOL % MBC * (1.0 - sINI % r0) ![mg C/g soil],Active Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % ENZP = (/dINI(7),dINI(8)/)!!(/1.1d-3, 1.1d-3/) ![mg C/g soil],ENZyme for POC        !reshape((/DBLE(1.1d-3),DBLE(1.1d-3)/),shape(sINP%CPOOL%ENZP))   
    sINP % CPOOL % ENZM = dINI(9) !!1.4d-3 ![mg C/g soil],Enzyme for MAOC  
    sINP % CPOOL % CO2 = 0. ![mg C/g soil],CO2 in the soil
    sINP % CPOOL % SOC = sum(sINP % CPOOL % POC) + sINP % CPOOL % MOC + sINP % CPOOL % QOC
!!     END TYPE sMEND_INP
    do j = 1, const_nPOC
        sINP % CPOOLIFR % POC(j) = frISO
        sINP % CPOOLIFR % ENZP(j) = frISO
        sINP % CPOOLI % POC(j) = sINP % CPOOLIFR % POC(j) * sINP % CPOOL % POC(j)
        sINP % CPOOLI % ENZP(j) = sINP % CPOOLIFR % ENZP(j) * sINP % CPOOL % ENZP(j)
    end do

    sINP % CPOOLIFR % MOC = frISO
    sINP % CPOOLIFR % ENZM = frISO
    sINP % CPOOLIFR % QOC = frISO
    sINP % CPOOLIFR % DOC = frISO
    sINP % CPOOLIFR % MBC = frISO
    sINP % CPOOLIFR % MBCA = frISO
    sINP % CPOOLIFR % MBCD = frISO
    sINP % CPOOLIFR % CO2 = frISO

    sINP % CPOOLI % MOC = sINP % CPOOLIFR % MOC * sINP % CPOOL % MOC
    sINP % CPOOLI % ENZM = sINP % CPOOLIFR % ENZM * sINP % CPOOL % ENZM
    sINP % CPOOLI % QOC = sINP % CPOOLIFR % QOC * sINP % CPOOL % QOC
    sINP % CPOOLI % DOC = sINP % CPOOLIFR % DOC * sINP % CPOOL % DOC
    sINP % CPOOLI % MBC = sINP % CPOOLIFR % MBC * sINP % CPOOL % MBC
    sINP % CPOOLI % MBCA = sINP % CPOOLIFR % MBCA * sINP % CPOOL % MBCA
    sINP % CPOOLI % MBCD = sINP % CPOOLIFR % MBCD * sINP % CPOOL % MBCD
    sINP % CPOOLI % CO2 = sINP % CPOOLIFR % CO2 * sINP % CPOOL % CO2

    do j = 1, const_nISO
        sINP % CPOOLI(j) % SOC  = sum(sINP % CPOOLI(j) % POC) &
        &                       + sINP % CPOOLI(j) % MOC + sINP % CPOOLI(j) % QOC
    end do
    sINP % CPOOLIFR % SOC = sINP % CPOOLI % SOC/sINP % CPOOL % SOC

    do j = 1, const_nISO - 1
        do k = 1, const_nPOC
            sINP % CPOOLI_SIG(j) % POC(k) = &
            & fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % POC(k), sINP % CPOOLI(j + 1) % POC(k))
            sINP % CPOOLI_SIG(j) % ENZP(k) = &
            & fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % ENZP(k), sINP % CPOOLI(j + 1) % ENZP(k))
        end do

        sINP % CPOOLI_SIG(j) % MOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MOC, sINP % CPOOLI(j + 1) % MOC)
        sINP % CPOOLI_SIG(j) % QOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % QOC, sINP % CPOOLI(j + 1) % QOC)
        sINP % CPOOLI_SIG(j) % MBC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBC, sINP % CPOOLI(j + 1) % MBC)
        sINP % CPOOLI_SIG(j) % MBCA = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBCA, sINP % CPOOLI(j + 1) % MBCA)
        sINP % CPOOLI_SIG(j) % MBCD = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBCD, sINP % CPOOLI(j + 1) % MBCD)
        sINP % CPOOLI_SIG(j) % DOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % DOC, sINP % CPOOLI(j + 1) % DOC)
        sINP % CPOOLI_SIG(j) % ENZM = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % ENZM, sINP % CPOOLI(j + 1) % ENZM)
        sINP % CPOOLI_SIG(j) % CO2 = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % CO2, sINP % CPOOLI(j + 1) % CO2)
        sINP % CPOOLI_SIG(j) % SOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % SOC, sINP % CPOOLI(j + 1) % SOC)
    end do !j = 1, const_nISO - 1

    sINI % sINP = sINP

END SUBROUTINE subMEND_INI

!-----------------------------------------------------------------------------
SUBROUTINE subMEND_output_rate(sINI, sPAR, sINP, sOUT)
!    USE MOD_MEND
!!Hourly output for derived parameters
!!"Hour","kPOC1","kPOC2","kMOC","kDOC","kMBC","phi","rMBA"
    !!ARGUMENTS:
    TYPE(sMEND_INI), intent(in) :: sINI
    TYPE(sMEND_PAR), intent(in) :: sPAR
    TYPE(sMEND_INP), intent(in) :: sINP
    TYPE(sMEND_OUT), intent(in) :: sOUT   
    !CHARACTER(LEN=8),intent(in) :: sDate
    !INTEGER,         intent(in) :: ihr
    
    !!LOCAL VARIABLES 
    REAL(r8) kPOC1           !!Equivalent 1st-order decomposition rate; k=Vd*ENZP/(POC + Ks)
    REAL(r8) kPOC2           !!Equivalent 1st-order decomposition rate; k=Vd*ENZP/(POC + Ks)
    REAL(r8) kMOC            !!Equivalent 1st-order decomposition rate; k=Vd*ENZM/(MOC + Ks)
    REAL(r8) kDOC            !!Equivalent 1st-order turnover rate of DOC; k=[(Vg+Vm)/Yg]*MBC/(DOC + Ks)
    REAL(r8) kMBa            !!Equivalent 1st-order turnover rate of Active Mcirobes; k=respiration_rate + mortality_rate+Enzyme_production_rate
    REAL(r8) kMBa_in         !!Equivalent 1st-order assimilation rate of active microbes 
    REAL(r8) kMBd            !!Equivalent 1st-order turnover rate of dormant microbes
    REAL(r8) kMBd_in         !!Equivalent 1st-order microbial dormancy rate 
    REAL(r8) kMB             !!Equivalent 1st-order turnover rate of microbes
    REAL(r8) kMB_in          !!Equivalent 1st-order assimilation rate of total microbes
    REAL(r8) phi             !!DOC saturation level; phi=DOC/(DOC+Ks)
    REAL(r8) rMBa            !!Active Fraction of Microbes, r=MBCA/MBC
    REAL(r8) CUE             !!Apparent Microbial Carbon Use Efficiency, CUE=[DOM_to_MBA - CO2_gmo - death]/DOM_to_MBA
    
    phi = sINP%CPOOL%DOC/(sPAR%KsDOC+sINP%CPOOL%DOC)
    rMBa = sINP%CPOOL%MBCA/sINP%CPOOL%MBC
    kPOC1 = sPAR%VdPOC(1)*sINP%CPOOL%ENZP(1)/(sPAR%KsPOC(1)+sINP%CPOOL%POC(1)) 
    kPOC2 = sPAR%VdPOC(2)*sINP%CPOOL%ENZP(2)/(sPAR%KsPOC(2)+sINP%CPOOL%POC(2)) 
    kMOC = sPAR%VdMOC*sINP%CPOOL%ENZM/(sPAR%KsMOC+sINP%CPOOL%MOC)
    kDOC = (sPAR%Vg + sPAR%Vm)/sPAR%Yg*sINP%CPOOL%MBCA/(sPAR%KsDOC+sINP%CPOOL%DOC)
    kMBa = (sPAR%Vg + sPAR%Vm)*(1.0/sPAR%Yg - 1.0)*phi &
            + sPAR%rMORT + (sPAR%pENZP+sPAR%pENZM) * sPAR%Vm 
    kMBa_in = (sPAR%Vg + sPAR%Vm)/sPAR%Yg*phi + sOUT%CFLUX%MBCD_to_MBCA/sINP%CPOOL%MBCA
    kMBd_in = sOUT%CFLUX%MBCA_to_MBCD/sINP%CPOOL%MBCD
    kMBd = (sOUT%CFLUX%MBCD_to_MBCA + sOUT%CFLUX%CO2_maintn_dorm)/sINP%CPOOL%MBCD
    kMB = (sOUT%CFLUX%CO2_gm + sOUT%CFLUX%MBC_PM)/sINP%CPOOL%MBC
    kMB_in = sOUT%CFLUX%DOC_to_MBC/sINP%CPOOL%MBC       
    if (sOUT%CFLUX%DOC_to_MBC > 0.0) then
        CUE = (sOUT%CFLUX%DOC_to_MBC - sOUT%CFLUX%CO2_gm - sOUT%CFLUX%MBC_PM)/sOUT%CFLUX%DOC_to_MBC
    else
        CUE = 0.
    end if

END SUBROUTINE subMEND_output_rate
!-----------------------------------------------------------------------------
SUBROUTINE subMEND_RUN(xx, sPAR, sINI, sOUT)
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(inout) :: sPAR
    TYPE(sMEND_INI), intent(inout) :: sINI
    TYPE(sMEND_OUT), intent(inout) :: sOUT
    REAL(r8)        , intent(in)    :: xx(sINI%nPar)
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_INP) sINP
    INTEGER j  
!!    INTEGER*4 i, k 
 !   INTEGER j, lp, jbeg,jend, tstep, ibeg, iend, iday
 !   INTEGER iFunc, iObs_time, iHour
    REAL(r8) frISOadd(const_nISO)
!    INTEGER nHour,nday,nmon
!    REAL(8), ALLOCATABLE:: dSIM_d(:,:)  !!daily
!    REAL(8), ALLOCATABLE:: dSIM_m(:,:)  !!monthly
!    INTEGER nOBS
    
!    CHARACTER(len=6) sYM
!    CHARACTER(len=8) sDate
!    INTEGER iYYYYMM, iyr,imo,ida, ihr
!    REAL(8) dSIM_h(sINI%nHour_sim,sINI%nVARopt)  !!store hourly output
!    REAL(8) dSIM(24,sINI%nVARopt)  !!store hourly output
    
!    character(len=200) sFile_inp
!    character(len=200) sRead
!    INTEGER iRead
!    REAL(8) rRead(5)

!!ATTENTION: use ASSOCIATE will result in subroutine invisible in Navigator
!    ASSOCIATE(sINP => sINI%sINP) 
    
    !!replace sDate_beg_all & sDate_end_all with sDate_beg_sim & sDate_end_sim 
!    nday = nDaysbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
!    nmon = nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
!    ALLOCATE(dSIM_d(nday,sINI%nVARopt*2)) !!mean & sd
!    ALLOCATE(dSIM_m(nmon,sINI%nVARopt*2)) !!mean & sd
    
    sINI%LCI0       = xx(1) !initial LCI
    sINI%r0         = xx(2) !fraction of active biomass
    
    frISOadd(2) = 1.0/(1.0 + sINI%SIN_C12_C14) !C14
    frISOadd(1) = 1.0 - frISOadd(2)            !C12
    
    sINP = sINI%sINP
    
    !!Initial Step:
    sOUT % CPOOL    = sINP % CPOOL
    sOUT % CPOOLI   = sINP % CPOOLI
    sOUT % CPOOLIFR = sINP % CPOOLIFR
    
    sINP % CPOOL = sOUT % CPOOL
    sINP % CPOOLI = sOUT % CPOOLI
    sINP % CPOOLIFR = sOUT % CPOOLIFR                
                
    sINP%SIN = sINI%SIN
    sINP%tmp = sINI%STP
    sINP%SWP = sINI%SWP
    sINP%pH  = sINI%SpH

    !jbeg = 1
    !jend = 86400            
    !!External INPUT: BEGIN
    sINP % CADD % POCadd(1) = sINP%SIN * sINI%SIN_frac(1) !+ sINI%SIN_other(1,1)  !![mg POC/g soil/h], inputs to POC
    sINP % CADD % POCadd(2) = sINP%SIN * sINI%SIN_frac(2) !+ sINI%SIN_other(1,2)
    sINP % CADD % DOCadd    = sINP%SIN * sINI%SIN_frac(3) !+ sINI%SIN_other(1,3)  !![mg DOC/g soil/h], inputs to DOC

    !if(i.ge.jbeg.and.i.le.jend) then
    !    sINP%CADD%POCadd(1) = sINP%CADD%POCadd(1)+sINI%SIN_other(2,1)/DBLE(jend-jbeg+1)
    !    sINP%CADD%POCadd(2) = sINP%CADD%POCadd(2)+sINI%SIN_other(2,2)/DBLE(jend-jbeg+1)
    !    sINP%CADD%DOCadd    = sINP%CADD%DOCadd   +sINI%SIN_other(2,3)/DBLE(jend-jbeg+1)
    !end if
    !        write(*,*)i,sINP % CADD % POCadd(1:2),sINP % CADD % DOCadd

    do j = 1, const_nPOC
        sINP % CADDI % POCadd(j) = sINP%CADD%POCadd(j)*frISOadd
    end do
    sINP % CADDI % DOCadd = sINP%CADD%DOCadd*frISOadd
    !!External INPUT: END

    sINI%sINP = sINP
    call subMEND_PAR(xx, sPAR, sINI)   !!modify parameter values by soil temperature, water potential, pH
    call subMEND(sPAR, sINI, sOUT) !!MEND model 
    !output results for model simulation
    call subMEND_output_rate(sINI, sPAR, sINP, sOUT)

END SUBROUTINE subMEND_RUN
!-----------------------------------------------------------------------------

SUBROUTINE subMEND(sPAR, sINI, sOUT)
    !model components 
!    USE MOD_MEND
!    IMPLICIT NONE
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(in) :: sPAR
    TYPE(sMEND_INI), intent(in) :: sINI
    TYPE(sMEND_OUT), intent(out):: sOUT
!    REAL(8)        , intent(in) :: xx(sPAR % nPar)  !!used in subMEND_PAR
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_INP) sINP
    TYPE(sMM_PAR) smmPAR
    TYPE(sMM_INP) smmINP

    TYPE(sSORP_PAR) sSorpPAR
    TYPE(sSORP_INP) sSorpINP
    TYPE(sSORP_OUT) sSorpOUT

    INTEGER i, j, k, iFunc
    REAL(r8) OC1, OC2, sumPOC, sumENZP, phi
    REAL(r8) frPOC(const_nPOC), frENZP(const_nPOC), POCdec(const_nPOC), ENZP_loss(const_nPOC)
    
!    ASSOCIATE(sINP => sINI%sINP) 
    sINP = sINI%sINP

!    write(*,'(a15,50d10.3)')"sINP%CPOOL=",sINP%CPOOL
    !compute initial total mass of the system 
    sOUT%TOCbeg = sINP % CPOOL % SOC + sINP % CPOOL % DOC + sINP % CPOOL % MBC &
                + sum(sINP % CPOOL % ENZP) + sINP % CPOOL % ENZM
    sOUT%TOCinp = sum(sINP % CADD % POCadd) + sINP % CADD % DOCadd
    !Flux 2: POC decomposition
    DO i = 1, const_nPOC
        smmINP % substrate = sINP % CPOOL % POC(i)
        smmINP % enzyme = sINP % CPOOL % ENZP(i)
        smmPAR % vm = sPAR % VdPOC(i)
        smmPAR % km = sPAR % KsPOC(i)
        POCdec(i) = fV_MM(sINI%iKinetics, smmPAR, smmINP)
        sOUT % CFLUX % POCdec(i) = min(POCdec(i), sINP % CPOOL % POC(i))
        sOUT%CFLUX%POCdec_to_DOC(i) = sPAR%frPOC2DOC*sOUT % CFLUX % POCdec(i) 
        sOUT%CFLUX%POCdec_to_MOC(i) = (1.0 - sPAR%frPOC2DOC) * sOUT%CFLUX%POCdec(i) 
    end DO!for i = 1:sINP%nPOC

    !Flux 3: MOC decomposition
    smmINP % substrate = sINP % CPOOL % MOC
    smmINP % enzyme = sINP % CPOOL % ENZM
    smmPAR % vm = sPAR % VdMOC
    smmPAR % km = sPAR % KsMOC
    sOUT % CFLUX % MOC_to_DOC = fV_MM(sINI%iKinetics, smmPAR, smmINP) !MOC decomposition
    sOUT % CFLUX % MOC_to_DOC = min(sOUT % CFLUX % MOC_to_DOC, sINP % CPOOL % MOC)

   !Flux 1: DOC uptake by MB
    smmINP % substrate = sINP % CPOOL % DOC
    smmINP % enzyme = sINP % CPOOL % MBCA
    smmPAR % vm = (sPAR % Vg + sPAR % Vm)/sPAR % Yg !GROWTH + MAINTENANCE    
    smmPAR % km = sPAR % KsDOC
    sOUT % CFLUX % DOC_to_MBC = fV_MM(0, smmPAR, smmINP)
    
    if(sOUT % CFLUX % DOC_to_MBC.gt.sINP % CPOOL % DOC) then
!        print*,"DEBUG:",sOUT % CFLUX % DOC_to_MBC,sINP % CPOOL % DOC,min(sOUT % CFLUX % DOC_to_MBC,sINP % CPOOL % DOC)
        sOUT % CFLUX % DOC_to_MBC = sINP % CPOOL % DOC       
    end if
    !!FIRST UPDATE, used for sorption-desorption calculation
    sOUT % CPOOL % DOC = sINP % CPOOL % DOC - sOUT % CFLUX % DOC_to_MBC  
    
    !Flux 4 & 5: Growth and Maintenance Respiration
    !emission of CO2 from soil to air NOT considered  
    sOUT % CFLUX % CO2_maintn_dorm = sPAR%VmD * sINP%CPOOL%MBCD  !!Microbial Maintenance of Dormant Microbes 
    
!    sOUT % CFLUX % CO2_growth  = sOUT % CFLUX % DOC_to_MBC * (1.D0 - sPAR % Yg) * sPAR % Vg/(sPAR % Vg + sPAR % Vm)
!    sOUT % CFLUX % CO2_maintn  = sOUT % CFLUX % DOC_to_MBC * (1.D0 - sPAR % Yg) - sOUT % CFLUX % CO2_growth
    
    smmINP % substrate = sINP % CPOOL % DOC
    smmINP % enzyme = sINP % CPOOL % MBCA
    smmPAR % vm = sPAR % Vg * (1.0/sPAR % Yg - 1.0)
    smmPAR % km = sPAR % KsDOC
    sOUT % CFLUX % CO2_growth = fV_MM(0, smmPAR, smmINP)

    smmPAR % vm = sPAR % Vm * (1.0/sPAR % Yg - 1.0)
    sOUT % CFLUX % CO2_maintn = fV_MM(0, smmPAR, smmINP) 
                                    
    sOUT % CFLUX % CO2_gm = sOUT%CFLUX%CO2_growth + sOUT%CFLUX%CO2_maintn + sOUT%CFLUX%CO2_maintn_dorm

    !Flux 6 & 7: DOC adsorption-desorption
    sSorpINP % adsorbate = sOUT % CPOOL % DOC !!DOC has been update;  = sINP%CPOOL%DOC - sOUT%CFLUX%DOC_to_MBC
    sSorpINP % adsorbent = sINP % CPOOL % QOC
    sSorpPAR % Qmax = sPAR % Qmax
    sSorpPAR % Kads = sPAR % Kads
    sSorpPAR % Kdes = sPAR % Kdes
    iFunc = fAds(sSorpPAR, sSorpINP, sSorpOUT)
    sOUT % CFLUX % QOC_to_DOC = sSorpOUT % des
    sOUT % CFLUX % DOC_to_QOC = sSorpOUT % ads
    sOUT % CFLUX % DOC_to_QOC_net = sSorpOUT % ads_net

    !Flux 8: microbial mortality
!    sOUT % CFLUX % MBC_PM = sPAR % Vm * sINP % CPOOL % MBCA
!    sOUT % CFLUX % MBC_mortality = (1.d0 - sPAR % pENZP - sPAR % pENZM) * sOUT % CFLUX % MBC_PM
    sOUT % CFLUX % MBC_mortality = sPAR%rMORT * sINP%CPOOL%MBCA
    sOUT % CFLUX % MBC_to_DOC = sPAR % frMBC2DOC * sOUT % CFLUX % MBC_mortality
    sOUT % CFLUX % MBC_to_POC = (1.0 - sPAR % frMBC2DOC) * sOUT % CFLUX % MBC_mortality
    !Flux 9: enzyme synthesis
    sumPOC = sum(sINP % CPOOL % POC)
    frPOC = sINP % CPOOL % POC/sumPOC !compute the fraction of each type of POC, e%g%, lignin, cellulose 
!    sOUT % CFLUX % MBC_to_ENZP = frPOC * sPAR % pENZP * sOUT % CFLUX % MBC_PM
!    sOUT % CFLUX % MBC_to_ENZM = sPAR % pENZM * sOUT % CFLUX % MBC_PM
    sOUT % CFLUX % MBC_to_ENZP = frPOC * sPAR % pENZP * sPAR % Vm * sINP % CPOOL % MBCA
    sOUT % CFLUX % MBC_to_ENZM =         sPAR % pENZM * sPAR % Vm * sINP % CPOOL % MBCA
    
    sOUT % CFLUX % MBC_PM = sOUT % CFLUX % MBC_mortality + sum(sOUT % CFLUX % MBC_to_ENZP) + sOUT % CFLUX % MBC_to_ENZM

    !Flux 10: enzyme turnover
    sOUT % CFLUX % ENZP_to_DOC = sPAR % rENZP * sINP % CPOOL % ENZP
    sOUT % CFLUX % ENZM_to_DOC = sPAR % rENZM * sINP % CPOOL % ENZM
    
    !Internal Flux: microbial dormancy and reactivation
    phi = sINP % CPOOL % DOC/(sINP % CPOOL % DOC + sPAR % KsDOC)
    sOUT % CFLUX % MBCA_to_MBCD = (1.0 - phi) * sPAR % VmA2D * sINP % CPOOL % MBCA
    sOUT % CFLUX % MBCD_to_MBCA = phi * sPAR % VmD2A * sINP % CPOOL % MBCD

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !Mass Balance
    !Total Carbon
!    write(1001,*)"fMEND = ",sum(sINP % CADD % POCadd)+sINP % CADD % DOCadd, sOUT % CFLUX % CO2_gm
    
    sOUT % CPOOL % POC      = sINP % CPOOL % POC + sINP % CADD % POCadd - sOUT % CFLUX % POCdec !Lignin and Cellulose
    sOUT % CPOOL % POC(1)   = sOUT % CPOOL % POC(1) + sOUT % CFLUX % MBC_to_POC !ONLY Lignin!     
    sOUT % CPOOL % MOC      = sINP % CPOOL % MOC + sum(sOUT % CFLUX % POCdec_to_MOC) - sOUT % CFLUX % MOC_to_DOC
    sOUT % CPOOL % QOC      = sINP % CPOOL % QOC + sOUT % CFLUX % DOC_to_QOC_net !QOC mass balance
    sOUT % CPOOL % MBC      = sINP % CPOOL % MBC + sOUT % CFLUX % DOC_to_MBC &
    &                       - sOUT % CFLUX % CO2_gm - sOUT % CFLUX % MBC_PM
    sOUT % CPOOL % MBCA     = sINP % CPOOL % MBCA + sOUT % CFLUX % DOC_to_MBC &
    &                       - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) &
    &                       - sOUT % CFLUX % MBC_PM &
    &                       - sOUT % CFLUX % MBCA_to_MBCD + sOUT % CFLUX % MBCD_to_MBCA
    sOUT % CPOOL % MBCD     = sINP % CPOOL % MBCD - sOUT % CFLUX % CO2_maintn_dorm &
    &                       + sOUT % CFLUX % MBCA_to_MBCD - sOUT % CFLUX % MBCD_to_MBCA
    sOUT % CPOOL % DOC      = sINP % CPOOL % DOC + sINP % CADD % DOCadd &
    &                       + sum(sOUT % CFLUX % POCdec_to_DOC) &
    &                       + sOUT % CFLUX % MBC_to_DOC &
    &                       + sOUT % CFLUX % MOC_to_DOC &
    &                       + sum(sOUT % CFLUX % ENZP_to_DOC) &
    &                       + sOUT % CFLUX % ENZM_to_DOC &
    &                       - sOUT % CFLUX % DOC_to_MBC &
    &                       - sOUT % CFLUX % DOC_to_QOC_net
    sOUT % CPOOL % ENZP     = sINP % CPOOL % ENZP + sOUT % CFLUX % MBC_to_ENZP - sOUT % CFLUX % ENZP_to_DOC
    sOUT % CPOOL % ENZM     = sINP % CPOOL % ENZM + sOUT % CFLUX % MBC_to_ENZM - sOUT % CFLUX % ENZM_to_DOC
    sOUT % CPOOL % CO2      = sINP % CPOOL % CO2 + sOUT % CFLUX % CO2_gm
    sOUT % CPOOL % SOC      = sum(sOUT % CPOOL % POC) + sOUT % CPOOL % MOC + sOUT % CPOOL % QOC

    do j = 2, const_nISO !j = 1 for lightest isotope, e.g., C12
        !Isotopic components

        sOUT % CPOOLI(j) % POC      = sINP % CPOOLI(j) % POC + sINP % CADDI(j) % POCadd &
        &                           - sOUT % CFLUX % POCdec * sINP % CPOOLIFR(j) % POC !Lignin and Cellulose                             
        sOUT % CPOOLI(j) % POC(1)   = sOUT % CPOOLI(j) % POC(1) & 
        &                           + sOUT % CFLUX % MBC_to_POC * sINP % CPOOLIFR(j) % MBCA !ONLY Lignin!    
        sOUT % CPOOLI(j) % MOC      = sINP % CPOOLI(j) % MOC &
        &                           + sum(sOUT % CFLUX % POCdec_to_MOC * sINP % CPOOLIFR(j) % POC) &
        &                           - sOUT % CFLUX % MOC_to_DOC * sINP % CPOOLIFR(j) % MOC
        sOUT % CPOOLI(j) % QOC      = sINP % CPOOLI(j) % QOC &
        &                           + sOUT % CFLUX % DOC_to_QOC * sINP % CPOOLIFR(j) % DOC &
        &                           - sOUT % CFLUX % QOC_to_DOC * sINP % CPOOLIFR(j) % QOC!QOC mass balance
        sOUT % CPOOLI(j) % MBC      = sINP % CPOOLI(j) % MBC &
        &                           + sOUT % CFLUX % DOC_to_MBC * sINP % CPOOLIFR(j) % DOC &
        &                           - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn + sOUT % CFLUX % MBC_PM) &
        &                           * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBCD
        sOUT % CPOOLI(j) % MBCA     = sINP % CPOOLI(j) % MBCA &
        &                           + sOUT % CFLUX % DOC_to_MBC * sINP % CPOOLIFR(j) % DOC &
        &                           - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn + sOUT % CFLUX % MBC_PM) &
        &                           * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % MBCA_to_MBCD * sINP % CPOOLIFR(j) % MBCA &
        &                           + sOUT % CFLUX % MBCD_to_MBCA * sINP % CPOOLIFR(j) % MBCD
        sOUT % CPOOLI(j) % MBCD     = sINP % CPOOLI(j) % MBCD &
        &                           - sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBCD &
        &                           + sOUT % CFLUX % MBCA_to_MBCD * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % MBCD_to_MBCA * sINP % CPOOLIFR(j) % MBCD
        sOUT % CPOOLI(j) % DOC      = sINP % CPOOLI(j) % DOC + sINP % CADDI(j) % DOCadd &
        &                           + sum(sOUT % CFLUX % POCdec_to_DOC * sINP % CPOOLIFR(j) % POC) &
        &                           + sOUT % CFLUX % MBC_to_DOC * sINP % CPOOLIFR(j) % MBCA &
        &                           + sOUT % CFLUX % MOC_to_DOC * sINP % CPOOLIFR(j) % MOC &
        &                           + sum(sOUT % CFLUX % ENZP_to_DOC * sINP % CPOOLIFR(j) % ENZP) &
        &                           + sOUT % CFLUX % ENZM_to_DOC * sINP % CPOOLIFR(j) % ENZM &
        &                           - sOUT % CFLUX % DOC_to_MBC * sINP % CPOOLIFR(j) % DOC &
        &                           - sOUT % CFLUX % DOC_to_QOC * sINP % CPOOLIFR(j) % DOC &
        &                           + sOUT % CFLUX % QOC_to_DOC * sINP % CPOOLIFR(j) % QOC
        sOUT % CPOOLI(j) % ENZP     = sINP % CPOOLI(j) % ENZP &
        &                           + sOUT % CFLUX % MBC_to_ENZP * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % ENZP_to_DOC * sINP % CPOOLIFR(j) % ENZP
        sOUT % CPOOLI(j) % ENZM     = sINP % CPOOLI(j) % ENZM &
        &                           + sOUT % CFLUX % MBC_to_ENZM * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % ENZM_to_DOC * sINP % CPOOLIFR(j) % ENZM
        sOUT % CPOOLI(j) % CO2      = sINP % CPOOLI(j) % CO2 &
                                    + (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) * sINP % CPOOLIFR(j) % MBCA &
        &                           + sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBCD

        !!CO2_C14 or C13 flux
        sOUT%CO2_gm_iso             = (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) * sINP%CPOOLIFR(2)%MBCA &
        &                           + sOUT % CFLUX % CO2_maintn_dorm * sINP%CPOOLIFR(2)%MBCD
        
        sOUT % CPOOLI(j) % SOC      = sum(sOUT % CPOOLI(j) % POC) &
        &                           + sOUT % CPOOLI(j) % MOC + sOUT % CPOOLI(j) % QOC

    end do !j = 1, const_nISO

    !j = 1, compute the lightest isotope
    DO i = 1, const_nPOC
        sOUT % CPOOLI(1) % POC(i) = sOUT % CPOOL % POC(i) - sum(sOUT % CPOOLI(2:const_nISO) % POC(i)) !Lignin and Cellulose
        sOUT % CPOOLI(1) % ENZP(i) = sOUT % CPOOL % ENZP(i) - sum(sOUT % CPOOLI(2:const_nISO) % ENZP(i))
    END DO
    sOUT % CPOOLI(1) % MOC = sOUT % CPOOL % MOC - sum(sOUT % CPOOLI(2:const_nISO) % MOC)
    sOUT % CPOOLI(1) % QOC = sOUT % CPOOL % QOC - sum(sOUT % CPOOLI(2:const_nISO) % QOC)
    sOUT % CPOOLI(1) % MBC = sOUT % CPOOL % MBC - sum(sOUT % CPOOLI(2:const_nISO) % MBC)
    sOUT % CPOOLI(1) % MBCA = sOUT % CPOOL % MBCA - sum(sOUT % CPOOLI(2:const_nISO) % MBCA)
    sOUT % CPOOLI(1) % MBCD = sOUT % CPOOL % MBCD - sum(sOUT % CPOOLI(2:const_nISO) % MBCD)
    sOUT % CPOOLI(1) % DOC = sOUT % CPOOL % DOC - sum(sOUT % CPOOLI(2:const_nISO) % DOC)
    sOUT % CPOOLI(1) % ENZM = sOUT % CPOOL % ENZM - sum(sOUT % CPOOLI(2:const_nISO) % ENZM)
    sOUT % CPOOLI(1) % CO2 = sOUT % CPOOL % CO2 - sum(sOUT % CPOOLI(2:const_nISO) % CO2)
    sOUT % CPOOLI(1) % SOC = sOUT % CPOOL % SOC - sum(sOUT % CPOOLI(2:const_nISO) % SOC)

    do j = 1, const_nISO
        !Isotopic fractions
        sOUT % CPOOLIFR(j) % POC = sOUT % CPOOLI(j) % POC/sOUT % CPOOL % POC
        sOUT % CPOOLIFR(j) % MOC = sOUT % CPOOLI(j) % MOC/sOUT % CPOOL % MOC
        sOUT % CPOOLIFR(j) % QOC = sOUT % CPOOLI(j) % QOC/sOUT % CPOOL % QOC
        sOUT % CPOOLIFR(j) % MBC = sOUT % CPOOLI(j) % MBC/sOUT % CPOOL % MBC
        sOUT % CPOOLIFR(j) % MBCA = sOUT % CPOOLI(j) % MBCA/sOUT % CPOOL % MBCA
        sOUT % CPOOLIFR(j) % MBCD = sOUT % CPOOLI(j) % MBCD/sOUT % CPOOL % MBCD
        sOUT % CPOOLIFR(j) % DOC = sOUT % CPOOLI(j) % DOC/sOUT % CPOOL % DOC
        sOUT % CPOOLIFR(j) % ENZP = sOUT % CPOOLI(j) % ENZP/sOUT % CPOOL % ENZP
        sOUT % CPOOLIFR(j) % ENZM = sOUT % CPOOLI(j) % ENZM/sOUT % CPOOL % ENZM
        sOUT % CPOOLIFR(j) % CO2 = sOUT % CPOOLI(j) % CO2/sOUT % CPOOL % CO2
        sOUT % CPOOLIFR(j) % SOC = sOUT % CPOOLI(j) % SOC/sOUT % CPOOL % SOC

    end do !j = 1, const_nISO        

    !    print*, sum(sOUT%CPOOLIFR%MBC)!sum(sOUT%CPOOLI(2:const_nISO)%MBC),sOUT%CPOOLI(1:const_nISO)%MBC

    !Convert Isotope Concentration to Signatures [Permil]
    do j = 1, const_nISO - 1

        do k = 1, const_nPOC
            sOUT % CPOOLI_SIG(j) % POC(k) = &
            & fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % POC(k), sOUT % CPOOLI(j + 1) % POC(k))
            sOUT % CPOOLI_SIG(j) % ENZP(k) = &
            & fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % ENZP(k), sOUT % CPOOLI(j + 1) % ENZP(k))
        end do

        sOUT % CPOOLI_SIG(j) % MOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MOC, sOUT % CPOOLI(j + 1) % MOC)
        sOUT % CPOOLI_SIG(j) % QOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % QOC, sOUT % CPOOLI(j + 1) % QOC)
        sOUT % CPOOLI_SIG(j) % MBC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBC, sOUT % CPOOLI(j + 1) % MBC)
        sOUT % CPOOLI_SIG(j) % MBCA = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBCA, sOUT % CPOOLI(j + 1) % MBCA)
        sOUT % CPOOLI_SIG(j) % MBCD = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBCD, sOUT % CPOOLI(j + 1) % MBCD)
        sOUT % CPOOLI_SIG(j) % DOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % DOC, sOUT % CPOOLI(j + 1) % DOC)
        sOUT % CPOOLI_SIG(j) % ENZM = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % ENZM, sOUT % CPOOLI(j + 1) % ENZM)
        sOUT % CPOOLI_SIG(j) % CO2 = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % CO2, sOUT % CPOOLI(j + 1) % CO2)
        sOUT % CPOOLI_SIG(j) % SOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % SOC, sOUT % CPOOLI(j + 1) % SOC)
    end do !j = 1, const_nISO - 1

    !mass balance check for the whole system
    sOUT%TOCend = sOUT % CPOOL % SOC + sOUT % CPOOL % DOC + sOUT % CPOOL % MBC &
                + sum(sOUT % CPOOL % ENZP) + sOUT % CPOOL % ENZM 
    sOUT%TOCout = sOUT % CFLUX % CO2_gm
    sOUT%RE = (sOUT%TOCend - sOUT%TOCbeg) - (sOUT%TOCinp - sOUT%TOCout)
    if(isnan(sOUT % CPOOL%MBCA)) then
!        write(*,*)sPAR
        write(*,'(a15,50d10.3)')"sOUT%CPOOL=",sOUT%CPOOL
        write(*,'(a15,50d10.3)')"sOUT%CFLUX=",sOUT%CFLUX
        stop
    end if
    
    if(abs(sOUT % RE).gt.1.e-8) then
        print*,"Balance Check ERROR: sOUT%RE=",sOUT % RE
    !else
!        print*,"Balance Check: ",sOUT % RE
    end if
    
!    END ASSOCIATE
END SUBROUTINE subMEND!MEND model

!-----------------------------------------------------------------------------
SUBROUTINE subMEND_PAR(xx, sPAR, sINI)
!    MEND MODEL CODE 
!    USE MOD_MEND
!    IMPLICIT NONE
!    REAL(8), PARAMETER :: Tref = 20.d0  !!moved to MOD_MEND_TYPE
!    REAL(8) fTArh, fT_CUE !!FUNCTIONS
!    REAL(8) fSWP, fSWP_OPT, fSWP_A2D, fSWP_D2A !!FUNCTIONS
!    REAL(8) fpH  !!FUNCTIONS
!    ATTENTION: currently the following parameters are NOT used: xx(6): KPC (=KPL), xx(16): pEM (=pEP)
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(inout) :: sPAR
    TYPE(sMEND_INI), intent(inout) :: sINI
    REAL(r8)        , intent(in)    :: xx(sINI%nPar)
    
    !!LOCAL VARIABLES:
    REAL(r8) CUE_slope,CUE_ref, SWPmin !!CUE_slope = (-0.005, -0.016)  
    REAL(r8) Vm0  !!specific maintenance rate without any modification
    
    REAL(r8) tp_scalar       !!temperature scalar
    REAL(r8) wp_scalar       !!INCrease with increasing wp 
    REAL(r8) wp_scalar_low   !!DECrease with increasing wp 
    REAL(r8) wp_scalar_opt   !!with increasing wp, INcrease first to OPTimal condition, then decrease
    REAL(r8) pH_scalar   !!pH scalar
    
    REAL(r8) tp
    REAL(r8) wp
    REAL(r8) pH
    
    CHARACTER(len=3) BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
    CHARACTER(len=3) SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    CHARACTER(len=10) sCase
    
    SWPmin = -13.86  !!for Microbial Mortality
    
!    sPAR%iKinetics = sINI%iKinetics
    BIOME = trim(sINI%BIOME)
    SOM   = trim(sINI%SOM)
    
    tp = sINI%sINP%tmp
    wp = sINI%sINP%SWP
    pH = sINI%sINP%pH
    
    !![1] DECOMPOSITION of POC1 & POC2
    sCase = "LIG"
    tp_scalar       = fTArh(sCase, tp, const_Tref)
    wp_scalar_opt   = fSWP_OPT(wp)
    pH_scalar       = fpH(sCase,pH)
    sPAR % VdPOC(1)     = xx(3)*tp_scalar*wp_scalar_opt*pH_scalar !!(/xx(1), xx(2)/) ![mg POC/mg ENZP/h], maximum reaction rate for conversion of POC by ENZP
    
    sCase = "CEL"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    wp_scalar = fSWP(wp, BIOME, SOM)
    pH_scalar = fpH(sCase,pH)
!    sPAR % VdPOC(2)     = xx(4)*tp_scalar*wp_scalar*pH_scalar
    sPAR % VdPOC(2)     = xx(3)*tp_scalar*wp_scalar*pH_scalar  !!set Vd_cel = Vd_lig
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % KsPOC        = (/xx(6), xx(7)/)*tp_scalar ![mg POC/cm3], half-saturation constant for conversion of POC by ENZP

    !![2] DECOMPOSITION of MOC
    sCase = "ENZ"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    wp_scalar = fSWP(wp, BIOME, SOM)
    pH_scalar = fpH(sCase,pH)
!    sPAR % VdMOC        = xx(5)*tp_scalar*wp_scalar*pH_scalar ![mg MOC/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC
    sPAR % VdMOC        = xx(3)*tp_scalar*wp_scalar*pH_scalar 
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % KsMOC        = xx(8)*tp_scalar ![mg MOC/cm3], half-saturation constant for conversion of MAOC by ENZMAOC

    !![3] ADSORPTION-DESORPTION
    sPAR % Qmax         = xx(9) ![mg C/g soil], adsorption capacity
    sPAR % Kba          = xx(10) ![mg C/g soil/h], binding affinity
    
    sCase = "Kdes"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % Kdes         = xx(11)*tp_scalar ![mg DOC/h],desorption rate constant
    
    sCase = "Kads"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % Kads         = (xx(11) * sPAR % Kba)*tp_scalar ![mg DOC/mg DOC/h], adsorption rate constant = Kba*Kdes
    
    !![4] ENZYME TURNOVER & PRODUCTION
    sPAR % rENZP        = (/xx(12), xx(12)/) !![1/h],turnover rate of ENZP
    sPAR % rENZM        = xx(12) !![1/h], turnover rate of ENZMAOC
    sPAR % pENZP        = xx(13) ![mg ENZP/mg MBC/h], production rate of ENZP 
    sPAR % pENZM        = xx(13)*xx(14) !!wgs:6/22/2015: set pENZM = xx(16)*pENZP, xx(16) ![mg ENZM/mg MBC/h], production rate of ENZMAOC
    
    !![5] ALLOCATION to DOC
    sPAR % frPOC2DOC    = xx(15) ![-], fraction of decomposed POC allocated to DOC
    sPAR % frMBC2DOC    = xx(16) ![-], fraction of dead MBC allocated to DOC
   
    !![6] MICROBIAL UPTAKE of DOM
    sCase = "DOM"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % Vg        = xx(17)*tp_scalar !![mg DOC/mg MBC/h], maximum uptake rate of DOC by MB
    sPAR % alpha     = xx(18) ![-], alpha = Vm/(Vm + Vg) <= 0.5 
    
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    Vm0 = xx(17) * sPAR%alpha/(1.0 - sPAR%alpha)
    sPAR % Vm         = Vm0*tp_scalar !!*wp_scalar ![1/h], specific microbial maintenance rate
     
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)  
    sPAR % KsDOC        = xx(19)*tp_scalar ![mg DOC/cm3], half-saturation constant for uptake of DOC by MB
    
    sCase = "CUE"  !!label, no-use
    CUE_ref   = xx(20)
    CUE_slope = -1.0*xx(21)
    sPAR % Yg = fT_CUE(tp, const_Tref, CUE_slope, CUE_ref) ![-], carbon use efficiency in uptake of DOC by MB 
    
    !![7] MICROBIAL MORTALITY
!    sPAR%wdie = xx(22)
    sPAR%gamma= xx(23) 
!    wp_scalar_low =  fSWP_Death(wp,SWPmin,sPAR%wdie)  
!    sPAR % rMORT = (1.D0 - sPAR%pENZP - sPAR%pENZM)*Vm0*wp_scalar_low !!wgs: 6/19/2015,  !!tp_scalar
!    wp_scalar_low = fSWP_A2D(wp, sPAR % SWP_A2D, sPAR%wdorm)  !!negative effect  !!wdie
       
    !![8] MICROBIAL DORMANCY & RESUSCITATION
    sPAR % beta       = xx(24) !beta = Vm_dormant/Vm
    sPAR % VmD        = sPAR % beta*Vm0 !!*wp_scalar_low  !!<5/15/2015>
    sPAR % SWP_A2D      = xx(25)
    sPAR % tau          = xx(26)
    sPAR % SWP_D2A      = sPAR % tau * sPAR % SWP_A2D 
    sPAR % wdorm        = xx(27)
    
    wp_scalar     = fSWP_D2A(wp, sPAR % SWP_D2A, sPAR % wdorm)  !!positive effect, increase with increasing SWC/SWP
    wp_scalar_low = fSWP_A2D(wp, sPAR % SWP_A2D, sPAR % wdorm)  !!negative effect
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % VmA2D = Vm0*tp_scalar*wp_scalar_low 
    sPAR % VmD2A = Vm0*tp_scalar*wp_scalar
    
    sPAR % rMORT = (sPAR%gamma*Vm0)*wp_scalar_low 

END SUBROUTINE subMEND_PAR
!-----------------------------------------------------------------------------

!==================================================================
!! MODEL KINETICS
!------------------------------------------------------------------  
REAL(r8) FUNCTION fV_MM(iType, sPAR, sINP)
    !Michaelis-Menten Kinetics (Equation)
    !!ARGUMENTS:
    INTEGER,       INTENT(IN) :: iType
    TYPE(sMM_INP), INTENT(IN) :: sINP
    TYPE(sMM_PAR), INTENT(IN) :: sPAR
    
    SELECT CASE (iType)
        CASE (1)  !!First Order
            fV_MM = sPAR%vm * sINP%substrate
        CASE (2)  !!Second Order
            fV_MM = sPAR%vm * sINP%enzyme * sINP%substrate
        CASE DEFAULT !!M-M Kinetics, iType = 0
            fV_MM = sPAR % vm * sINP % substrate * sINP % enzyme/(sPAR % km + sINP % substrate)
            
    END SELECT
    fV_MM = min(fV_MM, sINP % substrate)
    return
END FUNCTION fV_MM
!------------------------------------------------------------------       
INTEGER FUNCTION fAds(sPAR, sINP, sOUT)
    !Adsorption and Desorption of DOC
    !eg: adsorbent: QOC, adsorbate: DOC
    !!ARGUMENTS:
    TYPE(sSORP_INP), INTENT(IN) :: sINP
    TYPE(sSORP_PAR), INTENT(IN) :: sPAR
    TYPE(sSORP_OUT), INTENT(OUT) :: sOUT
    
    sOUT % ads = sPAR % Kads * sINP % adsorbate * (1.0 - sINP % adsorbent/sPAR % Qmax)
    sOUT % des = sPAR % Kdes * sINP % adsorbent/sPAR % Qmax
    if (sOUT % des > (sINP % adsorbent + sOUT % ads)) then
        sOUT % des = sINP % adsorbent + sOUT % ads
    elseif(sOUT % ads > (sINP % adsorbate + sOUT % des)) then
        sOUT % ads = sINP % adsorbate + sOUT % des
    end if
    sOUT % ads_net = sOUT % ads - sOUT % des
    return
END FUNCTION fAds
!------------------------------------------------------------------

!!=================================================================
!! TEMPERATURE SCALAR: BEGIN
!------------------------------------------------------------------
REAL(r8) function fTArrhenius(tmp, Ea)
    !!ARGUMENTS:
    REAL(r8) tmp,Ea
    !Ea [kJ/mol]: Arrhenius activation energy
    !temp (C): temperature
    fTArrhenius = exp(-Ea * 1.e3/(const_R * (tmp + const_tmp_C2K)))
    return
END function fTArrhenius
!------------------------------------------------------------------
REAL(r8) function fTArh0(T, Tref, Ea)
!    USE MOD_MEND
    !fTvm,[-] temperature-adjusted factor for maximum reaction rate in M-M kinetics 
    !fT = 1 at T = Tref
    !Ea [KJ/mol]: Arrhenius activation energy
    !temp (C): temperature
    !!ARGUMENTS:
    REAL(r8), intent(in) :: T, Tref, Ea
    REAL(r8) TKref, TK
    TKref = Tref + const_tmp_C2K
    TK = T + const_tmp_C2K
    fTArh0 = exp(Ea*1.e3/const_R * (1.0/TKref - 1.0/TK))
    return
END function fTArh0

!------------------------------------------------------------------
REAL(r8) function fTArh(sCase, T, Tref)
!    USE MOD_MEND
    !fTvm,[-] temperature-adjusted factor for maximum reaction rate in M-M kinetics 
    !fT = 1 at T = Tref
    !Ea [J/mol]: Arrhenius activation energy
    !temp (C): temperature
    !!ARGUMENTS:  
    CHARACTER(*), intent(inout) :: sCase
    REAL(r8)          , intent(in) :: T, Tref
    
    !!LOCAL VARIABLES:
    REAL(r8) TKref, TK, Ea
    
!    sCase = trim(sCase)
    SELECT CASE (trim(sCase))
        CASE ("BG")  !! Beta-glucosidase
            Ea = 42.2
        CASE ("CBH") !! Cellobiohydrolase
            Ea = 32.2
        CASE ("EG")  !! Endo-glucanase
            Ea = 34.4
        CASE ("PER") !! Peroxidase
            Ea = 52.9
        CASE ("POX") !! Phenol oxidase
            Ea = 53.1
        CASE ("LIG") !! Ligninases
            Ea = 53.0
        CASE ("CEL") !! Cellulases
            Ea = 36.3
        CASE ("NAG") !! N-acetylglutamate synthase
            Ea = 52.9
        CASE ("PAC") !! Acid phosphatases
            Ea = 36.3
        CASE ("PAL") !! Alkaline phosphatases
            Ea = 23.6
        CASE ("PHO") !! PHOSPHATASES
            Ea = 34.4
        CASE ("Km")  !! half-saturation constant
            Ea = 30.
        CASE ("MR")  !! Microbial Maintenance
            Ea = 20.
        CASE ("Kads")!! adsorption
            Ea = 5.
        CASE ("Kdes")!! desorption
            Ea = 20.
        CASE ("DOM") !! DOM uptake
            Ea = 47.
        CASE DEFAULT
            Ea = 47.
    END SELECT
    fTArh = fTArh0(T, Tref, Ea)
    return
END function fTArh
!------------------------------------------------------------------
REAL(r8) function fT_Linear(T, Tref, slope, intercept)
    !fKmT: [mg/m3], half-saturation constant in M-M kinetics (Km)modified by Temperature
    !slope: [mg/m3/C]
    !intercept: [mg/m3]
    !!ARGUMENTS:
    REAL(r8), intent(in) :: T, Tref, slope, intercept
    
    fT_Linear = slope * (T - Tref) + intercept
    return
END function fT_Linear
!------------------------------------------------------------------
REAL(r8) function fT_CUE(T,Tref,slope,CUEref)
    !fT_CUE: [-], half-saturation constant in M-M kinetics (Km)modified by Temperature
    !slope: [1/degree C]
    !intercept: [-]
    !! parameter values: Wang et al. (2015), ISMEJ 
    !!ARGUMENTS:
    REAL(r8), intent(in):: T, Tref,slope,CUEref
!    REAL(8), PARAMETER :: Tref      = 0    !![degree C]
!    REAL(8), PARAMETER :: slope     = -0.01
!    REAL(8), PARAMETER :: intercept = 0.56  !! []
    
    !!LOCAL VARIABLES:
    REAL(r8), PARAMETER:: CUEmax = 0.9
    REAL(r8), PARAMETER:: CUEmin = 0.1e-1
    
    fT_CUE = fT_Linear(T, Tref, slope, CUEref)
    if(fT_CUE.gt.CUEmax) then
        fT_CUE = CUEmax
    elseif(fT_CUE.lt.CUEmin) then
        fT_CUE = CUEmin
    end if
    
    return
END function fT_CUE
!------------------------------------------------------------------
REAL(r8) function fTAG(tmp)
    !Arrhenius Equation used in ECOSYS
    !Grant, et al., 1993. SBB 25, 1317-1329.
    !!ARGUMENTS:
    REAL(r8) tmp
    
    !!LOCAL VARIABLES:
    REAL(r8) A, S, Ha, Hdl, Hdh, b1, b2, a1, a2, a3
    A = 17.1124 !a parameter selected such that fmt=1 at temp = 30 C
    S = 710 ![J/mol/K], the change in entropy
    Ha = 57500 ![J/mol], energy of activation
    Hdl = 192500 ![J/mol], energy of low temperature deactivation
    Hdh = 222500 ![J/mol], energy of high temperature deactivation
    b1 = const_tmp_C2K + tmp
    b2 = const_R * b1
    a1 = exp(A - Ha/b2)
    a2 = exp((Hdl - S * b1)/b2)
    a3 = exp((-Hdh + S * b1)/b2)
    fTAG = b1 * a1/(1.0 + a2 + a3)
    return
END function fTAG
!!-----------------------------------------------------------------
!! TEMPERATURE SCALAR: END
!!=================================================================

!!=================================================================
!! Soil Water Scalar: BEGIN
!!-----------------------------------------------------------------
REAL(r8) function fSWP2SWC(SWP,SWP_units,SWCres,SWCsat,alpha,n)
!!van-Genuchten equation
!!convert SWP(cm) to SWC (0-1)
!!SWP will be converted to cm
    !!ARGUMENTS:
    CHARACTER(len=*) SWP_units  
    REAL(r8),intent(in):: SWP,SWCres,SWCsat,alpha,n
    
    !!LOCAL VARIABLES:
    REAL(r8) SWP_cm, m, eff_sat  !!effective saturation

    if (trim(SWP_units).eq."MPa") then
        SWP_cm = SWP/const_cm2MPa
    else if (trim(SWP_units).eq."bar") then  !!1 bar = 100 kPa = 0.1 MPa
        SWP_cm = SWP*0.1/const_cm2MPa
    else if (trim(SWP_units).eq."kPa") then  !!1 kPa = 1d-3 MPa
        SWP_cm = SWP*1e-3/const_cm2MPa
    else if (trim(SWP_units).eq."Pa") then  !!1 kPa = 1d-6 MPa
        SWP_cm = SWP*1e-6/const_cm2MPa
    else if (trim(SWP_units).eq."mm") then
        SWP_cm = SWP*0.1
    else if (trim(SWP_units).eq."m") then
        SWP_cm = SWP*100.
    end if
    
    m = 1.0-1.0/n
    
    eff_sat = (1/(1+(alpha*abs(SWP_cm))**n))**m
    fSWP2SWC = SWCres + (SWCsat - SWCres)*eff_sat
    return
   
END function fSWP2SWC
!!-----------------------------------------------------------------
REAL(r8) function fSWC2SWP(SWC0,SWCres,SWCsat,alpha,n,SWPmin)
!!van-Genuchten equation, SWP in units of [cm], SWC if fraction (0-1)
!!return fSWC2SWP), actually fSWC2SWP<0
!!convert SWC(0-1) to SWP(MPa)
!!SWC needs be converted to fraction first
!    USE MOD_MEND
!    IMPLICIT NONE
    REAL(r8),PARAMETER::rlim = 1.01
    !!ARGUMENTS:
    REAL(r8),intent(in):: SWC0
    REAL(r8) SWCres,SWCsat,alpha,n
    
    !!LOCAL VARIABLES:
    REAL(r8) SWC,SWPmin  !!min SWP <0
    REAL(r8) m, eff_sat  !!effective saturation
    
    m = 1.0-1.0/n

    if (SWC0.le.SWCres*rlim) then
!        fSWC2SWP = SWPmin
        SWC = SWCres*rlim
    else
        SWC = SWC0
    end if

    if (SWC.lt.SWCsat) then
        eff_sat = (SWC - SWCres)/(SWCsat - SWCres)
        fSWC2SWP = (1.0/(eff_sat**(1.0/m)) - 1.0)**(1.0/n)/alpha
        fSWC2SWP = -1.0*fSWC2SWP*const_cm2MPa
    else
        fSWC2SWP = 0
    end if

    return
END function fSWC2SWP
!------------------------------------------------------------------
REAL(r8) function fSWP0(SWP,SWPmin,w)
    !Rate Scalar for Soil Water Potential (SWP)
    !Manzoni et al (2012), Ecology, 93: 930-938
!    USE MOD_MEND
!    IMPLICIT NONE
    
    REAL(r8), PARAMETER :: SWP_FC = -0.033 ![MPa], field capacity SWP 
    !!ARGUMENTS:
    REAL(r8) SWP, SWPmin ![MPa]
    REAL(r8) w
    
    if (SWP.lt.SWPmin) then
        fSWP0 = 0.0  
    else if (SWP.lt.SWP_FC) then
        fSWP0 = 1-(log(SWP/SWP_FC)/log(SWPmin/SWP_FC))**w
    else 
        fSWP0 = 1.0
    end if
    return
END function fSWP0

!------------------------------------------------------------------
REAL(r8) function fSWP(SWP,BIOME,SOM)
    !SWP Scalar for SOM (Cellulose) decomposition
    !Manzoni et al (2012), Ecology, 93: 930-938
    !!fSWP increases with increasing SWP (wetter condition)
    
    !REAL(8), PARAMETER :: SWP_FC = -0.033 ![MPa], field capacity SWP 

    !!ARGUMENTS:
    REAL(r8) SWP ![MPa]
    CHARACTER(len=3) BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
    CHARACTER(len=3) SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    
    !!LOCAL VARIABLES:
    REAL(r8) SWPmin, w ![MPa]
    
    if (trim(SOM).eq."SOD") then
        SWPmin = -1.71
        w      = 1.43
    else if (trim(SOM).eq."SOL") then
        SWPmin = -13.86
        w      = 1.20
    else if (trim(SOM).eq."LIT") then
        SWPmin = -36.49
        w      = 1.04
    end if
    
    SELECT CASE (trim(BIOME)) 
        CASE ("ASM")  !!Arid, Semiarid, & Mediterranean
            if (trim(SOM).eq."SOL") then
                SWPmin = -10.95
                w      = 1.26
            end if
        CASE ("MGC")  !!Mesic Grassland & Cropland
            if (trim(SOM).eq."SOD") then
                SWPmin = -1.71
                w      = 1.43
            else if (trim(SOM).eq."SOL") then
                SWPmin = -22.61
                w      = 1.11
            else if (trim(SOM).eq."LIT") then
                SWPmin = -39.73
                w      = 0.89
            end if
        CASE ("MDF")    !!Mesic Deciduous Forest
            if (trim(SOM).eq."SOL") then
                SWPmin = -4.97
                w      = 1.07
            else if (trim(SOM).eq."LIT") then
                SWPmin = -29.00
                w      = 1.27
            end if
        CASE ("MCF")    !!MEsic Conifer Forest
            if (trim(SOM).eq."SOL") then
                SWPmin = -8.24
                w      = 1.40
            else if (trim(SOM).eq."LIT") then
                SWPmin = -39.85
                w      = 1.06
            end if
        CASE DEFAULT    !!All Biome Average
            if (trim(SOM).eq."SOD") then
                SWPmin = -1.71
                w      = 1.43
            else if (trim(SOM).eq."SOL") then
                SWPmin = -13.86
                w      = 1.20
            else if (trim(SOM).eq."LIT") then
                SWPmin = -36.49
                w      = 1.04
            end if
    END SELECT
    fSWP = fSWP0(SWP,SWPmin,w)
    return
END function fSWP
!------------------------------------------------------------------
REAL(r8) FUNCTION fSWP_Death(SWP,SWPmin,w)
    !!SWP Scalar for Microbial Mortality
    !!fSWP_Death DEcreases with increasing SWP (wetter condition)

    !!ARGUMENTS:
    REAL(r8) SWP ![MPa]
    !!ARGUMENTS:
    REAL(r8) SWPmin ![MPa]
    REAL(r8) w
!    CHARACTER(len=3) BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
!    CHARACTER(len=3) SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    
    fSWP_Death = 1.0 - fSWP0(SWP,SWPmin,w)
!    fSWP_Death = 1.d0 - fSWP(SWP,BIOME,SOM)
    
    return
END FUNCTION fSWP_Death
!------------------------------------------------------------------
REAL(r8) function fSWP_OPT(SWP)
    !SWP Scalar for SOM (lignin) decomposition
    !Hansen et al (1990), DAISY Model, page 105, Eq (6-16)
    
    REAL(r8), PARAMETER :: SWPmin = -exp(2.5*log(10.))   !![MPa]
    REAL(r8), PARAMETER :: SWPlow = -exp(-1.5*log(10.))
    REAL(r8), PARAMETER :: SWPhigh= -exp(-2.5*log(10.))
    REAL(r8), PARAMETER :: SWPmax = -exp(-4.0*log(10.))
    !!ARGUMENTS:
    REAL(r8) SWP ![MPa]
    
    if (SWP.le.SWPmin) then
        fSWP_OPT = 0.0  
    else if (SWP.le.SWPlow) then
        fSWP_OPT = 0.625-0.25*log10(abs(SWP))
    else if (SWP.le.SWPhigh) then
        fSWP_OPT = 1.0
    else if (SWP.le.SWPmax) then
        fSWP_OPT = (2.5+0.4*log10(abs(SWP)))/1.5
    else
        fSWP_OPT = 0.6
    end if
    return
END function fSWP_OPT
!------------------------------------------------------------------
REAL(r8) function fSWP_CLM(SWP,SWPsat)
    !Rate Scalar for Soil Water Potential (SWP)
    !CLM4.5 Technical Note, page 285, Eq (15.6)
    !Andren and Paustain (1987), Orchard and Cook (1983).
    
    REAL(r8), PARAMETER :: SWPmin = -10. ![MPa], lower limit for SWP control on decomposition rate
    !!ARGUMENTS:
    REAL(r8) SWP, SWPsat ![MPa], saturated SWP
    
    if (SWP.lt.SWPmin) then
        fSWP_CLM = 0.0
    else if (SWP.le.SWPsat) then
        fSWP_CLM = log(SWPmin/SWP)/log(SWPmin/SWPsat)
    else
        fSWP_CLM = 1.0
    end if
    return
END function fSWP_CLM

!------------------------------------------------------------------
REAL(r8) function fSWPsat(p_sand,p_clay)
    !Saturated Soil Water Potential (SWP) [MPa]
    !CLM4.5 Technical Note, page 285, Eq (15.7)
    !Cosby et al. (1984)
    !p_sand (0-100): volume % of sand
    !p_clay (0-100): volume % of clay

    
    REAL(r8), PARAMETER :: SWP_base = -9.8e-5 ![MPa], base SWP
    REAL(r8), PARAMETER :: a = 1.54
    REAL(r8), PARAMETER :: b = 9.5e-3
    REAL(r8), PARAMETER :: c = 6.3e-3
    
    !!ARGUMENTS:
    REAL(r8) p_sand, p_clay, p_silt
    
    p_silt = 100. - p_sand - p_clay
    fSWPsat = SWP_base*exp((a-b*p_sand+c*p_silt)*log(10.))
    return
END function fSWPsat

!------------------------------------------------------------------
REAL(r8) function fSWP_A2D(SWP, SWP_A2D, w)
    !Soil Water Scalar for Microbial Dormancy
    !! Manzoni (2014) SBB, 73: 69-83

    !!ARGUMENTS:
    REAL(r8), intent(in):: SWP_A2D !! = -0.4 [MPa], 
    REAL(r8), intent(in):: SWP
    REAL(r8), intent(in):: w       !! = 4 
    
    fSWP_A2D = abs(SWP)**w/(abs(SWP)**w + abs(SWP_A2D)**w)
!    if(ISNAN(fSWP_A2D)) then
!        print*,"wp_scalar_low=",fSWP_A2D
!    end if
    return
END function fSWP_A2D

!------------------------------------------------------------------
REAL(r8) function fSWP_D2A(SWP, SWP_D2A, w)
    !!Soil Water Scalar for Microbial reactivation
    !!ARGUMENTS:
    REAL(r8), intent(in):: SWP_D2A !! = 1/4*SWP_A2D [MPa]
    REAL(r8), intent(in):: SWP
    REAL(r8), intent(in):: w       !! = 4 
    
    fSWP_D2A = abs(SWP_D2A)**w/(abs(SWP)**w + abs(SWP_D2A)**w)
    return
END function fSWP_D2A
!!-----------------------------------------------------------------
!! Soil Water Scalar: END
!!=================================================================


!!=================================================================
!! pH Scalar: BEGIN
!------------------------------------------------------------------
REAL(r8) FUNCTION fpH(sCase,pH)
!    REAL(8) fpH0 !! function
    !!ARGUMENTS:
    CHARACTER(len=*) , intent(in) :: sCase
    REAL(r8)          , intent(in) :: pH !pH value
    
    !!LOCAL VARIABLES
    REAL(r8) pHopt !optimum pH
    REAL(r8) pHsen !pH sensitivity
    
!    sCase = trim(sCase)
    SELECT CASE (trim(sCase))
        CASE ("BG")  !! Beta-glucosidase
            pHopt = 5.6
            pHsen = 1.7
        CASE ("CBH") !! Cellobiohydrolase
            pHopt = 5.1
            pHsen = 2.1
        CASE ("EG")  !! Endo-glucanase
            pHopt = 5.1
            pHsen = 1.6
        CASE ("PER") !! Peroxidase
            pHopt = 4.5
            pHsen = 1.5
        CASE ("POX") !! Phenol oxidase
            pHopt = 4.1
            pHsen = 1.4
        CASE ("LIG") !! Ligninases
            pHopt = 4.2
            pHsen = 1.4
        CASE ("CEL") !! Cellulases
            pHopt = 5.3
            pHsen = 1.7
        CASE ("PAC") !! Acid phosphatases
            pHopt = 5.2
            pHsen = 1.8
        CASE ("PAL") !! Alkaline phosphatases
            pHopt = 9.5
            pHsen = 2.6
        CASE ("PHO") !! PHOSPHATASES
            pHopt = 6.0
            pHsen = 2.0
        CASE ("ENZ") !! ENZ for MOM (Mineral-Associated Organic Matter)
            pHopt = 4.8
            pHsen = 1.6
        CASE DEFAULT
            pHopt = 6.0 !!mean pH of 763 soil samples
            pHsen = 2.0
    END SELECT
    
    fpH = fpH0(pH,pHopt,pHsen)
    return
END FUNCTION fpH
!------------------------------------------------------------------
REAL(r8) FUNCTION fpH0(pH, pHopt, pHsen)
    !!ARGUMENTS:
    REAL(r8), intent(in) :: pH !pH value
    REAL(r8), intent(in) :: pHopt !optimum pH
    REAL(r8), intent(in) :: pHsen !pH sensitivity
    fpH0 = exp(-1.0 * ((pH - pHopt)/pHsen)**2.0)
    return
END FUNCTION fpH0
!!-----------------------------------------------------------------
!! pH Scalar: END
!!=================================================================

!!=================================================================
!! ISOTOPES: BEGIN
!!-----------------------------------------------------------------
REAL(r8) FUNCTION fPermil(iOpt, Rstd, iso1, iso2)
    !Convert isotope concentration to signature/abundance [‰]
    !Rastetter et al. 2005. Ecological Applications 15, 1772-1782.
    !!ARGUMENTS:
    INTEGER, intent(in) :: iOpt !option, 0-input concentration of both iso1 and iso2, otherwise, iso2 = ratio of iso2/iso1
    REAL(r8), intent(in) :: Rstd !standard ratio of iso2/iso1, e.g., C14/C12 = 1d-12, C13/C12 = 0.0112372
    REAL(r8), intent(in) :: iso1 !concentration of iso1
    REAL(r8), intent(in) :: iso2 !concentration of iso2, or iso2/iso1
    if (iOpt .eq. 0) then
        fPermil = ((iso2/iso1)/Rstd - 1.0) * 1000.0
    else
        fPermil = (iso2/Rstd - 1.0) * 1000.0
    end if
    return
END FUNCTION fPermil
!!-----------------------------------------------------------------
!! ISOTOPES: END
!!=================================================================
        
END MODULE MOD_MEND
