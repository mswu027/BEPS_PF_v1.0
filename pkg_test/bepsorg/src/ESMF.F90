

! TBH:  This version is for use with the ESMF library embedded in the WRF 
! TBH:  distribution.  
MODULE ESMF
   USE ESMF_AlarmMod
   USE ESMF_BaseMod
   USE ESMF_BaseTimeMod
   USE ESMF_CalendarMod
   USE ESMF_ClockMod
   USE ESMF_FractionMod
   USE ESMF_TimeIntervalMod
   USE ESMF_TimeMod
   USE ESMF_ShrTimeMod
   USE ESMF_AlarmClockMod
   USE ESMF_Stubs   ! add new dummy interfaces and typedefs here as needed
   USE MeatMod
























! Note that MAX_ALARMS must match MAX_WRF_ALARMS defined in 
! ../../frame/module_domain.F !!!  Eliminate this dependence with 
! grow-as-you-go AlarmList in ESMF_Clock...  


   INTEGER, PARAMETER :: ESMF_MAX_ALARMS=60
!
END MODULE ESMF
