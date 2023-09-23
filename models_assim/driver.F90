!******************************************************
!! main program for BEPS
!! Editted : J.Wang
!! Date    : June2017
!******************************************************

program main
use shr_kind_mod, only: r8 =>shr_kind_r8
use controlInput_mod
use beps_time_manager
use bepstype
use bepstypeInit
use ecoRespMod
use mid_results
!use mpi_mod
use beps_soilMod
use beps_cropMod
use beps_par
use outputMod
use PF_run
use restart
use esmf
implicit none

type(climatedata)    :: meteo
type(results)        :: mid_res
real(r8)             :: CosZs,hr_loc,hr_arc             !! solar zenith angle, local time, local time arc
integer              :: i,j,k,llll,jj,ii,kk
type(soil)           :: soilp                !! at single point
real(r8)             :: Ccd(0:4),Cssd(0:4),Csmd(0:4),Cfsd(0:4),Cfmd(0:4),Csm(0:4),Cm(0:4),Cs(0:4),Cp(0:4)
real(r8)             :: param(0:49),var_o(0:40),var_n(0:40),coef(0:49)
real(r8)             :: inparticles(200,20) ! 19 parameters for optimization, 20th is the pf weight
real(r8)             :: pfweightnormalize(200)
real(r8)             :: lai
type(para),pointer   :: ppar                  !! initial parameters
type(surf),pointer   :: bfields               !! boundary fields
type(forc),pointer   :: climate               !! climate fieldsii =
type(soils),pointer  :: psoil                 !! GLobally
type(res),pointer    :: pp                    !! for output
type(PF_para),pointer   :: PF_ppar            !! Particle Fliter parameters
type(PF_obs0),pointer   :: PF_obsva           !! Particle Fliter observation at nd hour
type(PF_resample0),pointer   :: PF_resa       !! Particle Fliter resample
type(param_vars)      :: p_param      !! C4 Crop from Xiaorong ---->Xiuli
type(Phen)            :: p_phenp            !! C4 Crop from Xiaorong ---->Xiuli
type(param_gdd),pointer   :: p_pgdd        !! C4 Crop from Xiaorong ---->Xiuli


real(r8)             :: ratio_cloud,shortRad_df,shortRad_dir,PFweightupdatesum,SS,PFweightupdatesum_resample
!real(r8)             :: NPP_yr_acc(npoints,PFT)            !! for storing yearly accumulated NPP, Mh/ha/s
! real(r8)             :: agb2vod
! real(r8)             :: D0(1:9)
! real(r8)             :: taweff(1:9)
integer              :: kount,rst_nstep,p,run_pf,p1,p2,p3,p4,p5
integer              :: yr,mn,dy,tod,caldy,n_meteo,n_lai, n_vod_pf !n_gpp_pf
!-- iLab::converted from real(r8) to integer
integer              :: doys
integer              :: ierr
real(r8)             :: daylen
!-- iLab::for revised interface to 'av_output'
character(len=len('YYYY-MM-DDTHH:MM:SS')), save :: ref_date = ''
integer :: yr_ref, mn_ref, dy_ref, tod_ref
real(r8)             :: secs_elapsed,secs_meteo,days_lai, secs_vod_pf ! secs_gpp_pf
! .. Parameters for model-internal use
real(r8)             :: pio180 = PI/180.
! variables for daily input, used in climin and getmonth
real(r8)             :: spds, cpds
! .. Local Arrays
real(r8)             :: atmean, atrange
! .. Local Scalars ..^M
real(r8)             :: r
real(r8)             :: rdaymid, delta, arg, h0, h1, sd, sd1, dhour, tmin, tmp1
real(r8)             :: a, b, sunset_arc
integer              :: nd                   !! for counting the time-step number, i.e. ith step

! ..................... related to crop module .........................................................
real(r8)             :: temp_gpp,outGPP   		!crop module
real(r8)             :: temp_npp,outNPP   		!crop module
real(r8)             :: temp_accu_gpp  !crop module    accumulate  GPP(kg m-2 d-1)
real(r8)             :: temp_accu_npp  !crop module    accumulate  NPP(kg m-2 d-1)
real(r8)             :: temp_accu_temp !crop module    accumulate  temperature(degree d-1)
logical :: is_end_day
logical :: is_end_week
!.........................................................................................

! .. Intrinsic Functions ..
intrinsic ACOS, COS, SIN, MOD,atan, REAL,int
! parameters used for calculating VOD,@MOUSONG.WU,2019-11
!NPP_yr_acc(:,:) = 0.
!agb2vod = 0.9517
!D0 = (/0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05/)
!taweff = (/0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006/)

!! setting up MPI enviroments with the namelists
!   call Initmpi()

call rdnamelist()
if(nscale == 1) then
  nlp = n_site
  npoints = nlp
  write(*,*) 'site points check', npoints
end if

!! Initialize the beps types
call Initbepstype()

!! Initialize output
call Init_output
call Init_output_PF

!--------------------------------------------------------------------------
! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
! because it is needed for the time manager, even if the ESMF_INTERFACE
! is not used.
!--------------------------------------------------------------------------
call ESMF_Initialize()

!! setting time manager
if(nsrest == nsrStartup) then
   !!! calling time_manager set init
    call set_timemgr_init(calendar_in  = calendar,&
                         start_ymd_in = icdate  ,&
                         start_tod_in = icsec   ,&
                         nelapse_in   = sim_duration,&
                         dtime_in     = step)

else if(nsrest == nsrContinue) then
    call restart_io("read")

    call set_timemgr_init(calendar_in  = calendar,   &
                            start_ymd_in = rst_icdate, &
                            start_tod_in = rst_icsec , &
                            nelapse_in   = sim_duration,&
                            dtime_in     = step)
end if

call timemgr_init()

if (nscale == 0) then     ! nscale = 0 for global simulation, 1 for site simulation
   print *, 'BEPS run at global scale!'
    !! Reading boundary fields, yearly data,and soil Cpools for BEPS @J.Wang (note: if yearly and C pools data fields will change year by year, these datasets should be read in the time looping
   call read_boundary()
   call read_yrdata()
   call read_cpools()
else
    print *, 'BEPS run at site scale!'
    call read_boundary_site()   ! read site data, including yrdata, boundary data, and carbon pools
    print *, 'read boundary site successfully!'
end if

bfields  => bound
climate  => clim
psoil    => soilstat
pp       => output
ppar     => assim
PF_ppar  => PF
PF_obsva => PF_obs
PF_resa  => PF_resample
p_pgdd   => pgdd

! 0 for forwar model; 1 for praticle filter
run_pf=1
PFweightupdatesum=0.
PFweightupdatesum_resample=0.

temp_gpp=0.
outGPP=0.
temp_npp=0.
outNPP=0.
temp_accu_gpp=0.
temp_accu_npp=0.
temp_accu_temp=0.
!..................................................

 !**********************************************************************************run BEPS with default parameters parameters**************************************************
if (run_pf==0) then

 call read_prior_para()       ! put the parameters to be optimized in a NETCDF
                                !file and read them as well their
                                !uncertainties,@MOUSONG.WU,2019-11
 print *, 'read prior para successfully!'

! do p =1,10 !! start parameter loop
! do p =1,10!!1,nparameters !! start parameter loop
!  write(*,*) 'p=',p
  !! setting time manager
  !!! calling time_manager set init
!  call set_timemgr_init(calendar_in  = calendar,&
!                         start_ymd_in = icdate  ,&
!                         start_tod_in = icsec   ,&
!                         nelapse_in   = sim_duration,&
!                         dtime_in     = step)

!  call timemgr_init()

  nd = 0
!  yr=2013
!  mn=1
!  dy=1
!  tod=0
  do     !! start time looping
    nd = nd + 1
    !-- iLab::inserted in order to pass 'calday' downstream
    caldy = get_curr_calday()
    call get_curr_date(yr,mn,dy,tod)
    doys = get_doys(yr)
    write(*,*) 'nd=', nd
    !-- iLab::inserted in order to pass ref_date downstream
    if( nd.eq.1 ) then
        yr_ref = yr
        mn_ref = mn
        dy_ref = dy
        tod_ref = tod
        write(ref_date(1:4),   '(i4.4)') yr_ref
        write(ref_date(5:5),   '(a)')    '-'
        write(ref_date(6:7),   '(i2.2)') mn_ref
        write(ref_date(8:8),   '(a)')    '-'
        write(ref_date(9:10),  '(i2.2)') dy_ref
        write(ref_date(11:19), '(a)')    'T00:00:00'
    !else
    !   call get_curr_date(yr,mn,dy,tod)
    !   !write(*,*) 'yr1=', yr
    !   !write(*,*) 'mn1=', mn
    !   !write(*,*) 'dy1=', dy
    !   !write(*,*) 'tod1=', tod
    end if
    !write(*,*) 'date = ',yr,mn,dy,tod
    !-- iLab::determine seconds elapsed since reference time
    call timemgr_diff_secs(yr_ref*10000+mn_ref*100+dy_ref, tod_ref, yr*10000+mn*100+dy, tod,&
            secs_elapsed)
    call get_CO2_concentration(yr,CO2_air)
    call get_COS_concentration(yr,COS_air)
    !! change hourly input into daily input for further using model for long-term simulations, @MOUSONG.WU, 201905
    if (meteo_input >= 0) then  ! call hourly meteo. input
      if (nscale == 0) then
          call read_meteo_hourly(yr, mn, dy, tod)
      else
          call timemgr_diff_secs(2001*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            secs_meteo)
          n_meteo = int(secs_meteo/3600 + 1)
          call read_meteo_site(n_meteo)
          print *, 'read site meteo successfully!'
      end if
    else
      if(is_first_step() .or. is_end_curr_day()) then
        call read_meteo_daily(yr, mn, dy, tod)
      end if
    end if

    if (lai_input >=0) then
	    call timemgr_datediff(2001*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            days_lai)
        n_lai = int(days_lai)
        if (is_first_step()) print *, 'lai is input!'
        if (is_first_step() .or. is_end_curr_day()) then
          if (nscale == 0) then
             call read_lai(yr, mn, dy, tod, n_lai)
          else
             call read_lai_site(n_lai)
             print *, 'read site lai successfully!'
          end if
        end if
    else
        if (is_first_step()) then
           print *, 'lai is simulated with phenology scheme!'
        end if
    end if

    !call mpi_barrier(mpi_comm_world,ierr)

       do i = 1,npoints    !! spatial iteration

          !! calculate the solar zenith
          call s_coszs(yr, mn, dy, tod, caldy, doys, &
               bfields%latitude(i),bfields%longitude(i),CosZs,hr_loc,hr_arc)
          !!if(myid == 0) write(*,*) "hr_loc=",hr_loc
          !! retrieve meteo for this point

          meteo%Srad                   = climate%Srad(i)
          meteo%wind                   = climate%Wind(i)
          meteo%rainfall               = climate%Rain(i)
          meteo%snow                   = climate%Snow(i)
          meteo%rh                     = climate%Rh(i)

          ! .. compute daily course of temperature and daylength
          rdaymid = REAL (sim_duration+1) / 2.
          delta = -23.4*COS(2.*PI*(rdaymid+10.)/365.)
          spds = SIN(bfields%latitude(i)*pio180)*SIN(delta*pio180)
          cpds = COS(bfields%latitude(i)*pio180)*COS(delta*pio180)
          arg = -spds/cpds
          IF (arg>1.) THEN
              !polar night:
              daylen = 0.
          ELSE IF (arg<-1) THEN
              !polar day:
              daylen = 24.
          ELSE
              !normal day / night:
              daylen = ACOS(arg)/PI*24.
          END IF

          ! Calculate cloud fraction, separate shortwave radiation
          if(CosZs < 0.001) then
              ratio_cloud=0.
          else
              ratio_cloud=meteo%Srad/(1367.*CosZs)
          end if

          if (ratio_cloud > 0.8) then
              shortRad_df  = 0.13*meteo%Srad
          else
              shortRad_df  = (0.943+0.734*ratio_cloud-4.9*ratio_cloud**2+1.796*ratio_cloud**3+2.058*ratio_cloud**4)*&
                                   & meteo%Srad
          end if
          shortRad_df   = min(shortRad_df,meteo%Srad)
          shortRad_df   = max(shortRad_df,0.)

          shortRad_dir  = meteo%Srad - shortRad_df

          meteo%S_dff   = shortRad_df
          meteo%S_dir   = shortRad_dir
          climate%Swdr(i)   = shortRad_dir
          climate%Swdf(i)   = shortRad_df

          !write(*,*) 'Before start_PFT_loop!!!!!'
          do j = 11,11!1,PFT    !! PFT iteration  ! 2023/06/30 beps land cover id
               !write(*,*) 'start_PFT_loop!!!!'
               if(bfields%lcno(i,j) > 0 .and. bfields%sw(i) >0. .and. bfields%stext(i) >0 .and. bfields%clumping(i) > 0.5) then
                    call readparam(bfields%lcno(i,j),param)

                    if(lai_input >= 0) then
                       lai = bfields%lai(i,j)
                    else
                       if (is_first_step()) then
                          lai = bfields%laiyr(i,j)
                          mid_res%lai_old = lai
                       else
                          lai = mid_res%lai_new
                       end if
                    end if

                    lai = lai*param(2)/bfields%clumping(i)
                    call readcoef(bfields%lcno(i,j),bfields%stext(i),coef)

                    if(nsrest == nsrStartup .and. is_first_step()) then
                        !write(*,*) 'Startup!!!!'
		        ! 2023/06/30 ????????????
                        call Init_soil_parameters(bfields%lcno(i,j),bfields%laiyr(i,j),bfields%stext(i), &
                                            ppar%p_f_decay(j,i), param(27),soilp) ! 2023/06/30
                        soilp%r_drainage = param(26)
                        !soilp%r_drainage = ppar%p_drainage(j)            ! read this
                                                                    !para from NC file,@MOUSONG.WU,2019-11
                        ii = bfields%stext(i)

                        !write(*,*) 'Ksat_old = ',soilp%Ksat(0)
                        do kk = 0,4
                            soilp%Ksat(kk) = ppar%p_Ksat_scalar(ii,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = ppar%p_b_scalar(ii,i)*soilp%b(kk)
                        end do

                        ! replace these three para. above with values from nc
                        ! file,@MOUSONG.WU,2019-11
                        call Init_soil_status(soilp,bfields%st(i),climate%Temp(i),bfields%sw(i),bfields%sdp(i))

                        do k = 0,40
                            var_o(k)   = 0.
                        end do

                        do k = 3,8
                           var_o(k)   = climate%Temp(i)
                        end do

                        do k = 9,14
                           var_o(k)   = soilp%temp_soil_p(k-9)
                        end do

                        do k = 21,26
                           var_o(k)   = soilp%thetam_prev(k-21)
                        end do

                        do k = 27,32
                           var_o(k)    = soilp%ice_ratio(k-27)
                        end do

                    else if(nsrest == nsrContinue .and. is_first_step()) then
                         !write(*,*) 'Continue!!!!'
                        call Init_soil_parameters(bfields%lcno(i,j),bfields%laiyr(i,j),bfields%stext(i), &
                                      ppar%p_f_decay(j,i), param(27),soilp) ! 2023/06/30

                        soilp%r_drainage = param(26)
                        !soilp%r_drainage = ppar%p_drainage(j)  ! read from nc
                                                          ! file,MOUSONG.WU@2019-11
                        ii = bfields%stext(i)
                        do kk = 0,4
                            soilp%Ksat(kk) = ppar%p_Ksat_scalar(ii,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = ppar%p_b_scalar(ii,i)*soilp%b(kk)
                        end do
                        !Replace with NC values,for optimization
                        !purpose,@MOUSONG.WU,2019-11
                        var_o(:)         = v2last(i,:,j)
                        do k= 9,14
                           soilp%temp_soil_c(k-9)  = var_o(k)
                        end do
                        do k= 21,26
                           soilp%thetam(k-21)      = var_o(k)
                        end do
                        !do k = 27,32
                        !soilp%ice_ratio(k-27)   = var_o(k)
                        !end do

                    else
                        !write(*,*) 'others!!!!'
                        var_o(:)         = v2last(i,:,j)
                        call retrive_soilp(soilp,i,j,0)
                    end if
                    !write(*,*) 'layer = ', soilp%n_layer

                    call inter_prg(yr, mn, dy, tod, &
                         lai,lai_input,bfields%lcno(i,j),bfields%clumping(i),ppar%p_Vcmax(j,i),ppar%p_VJ_slope(j,i),&
                         ppar%p_VN_slope(j,i),ppar%p_b_h2o(j,i),ppar%p_m_h2o(j,i),ppar%p_f_leaf(j,i),&
                         ppar%p_kc25(j,i),ppar%p_ko25(j,i),ppar%p_tau25(j,i),ppar%p_sif_alpha(j,i),ppar%p_sif_beta(j,i),&
                         ppar%p_a(j,i),ppar%p_b(j,i),ppar%p_c(j,i),&
                         param,meteo,CosZs,var_o,var_n,soilp,mid_res,daylen)

                    !print *, 'end of inter_prg'

                    v2last(i,:,j)  = var_n(:)
                    call retrive_soilp(soilp,i,j,1)
                    !print *, 'end of retrive_soilp'
                    !!! simluating Ra
                    call plant_resp(ppar%p_q10(j,i),bfields%lcno(i,j),mid_res,bfields%laiyr(i,j),lai,meteo%temp,&
                               soilp%temp_soil_c(1),CosZs)

                    !print *, 'end of plant_resp'
                    !USE p_q10 here to adjust q10,p_q10 is read from initial para. NC file, for
                    !optimization purpose,@MOUSONG.WU,2019-11
                    !!! simulating Rh
                    Ccd(0)       = bfields%Ccd(i,j)
                    Cssd(0)      = bfields%Cssd(i,j)
                    Csmd(0)      = bfields%Csmd(i,j)
                    Cfsd(0)      = bfields%Cfsd(i,j)
                    Cfmd(0)      = bfields%Cfmd(i,j)
                    Csm(0)       = bfields%Csm(i,j)
                    Cm(0)        = bfields%Cm(i,j)
                    Cs(0)        = bfields%Cm(i,j)
                    Cp(0)        = bfields%Cp(i,j)
                    ! to get soil texture for this point,@MOUSONG.WU,2019-11
                    jj = bfields%stext(i)

                    call soil_resp(ppar%p_f_resp(j,i),Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,&
                                 Cp,bfields%nppyr(i,j),coef,bfields%stext(i),soilp,mid_res)

                    !print *, 'end of soil_resp'

                    !! for output variables
                    pp%GPPpft(i,j)   = mid_res%GPP*bfields%PCT_PFT(i,j)/100.
                    pp%SIFpft(i,j)   = mid_res%SIF*bfields%PCT_PFT(i,j)/100.
                    pp%NPPpft(i,j)   = mid_res%NPP*bfields%PCT_PFT(i,j)/100.
                    pp%NEPpft(i,j)   = mid_res%NEP*bfields%PCT_PFT(i,j)/100.
                    pp%SHpft(i,j)    = mid_res%SH*bfields%PCT_PFT(i,j)/100.
                    pp%LHpft(i,j)    = mid_res%LH*bfields%PCT_PFT(i,j)/100.
                    pp%Transpft(i,j) = mid_res%Trans*bfields%PCT_PFT(i,j)/100.
                    pp%Evappft(i,j)  = mid_res%Evap*bfields%PCT_PFT(i,j)/100.
                    pp%Net_Radpft(i,j)  = mid_res%Net_Rad*bfields%PCT_PFT(i,j)/100.
                    pp%LAIpft(i,j)   = lai*bfields%PCT_PFT(i,j) /100.
                    pp%Thetampft(i,j)   = mid_res%thetam_surf*bfields%PCT_PFT(i,j)/100.
                    pp%fAPARpft(i,j) = mid_res%fAPAR*bfields%PCT_PFT(i,j)/100.
                    pp%COS_fluxpft(i,j) = (mid_res%COS_plant+mid_res%COS_grnd)*bfields%PCT_PFT(i,j)/100.
                    pp%PWSpft(i,j) =  mid_res%PWS*bfields%PCT_PFT(i,j)/100.
                    pp%ETapft(i,j) = mid_res%ETa*bfields%PCT_PFT(i,j)/100.
                    pp%VODpft(i,j) = mid_res%VOD*bfields%PCT_PFT(i,j)/100.
                    pp%fei_leafpft(i,j) = mid_res%fei_leaf*bfields%PCT_PFT(i,j)/100.
                    !write(*,*) 'pp%VODpft=', pp%VODpft(i,j)
                    if(hr_loc >= 13. .and. hr_loc <14.) then
                        pp%SIFpft_sat(i,j)  = mid_res%SIF*bfields%PCT_PFT(i,j)/100.
                    else
                        pp%SIFpft_sat(i,j)  = 0.
                    end if
                    !pp%SIFpft_sat(i,j) = max(pp%SIFpft_sat(i,j),0.)
                    !write(*,*) 'SIFpft_sat = ', pp%SIFpft_sat(i,j)
                end if
          end do   !! end PFT loop

          pp%GPP(i)   = sum(pp%GPPpft(i,:))
          pp%SIF(i)   = sum(pp%SIFpft(i,:))
          pp%SIF_sat(i)  = sum(pp%SIFpft_sat(i,:))
          pp%NPP(i)   = sum(pp%NPPpft(i,:))
          pp%NEP(i)   = sum(pp%NEPpft(i,:))
          pp%SH(i)    = sum(pp%SHpft(i,:))
          pp%LH(i)    = sum(pp%LHpft(i,:))
          pp%Trans(i) = sum(pp%Transpft(i,:))
          pp%Evap(i)  = sum(pp%Evappft(i,:))
          pp%Net_Rad(i) = sum(pp%Net_Radpft(i,:))
          pp%LAI(i)     = sum(pp%LAIpft(i,:))
          pp%Thetam(i)  = sum(pp%Thetampft(i,:))
          pp%fAPAR(i)   = sum(pp%fAPARpft(i,:))
          pp%VOD(i)     = sum(pp%VODpft(i,:))
          pp%COS_flux(i)     = sum(pp%COS_fluxpft(i,:))
          pp%PWS(i)     = sum(pp%PWSpft(i,:))
          pp%ETa(i)     = sum(pp%ETapft(i,:))
          pp%fei_leaf(i)     = sum(pp%fei_leafpft(i,:))

       end do      !! end spatial loop
       !          call mpi_barrier(mpi_comm_world,ierr)
       !! advance time
       call advance_timestep()

       !!! write oout data and restart fields
       !call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed,p)
       call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed)
       !-- iLab::restart requires next time-step
       call get_curr_date(yr,mn,dy,tod)

       if(is_last_step()) exit
  end do   !! end time loop


 !*************************************************************************************Particle Fliter ****************************************************************************
else if (run_pf==1) then

 write(*,*) 'run_pf=',run_pf

 !call read_prior_para()       ! put the parameters to be optimized in a NETCDF
                                !file and read them as well their
                                !uncertainties,@MOUSONG.WU,2019-11

 call Create_PF_para(11,8,parloop) !"1" is PFT,"2" is soil texture !!create particles. xiuli
      !***************************************************************************************Main program***********************************************************
 nd = 0
 do     !! start time looping
    nd = nd + 1
    write(*,*) 'nd=',nd
    !-- iLab::inserted in order to pass 'calday' downstream
    caldy = get_curr_calday()
    call get_curr_date(yr,mn,dy,tod)
    doys = get_doys(yr)
    !-- iLab::inserted in order to pass ref_date downstream
    if( nd.eq.1 ) then
        yr_ref = yr
        mn_ref = mn
        dy_ref = dy
        tod_ref = tod
        write(ref_date(1:4),   '(i4.4)') yr_ref
        write(ref_date(5:5),   '(a)')    '-'
        write(ref_date(6:7),   '(i2.2)') mn_ref
        write(ref_date(8:8),   '(a)')    '-'
        write(ref_date(9:10),  '(i2.2)') dy_ref
        write(ref_date(11:19), '(a)')    'T00:00:00'
    end if
    !write(*,*) "yr=" , yr
    !write(*,*) "mn=" , mn
    !write(*,*) "dy=" , dy
    !write(*,*) "tod=" , tod
    !write(*,*) "caldy=" , caldy
    is_end_day = ((yr>yr_ref .or. mn>mn_ref .or. dy>dy_ref) .and. (tod == 0)) !--iLab: taken from BEPS time manager
    is_end_week = ((yr>yr_ref .or. mn>mn_ref .or. dy>dy_ref) .and. (MOD((dy-dy_ref+1),7)==0) .and. (tod == 82800))
    !-- iLab::determine seconds elapsed since reference time
    call timemgr_diff_secs(yr_ref*10000+mn_ref*100+dy_ref, tod_ref, yr*10000+mn*100+dy, tod,&
            secs_elapsed)
    call get_CO2_concentration(yr,CO2_air)
    call get_COS_concentration(yr,COS_air)
    !! change hourly input into daily input for further using model for long-term simulations, @MOUSONG.WU, 201905
    if (meteo_input >= 0) then  ! call hourly meteo. input
      if (nscale == 0) then
          call read_meteo_hourly(yr, mn, dy, tod)
      else
          call timemgr_diff_secs(2001*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            secs_meteo)
          n_meteo = int(secs_meteo/3600 + 1)
          call read_meteo_site(n_meteo)
          print *, 'read site meteo successfully!'
      end if
    else
      if(is_first_step() .or. is_end_curr_day()) then
        call read_meteo_daily(yr, mn, dy, tod)
      end if
    end if

    if (lai_input >=0) then
        call timemgr_datediff(2001*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            days_lai)
        n_lai = int(days_lai)
        if (is_first_step()) print *, 'lai is input!'
        if (is_first_step() .or. is_end_curr_day()) then
          if (nscale == 0) then
             call read_lai(yr, mn, dy, tod, n_lai)
          else
             write(*,*) "n_lai=", n_lai
             call read_lai_site(n_lai)
             print *, 'read site lai successfully!'
          end if
        end if
    else
        if (is_first_step()) then
           print *, 'lai is simulated with phenology scheme!'
        end if
    end if


    !call mpi_barrier(mpi_comm_world,ierr)
    call timemgr_diff_secs(2001*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            secs_vod_pf)
    ! n_gpp_pf = int(secs_gpp_pf/3600 + 1)
    n_vod_pf = int(secs_vod_pf/3600 + 1)
    call read_PF_obs(n_vod_pf) !xiuli

    !if (PF_obsva%obs_GPP(1) == PF_obsva%obs_GPP(1)) then !if observation exsits xiuli
    if (PF_obsva%obs_VOD(1) == PF_obsva%obs_VOD(1)) then
     do p1 =1,parloop !! start particles loop xiuli
     !do p =1,10!!1,nparameters !! start parameter loop

       write(*,*) 'p1=',p1
       do i = 1,npoints    !! spatial iteration

          !! calculate the solar zenith
          call s_coszs(yr, mn, dy, tod, caldy, doys, &
               bfields%latitude(i),bfields%longitude(i),CosZs,hr_loc,hr_arc)
          !!if(myid == 0) write(*,*) "hr_loc=",hr_loc
          !! retrieve meteo for this point

          meteo%Srad                   = climate%Srad(i)
          meteo%wind                   = climate%Wind(i)
          meteo%rainfall               = climate%Rain(i)
          meteo%snow                   = climate%Snow(i)
          meteo%rh                     = climate%Rh(i)


          ! .. compute daily course of temperature and daylength
          rdaymid = REAL (sim_duration+1) / 2.
          delta = -23.4*COS(2.*PI*(rdaymid+10.)/365.)
          spds = SIN(bfields%latitude(i)*pio180)*SIN(delta*pio180)
          cpds = COS(bfields%latitude(i)*pio180)*COS(delta*pio180)
          arg = -spds/cpds
          IF (arg>1.) THEN
              !polar night:
              daylen = 0.
          ELSE IF (arg<-1) THEN
              !polar day:
              daylen = 24.
          ELSE
              !normal day / night:
              daylen = ACOS(arg)/PI*24.
          END IF

                   ! Calculate cloud fraction, separate shortwave radiation
          if(CosZs < 0.001) then
              ratio_cloud=0.
          else
              ratio_cloud=meteo%Srad/(1367.*CosZs)
          end if

          if (ratio_cloud > 0.8) then
              shortRad_df  = 0.13*meteo%Srad
          else
              shortRad_df  = (0.943+0.734*ratio_cloud-4.9*ratio_cloud**2+1.796*ratio_cloud**3+2.058*ratio_cloud**4)*&
                                   & meteo%Srad
          end if
          shortRad_df   = min(shortRad_df,meteo%Srad)
          shortRad_df   = max(shortRad_df,0.)

          shortRad_dir  = meteo%Srad - shortRad_df

          meteo%S_dff   = shortRad_df
          meteo%S_dir   = shortRad_dir
          climate%Swdr(i)   = shortRad_dir
          climate%Swdf(i)   = shortRad_df

          !write(*,*) 'Before start_PFT_loop!!!!!'
          do j = 11,11!1,PFT    !! PFT iteration  ! 2023/06/30 beps land cover id
               !write(*,*) 'start_PFT_loop!!!!'
               if(bfields%lcno(i,j) > 0 .and. bfields%sw(i) >0. .and. bfields%stext(i) >0 .and. bfields%clumping(i) > 0.5) then
                    call readparam(bfields%lcno(i,j),param)

                    if(lai_input >= 0) then
                       lai = bfields%lai(i,j)
                    else
                       if (is_first_step()) then
                          lai = bfields%laiyr(i,j)
                          mid_res%lai_old = lai
                       else
                          lai = mid_res%lai_new
                       end if
                    end if

                    lai = lai*param(2)/bfields%clumping(i)
                    call readcoef(bfields%lcno(i,j),bfields%stext(i),coef)

                    if(nsrest == nsrStartup .and. is_first_step()) then
                        !write(*,*) 'Startup!!!!'
		        ! 2023/06/30 ????????????
                        call Init_soil_parameters(bfields%lcno(i,j),bfields%laiyr(i,j),bfields%stext(i), &
                                            PF_ppar%f_decay(p1,i), param(27),soilp) ! 2023/06/30
                        soilp%r_drainage = param(26)
                        !soilp%r_drainage = ppar%p_drainage(j)            ! read this
                                                                    !para from NC file,@MOUSONG.WU,2019-11
                        ii = bfields%stext(i)

                        !write(*,*) 'Ksat_old = ',soilp%Ksat(0)
                        do kk = 0,4
                            soilp%Ksat(kk) = PF_ppar%Ksat_scalar(p1,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = PF_ppar%b_scalar(p1,i)*soilp%b(kk)
                        end do

                        ! replace these three para. above with values from nc
                        ! file,@MOUSONG.WU,2019-11
                        call Init_soil_status(soilp,bfields%st(i),climate%Temp(i),bfields%sw(i),bfields%sdp(i))

                        do k = 0,40
                            var_o(k)   = 0.
                        end do

                        do k = 3,8
                           var_o(k)   = climate%Temp(i)
                        end do

                        do k = 9,14
                           var_o(k)   = soilp%temp_soil_p(k-9)
                        end do

                        do k = 21,26
                           var_o(k)   = soilp%thetam_prev(k-21)
                        end do

                        do k = 27,32
                           var_o(k)    = soilp%ice_ratio(k-27)
                        end do

                    else if(nsrest == nsrContinue .and. is_first_step()) then
                         !write(*,*) 'Continue!!!!'
                        call Init_soil_parameters(bfields%lcno(i,j),bfields%laiyr(i,j),bfields%stext(i), &
                                      PF_ppar%f_decay(p1,i), param(27),soilp) ! 2023/06/30

                        soilp%r_drainage = param(26)
                        !soilp%r_drainage = ppar%p_drainage(j)  ! read from nc
                                                          ! file,MOUSONG.WU@2019-11
                        ii = bfields%stext(i)
                        do kk = 0,4
                            soilp%Ksat(kk) = PF_ppar%Ksat_scalar(p1,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = PF_ppar%b_scalar(p1,i)*soilp%b(kk)
                        end do
                        !Replace with NC values,for optimization
                        !purpose,@MOUSONG.WU,2019-11
                        var_o(:)         = v2last(i,:,j)
                        do k= 9,14
                           soilp%temp_soil_c(k-9)  = var_o(k)
                        end do
                        do k= 21,26
                           soilp%thetam(k-21)      = var_o(k)
                        end do
                        !do k = 27,32
                        !soilp%ice_ratio(k-27)   = var_o(k)
                        !end do

                    else
                        !write(*,*) 'others!!!!'
                        var_o(:)         = v2last(i,:,j)
                        call retrive_soilp(soilp,i,j,0)
                    end if
                    !write(*,*) 'layer = ', soilp%n_layer

                    call inter_prg(yr, mn, dy, tod, &
                         lai,lai_input,bfields%lcno(i,j),bfields%clumping(i),PF_ppar%Vcmax(p1,i),PF_ppar%VJ_slope(p1,i),&
                         PF_ppar%VN_slope(p1,i),PF_ppar%b_h2o(p1,i),PF_ppar%m_h2o(p1,i),PF_ppar%f_leaf(p1,i),&
                         PF_ppar%kc25(p1,i),PF_ppar%ko25(p1,i),PF_ppar%tau25(p1,i),PF_ppar%sif_alpha(p1,i),PF_ppar%sif_beta(p1,i),&
                         PF_ppar%a(p1,i),PF_ppar%b(p1,i),PF_ppar%c(p1,i),param,meteo,CosZs,var_o,var_n,soilp,mid_res,daylen)
                    !print *, 'end of inter_prg'
                    ! CHANGE Vcmax read from Vcmax file with the Vcmax read from initial para. NC
                    ! file, for optimization purpose,@MOUSONG.WU,2019-11
                    !              do llll =0,40
                    !                write(*,*)  "DG004: Var_n = ",llll,var_n(llll)
                    !             end do
                    v2last(i,:,j)  = var_n(:)
                    call retrive_soilp(soilp,i,j,1)
                    !write(*,*) 'porosity = ', soilp%fei(0)
                    !print *, 'end of retrive_soilp'
                    !!! simluating Ra
                    call plant_resp(PF_ppar%q10(p1,i),bfields%lcno(i,j),mid_res,bfields%laiyr(i,j),lai,meteo%temp,&
                               soilp%temp_soil_c(1),CosZs)
                    !print *, 'end of plant_resp'
                    !USE p_q10 here to adjust q10,p_q10 is read from initial para. NC file, for
                    !optimization purpose,@MOUSONG.WU,2019-11
                    !!! simulating Rh
                    Ccd(0)       = bfields%Ccd(i,j)
                    Cssd(0)      = bfields%Cssd(i,j)
                    Csmd(0)      = bfields%Csmd(i,j)
                    Cfsd(0)      = bfields%Cfsd(i,j)
                    Cfmd(0)      = bfields%Cfmd(i,j)
                    Csm(0)       = bfields%Csm(i,j)
                    Cm(0)        = bfields%Cm(i,j)
                    Cs(0)        = bfields%Cm(i,j)
                    Cp(0)        = bfields%Cp(i,j)
                    ! to get soil texture for this point,@MOUSONG.WU,2019-11
                    jj = bfields%stext(i)

                    call soil_resp(PF_ppar%f_resp(p1,i),Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,&
                                 Cp,bfields%nppyr(i,j),coef,bfields%stext(i),soilp,mid_res)
                    !print *, 'end of soil_resp'

                    if (j>10) then
                        if (mid_res%GPP*1000.*3600. >=4) then
                            outGPP=temp_gpp
                            !write(6,*)   "outGPP=temp_gpp"
                        else
                            outGPP=mid_res%GPP*1000.*3600.
                            !write(6,*)   "outGPP=mid_res%GPP*1000.*3600."
                        end if

                        if (mid_res%NPP*1000.*3600. >=4) then
                            outNPP=temp_npp
                            !write(6,*)   "outNPP=temp_npp"
                        else
                            outNPP=mid_res%NPP*1000.*3600.
                            !write(6,*)   "outNPP=mid_res%NPP*1000.*3600."
                        end if

                        temp_gpp=outGPP
                        !write(6,*)   "temp_gpp=outGPP"
                        temp_npp=outNPP
                        !write(6,*)   "temp_npp=outNPP"

                        temp_accu_gpp=temp_accu_gpp+outGPP
                        !write(6,*)   "temp_accu_gpp=temp_accu_gpp+outGPP"
                        temp_accu_npp=temp_accu_npp+outNPP
                        !write(6,*)   "temp_accu_npp=temp_accu_npp+outNPP"
                        temp_accu_temp=temp_accu_temp+(meteo%temp+273.15)
                        !write(6,*)   "temp_accu_temp=temp_accu_temp+(meteo%temp+273.15)"
                        !write(*,*)   "temp_accu_temp=", temp_accu_temp

                        if (is_end_day) then
                            p_phenp%temp_daily = temp_accu_temp/24.     !crop module	temperature(K d-1)
                            !write(6,*)   "p_phenp%temp_daily = temp_accu_temp/24. "
                            p_phenp%gpp_ft_acc=temp_accu_gpp/1000.		!crop module	GPP(kg m-2 d-1)
                            !write(6,*)   "p_phenp%gpp_ft_acc=temp_accu_gpp/1000."
			                p_phenp%npp_ft_acc=temp_accu_npp/1000.		!crop module	NPP(kg m-2 d-1)
                            !write(6,*)   "p_phenp%npp_ft_acc=temp_accu_npp/1000."

                            if (secs_elapsed/86400.>=p_pgdd%emer_doy(i,j) .and. secs_elapsed/86400.<=p_pgdd%har_doy(i,j)) then

			                    call beps_crop_Development(i,j,p_param,p_phenp) 		!crop module
                                !write(6,*)   "call beps_crop_Development(i,j,p_param,p_phenp)"
			                    call beps_crop_pool_alloc(p_param,p_phenp)                !crop module
                                !write(6,*)   "call beps_crop_pool_alloc(p_param,p_phenp)"

                                temp_accu_temp=0.		!crop module
		                        temp_accu_gpp=0.		!crop module
		                        temp_accu_npp=0.		!crop module
                            end if
                        end if
                    end if


                    !! for output variables
                    pp%GPPpft(i,j)   = mid_res%GPP*bfields%PCT_PFT(i,j)/100.
                    pp%SIFpft(i,j)   = mid_res%SIF*bfields%PCT_PFT(i,j)/100.
                    pp%NPPpft(i,j)   = mid_res%NPP*bfields%PCT_PFT(i,j)/100.
                    pp%NEPpft(i,j)   = mid_res%NEP*bfields%PCT_PFT(i,j)/100.
                    pp%SHpft(i,j)    = mid_res%SH*bfields%PCT_PFT(i,j)/100.
                    pp%LHpft(i,j)    = mid_res%LH*bfields%PCT_PFT(i,j)/100.
                    pp%Transpft(i,j) = mid_res%Trans*bfields%PCT_PFT(i,j)/100.
                    pp%Evappft(i,j)  = mid_res%Evap*bfields%PCT_PFT(i,j)/100.
                    pp%Net_Radpft(i,j)  = mid_res%Net_Rad*bfields%PCT_PFT(i,j)/100.
                    pp%LAIpft(i,j)   = lai*bfields%PCT_PFT(i,j) /100.
                    pp%Thetampft(i,j)   = mid_res%thetam_surf*bfields%PCT_PFT(i,j)/100.
                    pp%fAPARpft(i,j) = mid_res%fAPAR*bfields%PCT_PFT(i,j)/100.
                    pp%COS_fluxpft(i,j) = (mid_res%COS_plant+mid_res%COS_grnd)*bfields%PCT_PFT(i,j)/100.

                    pp%PWSpft(i,j) =  mid_res%PWS*bfields%PCT_PFT(i,j)/100.
                    pp%ETapft(i,j) = mid_res%ETa*bfields%PCT_PFT(i,j)/100.
                    pp%VODpft(i,j) = mid_res%VOD*bfields%PCT_PFT(i,j)/100.
                    pp%fei_leafpft(i,j) = mid_res%fei_leaf*bfields%PCT_PFT(i,j)/100.

                    !pp%NPP_yr_acc(i,j) = pp%NPP_yr_acc(i,j) + mid_res%NPP*bfields%PCT_PFT(i,j)/100.*1.e-2*step       ! convert NPP to Mg/ha for calculation of VOD

                    !if (is_end_curr_year()) then
                    ! calculate VOD (vegetation optical depth) with results derived from SMOS-IC product, @Mousong.Wu, 201905,taweff is a PFT specific parameter
                      ! pp%VODpft(i,j)   = PF_ppar%agb2vod(p1,i)*atan(PF_ppar%taweff(p1,i)*pp%NPP_yr_acc(i,j)) + PF_ppar%D0(p1,i)*lai
                      ! pp%NPP_yr_acc(i,j) = 0.
                    !                  write(*,*) 'VOD= ',pp%VODpft(i,j)
                    !else
                      ! pp%VODpft(i,j) = 0.
                    !end if

                    ! write(*,*) 'hr_loc = ', hr_loc
                    !! calculate the OCO-2 SIF   across at 1:30pm
                    !write(*,*) 'SIFpft= ',mid_res%SIF*bfields%PCT_PFT(i,j)
                    if(hr_loc >= 13. .and. hr_loc <14.) then
                        pp%SIFpft_sat(i,j)  = mid_res%SIF*bfields%PCT_PFT(i,j)/100.
                    else
                        pp%SIFpft_sat(i,j)  = 0.
                    end if
                    !pp%SIFpft_sat(i,j) = max(pp%SIFpft_sat(i,j),0.)
                    !write(*,*) 'SIFpft_sat = ', pp%SIFpft_sat(i,j)
                end if

          end do   !! end PFT loop

          pp%GPP(i)   = sum(pp%GPPpft(i,:))
          pp%SIF(i)   = sum(pp%SIFpft(i,:))
          pp%SIF_sat(i)  = sum(pp%SIFpft_sat(i,:))
          pp%NPP(i)   = sum(pp%NPPpft(i,:))
          pp%NEP(i)   = sum(pp%NEPpft(i,:))
          pp%SH(i)    = sum(pp%SHpft(i,:))
          pp%LH(i)    = sum(pp%LHpft(i,:))
          pp%Trans(i) = sum(pp%Transpft(i,:))
          pp%Evap(i)  = sum(pp%Evappft(i,:))
          pp%Net_Rad(i) = sum(pp%Net_Radpft(i,:))
          pp%LAI(i)     = sum(pp%LAIpft(i,:))
          pp%Thetam(i)  = sum(pp%Thetampft(i,:))
          pp%fAPAR(i)   = sum(pp%fAPARpft(i,:))
          pp%VOD(i)     = sum(pp%VODpft(i,:))
          pp%COS_flux(i)     = sum(pp%COS_fluxpft(i,:))

          pp%PWS(i)     = sum(pp%PWSpft(i,:))
          pp%ETa(i)     = sum(pp%ETapft(i,:))
          pp%fei_leaf(i)     = sum(pp%fei_leafpft(i,:))

        end do      !! end spatial loop

        !weight calculation
        call PF_weight(p1,yr, mn, dy, tod,get_nstep(), is_end_curr_month(),is_end_week)

     end do   !! end particles loop
     PFweightupdatesum=sum(PF_ppar%pfweightupdate(:,1)) !do 循环结束后自动增加1，这个时候还是算的i=1的时候的weight，但i已经自动成2了
     write(*,*) 'PFweightupdatesum=' , PFweightupdatesum

     do p2=1,parloop
          !write(*,*) "p2=", p2
          pfweightnormalize(p2)=PF_ppar%pfweightupdate(p2,1)/PFweightupdatesum
          if (pfweightnormalize(p2)>0. .and. pfweightnormalize(p2)<1.e-50) pfweightnormalize(p2)=0.
          !write(*,*) "pfweightnormalize(p2)=", pfweightnormalize(p2)
     end do
     !print *, ' pfweightnormalize(p2)!'

     SS=0
     do p3=1,parloop
         !write(*,*) "p3=", p3
         !write(*,*) "pfweightnormalize(p3)=", pfweightnormalize(p3)
         SS=SS+(pfweightnormalize(p3))**2

         !write(*,*) "SS=", SS
     end do
     !print *, ' SS=SS+(pfweightnormalize(p3))**2!'
     SS=1.0/SS
     !print *, ' SS=1.0/SS!'
     if (SS<parloop*0.95) then
        print *, "resample"
        inparticles(:,1)=PF_ppar%Vcmax(:,1)
        inparticles(:,2)=PF_ppar%q10(:,1)
        inparticles(:,3)=PF_ppar%VJ_slope(:,1)
        ! inparticles(:,4)=PF_ppar%N_leaf(:,1)
        inparticles(:,4)=PF_ppar%VN_slope(:,1)
        !inparticles(:,5)=PF_ppar%r_decay(:,1) ! use f_decay not r_decay
        inparticles(:,5)=PF_ppar%b_h2o(:,1)
        inparticles(:,6)=PF_ppar%sif_alpha(:,1)
        inparticles(:,7)=PF_ppar%sif_beta(:,1)
        !inparticles(:,9)=PF_ppar%taweff(:,1)
        !inparticles(:,10)=PF_ppar%D0(:,1)
        inparticles(:,8)=PF_ppar%Ksat_scalar(:,1)
        inparticles(:,9)=PF_ppar%b_scalar(:,1)
        inparticles(:,10)=PF_ppar%m_h2o(:,1)
        inparticles(:,11)=PF_ppar%f_leaf(:,1)
        inparticles(:,12)=PF_ppar%kc25(:,1)
        inparticles(:,13)=PF_ppar%ko25(:,1)
        inparticles(:,14)=PF_ppar%tau25(:,1)
        !inparticles(:,18)=PF_ppar%agb2vod(:,1)
        inparticles(:,15)=PF_ppar%f_resp(:,1)
        inparticles(:,16)=PF_ppar%f_decay(:,1)
        inparticles(:,17)=PF_ppar%a(:,1)
        inparticles(:,18)=PF_ppar%b(:,1)
        inparticles(:,19)=PF_ppar%c(:,1)

        inparticles(:,20)=PF_ppar%pfweight(:,1)

        call resample(inparticles,PF_ppar%pfweight(:,1))

        PF_ppar%Vcmax(:,1)= PF_resa%outparticles(:,1)
        PF_ppar%q10(:,1)=PF_resa%outparticles(:,2)
        PF_ppar%VJ_slope(:,1)=PF_resa%outparticles(:,3)
        !PF_ppar%N_leaf(:,1)=PF_resa%outparticles(:,4)
        PF_ppar%VN_slope(:,1)=PF_resa%outparticles(:,4)
        !PF_ppar%r_decay(:,1)=PF_resa%outparticles(:,5)
        PF_ppar%b_h2o(:,1)=PF_resa%outparticles(:,5)
        PF_ppar%sif_alpha(:,1)=PF_resa%outparticles(:,6)
        PF_ppar%sif_beta(:,1)=PF_resa%outparticles(:,7)
        !PF_ppar%taweff(:,1)=PF_resa%outparticles(:,9)
        !PF_ppar%D0(:,1)=PF_resa%outparticles(:,10)
        PF_ppar%Ksat_scalar(:,1)=PF_resa%outparticles(:,8)
        PF_ppar%b_scalar(:,1)=PF_resa%outparticles(:,9)
        PF_ppar%m_h2o(:,1)=PF_resa%outparticles(:,10)
        PF_ppar%f_leaf(:,1)=PF_resa%outparticles(:,11)
        PF_ppar%kc25(:,1)=PF_resa%outparticles(:,12)
        PF_ppar%ko25(:,1)=PF_resa%outparticles(:,13)
        PF_ppar%tau25(:,1)=PF_resa%outparticles(:,14)
        !PF_ppar%agb2vod(:,1)=PF_resa%outparticles(:,18)
        PF_ppar%f_resp(:,1)=PF_resa%outparticles(:,15)
        PF_ppar%f_decay(:,1)=PF_resa%outparticles(:,16)
        PF_ppar%a(:,1)=PF_resa%outparticles(:,17)
        PF_ppar%b(:,1)=PF_resa%outparticles(:,18)
        PF_ppar%c(:,1)=PF_resa%outparticles(:,19)
        PF_resa%resample_weight(:)=PF_resa%outparticles(:,20)

        do p4=1,parloop
           call PF_weight_update_resample(p4,PF_resa%resample_weight(p4))
        end do
        PFweightupdatesum_resample=sum(PF_resa%resample_weight_update(:))

        do p5=1,parloop
          PF_ppar%pfweight(p5,1)=PF_resa%resample_weight_update(p5)/PFweightupdatesum_resample
        end do

     else
        PF_ppar%Vcmax(:,1)= PF_ppar%Vcmax(:,1)
        PF_ppar%q10(:,1)=PF_ppar%q10(:,1)
        PF_ppar%VJ_slope(:,1)=PF_ppar%VJ_slope(:,1)
        !PF_ppar%N_leaf(:,1)=PF_ppar%N_leaf(:,1)
        PF_ppar%VN_slope(:,1)=PF_ppar%VN_slope(:,1)
        !PF_ppar%r_decay(:,1)=PF_ppar%r_decay(:,1)
        PF_ppar%b_h2o(:,1)=PF_ppar%b_h2o(:,1)
        PF_ppar%sif_alpha(:,1)=PF_ppar%sif_alpha(:,1)
        PF_ppar%sif_beta(:,1)=PF_ppar%sif_beta(:,1)
        !PF_ppar%taweff(:,1)=PF_ppar%taweff(:,1)
        !PF_ppar%D0(:,1)=PF_ppar%D0(:,1)
        PF_ppar%Ksat_scalar(:,1)=PF_ppar%Ksat_scalar(:,1)
        PF_ppar%b_scalar(:,1)=PF_ppar%b_scalar(:,1)
        PF_ppar%m_h2o(:,1)=PF_ppar%m_h2o(:,1)
        PF_ppar%f_leaf(:,1)=PF_ppar%f_leaf(:,1)
        PF_ppar%kc25(:,1)=PF_ppar%kc25(:,1)
        PF_ppar%ko25(:,1)=PF_ppar%ko25(:,1)
        PF_ppar%tau25(:,1)=PF_ppar%tau25(:,1)
        !PF_ppar%agb2vod(:,1)=PF_ppar%agb2vod(:,1)
        PF_ppar%f_resp(:,1)=PF_ppar%f_resp(:,1)
        PF_ppar%f_decay(:,1)=PF_ppar%f_decay(:,1)
        PF_ppar%a(:,1)=PF_ppar%a(:,1)
        PF_ppar%b(:,1)=PF_ppar%b(:,1)
        PF_ppar%c(:,1)=PF_ppar%c(:,1)
        do p5=1,parloop
          PF_ppar%pfweight(p5,1)=PF_ppar%pfweightupdate(p5,1)/PFweightupdatesum
        end do
     end if

     if (is_end_week) then
        call Create_PF_para(11,8,parloop)
        write(6,*)   "call Create_PF_para(14,2,parloop)"
        write(*,*)   "day=", dy
     end if
     !*****************************************************************end PF**********************************************************
     !call mpi_barrier(mpi_comm_world,ierr)
     !! advance time
     call advance_timestep()

     !!! write oout data and restart fields
     !call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed,p1)
     !call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed)

     !-- iLab::restart requires next time-step
     call get_curr_date(yr,mn,dy,tod)
     if(restart_frq <0) then   !! <0 ndays
        kount = get_nstep()
        rst_nstep  = -restart_frq*86400/step
        if(mod(kount,rst_nstep) ==0) call restart_io("write",yr,mn,dy,tod)
     else if(restart_frq ==0) then
        if(is_end_curr_month()) call restart_io("write",yr,mn,dy,tod)
     end if

     !!! determine whether it's the last step
     if(is_last_step()) exit

    else

     call advance_timestep()

     !!! write oout data and restart fields
     !call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed)

     !-- iLab::restart requires next time-step
     call get_curr_date(yr,mn,dy,tod)
     if(restart_frq <0) then   !! <0 ndays
        kount = get_nstep()
        rst_nstep  = -restart_frq*86400/step
        if(mod(kount,rst_nstep) ==0) call restart_io("write",yr,mn,dy,tod)
     else if(restart_frq ==0) then
        if(is_end_curr_month()) call restart_io("write",yr,mn,dy,tod)
     end if

     !!! determine whether it's the last step
     if(is_last_step()) exit

    endif

 end do !! end time loop

end if

!! For clean up
call ESMF_Finalize()

!if(myid==0)
write(6,*)   "BEPS run finished successfully!"
!call MPI_Finalize(ierr)
!stop
end program
