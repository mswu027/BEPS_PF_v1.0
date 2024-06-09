!******************************************************
!! main program for BEPS
!! Editted : J.Wang
!! Date    : June2017
!! add particle filter(PA) by @Xiuli Xing, 2023/04
!! add Monte Carlo (MC) and update PA by @Lu Hu,  2023/08
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
use MC_run
use restart
use esmf
implicit none

type(climatedata)    :: meteo
type(results)        :: mid_res
real(r8)             :: CosZs,hr_loc,hr_arc             !! solar zenith angle, local time, local time arc
integer              :: i,j,k,llll,jj,ii,kk,pt_ind,lc_ind,st_ind, mc_ind
type(soil)           :: soilp                !! at single point
real(r8)             :: Ccd(0:4),Cssd(0:4),Csmd(0:4),Cfsd(0:4),Cfmd(0:4),Csm(0:4),Cm(0:4),Cs(0:4),Cp(0:4)
real(r8)             :: param(0:49),var_o(0:40),var_n(0:40),coef(0:49)
real(r8)             :: inparticles(200,46) ! 45 parameters for optimization, 46th is the pf weight
real(r8)             :: pfweightnormalize(200)
real(r8)             :: lai
type(para),pointer   :: ppar                  !! initial parameters
type(surf),pointer   :: bfields               !! boundary fields
type(forc),pointer   :: climate               !! climate fieldsii =
type(soils),pointer  :: psoil                 !! GLobally
type(res),pointer    :: pp                    !! for output
type(mc_res),pointer    :: pp_mc              !! for MC temp output
type(PF_para),pointer   :: PF_ppar            !! Particle Fliter parameters
type(PF_obs0),pointer   :: PF_obsva         !! Particle Fliter observation at nd hour
type(MC_para),pointer   :: MC_ppar
type(MC_obs0),pointer   :: MC_obsva
type(PF_resample0),pointer   :: PF_resa       !! Particle Fliter resample
type(param_vars)      :: p_param      !! C4 Crop from Xiaorong ---->Xiuli
type(Phen)            :: p_phenp            !! C4 Crop from Xiaorong ---->Xiuli
type(param_gdd),pointer   :: p_pgdd        !! C4 Crop from Xiaorong ---->Xiuli
! 2023/12/15
type(Phydraulic)        :: PHydra !! parameters related to plant hydraulics stored in struct  PH_para

real(r8)             :: ratio_cloud,shortRad_df,shortRad_dir,PFweightupdatesum,SS,PFweightupdatesum_resample
!real(r8)             :: NPP_yr_acc(npoints,PFT)            !! for storing yearly accumulated NPP, Mh/ha/s
! real(r8)             :: agb2vod
! real(r8)             :: D0(1:9)
! real(r8)             :: taweff(1:9)
integer              :: kount,rst_nstep,p,run_pf,p1,p2,p3,p4,p5
integer              :: start_year
integer              :: yr,mn,dy,tod,caldy,n_meteo,n_lai, n_obs_pf, n_obs_mc
!-- iLab::converted from real(r8) to integer
integer              :: doys
integer              :: ierr
real(r8)             :: daylen
!-- iLab::for revised interface to 'av_output'
character(len=len('YYYY-MM-DDTHH:MM:SS')), save :: ref_date = ''
integer :: yr_ref, mn_ref, dy_ref, tod_ref
real(r8)             :: secs_elapsed,secs_meteo,days_lai, secs_obs_pf, secs_obs_mc
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
! 2024/01/07
integer :: fw_flag
! 2024/03/12, flag for choosing the form of water stress on stomatal conductance or Vcmax
! 0: fws = 1.0, do not consider the water stress from soil or leaf water potential
! 1: fws = f_soilwater, soil moisture stress on BWB slope
! 2: fws = 1.0, vcmax = vcmax*f_feileaf,leaf water potential stress on Vcmax
! 3: fws = f_feileaf,leaf water potential stress on BWB slope
! 4: fws = 1.0, but f_feileaf only works on Etp, Eta = Etp * f_feileaf
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
!mc_length = sim_duration * 24 ! data length of obs, 'mc_length' in beps_par.F90
! nlp  regional total numbers of pixels
! npoints, numbers of pixels for each sub task in MPI
! as for site, nlp is the total numbers of sites, do not consider the MPI, npoints=nlp
if(nscale == 1) then
  nlp = n_site
  npoints = nlp   !
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

! add a time display flag by @Lu Hu in 'timemgr_init' subroutine in beps_time_manager.F90
! '0' means printing the time message, ~0 means not displaying
call timemgr_init(0)

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
! required pointer for forward model
pp       => output
ppar     => assim
! required pointer for Monte Carlo
pp_mc    => mc_output
MC_ppar  => MC
MC_obsva => MC_obs
! required pointer for particle filter
PF_ppar  => PF
PF_obsva => PF_obs
PF_resa  => PF_resample
!
p_pgdd   => pgdd


! 0 for forward model; 1 for particle filter; <0 for Monte Carlo
! US_MOz, data start year defined by users!
start_year=2004
run_pf = 0
PFweightupdatesum=0.
PFweightupdatesum_resample=0.


fw_flag = 2
! 2024/03/12, flag for choosing the form of water stress on stomatal conductance or Vcmax
! 0: fws = 1.0, do not consider the water stress from soil or leaf water potential
! 1: fws = f_soilwater, soil moisture stress on BWB slope
! 2: fws = 1.0, vcmax = vcmax*f_feileaf,leaf water potential stress on Vcmax
! 3: fws = f_feileaf,leaf water potential stress on BWB slope
! 4: fws = 1.0, but f_feileaf only works on Etp, Eta = Etp * f_feileaf
! 5: fws = 1.0, vcmax = vcmax*f_soilwater,soil moisture stress on Vcmax
! 6: fws = 1.0, vcmax = vcmax*min(f_soilwater,f_feileaf)

temp_gpp=0.
outGPP=0.
temp_npp=0.
outNPP=0.
temp_accu_gpp=0.
temp_accu_npp=0.
temp_accu_temp=0.

!..................................................

!**************************************run BEPS with Monte Carlo****************************************
if (run_pf<0) then

 write(*,*) '****run BEPS with Monte Carlo******'

 do mc_ind =1, nparameters ! nparameters defined in bepspar.f90

    write(*,*) "Monte Carlo times=", mc_ind

    call Create_MC_para(3,5,1)  ! 1st means the land cover index; 2nd means the soil texture index

    call set_timemgr_init(calendar_in = calendar, start_ymd_in = icdate, &
                           start_tod_in = icsec, nelapse_in = sim_duration,&
                           dtime_in = step)
    call timemgr_init(1)

    nd = 0

    do  !! start time looping
        nd = nd + 1
        ! write(*,*) 'nd=',nd
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
            !write(ref_date(1:4),   '(i4.4)') yr_ref
            !write(ref_date(5:5),   '(a)')    '-'
            !write(ref_date(6:7),   '(i2.2)') mn_ref
            !write(ref_date(8:8),   '(a)')    '-'
            !write(ref_date(9:10),  '(i2.2)') dy_ref
            !write(ref_date(11:19), '(a)')    'T00:00:00'
        end if

        call timemgr_diff_secs(yr_ref*10000+mn_ref*100+dy_ref, tod_ref, yr*10000+mn*100+dy, tod,&
            secs_elapsed)
        call get_CO2_concentration(yr,CO2_air)
        call get_COS_concentration(yr,COS_air)
         !! change hourly input into daily input for further using model for long-term simulations, @MOUSONG.WU, 201905
        if (meteo_input >= 0) then  ! call hourly meteo. input
          if (nscale == 0) then
              call read_meteo_hourly(yr, mn, dy, tod)
          else
             ! Note: start time should be different for each site!!
              call timemgr_diff_secs(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
                   secs_meteo)
              n_meteo = int(secs_meteo/3600 + 1)
              call read_meteo_site(n_meteo)
              if (is_first_step()) print *, 'read site meteo successfully!'
              ! print *, 'read site meteo successfully!'
          end if
        else
          if(is_first_step() .or. is_end_curr_day()) then
            call read_meteo_daily(yr, mn, dy, tod)
          end if
        end if

        if (lai_input >=0) then
            ! Note: start time should be different for each site!!
            call timemgr_datediff(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
                 days_lai)
            n_lai = int(days_lai)
            if (is_first_step()) print *, 'lai is input!'
            if (is_first_step() .or. is_end_curr_day()) then
                if (nscale == 0) then
                   call read_lai(yr, mn, dy, tod, n_lai)
                else
                   !write(*,*) "n_lai=", n_lai
                   call read_lai_site(n_lai)
                   !print *, 'read site lai successfully!'
                end if
            end if
        else
            if (is_first_step()) then
               print *, 'lai is simulated with phenology scheme!'
            end if
        end if

        ! call mpi_barrier(mpi_comm_world,ierr)
        ! Note: start time should be different for each site!!
        call timemgr_diff_secs(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
                  secs_obs_mc)


        ! n_gpp_pf = int(secs_gpp_pf/3600 + 1)
        ! n_vod_mc is the time index for read observed values
        n_obs_mc = int(secs_obs_mc/3600 + 1)

        !write(*,*) 'n_obs_mc=', n_obs_mc
        call read_MC_obs(n_obs_mc) !xiuli

        ! spatial iteration
        do pt_ind =1, npoints
           !! calculate the solar zenith
           call s_coszs(yr, mn, dy, tod, caldy, doys, &
               bfields%latitude(pt_ind),bfields%longitude(pt_ind),CosZs,hr_loc,hr_arc)
           !!if(myid == 0) write(*,*) "hr_loc=",hr_loc
           !! retrieve meteo for this point
#ifdef COUP_CSM
          meteo%LR                     = climate%Lwdn(pt_ind)
          meteo%rainfall               = climate%Rain(pt_ind)
          meteo%snow                   = climate%Snow(pt_ind)
          meteo%S_dff                  = climate%Swdf(pt_ind)
          meteo%S_dir                  = climate%Swdr(pt_ind)
          meteo%Srad                   = meteo%S_dff + meteo%S_dir
          meteo%wind                   = climate%Wind(pt_ind)
          !meteo%rh                     = ...
#else
          meteo%Srad                   = climate%Srad(pt_ind)
          meteo%wind                   = climate%Wind(pt_ind)
          meteo%rainfall               = climate%Rain(pt_ind)
          meteo%snow                   = climate%Snow(pt_ind)
          meteo%rh                     = climate%Rh(pt_ind)
          if (meteo_input < 0) then
            meteo%tempmx                = climate%Tempmx(pt_ind)    ! read daily max and min temperatures instead
            meteo%tempmn                = climate%Tempmn(pt_ind)
          else
            meteo%temp                  = climate%Temp(pt_ind)
          end if

          ! .. compute daily course of temperature and daylength
           rdaymid = REAL (sim_duration+1) / 2.
           delta = -23.4*COS(2.*PI*(rdaymid+10.)/365.)
           spds = SIN(bfields%latitude(pt_ind)*pio180)*SIN(delta*pio180)
           cpds = COS(bfields%latitude(pt_ind)*pio180)*COS(delta*pio180)
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

          !###########Compute subdaily temperature based on daily input,@MOUSONG.WU,201905#####################

          if(meteo_input < 0) then
              ! .. compute average conditions
              atmean = (meteo%tempmx + meteo%tempmn)/ 2.
              atrange = meteo%tempmx - meteo%tempmn

              !hour angle at sunset, added by MOUSONG.WU, 201905
              sunset_arc   = (daylen/2.)*2.0*PI/24.0

              IF (daylen>=4. .AND. daylen<=20.) THEN
                  !sunrise
                  h0 = 12. - daylen/2.
                  !sundown
                  h1 = 12. + daylen/2.
                  !at sundown:
                  sd1 = SIN(PI*(2.*h1+(daylen-52.)/2.)/(daylen+4.))

                  !! unroll zum vektorisieren
                  IF (hr_loc>h0 .AND. hr_loc<h1) THEN
                      sd = SIN(PI*(2.*hr_loc+(daylen-52.)/2.)/(daylen+4.))
                      meteo%temp = atmean + atrange/2.*sd
                  ELSE
                      ! temperature at sundown
                      tmp1 = atmean + atrange/2.*sd1
                      ! hours since sundown
                      dhour = MOD(hr_loc-h1+24.,24.)
                      tmin = atmean - atrange/2.
                      meteo%temp = tmp1 - (tmp1-tmin)*(dhour/(24.-daylen))
                  END IF
              ELSEIF (daylen>20.) THEN
                  sd = COS(PI*(hr_loc-14.)/(daylen/2.+2.))
                  meteo%temp = atmean + atrange/2.*sd
              ELSE
                  meteo%temp = atmean
              END IF
              climate%Temp(pt_ind) = meteo%temp
              ! convert daily solar radiation into hourly using the method by M. Collares-Pereira and A. Rabl,
              ! “The average distribution of solar radiation-correlations between diffuse and hemispherical
              !and between daily and hourly insolation values,” Solar Energy,vol. 22, no. 2, pp. 155–164, 1979.
              a = 0.409 + 0.5016 * SIN(sunset_arc - 60.)
              b = 0.6609 - 0.4767 * SIN(sunset_arc - 60.)
              meteo%Srad = meteo%Srad*(a+b*COS(hr_arc))*(PI/24.)*(COS(hr_arc)-COS(sunset_arc))/ &
                                &(SIN(sunset_arc)-(2*PI*sunset_arc/360.)*COS(sunset_arc))
          end if

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
           climate%Swdr(pt_ind)   = shortRad_dir
           climate%Swdf(pt_ind)   = shortRad_df
#endif
           do lc_ind =3, 3 ! !! PFT iteration  !

              if (bfields%lcno(pt_ind,lc_ind) > 0 .and. bfields%sw(pt_ind) >0. .and. &
                bfields%stext(pt_ind) >0 .and. bfields%clumping(pt_ind) > 0.5) then

                 call readparam(bfields%lcno(pt_ind,lc_ind),param)

                 if(lai_input >= 0) then
                    lai = bfields%lai(pt_ind,lc_ind)
                 else
                    if (is_first_step()) then
                       lai = bfields%laiyr(pt_ind,lc_ind)
                       mid_res%lai_old = lai
                    else
                       lai = mid_res%lai_new
                    end if
                 end if

                 lai = lai*param(2)/bfields%clumping(pt_ind)
                 call readcoef(bfields%lcno(pt_ind,lc_ind),bfields%stext(pt_ind),coef)

                 if (nsrest == nsrStartup .and. is_first_step()) then

                    call Init_soil_parameters(bfields%lcno(pt_ind,lc_ind),bfields%laiyr(pt_ind,lc_ind),&
                              bfields%stext(pt_ind),MC_ppar%f_decay(1,pt_ind), param(27),soilp)

                    soilp%r_drainage = param(26)*MC_ppar%drainage_scalar(1,pt_ind)

                    ! index of soil texture
                    st_ind = bfields%stext(pt_ind)

                    do kk =0, 4
                       soilp%Ksat(kk) = MC_ppar%Ksat_scalar(1,pt_ind)*soilp%Ksat(kk)
                       soilp%b(kk)    = MC_ppar%b_scalar(1,pt_ind)*soilp%b(kk)
                       soilp%fei(kk)    = MC_ppar%porosity_scalar(1,pt_ind)*soilp%fei(kk)
                       soilp%theta_vfc(kk)    = MC_ppar%vfc_scalar(1,pt_ind)*soilp%theta_vfc(kk)
                       soilp%theta_vwp(kk)    = MC_ppar%vwp_scalar(1,pt_ind)*soilp%theta_vwp(kk)
                       soilp%psi_sat(kk)    = MC_ppar%psisat_scalar(1,pt_ind)*soilp%psi_sat(kk)
                    end do

                    call Init_soil_status(soilp,bfields%st(pt_ind),climate%Temp(pt_ind),bfields%sw(pt_ind),bfields%sdp(pt_ind))

                    do k = 0,40
                       var_o(k)   = 0.
                    end do

                    do k = 3,8
                       var_o(k)   = climate%Temp(pt_ind)
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

                 else if (nsrest == nsrContinue .and. is_first_step()) then

                      call Init_soil_parameters(bfields%lcno(pt_ind,lc_ind),bfields%laiyr(pt_ind,lc_ind),&
                              bfields%stext(pt_ind),MC_ppar%f_decay(1,pt_ind), param(27),soilp)

                      soilp%r_drainage = param(26)*MC_ppar%drainage_scalar(1,pt_ind)

                      st_ind = bfields%stext(pt_ind)

                      do kk =0, 4
                         soilp%Ksat(kk) = MC_ppar%Ksat_scalar(1,pt_ind)*soilp%Ksat(kk)
                         soilp%b(kk)    = MC_ppar%b_scalar(1,pt_ind)*soilp%b(kk)
                         soilp%fei(kk)    = MC_ppar%porosity_scalar(1,pt_ind)*soilp%fei(kk)
                         soilp%theta_vfc(kk)    = MC_ppar%vfc_scalar(1,pt_ind)*soilp%theta_vfc(kk)
                         soilp%theta_vwp(kk)    = MC_ppar%vwp_scalar(1,pt_ind)*soilp%theta_vwp(kk)
                         soilp%psi_sat(kk)    = MC_ppar%psisat_scalar(1,pt_ind)*soilp%psi_sat(kk)
                      end do

                      var_o(:)         = v2last(pt_ind,:,lc_ind)

                      do k= 9, 14
                         soilp%temp_soil_c(k-9)  = var_o(k)
                      end do

                      do k= 21, 26
                         soilp%thetam(k-21)      = var_o(k)
                      end do

                 else
                      var_o(:)         = v2last(pt_ind,:,lc_ind)
                      call retrive_soilp(soilp,pt_ind,lc_ind,0)
                 end if

                 call Init_planthydraulics_para(PHydra)
                 ! plant hydraulics parameters
                 PHydra%theta_Amin = MC_ppar%theta_Amin(1,pt_ind)
                 PHydra%pox = MC_ppar%pox(1,pt_ind)
                 PHydra%fei_c = MC_ppar%fei_c(1,pt_ind)
                 PHydra%spac_p1 = MC_ppar%spac_p1(1,pt_ind)
                 PHydra%spac_p2 = MC_ppar%spac_p2(1,pt_ind)
                 PHydra%tWA = MC_ppar%tWA(1,pt_ind)
                 PHydra%tWB = MC_ppar%tWB(1,pt_ind)
                 PHydra%Ttrig = MC_ppar%Ttrig(1,pt_ind)
                 PHydra%r_xylem = MC_ppar%r_xylem(1,pt_ind)
                 PHydra%r_r = MC_ppar%r_r(1,pt_ind)
                 PHydra%Lr = MC_ppar%Lr(1,pt_ind)
                 PHydra%deltal_min = MC_ppar%deltal_min(1,pt_ind)
                 PHydra%deltal_max = MC_ppar%deltal_max(1,pt_ind)
                 PHydra%p_delta = MC_ppar%p_delta(1,pt_ind)
                 PHydra%ppslh = MC_ppar%ppslh(1,pt_ind)
                 PHydra%fei_min = MC_ppar%fei_min(1,pt_ind)
                 PHydra%fei_th = MC_ppar%fei_th(1,pt_ind)
                 PHydra%p_excess = MC_ppar%p_excess(1,pt_ind)
                 PHydra%Tleaf_H = MC_ppar%Tleaf_H(1,pt_ind)
                 PHydra%Tleaf_L = MC_ppar%Tleaf_L(1,pt_ind)
                 PHydra%Tleaf_O = MC_ppar%Tleaf_O(1,pt_ind)

                 call inter_prg(yr, mn, dy, tod, lai,lai_input,bfields%lcno(pt_ind,lc_ind),bfields%clumping(pt_ind), &
                                MC_ppar%Vcmax(1,pt_ind),MC_ppar%VJ_slope(1,pt_ind),MC_ppar%VN_slope(1,pt_ind), &
                                MC_ppar%b_h2o(1,pt_ind),MC_ppar%m_h2o(1,pt_ind),MC_ppar%f_leaf(1,pt_ind),&
                                MC_ppar%kc25(1,pt_ind),MC_ppar%ko25(1,pt_ind),MC_ppar%tau25(1,pt_ind),&
                                MC_ppar%sif_alpha(1,pt_ind),MC_ppar%sif_beta(1,pt_ind),&
                                MC_ppar%vod_a(1,pt_ind),MC_ppar%vod_b(1,pt_ind),MC_ppar%vod_c(1,pt_ind),PHydra,&
                                param,meteo,CosZs,var_o,var_n,soilp,mid_res,daylen,fw_flag)

                 v2last(pt_ind,:,lc_ind)  = var_n(:)
                 call retrive_soilp(soilp,pt_ind,lc_ind,1)
                 !!! simluating Ra
                 call plant_resp(MC_ppar%q10(1,pt_ind),bfields%lcno(pt_ind,lc_ind),mid_res,bfields%laiyr(pt_ind,lc_ind), &
                                lai,meteo%temp,soilp%temp_soil_c(1),CosZs)

                 !!! simulating Rh
                 Ccd(0)       = bfields%Ccd(pt_ind,lc_ind)
                 Cssd(0)      = bfields%Cssd(pt_ind,lc_ind)
                 Csmd(0)      = bfields%Csmd(pt_ind,lc_ind)
                 Cfsd(0)      = bfields%Cfsd(pt_ind,lc_ind)
                 Cfmd(0)      = bfields%Cfmd(pt_ind,lc_ind)
                 Csm(0)       = bfields%Csm(pt_ind,lc_ind)
                 Cm(0)        = bfields%Cm(pt_ind,lc_ind)
                 Cs(0)        = bfields%Cm(pt_ind,lc_ind)
                 Cp(0)        = bfields%Cp(pt_ind,lc_ind)

                 st_ind = bfields%stext(pt_ind)

                 call soil_resp(MC_ppar%f_resp(1,pt_ind),Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,&
                                 Cp,bfields%nppyr(pt_ind,lc_ind),coef,bfields%stext(pt_ind),soilp,mid_res)
                 !write(*,*), 'pt_ind=', pt_ind
                 !write(*,*), 'lc_ind=', lc_ind
                 !write(*,*), 'vod*pct_pft=', mid_res%VOD*bfields%PCT_PFT(pt_ind,lc_ind)/100
                 !pp_mc%SMpft(pt_ind,lc_ind) = mid_res%thetam_surf*bfields%PCT_PFT(pt_ind,lc_ind)/100.
                 !2023/12/17 moflux sm depth 0.38 m, using beps layer 3 assimilation(0.35 m)
                 pp_mc%SMpft(pt_ind,lc_ind) = mid_res%thetam3*bfields%PCT_PFT(pt_ind,lc_ind)/100.
                 pp_mc%LWPpft(pt_ind,lc_ind) = mid_res%fei_leaf*bfields%PCT_PFT(pt_ind,lc_ind)/100.

            end if



          end do  !! end PFT loop

          pp_mc%SM(nd,pt_ind)     = sum(pp_mc%SMpft(pt_ind,:))
          pp_mc%tempSM(pt_ind) = sum(pp_mc%SMpft(pt_ind,:))
          pp_mc%obsSM(nd,pt_ind) = MC_obsva%obs_SM(pt_ind)
          pp_mc%LWP(nd,pt_ind) = sum(pp_mc%LWPpft(pt_ind,:))
          pp_mc%obsLWP(nd,pt_ind) = MC_obsva%obs_LWP(pt_ind)

        end do  !! end spatial loop

        !! advance time
        call advance_timestep()

       !!! write MC result
        !call write_MCoutput_site(yr,mn,dy,tod,ref_date,secs_elapsed,mc_ind)
       !call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed,p)
       !call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed)

       !-- iLab::restart requires next time-step
        call get_curr_date(yr,mn,dy,tod)

        if(is_last_step()) exit


    end do !! end time loop

    call Performance_calc()

    call LWP_Performance_calc()
    !print *, 'end of LWP_Performance_calc'

    call write_MC_para(mc_ind)

 end do !! end MC loop




!**************************************run BEPS with default parameters parameters**********************
else if (run_pf==0) then

 call read_prior_para()       ! put the parameters to be optimized in a NETCDF
                                !file and read them as well their
                                !uncertainties,@MOUSONG.WU,2019-11
 print *, 'read prior para successfully!'

  !! setting time manager
  !!! calling time_manager set init
!  call set_timemgr_init(calendar_in  = calendar,&
!                         start_ymd_in = icdate  ,&
!                         start_tod_in = icsec   ,&
!                         nelapse_in   = sim_duration,&
!                         dtime_in     = step)

!  call timemgr_init()

  nd = 0

  do     !! start time looping
    nd = nd + 1
    !-- iLab::inserted in order to pass 'calday' downstream
    caldy = get_curr_calday()
    call get_curr_date(yr,mn,dy,tod)
    doys = get_doys(yr)
    !write(*,*) 'nd=', nd
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
          ! each site should be different!
          call timemgr_diff_secs(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            secs_meteo)
          n_meteo = int(secs_meteo/3600 + 1)
          call read_meteo_site(n_meteo)
          !print *, 'read site meteo successfully!'
      end if
    else
      if(is_first_step() .or. is_end_curr_day()) then
        call read_meteo_daily(yr, mn, dy, tod)
      end if
    end if

    if (lai_input >=0) then
	    call timemgr_datediff(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            days_lai)
        n_lai = int(days_lai)
        if (is_first_step()) print *, 'lai is input!'
        if (is_first_step() .or. is_end_curr_day()) then
          if (nscale == 0) then
             call read_lai(yr, mn, dy, tod, n_lai)
          else
             call read_lai_site(n_lai)
             !print *, 'read site lai successfully!'
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
#ifdef COUP_CSM
          meteo%LR                     = climate%Lwdn(i)
          meteo%rainfall               = climate%Rain(i)
          meteo%snow                   = climate%Snow(i)
          meteo%S_dff                  = climate%Swdf(i)
          meteo%S_dir                  = climate%Swdr(i)
          meteo%Srad                   = meteo%S_dff + meteo%S_dir
          meteo%wind                   = climate%Wind(i)
          !meteo%rh                     = ...
#else
          meteo%Srad                   = climate%Srad(i)
          meteo%wind                   = climate%Wind(i)
          meteo%rainfall               = climate%Rain(i)
          meteo%snow                   = climate%Snow(i)
          meteo%rh                     = climate%Rh(i)
          if (meteo_input < 0) then
            meteo%tempmx                = climate%Tempmx(i)    ! read daily max and min temperatures instead
            meteo%tempmn                = climate%Tempmn(i)
          else
            meteo%temp                  = climate%Temp(i)
          end if

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

          !###########Compute subdaily temperature based on daily input,@MOUSONG.WU,201905#####################

          if(meteo_input < 0) then
              ! .. compute average conditions
              atmean = (meteo%tempmx + meteo%tempmn)/ 2.
              atrange = meteo%tempmx - meteo%tempmn

              !hour angle at sunset, added by MOUSONG.WU, 201905
              sunset_arc   = (daylen/2.)*2.0*PI/24.0

              IF (daylen>=4. .AND. daylen<=20.) THEN
                  !sunrise
                  h0 = 12. - daylen/2.
                  !sundown
                  h1 = 12. + daylen/2.
                  !at sundown:
                  sd1 = SIN(PI*(2.*h1+(daylen-52.)/2.)/(daylen+4.))

                  !! unroll zum vektorisieren
                  IF (hr_loc>h0 .AND. hr_loc<h1) THEN
                      sd = SIN(PI*(2.*hr_loc+(daylen-52.)/2.)/(daylen+4.))
                      meteo%temp = atmean + atrange/2.*sd
                  ELSE
                      ! temperature at sundown
                      tmp1 = atmean + atrange/2.*sd1
                      ! hours since sundown
                      dhour = MOD(hr_loc-h1+24.,24.)
                      tmin = atmean - atrange/2.
                      meteo%temp = tmp1 - (tmp1-tmin)*(dhour/(24.-daylen))
                  END IF
              ELSEIF (daylen>20.) THEN
                  sd = COS(PI*(hr_loc-14.)/(daylen/2.+2.))
                  meteo%temp = atmean + atrange/2.*sd
              ELSE
                  meteo%temp = atmean
              END IF
              climate%Temp(i) = meteo%temp
              ! convert daily solar radiation into hourly using the method by M. Collares-Pereira and A. Rabl,
              ! “The average distribution of solar radiation-correlations between diffuse and hemispherical
              !and between daily and hourly insolation values,” Solar Energy,vol. 22, no. 2, pp. 155–164, 1979.
              a = 0.409 + 0.5016 * SIN(sunset_arc - 60.)
              b = 0.6609 - 0.4767 * SIN(sunset_arc - 60.)
              meteo%Srad = meteo%Srad*(a+b*COS(hr_arc))*(PI/24.)*(COS(hr_arc)-COS(sunset_arc))/ &
                                &(SIN(sunset_arc)-(2*PI*sunset_arc/360.)*COS(sunset_arc))
          end if

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
#endif
          do j = 1,PFT !3,3!1,PFT    !! PFT iteration  ! 2023/06/30 beps land cover id
               ! bfields%clumping(i) > 0.5 should be improved
               if(bfields%lcno(i,j) > 0 .and. bfields%sw(i) >0. .and. bfields%stext(i) >0 .and. bfields%clumping(i) > 0.5) then
                    call readparam(bfields%lcno(i,j),param)

                    !param(29)=24.2 ! actual canopy height from US-Moz BASE file
                    !param(29)=27.0 ! actual canopy height from US-MMS BASE file
                    !param(29)=2.75 ! actual canopy height from US-Ses BASE file
                    param(29) = bfields%HeightC(i)
                    !param(30) = 0.5
                    !write(*,*), 'height=',param(29)

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

                        call Init_soil_parameters(bfields%lcno(i,j),bfields%laiyr(i,j),bfields%stext(i), &
                                            ppar%p_f_decay(j,i), param(27),soilp) ! 2023/06/30

                        !soilp%r_drainage = param(26)
                        soilp%r_drainage = param(26)*ppar%p_drainage_scalar(j,i)

                        !soilp%r_drainage = ppar%p_drainage(j)            ! read this
                                                                    !para from NC file,@MOUSONG.WU,2019-11
                        ii = bfields%stext(i)

                        !write(*,*) 'Ksat_old = ',soilp%Ksat(0)
                        do kk = 0,4
                            soilp%Ksat(kk) = ppar%p_Ksat_scalar(ii,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = ppar%p_b_scalar(ii,i)*soilp%b(kk)
                            soilp%fei(kk)    = ppar%p_porosity_scalar(ii,i)*soilp%fei(kk)
                            soilp%theta_vfc(kk)    = ppar%p_vfc_scalar(ii,i)*soilp%theta_vfc(kk)
                            soilp%theta_vwp(kk)    = ppar%p_vwp_scalar(ii,i)*soilp%theta_vwp(kk)
                            soilp%psi_sat(kk)    = ppar%p_psisat_scalar(ii,i)*soilp%psi_sat(kk)
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

                        !soilp%r_drainage = param(26)
                        soilp%r_drainage = param(26)*ppar%p_drainage_scalar(j,i)

                        !soilp%r_drainage = ppar%p_drainage(j)  ! read from nc
                                                          ! file,MOUSONG.WU@2019-11
                        ii = bfields%stext(i)

                        do kk = 0,4
                            soilp%Ksat(kk) = ppar%p_Ksat_scalar(ii,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = ppar%p_b_scalar(ii,i)*soilp%b(kk)
                            soilp%fei(kk)    = ppar%p_porosity_scalar(ii,i)*soilp%fei(kk)
                            soilp%theta_vfc(kk)    = ppar%p_vfc_scalar(ii,i)*soilp%theta_vfc(kk)
                            soilp%theta_vwp(kk)    = ppar%p_vwp_scalar(ii,i)*soilp%theta_vwp(kk)
                            soilp%psi_sat(kk)    = ppar%p_psisat_scalar(ii,i)*soilp%psi_sat(kk)
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
                    ! 2023/12/15
                    call Init_planthydraulics_para(PHydra)
                    ! plant hydraulics parameters
                    PHydra%theta_Amin = ppar%p_theta_Amin(j,i)
                    PHydra%pox = ppar%p_pox(j,i)
                    PHydra%fei_c = ppar%p_fei_c(j,i)
                    PHydra%spac_p1 = ppar%p_spac_p1(j,i)
                    PHydra%spac_p2 = ppar%p_spac_p2(j,i)
                    PHydra%tWA = ppar%p_tWA(j,i)
                    PHydra%tWB = ppar%p_tWB(j,i)
                    PHydra%Ttrig = ppar%p_Ttrig(j,i)
                    PHydra%r_xylem = ppar%p_r_xylem(j,i)
                    PHydra%r_r = ppar%p_r_r(j,i)
                    PHydra%Lr = ppar%p_Lr(j,i)
                    PHydra%deltal_min = ppar%p_deltal_min(j,i)
                    PHydra%deltal_max = ppar%p_deltal_max(j,i)
                    PHydra%p_delta = ppar%p_p_delta(j,i)
                    PHydra%ppslh = ppar%p_ppslh(j,i)
                    PHydra%fei_min = ppar%p_fei_min(j,i)
                    PHydra%fei_th = ppar%p_fei_th(j,i)
                    PHydra%p_excess = ppar%p_p_excess(j,i)
                    !PHydra%Tleaf_H = ppar%p_Tleaf_H(j,i)
                    !PHydra%Tleaf_L = ppar%p_Tleaf_L(j,i)
                    !PHydra%Tleaf_O = ppar%p_Tleaf_O(j,i)
                    ! for US-MOz
                    !PHydra%fei_min=350.
                    !PHydra%fei_th=10.
                    !ppar%p_m_h2o(j,i)=6.0
                    ! 2024/05/13
                    !ppar%p_Vcmax(j,i) = 50.0


                    ! for us-mms
                    !PHydra%fei_min=250.
                    !ppar%p_Vcmax(j,i) = 48.0
                    !ppar%p_m_h2o(j,i)=6.0

                    ! for us-ses
                    !PHydra%fei_min=600.
                    !ppar%p_Vcmax(j,i) = 48.0
                    !ppar%p_m_h2o(j,i)=6.0


                    !ppar%p_b_h2o(j,i)=0.01
                    !ppar%p_Vcmax(j,i) = 50.0

                    !write(*,*), 'ppar%p_m_h2o(j,i)=',ppar%p_m_h2o(j,i)
                    call inter_prg(yr, mn, dy, tod,lai,lai_input,bfields%lcno(i,j),bfields%clumping(i),&
                            ppar%p_Vcmax(j,i),ppar%p_VJ_slope(j,i),ppar%p_VN_slope(j,i),ppar%p_b_h2o(j,i),&
                            ppar%p_m_h2o(j,i),ppar%p_f_leaf(j,i),ppar%p_kc25(j,i),ppar%p_ko25(j,i),&
                            ppar%p_tau25(j,i),ppar%p_sif_alpha(j,i),ppar%p_sif_beta(j,i),&
                            ppar%p_vod_a(j,i),ppar%p_vod_b(j,i),ppar%p_vod_c(j,i),PHydra,&
                            param,meteo,CosZs,var_o,var_n,soilp,mid_res,daylen,fw_flag)

                    !print *, 'end of inter_prg'
                    !write(*,*) 'ppar%p_a(j,i) = ', ppar%p_a(j,i)
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

                    pp%VODpft(i,j) = mid_res%VOD*bfields%PCT_PFT(i,j)/100.

                    pp%fei_leafpft(i,j) = mid_res%fei_leaf*bfields%PCT_PFT(i,j)/100.
                    pp%ETapft(i,j) = mid_res%ETa*bfields%PCT_PFT(i,j)/100.
                    pp%Quptpft(i,j) = mid_res%Qupt*bfields%PCT_PFT(i,j)/100.


                    ! these output for check
                    pp%Thetam2_pft(i,j)=mid_res%thetam2*bfields%PCT_PFT(i,j)/100.
                    pp%Thetam3_pft(i,j)=mid_res%thetam3*bfields%PCT_PFT(i,j)/100.
                    pp%Thetam4_pft(i,j)=mid_res%thetam4*bfields%PCT_PFT(i,j)/100.
                    pp%Thetam5_pft(i,j)=mid_res%thetam5*bfields%PCT_PFT(i,j)/100.
                    pp%SWP1_pft(i,j)=mid_res%swp1*bfields%PCT_PFT(i,j)/100.
                    pp%SWP2_pft(i,j)=mid_res%swp2*bfields%PCT_PFT(i,j)/100.
                    pp%SWP3_pft(i,j)=mid_res%swp3*bfields%PCT_PFT(i,j)/100.
                    pp%SWP4_pft(i,j)=mid_res%swp4*bfields%PCT_PFT(i,j)/100.
                    pp%SWP5_pft(i,j)=mid_res%swp5*bfields%PCT_PFT(i,j)/100.

                    pp%TS1_pft(i,j)=mid_res%TS1*bfields%PCT_PFT(i,j)/100.
                    pp%TS2_pft(i,j)=mid_res%TS2*bfields%PCT_PFT(i,j)/100.
                    pp%TS3_pft(i,j)=mid_res%TS3*bfields%PCT_PFT(i,j)/100.
                    pp%TS4_pft(i,j)=mid_res%TS4*bfields%PCT_PFT(i,j)/100.
                    pp%TS5_pft(i,j)=mid_res%TS5*bfields%PCT_PFT(i,j)/100.

                    !pp%PondWater_pft(i,j)=mid_res%Zp*bfields%PCT_PFT(i,j)/100.
                    !pp%Rain_g_pft(i,j)=mid_res%rain_g*bfields%PCT_PFT(i,j)/100.
                    pp%f_soilwater_pft(i,j)=mid_res%f_soilwater*bfields%PCT_PFT(i,j)/100.
                    pp%f_feileaf_pft(i,j)=mid_res%f_feileaf*bfields%PCT_PFT(i,j)/100.

                    pp%LHapft(i,j)    = mid_res%LHa*bfields%PCT_PFT(i,j)/100.

                    pp%GPP_o_sunlit_pft(i,j) = mid_res%gpp_o_sunlit*bfields%PCT_PFT(i,j)/100.
                    pp%GPP_o_shaded_pft(i,j) = mid_res%gpp_o_shaded*bfields%PCT_PFT(i,j)/100.
                    pp%GPP_u_sunlit_pft(i,j) = mid_res%gpp_u_sunlit*bfields%PCT_PFT(i,j)/100.
                    pp%GPP_u_shaded_pft(i,j) = mid_res%gpp_u_shaded*bfields%PCT_PFT(i,j)/100.

                    pp%TR_o_sunlit_pft(i,j) = mid_res%TR_o_sunlit*bfields%PCT_PFT(i,j)/100.
                    pp%TR_o_shaded_pft(i,j) = mid_res%TR_o_shaded*bfields%PCT_PFT(i,j)/100.
                    pp%TR_u_sunlit_pft(i,j) = mid_res%TR_u_sunlit*bfields%PCT_PFT(i,j)/100.
                    pp%TR_u_shaded_pft(i,j) = mid_res%TR_u_shaded*bfields%PCT_PFT(i,j)/100.

                    pp%Eil_pft(i,j) = mid_res%Eil*bfields%PCT_PFT(i,j)/100.
                    pp%Evap_soil_pft(i,j) = mid_res%Evap_soil*bfields%PCT_PFT(i,j)/100.
                    pp%Evap_SW_pft(i,j) = mid_res%Evap_SW*bfields%PCT_PFT(i,j)/100.
                    pp%EiS_pft(i,j) = mid_res%EiS*bfields%PCT_PFT(i,j)/100.
                    pp%Evap_SS_pft(i,j) = mid_res%Evap_SS*bfields%PCT_PFT(i,j)/100.
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
          pp%fei_leaf(i)     = sum(pp%fei_leafpft(i,:))
          pp%ETa(i) = sum(pp%ETapft(i,:))
          pp%Qupt(i) = sum(pp%Quptpft(i,:))

          ! 2023/10/30
          pp%Thetam_layer2(i) = sum(pp%Thetam2_pft(i,:))
          pp%Thetam_layer3(i) = sum(pp%Thetam3_pft(i,:))
          pp%Thetam_layer4(i) = sum(pp%Thetam4_pft(i,:))
          pp%Thetam_layer5(i) = sum(pp%Thetam5_pft(i,:))
          pp%SWP_layer1(i) = sum(pp%SWP1_pft(i,:))
          pp%SWP_layer2(i) = sum(pp%SWP2_pft(i,:))
          pp%SWP_layer3(i) = sum(pp%SWP3_pft(i,:))
          pp%SWP_layer4(i) = sum(pp%SWP4_pft(i,:))
          pp%SWP_layer5(i) = sum(pp%SWP5_pft(i,:))

          pp%TS_layer1(i) = sum(pp%TS1_pft(i,:))
          pp%TS_layer2(i) = sum(pp%TS2_pft(i,:))
          pp%TS_layer3(i) = sum(pp%TS3_pft(i,:))
          pp%TS_layer4(i) = sum(pp%TS4_pft(i,:))
          pp%TS_layer5(i) = sum(pp%TS5_pft(i,:))
          !pp%PondWater(i) = sum(pp%PondWater_pft(i,:))
          !pp%Rain_g(i) = sum(pp%Rain_g_pft(i,:))
          pp%f_soilwater(i) = sum(pp%f_soilwater_pft(i,:))
          pp%f_feileaf(i) = sum(pp%f_feileaf_pft(i,:))
          !pp%f_Tleaf(i) = sum(pp%f_Tleaf_pft(i,:))

          !2024/03/12
          pp%LHa(i)    = sum(pp%LHapft(i,:))

          pp%GPP_o_sunlit(i) = sum(pp%GPP_o_sunlit_pft(i,:))
          pp%GPP_o_shaded(i) = sum(pp%GPP_o_shaded_pft(i,:))
          pp%GPP_u_sunlit(i) = sum(pp%GPP_u_sunlit_pft(i,:))
          pp%GPP_u_shaded(i) = sum(pp%GPP_u_shaded_pft(i,:))

          pp%TR_o_sunlit(i) = sum(pp%TR_o_sunlit_pft(i,:))
          pp%TR_o_shaded(i) = sum(pp%TR_o_shaded_pft(i,:))
          pp%TR_u_sunlit(i) = sum(pp%TR_u_sunlit_pft(i,:))
          pp%TR_u_shaded(i) = sum(pp%TR_u_shaded_pft(i,:))

          pp%Eil(i) = sum(pp%Eil_pft(i,:))
          pp%Evap_soil(i) = sum(pp%Evap_soil_pft(i,:))
          pp%Evap_SW(i) = sum(pp%Evap_SW_pft(i,:))
          pp%EiS(i) = sum(pp%EiS_pft(i,:))
          pp%Evap_SS(i) = sum(pp%Evap_SS_pft(i,:))
          !write (*,*),'pp%EiS =',pp%EiS(i)


       end do      !! end spatial loop
       !          call mpi_barrier(mpi_comm_world,ierr)
       !! advance time
       call advance_timestep()

       !!! write out data and restart fields
       !call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed,p)
       call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed)
       !-- iLab::restart requires next time-step
       call get_curr_date(yr,mn,dy,tod)

       if(is_last_step()) exit
  end do   !! end time loop


!************************************Particle Filter***********************************
else if (run_pf==1) then
 write(*,*) 'run_pf=',run_pf

 !call read_prior_para()       ! put the parameters to be optimized in a NETCDF
                                !file and read them as well their
                                !uncertainties,@MOUSONG.WU,2019-11

 !call Create_PF_para(11,8,parloop) !"1" is PFT,"2" is soil texture !!create particles. xiuli
 call Create_PF_para(3,5,parloop)
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
    !is_end_month = ((yr>yr_ref .or. mn>mn_ref .or. dy>dy_ref) .and. (MOD((dy-dy_ref+1),7)==0) .and. (tod == 82800))
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
        ! Note: start time should be different for each site!!
          call timemgr_diff_secs(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
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
        ! Note: start time should be different for each site!!
        call timemgr_datediff(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
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
    ! Note: start time should be different for each site!!
    call timemgr_diff_secs(start_year*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            secs_obs_pf)
    ! n_gpp_pf = int(secs_gpp_pf/3600 + 1)
    n_obs_pf = int(secs_obs_pf/3600 + 1)
    !write(*,*) 'n_vod_pf=', n_vod_pf
    call read_PF_obs(n_obs_pf) !xiuli

    !if (PF_obsva%obs_GPP(1) == PF_obsva%obs_GPP(1)) then !if observation exsits xiuli
    if (PF_obsva%obs_var(1) /= -99999._r8) then
     do p1 =1,parloop !! start particles loop xiuli
     !do p =1,10!!1,nparameters !! start parameter loop

       if (p1==parloop) then
          write(*,*) 'p1=',p1
       end if

       do i = 1,npoints    !! spatial iteration

          !! calculate the solar zenith
          call s_coszs(yr, mn, dy, tod, caldy, doys, &
               bfields%latitude(i),bfields%longitude(i),CosZs,hr_loc,hr_arc)
          !!if(myid == 0) write(*,*) "hr_loc=",hr_loc
          !! retrieve meteo for this point
#ifdef COUP_CSM
          meteo%LR                     = climate%Lwdn(i)
          meteo%rainfall               = climate%Rain(i)
          meteo%snow                   = climate%Snow(i)
          meteo%S_dff                  = climate%Swdf(i)
          meteo%S_dir                  = climate%Swdr(i)
          meteo%Srad                   = meteo%S_dff + meteo%S_dir
          meteo%wind                   = climate%Wind(i)
          !meteo%rh                     = ...
#else
          meteo%Srad                   = climate%Srad(i)
          meteo%wind                   = climate%Wind(i)
          meteo%rainfall               = climate%Rain(i)
          meteo%snow                   = climate%Snow(i)
          meteo%rh                     = climate%Rh(i)
          if (meteo_input < 0) then
            meteo%tempmx                = climate%Tempmx(i)    ! read daily max and min temperatures instead
            meteo%tempmn                = climate%Tempmn(i)
          else
            meteo%temp                  = climate%Temp(i)
          end if

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

          !###########Compute subdaily temperature based on daily input,@MOUSONG.WU,201905#####################

          if(meteo_input < 0) then
              ! .. compute average conditions
              atmean = (meteo%tempmx + meteo%tempmn)/ 2.
              atrange = meteo%tempmx - meteo%tempmn

              !hour angle at sunset, added by MOUSONG.WU, 201905
              sunset_arc   = (daylen/2.)*2.0*PI/24.0

              IF (daylen>=4. .AND. daylen<=20.) THEN
                  !sunrise
                  h0 = 12. - daylen/2.
                  !sundown
                  h1 = 12. + daylen/2.
                  !at sundown:
                  sd1 = SIN(PI*(2.*h1+(daylen-52.)/2.)/(daylen+4.))

                  !! unroll zum vektorisieren
                  IF (hr_loc>h0 .AND. hr_loc<h1) THEN
                      sd = SIN(PI*(2.*hr_loc+(daylen-52.)/2.)/(daylen+4.))
                      meteo%temp = atmean + atrange/2.*sd
                  ELSE
                      ! temperature at sundown
                      tmp1 = atmean + atrange/2.*sd1
                      ! hours since sundown
                      dhour = MOD(hr_loc-h1+24.,24.)
                      tmin = atmean - atrange/2.
                      meteo%temp = tmp1 - (tmp1-tmin)*(dhour/(24.-daylen))
                  END IF
              ELSEIF (daylen>20.) THEN
                  sd = COS(PI*(hr_loc-14.)/(daylen/2.+2.))
                  meteo%temp = atmean + atrange/2.*sd
              ELSE
                  meteo%temp = atmean
              END IF
              climate%Temp(i) = meteo%temp
              ! convert daily solar radiation into hourly using the method by M. Collares-Pereira and A. Rabl,
              ! “The average distribution of solar radiation-correlations between diffuse and hemispherical
              !and between daily and hourly insolation values,” Solar Energy,vol. 22, no. 2, pp. 155–164, 1979.
              a = 0.409 + 0.5016 * SIN(sunset_arc - 60.)
              b = 0.6609 - 0.4767 * SIN(sunset_arc - 60.)
              meteo%Srad = meteo%Srad*(a+b*COS(hr_arc))*(PI/24.)*(COS(hr_arc)-COS(sunset_arc))/ &
                                &(SIN(sunset_arc)-(2*PI*sunset_arc/360.)*COS(sunset_arc))
          end if

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
#endif
          !write(*,*) 'Before start_PFT_loop!!!!!'
          ! Note: for each site, PFT may be different!
          do j = 3,3!1,PFT    !! PFT iteration
               !write(*,*) 'start_PFT_loop!!!!'
               if(bfields%lcno(i,j) > 0 .and. bfields%sw(i) >0. .and. bfields%stext(i) >0 .and. bfields%clumping(i) > 0.5) then
                    call readparam(bfields%lcno(i,j),param)
                    param(29)=24.2

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

                        call Init_soil_parameters(bfields%lcno(i,j),bfields%laiyr(i,j),bfields%stext(i), &
                                            PF_ppar%f_decay(p1,i), param(27),soilp) ! 2023/06/30
                        !soilp%r_drainage = param(26)
                        soilp%r_drainage = param(26)*PF_ppar%drainage_scalar(p1,i)
                        !soilp%r_drainage = ppar%p_drainage(j)            ! read this
                                                                    !para from NC file,@MOUSONG.WU,2019-11

                        ii = bfields%stext(i)

                        !write(*,*) 'Ksat_old = ',soilp%Ksat(0)
                        do kk = 0,4
                            soilp%Ksat(kk) = PF_ppar%Ksat_scalar(p1,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = PF_ppar%b_scalar(p1,i)*soilp%b(kk)
                            soilp%fei(kk)    = PF_ppar%porosity_scalar(p1,i)*soilp%fei(kk)
                            soilp%theta_vfc(kk)    = PF_ppar%vfc_scalar(p1,i)*soilp%theta_vfc(kk)
                            soilp%theta_vwp(kk)    = PF_ppar%vwp_scalar(p1,i)*soilp%theta_vwp(kk)
                            soilp%psi_sat(kk)    = PF_ppar%psisat_scalar(p1,i)*soilp%psi_sat(kk)
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

                        !soilp%r_drainage = param(26)
                        soilp%r_drainage = param(26)*PF_ppar%drainage_scalar(p1,i)
                        !soilp%r_drainage = ppar%p_drainage(j)  ! read from nc
                                                          ! file,MOUSONG.WU@2019-11
                        ii = bfields%stext(i)
                        do kk = 0,4
                            soilp%Ksat(kk) = PF_ppar%Ksat_scalar(p1,i)*soilp%Ksat(kk)
                            soilp%b(kk)    = PF_ppar%b_scalar(p1,i)*soilp%b(kk)
                            soilp%fei(kk)    = PF_ppar%porosity_scalar(p1,i)*soilp%fei(kk)
                            soilp%theta_vfc(kk)    = PF_ppar%vfc_scalar(p1,i)*soilp%theta_vfc(kk)
                            soilp%theta_vwp(kk)    = PF_ppar%vwp_scalar(p1,i)*soilp%theta_vwp(kk)
                            soilp%psi_sat(kk)    = PF_ppar%psisat_scalar(p1,i)*soilp%psi_sat(kk)
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

                    call Init_planthydraulics_para(PHydra)
                    ! plant hydraulics parameters
                    PHydra%theta_Amin = PF_ppar%theta_Amin(p1,i)
                    PHydra%pox = PF_ppar%pox(p1,i)
                    PHydra%fei_c = PF_ppar%fei_c(p1,i)
                    PHydra%spac_p1 = PF_ppar%spac_p1(p1,i)
                    PHydra%spac_p2 = PF_ppar%spac_p2(p1,i)
                    PHydra%tWA = PF_ppar%tWA(p1,i)
                    PHydra%tWB = PF_ppar%tWB(p1,i)
                    PHydra%Ttrig = PF_ppar%Ttrig(p1,i)
                    PHydra%r_xylem = PF_ppar%r_xylem(p1,i)
                    PHydra%r_r = PF_ppar%r_r(p1,i)
                    PHydra%Lr = PF_ppar%Lr(p1,i)
                    PHydra%deltal_min = PF_ppar%deltal_min(p1,i)
                    PHydra%deltal_max = PF_ppar%deltal_max(p1,i)
                    PHydra%p_delta = PF_ppar%p_delta(p1,i)
                    PHydra%ppslh = PF_ppar%ppslh(p1,i)
                    PHydra%fei_min = PF_ppar%fei_min(p1,i)
                    PHydra%fei_th = PF_ppar%fei_th(p1,i)
                    PHydra%p_excess = PF_ppar%p_excess(p1,i)
                    PHydra%Tleaf_H = PF_ppar%Tleaf_H(p1,i)
                    PHydra%Tleaf_L = PF_ppar%Tleaf_L(p1,i)
                    PHydra%Tleaf_O = PF_ppar%Tleaf_O(p1,i)

                    call inter_prg(yr, mn, dy, tod,lai,lai_input,bfields%lcno(i,j),bfields%clumping(i),&
                            PF_ppar%Vcmax(p1,i),PF_ppar%VJ_slope(p1,i),PF_ppar%VN_slope(p1,i),&
                            PF_ppar%b_h2o(p1,i),PF_ppar%m_h2o(p1,i),PF_ppar%f_leaf(p1,i),&
                            PF_ppar%kc25(p1,i),PF_ppar%ko25(p1,i),PF_ppar%tau25(p1,i),PF_ppar%sif_alpha(p1,i),&
                            PF_ppar%sif_beta(p1,i),PF_ppar%vod_a(p1,i),PF_ppar%vod_b(p1,i),PF_ppar%vod_c(p1,i),&
                            PHydra,param,meteo,CosZs,var_o,var_n,soilp,mid_res,daylen,fw_flag)
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

                    pp%Thetam3_pft(i,j)=mid_res%thetam3*bfields%PCT_PFT(i,j)/100.

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
          pp%Thetam_layer3(i) = sum(pp%Thetam3_pft(i,:))

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
        inparticles(:,2)=PF_ppar%VJ_slope(:,1)
        inparticles(:,3)=PF_ppar%VN_slope(:,1)
        inparticles(:,4)=PF_ppar%b_h2o(:,1)
        inparticles(:,5)=PF_ppar%m_h2o(:,1)
        inparticles(:,6)=PF_ppar%f_leaf(:,1)
        inparticles(:,7)=PF_ppar%kc25(:,1)
        inparticles(:,8)=PF_ppar%ko25(:,1)
        inparticles(:,9)=PF_ppar%tau25(:,1)

        inparticles(:,10)=PF_ppar%sif_alpha(:,1)
        inparticles(:,11)=PF_ppar%sif_beta(:,1)

        inparticles(:,12)=PF_ppar%q10(:,1)
        inparticles(:,13)=PF_ppar%f_resp(:,1)

        !inparticles(:,5)=PF_ppar%r_decay(:,1) ! use f_decay not r_decay
        inparticles(:,14)=PF_ppar%f_decay(:,1)
        inparticles(:,15)=PF_ppar%Ksat_scalar(:,1)
        inparticles(:,16)=PF_ppar%b_scalar(:,1)
        inparticles(:,17)=PF_ppar%porosity_scalar(:,1)
        inparticles(:,18)=PF_ppar%vfc_scalar(:,1)
        inparticles(:,19)=PF_ppar%vwp_scalar(:,1)
        inparticles(:,20)=PF_ppar%psisat_scalar(:,1)
        inparticles(:,21)=PF_ppar%drainage_scalar(:,1)

        inparticles(:,22)=PF_ppar%vod_a(:,1)
        inparticles(:,23)=PF_ppar%vod_b(:,1)
        inparticles(:,24)=PF_ppar%vod_c(:,1)

        inparticles(:,25)=PF_ppar%theta_Amin(:,1)
        inparticles(:,26)=PF_ppar%pox(:,1)
        inparticles(:,27)=PF_ppar%fei_c(:,1)
        inparticles(:,28)=PF_ppar%spac_p1(:,1)
        inparticles(:,29)=PF_ppar%spac_p2(:,1)
        inparticles(:,30)=PF_ppar%tWA(:,1)
        inparticles(:,31)=PF_ppar%tWB(:,1)
        inparticles(:,32)=PF_ppar%Ttrig(:,1)
        inparticles(:,33)=PF_ppar%r_xylem(:,1)
        inparticles(:,34)=PF_ppar%r_r(:,1)
        inparticles(:,35)=PF_ppar%Lr(:,1)
        inparticles(:,36)=PF_ppar%deltal_min(:,1)
        inparticles(:,37)=PF_ppar%deltal_max(:,1)
        inparticles(:,38)=PF_ppar%p_delta(:,1)
        inparticles(:,39)=PF_ppar%ppslh(:,1)
        inparticles(:,40)=PF_ppar%fei_min(:,1)
        inparticles(:,41)=PF_ppar%fei_th(:,1)
        inparticles(:,42)=PF_ppar%p_excess(:,1)
        inparticles(:,43)=PF_ppar%Tleaf_H(:,1)
        inparticles(:,44)=PF_ppar%Tleaf_L(:,1)
        inparticles(:,45)=PF_ppar%Tleaf_O(:,1)

        inparticles(:,46)=PF_ppar%pfweight(:,1)

        call resample(inparticles,PF_ppar%pfweight(:,1))

        PF_ppar%Vcmax(:,1)= PF_resa%outparticles(:,1)
        PF_ppar%VJ_slope(:,1)=PF_resa%outparticles(:,2)
        PF_ppar%VN_slope(:,1)=PF_resa%outparticles(:,3)
        PF_ppar%b_h2o(:,1)=PF_resa%outparticles(:,4)
        PF_ppar%m_h2o(:,1)=PF_resa%outparticles(:,5)
        PF_ppar%f_leaf(:,1)=PF_resa%outparticles(:,6)
        PF_ppar%kc25(:,1)=PF_resa%outparticles(:,7)
        PF_ppar%ko25(:,1)=PF_resa%outparticles(:,8)
        PF_ppar%tau25(:,1)=PF_resa%outparticles(:,9)

        PF_ppar%sif_alpha(:,1)=PF_resa%outparticles(:,10)
        PF_ppar%sif_beta(:,1)=PF_resa%outparticles(:,11)

        PF_ppar%q10(:,1)=PF_resa%outparticles(:,12)
        PF_ppar%f_resp(:,1)=PF_resa%outparticles(:,13)

        PF_ppar%f_decay(:,1)=PF_resa%outparticles(:,14)
        PF_ppar%Ksat_scalar(:,1)=PF_resa%outparticles(:,15)
        PF_ppar%b_scalar(:,1)=PF_resa%outparticles(:,16)
        PF_ppar%porosity_scalar(:,1)=PF_resa%outparticles(:,17)
        PF_ppar%vfc_scalar(:,1)=PF_resa%outparticles(:,18)
        PF_ppar%vwp_scalar(:,1)=PF_resa%outparticles(:,19)
        PF_ppar%psisat_scalar(:,1)=PF_resa%outparticles(:,20)
        PF_ppar%drainage_scalar(:,1)=PF_resa%outparticles(:,21)

        PF_ppar%vod_a(:,1)=PF_resa%outparticles(:,22)
        PF_ppar%vod_b(:,1)=PF_resa%outparticles(:,23)
        PF_ppar%vod_c(:,1)=PF_resa%outparticles(:,24)

        PF_ppar%theta_Amin(:,1)=PF_resa%outparticles(:,25)
        PF_ppar%pox(:,1)=PF_resa%outparticles(:,26)
        PF_ppar%fei_c(:,1)=PF_resa%outparticles(:,27)
        PF_ppar%spac_p1(:,1)=PF_resa%outparticles(:,28)
        PF_ppar%spac_p2(:,1)=PF_resa%outparticles(:,29)
        PF_ppar%tWA(:,1)=PF_resa%outparticles(:,30)
        PF_ppar%tWB(:,1)=PF_resa%outparticles(:,31)
        PF_ppar%Ttrig(:,1)=PF_resa%outparticles(:,32)
        PF_ppar%r_xylem(:,1)=PF_resa%outparticles(:,33)
        PF_ppar%r_r(:,1)=PF_resa%outparticles(:,34)
        PF_ppar%Lr(:,1)=PF_resa%outparticles(:,35)
        PF_ppar%deltal_min(:,1)=PF_resa%outparticles(:,36)
        PF_ppar%deltal_max(:,1)=PF_resa%outparticles(:,37)
        PF_ppar%p_delta(:,1)=PF_resa%outparticles(:,38)
        PF_ppar%ppslh(:,1)=PF_resa%outparticles(:,39)
        PF_ppar%fei_min(:,1)=PF_resa%outparticles(:,40)
        PF_ppar%fei_th(:,1)=PF_resa%outparticles(:,41)
        PF_ppar%p_excess(:,1)=PF_resa%outparticles(:,42)
        PF_ppar%Tleaf_H(:,1)=PF_resa%outparticles(:,43)
        PF_ppar%Tleaf_L(:,1)=PF_resa%outparticles(:,44)
        PF_ppar%Tleaf_O(:,1)=PF_resa%outparticles(:,45)

        PF_resa%resample_weight(:)=PF_resa%outparticles(:,46)

        do p4=1,parloop
           call PF_weight_update_resample(p4,PF_resa%resample_weight(p4))
        end do
        PFweightupdatesum_resample=sum(PF_resa%resample_weight_update(:))

        do p5=1,parloop
          PF_ppar%pfweight(p5,1)=PF_resa%resample_weight_update(p5)/PFweightupdatesum_resample
        end do

     else
        PF_ppar%Vcmax(:,1)= PF_ppar%Vcmax(:,1)
        PF_ppar%VJ_slope(:,1)=PF_ppar%VJ_slope(:,1)
        PF_ppar%VN_slope(:,1)=PF_ppar%VN_slope(:,1)
        PF_ppar%b_h2o(:,1)=PF_ppar%b_h2o(:,1)
        PF_ppar%m_h2o(:,1)=PF_ppar%m_h2o(:,1)
        PF_ppar%f_leaf(:,1)=PF_ppar%f_leaf(:,1)
        PF_ppar%kc25(:,1)=PF_ppar%kc25(:,1)
        PF_ppar%ko25(:,1)=PF_ppar%ko25(:,1)
        PF_ppar%tau25(:,1)=PF_ppar%tau25(:,1)

        PF_ppar%sif_alpha(:,1)=PF_ppar%sif_alpha(:,1)
        PF_ppar%sif_beta(:,1)=PF_ppar%sif_beta(:,1)

        PF_ppar%q10(:,1)=PF_ppar%q10(:,1)
        PF_ppar%f_resp(:,1)=PF_ppar%f_resp(:,1)

        PF_ppar%f_decay(:,1)=PF_ppar%f_decay(:,1)
        PF_ppar%Ksat_scalar(:,1)=PF_ppar%Ksat_scalar(:,1)
        PF_ppar%b_scalar(:,1)=PF_ppar%b_scalar(:,1)
        PF_ppar%porosity_scalar(:,1)=PF_ppar%porosity_scalar(:,1)
        PF_ppar%vfc_scalar(:,1)=PF_ppar%vfc_scalar(:,1)
        PF_ppar%vwp_scalar(:,1)=PF_ppar%vwp_scalar(:,1)
        PF_ppar%psisat_scalar(:,1)=PF_ppar%psisat_scalar(:,1)
        PF_ppar%drainage_scalar(:,1)=PF_ppar%drainage_scalar(:,1)

        PF_ppar%vod_a(:,1)=PF_ppar%vod_a(:,1)
        PF_ppar%vod_b(:,1)=PF_ppar%vod_b(:,1)
        PF_ppar%vod_c(:,1)=PF_ppar%vod_c(:,1)

        PF_ppar%theta_Amin(:,1)=PF_ppar%theta_Amin(:,1)
        PF_ppar%pox(:,1)=PF_ppar%pox(:,1)
        PF_ppar%fei_c(:,1)=PF_ppar%fei_c(:,1)
        PF_ppar%spac_p1(:,1)=PF_ppar%spac_p1(:,1)
        PF_ppar%spac_p2(:,1)=PF_ppar%spac_p2(:,1)
        PF_ppar%tWA(:,1)=PF_ppar%tWA(:,1)
        PF_ppar%tWB(:,1)=PF_ppar%tWB(:,1)
        PF_ppar%Ttrig(:,1)=PF_ppar%Ttrig(:,1)
        PF_ppar%r_xylem(:,1)=PF_ppar%r_xylem(:,1)
        PF_ppar%r_r(:,1)=PF_ppar%r_r(:,1)
        PF_ppar%Lr(:,1)=PF_ppar%Lr(:,1)
        PF_ppar%deltal_min(:,1)=PF_ppar%deltal_min(:,1)
        PF_ppar%deltal_max(:,1)=PF_ppar%deltal_max(:,1)
        PF_ppar%p_delta(:,1)=PF_ppar%p_delta(:,1)
        PF_ppar%ppslh(:,1)=PF_ppar%ppslh(:,1)
        PF_ppar%fei_min(:,1)=PF_ppar%fei_min(:,1)
        PF_ppar%fei_th(:,1)=PF_ppar%fei_th(:,1)
        PF_ppar%p_excess(:,1)=PF_ppar%p_excess(:,1)
        PF_ppar%Tleaf_H(:,1)=PF_ppar%Tleaf_H(:,1)
        PF_ppar%Tleaf_L(:,1)=PF_ppar%Tleaf_L(:,1)
        PF_ppar%Tleaf_O(:,1)=PF_ppar%Tleaf_O(:,1)

        do p5=1,parloop
          PF_ppar%pfweight(p5,1)=PF_ppar%pfweightupdate(p5,1)/PFweightupdatesum
        end do
     end if

     if (is_end_week) then
        call Create_PF_para(3,5,parloop)
        write(6,*)   "call Create_PF_para(3,5,parloop)"
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
