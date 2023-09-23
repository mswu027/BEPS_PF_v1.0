!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file bepsfunc.F90
!> \brief defines a functional implementation for running BEPS model,
!>        and defines the interface for using BEPS within
!>        a sensitivity and/or optimisation framework.
!>        It defines the top-level interface 'evalf' that maps
!>        a one-dimensional control vector (normalised coordinates)
!>        to a one-dimensional simulation vector.
!>        Please note, that the implementation core of this routine
!>        was taken/transferrred from the original BEPS driver.
!>
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  January 2020 (with several updates applied thereafter)
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!***********************************************************
!     evalf
!
!> @brief function like BEPS model evaluation (or simulation)
!>        that maps a given one-dimensional control vector x ("independents")
!>        in normalised units to a one-dimensional simulation vector y ("dependents")
!>        The 1D control vector captures 3 BEPS quantities (SIF, Thetam, COSflux),
!>        ordering in the 1D vector is (varying slowest to fastest)
!>        time (hourly), site, BEPS-quantity.
!>
!> @param[in]   n   length of control vector
!> @param[in]   x   control vector (in normalised units)
!> @param[in]   m   length of simulation vector
!> @param[out]  y   simulation vector
!
!> \authors Michael Vossbeck, The Inversion Lab
!> \date    January 2020
!> \last    April 2021
!>
subroutine evalf(n, x, m, y)
  use shr_kind_mod, only: r8 =>shr_kind_r8
  use beps_con, only:PI
  use controlInput_mod, only:&
       lai_input, meteo_input, sim_duration, restart_frq, nscale, &
       read_meteo_site, read_meteo_daily, read_meteo_hourly, &
       read_lai_site, read_lai
  ! use controlInput_mod
  use bepstype, only:bound, clim, assim
  use ecoRespMod, only:plant_resp, soil_resp
  use mid_results
  use beps_soilMod, only:soil, init_soil_parameters, init_soil_status
  use beps_par
  use outputMod, only:output, av_output
  use restart, only:v2last, restart_io
  use mo_bepsfunc_ctl
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: y(m)
  ! local declarations
  character(len=*), parameter :: sub = 'evalf'
  type(climatedata)    :: meteo
  type(results)        :: mid_res
  real(r8)             :: CosZs,hr_loc,hr_arc  !! solar zenith angle, local time, local time arc
  integer              :: i,j,k,jj,ii,kk
  type(soil)           :: soilp                !! at single point
  real(r8)             :: Ccd(0:4),Cssd(0:4),Csmd(0:4),Cfsd(0:4),Cfmd(0:4),Csm(0:4),Cm(0:4),Cs(0:4),Cp(0:4)
  real(r8)             :: param(0:49),var_o(0:40),var_n(0:40),coef(0:49)
  real(r8)             :: lai

  real(r8)             :: ratio_cloud,shortRad_df,shortRad_dir
  !-- iLab::similar changes as in driver.F90 by MSWU@2020-09-21
  ! real(r8)             :: NPP_yr_acc           !! for storing yearly accumulated NPP, Mh/ha/s
  ! real(r8)             :: agb2vod
  ! real(r8)             :: D0(1:9)
  ! real(r8)             :: taweff(1:9)
  integer              :: kount,rst_nstep
  integer              :: yr,mn,dy,tod,caldy
  integer              :: yr_next, mn_next, dy_next, tod_next
  real(r8)             :: daylen
  ! .. Parameters for model-internal use
  real(r8)             :: pio180 = PI/180.
  ! variables for daily input, used in climin and getmonth
  real(r8)             :: spds, cpds
  ! .. Local Arrays
  real(r8)             :: atmean, atrange
  ! .. Local Scalars ..
  ! real(r8)             :: r
  real(r8)             :: rdaymid, delta, arg, h0, h1, sd, sd1, dhour, tmin, tmp1
  real(r8)             :: a, b, sunset_arc
  integer              :: nd                   !! for counting the time-step number, i.e. ith step
  !-- iLab::ported from MSWU changes () in driver.f90
  integer              :: n_meteo              !-- index being passed to read_meteo_site
  ! .. Intrinsic Functions ..
  intrinsic ACOS, COS, SIN, MOD, ATAN, REAL
  integer :: jcnt !-- counter for 1D simulation vector
  integer :: doys !-- number of days in (actual) year
  logical :: is_end_curr_year, is_end_curr_month, is_end_curr_day, is_first_step

  logical :: bepsf_debug = .false.

  !-- consistency
  if( nscale.ne.1 ) then
     write(*, '(a)') ' FATAL::'//sub//':setup was *not* set for site-level!'
     stop
  endif
  if( meteo_input.lt.0 ) then
     write(*, '(a)') ' FATAL::'//sub//':for site-level only hourly meteorological input'//&
          ' is supported.'
     stop
  endif
  if( lai_input.lt.0 ) then
     write(*, '(a)') ' FATAL::'//sub//':for site-level LAI is expected to be forcing input.'
     stop
  endif

  ! parameters used for calculating VOD,@MOUSONG.WU,2019-11
  !-- iLab::similar changes as in driver.F90 by MSWU@2020-09-21
  ! NPP_yr_acc = 0.
  ! agb2vod = 0.9517
  ! D0 = (/0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05/)
  ! taweff = (/15.0, 15.0, 0., 0., 0., 15., 0., 0., 15./)
  output%NPP_yr_acc = 0._r8

  !-- initialise simulation vector
  y = fillval_r8

  !-- map actual control vector to model parameter
  call x2beps(n,x)

  jcnt = 1
  !! start time looping
  timeloop:do nd=1,ntp
     yr = time_points(1,nd)
     mn = time_points(2,nd)
     dy = time_points(3,nd)
     tod = time_points(4,nd)
     caldy = time_points(5,nd)
     doys = time_points(6,nd)
     is_first_step = (nd.eq.1)
     !--iLab::equations taken from beps_time_manager
     is_end_curr_year  = (mn.eq.1 .and. dy.eq.1 .and. tod.eq.0)
     is_end_curr_month = (dy.eq.1 .and. tod.eq.0)
     is_end_curr_day   = (tod.eq.0)
     call get_CO2_concentration(yr,CO2_air)
     call get_COS_concentration(yr,COS_air)
     !! change hourly input into daily input for further using model for long-term simulations, @MOUSONG.WU, 201905
     if (meteo_input >= 0) then  ! call hourly meteo. input
        if (nscale == 0) then
           call read_meteo_hourly(yr, mn, dy, tod)
        else
           !-- iLab::adapted index (similar to changes in driver.f90 by MSWU,bepspkg_2020-09-21_essi)
           call get_nmeteo(yr, mn, dy, tod, n_meteo)
           call read_meteo_site(n_meteo)
        end if
     else
        if(is_first_step .or. is_end_curr_day) then
           call read_meteo_daily(yr,mn, dy, tod)
        end if
     end if

     if (lai_input >=0) then
        ! print *, 'lai is input!'
        if (is_first_step .or. is_end_curr_day) then
           if (nscale == 0) then
              call read_lai(yr,mn, dy, tod, caldy)
           else
              call read_lai_site(caldy)
           end if
        end if
     else
        ! print *, 'lai is simulated with phenology scheme!'
     end if

     ! Reading the Vcmax here for assimilation usage @J.Wang
     !if(is_first_step() .or. is_end_curr_month()) then
     !   if (nscale == 0) then
     !    call read_Vcmax(yr, mn, dy, tod)
     !   else
     !    call read_Vcmax_site(yr, mn, dy, tod)
     !   end if
     !end if

     !    call mpi_barrier(mpi_comm_world,ierr)
     !! spatial iteration
     pointloop:do i = 1,npoints

        !! calculate the solar zenith
        call s_coszs(yr, mn, dy, tod, caldy, doys, &
             bound%latitude(i),bound%longitude(i),CosZs,hr_loc,hr_arc)
        !!if(myid == 0) write(*,*) "hr_loc=",hr_loc

        !! retrieve meteo for this point
#ifdef COUP_CSM
        meteo%LR                     = clim%Lwdn(i)
        meteo%rainfall               = clim%Rain(i)
        meteo%snow                   = clim%Snow(i)
        meteo%S_dff                  = clim%Swdf(i)
        meteo%S_dir                  = clim%Swdr(i)
        meteo%Srad                   = meteo%S_dff + meteo%S_dir
        meteo%wind                   = clim%Wind(i)
        !              meteo%rh                     = ...
#else
        meteo%Srad                   = clim%Srad(i)
        meteo%wind                   = clim%Wind(i)
        meteo%rainfall               = clim%Rain(i)
        meteo%snow                   = clim%Snow(i)
        meteo%rh                     = clim%Rh(i)
        if (meteo_input < 0) then
           meteo%tempmx                = clim%Tempmx(i)    ! read daily max and min temperatures instead
           meteo%tempmn                = clim%Tempmn(i)
        else
           meteo%temp                  = clim%Temp(i)
        end if

        ! .. compute daily course of temperature and daylength
        rdaymid = REAL (sim_duration+1) / 2.
        delta = -23.4*COS(2.*PI*(rdaymid+10.)/365.)
        spds = SIN(bound%latitude(i)*pio180)*SIN(delta*pio180)
        cpds = COS(bound%latitude(i)*pio180)*COS(delta*pio180)
        arg = -spds/cpds
        IF (arg>1.) THEN
           !        polar night:
           daylen = 0.
        ELSE IF (arg<-1) THEN
           !        polar day:
           daylen = 24.
        ELSE
           !        normal day / night:
           daylen = ACOS(arg)/PI*24.
        END IF

        !###########Compute subdaily temperature based on daily input,@MOUSONG.WU,201905#####################

        if(meteo_input < 0) then
           ! .. compute average conditions
           atmean = (meteo%tempmx + meteo%tempmn)/ 2.
           atrange = meteo%tempmx - meteo%tempmn

           !  hour angle at sunset, added by MOUSONG.WU, 201905
           sunset_arc   = (daylen/2.)*2.0*PI/24.0

           IF (daylen>=4. .AND. daylen<=20.) THEN
              !        sunrise
              h0 = 12. - daylen/2.
              !        sundown
              h1 = 12. + daylen/2.
              !        at sundown:
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
           clim%Temp(i) = meteo%temp
           ! convert daily solar radiation into hourly using the method by M. Collares-Pereira and A. Rabl,
           ! The average distribution of solar radiation-correlations between diffuse and hemispherical
           !and between daily and hourly insolation values, Solar Energy,vol. 22, no. 2, pp. 155164, 1979.
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
        clim%Swdr(i)   = shortRad_dir
        clim%Swdf(i)   = shortRad_df
#endif
        !! PFT iteration
        pftloop: do j = 1,PFT
           if(bound%lcno(i,j) > 0 .and. bound%sw(i) >0. .and. bound%stext(i) >0 .and. bound%clumping(i) > 0.5) then
              call readparam(bound%lcno(i,j),param)
              if(lai_input >= 0) then
                 lai = bound%lai(i,j)
              else 
                 if (is_first_step) then
                    lai = bound%laiyr(i,j)
                    mid_res%lai_old = lai
                 else
                    lai = mid_res%lai_new
                 end if
              end if
              lai = lai*param(2)/bound%clumping(i)
              call readcoef(bound%lcno(i,j),bound%stext(i),coef)

              if(nsrest == nsrStartup .and. is_first_step) then
                 call Init_soil_parameters(bound%lcno(i,j),bound%stext(i), param(27),soilp)
                 !-- iLab::applied similar changes as in driver.F90 (MSWU@2020-09-21)
                 soilp%r_drainage = param(26)
                 ! soilp%r_drainage = assim%p_drainage(j)            ! read this
                 !para from NC file,@MOUSONG.WU,2019-11
                 ii = bound%stext(i)
                 ! write(*,*) 'p_Ksat = ',assim%p_Ksat(ii) 
                 ! write(*,*) 'Ksat_old = ',soilp%Ksat(0)
                 do kk = 0,4
                    soilp%Ksat(kk) = assim%p_Ksat_scalar(ii)*soilp%Ksat(kk)
                    soilp%b(kk)    = assim%p_b_scalar(ii)*soilp%b(kk)
                 end do
                 ! replace these three para. above with values from nc
                 ! file,@MOUSONG.WU,2019-11
                 call Init_soil_status(soilp,bound%st(i),clim%Temp(i),bound%sw(i),bound%sdp(i))
                 do k = 0,40
                    var_o(k)   = 0.
                 end do
                 do k = 3,8
                    var_o(k)   = clim%Temp(i)
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
              else if(nsrest == nsrContinue .and. is_first_step) then
                 call Init_soil_parameters(bound%lcno(i,j),bound%stext(i),param(27),soilp)
                 !-- iLab::applied similar changes as in driver.F90 (MSWU@2020-09-21)
                 soilp%r_drainage = param(26)
                 ! soilp%r_drainage = assim%p_drainage(j)  ! read from nc
                 !                                         ! file,MOUSONG.WU@2019-11
                 ii = bound%stext(i)
                 do kk = 0,4
                    soilp%Ksat(kk) = assim%p_Ksat_scalar(ii)*soilp%Ksat(kk)
                    soilp%b(kk)    = assim%p_b_scalar(ii)*soilp%b(kk)
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
                 !                 do k = 27,32
                 !                    soilp%ice_ratio(k-27)   = var_o(k)
                 !                 end do

              else
                 var_o(:)         = v2last(i,:,j)
                 call retrive_soilp(soilp,i,j,0)
              end if

!!! simulating photosynthesis
              !              do llll = 0,40
              !                 write(*,*)  "DG004: Var_o = ",llll, var_o(llll)
              !              end do

              call inter_prg(yr, mn, dy, tod, &
                   lai,lai_input,bound%lcno(i,j),bound%clumping(i), &
                   assim%p_Vcmax(j),assim%p_VJ_slope(j),&
                   assim%p_f_leaf,assim%p_kc25,assim%p_ko25,assim%p_tau25,&
                   assim%p_sif_alpha(j),assim%p_sif_beta(j),&
                   param,meteo,CosZs, &
                   var_o,var_n,soilp,mid_res,daylen)
              ! CHANGE Vcmax read from Vcmax file with the Vcmax read from initial para. NC
              ! file, for optimization purpose,@MOUSONG.WU,2019-11
              !              do llll =0,40
              !                write(*,*)  "DG004: Var_n = ",llll,var_n(llll)
              !             end do
              v2last(i,:,j)  = var_n(:)
              call retrive_soilp(soilp,i,j,1)
!!! simluating Ra
              call plant_resp(assim%p_q10(j),bound%lcno(i,j),mid_res,bound%laiyr(i,j),lai,meteo%temp,&
                   soilp%temp_soil_c(1),CosZs)
              !USE p_q10 here to adjust q10,p_q10 is read from initial para. NC file, for
              !optimization purpose,@MOUSONG.WU,2019-11
!!! simulating Rh
              Ccd(0)       = bound%Ccd(i,j)
              Cssd(0)      = bound%Cssd(i,j)
              Csmd(0)      = bound%Csmd(i,j)
              Cfsd(0)      = bound%Cfsd(i,j)
              Cfmd(0)      = bound%Cfmd(i,j)
              Csm(0)       = bound%Csm(i,j)
              Cm(0)        = bound%Cm(i,j)
              Cs(0)        = bound%Cm(i,j)
              Cp(0)        = bound%Cp(i,j)
              ! to get soil texture for this point,@MOUSONG.WU,2019-11
              jj = bound%stext(i)

              !-- iLab::similar changes as in driver.F90 (MSWU@2020-09-21),
              !         'beta' no longer parameter
              call soil_resp(Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,&
                   Cp,bound%nppyr(i,j),coef,bound%stext(i),soilp,mid_res)
              ! call soil_resp(Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,&
              !      Cp,bound%nppyr(i,j),coef,bound%stext(i),assim%p_beta(jj),soilp,mid_res)
              ! !! use p_beta read from initial para. NC file, to adjust slow
              ! !carbon pool, for optimization purpose,@MOUSONG.WU,2019-11
              ! mid_res%COS_grnd2 = mid_res%COS_grnd2 * (mid_res%NPP - mid_res%NEP)
              ! ! COS_flux_buf = mid_res%COS_plant + max(mid_res%COS_grnd1,mid_res%COS_grnd2)
              ! mid_res%COS_flux  = mid_res%COS_plant + max(mid_res%COS_grnd1,mid_res%COS_grnd2)

              !! for output variables
              output%GPPpft(i,j)   = mid_res%GPP*bound%PCT_PFT(i,j)/100. ! remove 100. for site since the PCT_PFT is fraction, @MOUSONG.WU,2019-11
              output%SIFpft(i,j)   = mid_res%SIF*bound%PCT_PFT(i,j)/100.
              output%NPPpft(i,j)   = mid_res%NPP*bound%PCT_PFT(i,j)/100.
              output%NEPpft(i,j)   = mid_res%NEP*bound%PCT_PFT(i,j)/100.
              output%SHpft(i,j)    = mid_res%SH*bound%PCT_PFT(i,j)/100.
              output%LHpft(i,j)    = mid_res%LH*bound%PCT_PFT(i,j)/100.
              output%Transpft(i,j) = mid_res%Trans*bound%PCT_PFT(i,j)/100.
              output%Evappft(i,j)  = mid_res%Evap*bound%PCT_PFT(i,j)/100.
              output%Net_Radpft(i,j)  = mid_res%Net_Rad*bound%PCT_PFT(i,j)/100.
              output%LAIpft(i,j)   = lai*bound%PCT_PFT(i,j)/100. 
              output%Thetampft(i,j)   = mid_res%thetam_surf*bound%PCT_PFT(i,j)/100.
              output%fAPARpft(i,j) = mid_res%fAPAR*bound%PCT_PFT(i,j)/100.
              !-- iLab::COS_flux computation now similar as in driver.F90 (MSWU@2020-09-21)
              output%COS_fluxpft(i,j) = (mid_res%COS_plant+mid_res%COS_grnd)*bound%PCT_PFT(i,j)/100.

              !-- iLab::annual NPP, VOD computation now similar as in driver.F90 (MSWU@2020-09-21)
              output%NPP_yr_acc(i,j) = output%NPP_yr_acc(i,j) + mid_res%NPP*bound%PCT_PFT(i,j)/100.*1.e-2*step       ! convert NPP to Mg/ha for calculation of VOD
               
              if( is_end_curr_year ) then
                 ! calculate VOD (vegetation optical depth) with results derived from SMOS-IC product, @Mousong.Wu, 201905,taweff is a PFT specific parameter
                 output%VODpft(i,j)   = assim%p_agb2vod*atan(assim%p_taweff(j)*output%NPP_yr_acc(i,j)) + assim%p_D0(j)*lai
                 output%NPP_yr_acc(i,j) = 0.
                  !                  write(*,*) 'VOD= ',output%VODpft(i,j)
              else
                 output%VODpft(i,j) = 0.
              end if

              ! write(*,*) 'hr_loc = ', hr_loc
              !! calculate the OCO-2 SIF   across at 1:30pm
              !write(*,*) 'SIFpft= ',mid_res%SIF*bound%PCT_PFT(i,j)/100. 
              if(hr_loc >= 13. .and. hr_loc <14.) then
                 output%SIFpft_sat(i,j)  = mid_res%SIF*bound%PCT_PFT(i,j)/100.
              else
                 output%SIFpft_sat(i,j)  = 0.
              end if
              !output%SIFpft_sat(i,j) = max(output%SIFpft_sat(i,j),0.)
              !write(*,*) 'SIFpft_sat = ', output%SIFpft_sat(i,j)
           end if
        end do pftloop !! end PFT loop

        output%GPP(i)   = sum(output%GPPpft(i,:))
        output%SIF(i)   = sum(output%SIFpft(i,:))
        output%SIF_sat(i)  = sum(output%SIFpft_sat(i,:))
        output%NPP(i)   = sum(output%NPPpft(i,:))
        output%NEP(i)   = sum(output%NEPpft(i,:))
        output%SH(i)    = sum(output%SHpft(i,:))
        output%LH(i)    = sum(output%LHpft(i,:))
        output%Trans(i) = sum(output%Transpft(i,:))
        output%Evap(i)  = sum(output%Evappft(i,:))
        output%Net_Rad(i) = sum(output%Net_Radpft(i,:))
        output%LAI(i)     = sum(output%LAIpft(i,:))
        output%Thetam(i)  = sum(output%Thetampft(i,:))
        output%fAPAR(i)   = sum(output%fAPARpft(i,:))
        output%VOD(i)     = sum(output%VODpft(i,:))
        output%COS_flux(i)     = sum(output%COS_fluxpft(i,:))

        !--iLab::mapping to target vector:SIF/Thetam/COSflux ->per time ->per point
        !        (see mo_bepsfunc_ctl.f90 for order of simulated BEPS quantities)
        y(jcnt)   = output%SIF(i)
        y(jcnt+1) = output%Thetam(i)
        y(jcnt+2) = output%COS_flux(i)
        jcnt = jcnt+3
     end do pointloop !! end spatial loop
     !          call mpi_barrier(mpi_comm_world,ierr)

     !! write out data and restart fields
     if( enable_netcdf_out ) then
        call av_output(yr, mn, dy, tod, nd+1, is_end_curr_month, ref_date, seconds_since_ref(nd))
     endif

     if(restart_frq <0) then   !! <0 ndays
        ! kount = get_nstep()
        kount = nd
        rst_nstep  = -restart_frq*86400/step
        yr_next  = time_points(1, kount+1)
        mn_next  = time_points(2, kount+1)
        dy_next  = time_points(3, kount+1)
        tod_next = time_points(4, kount+1)
        if(mod(kount,rst_nstep) ==0) call restart_io("write", yr_next, mn_next, dy_next, tod_next)
     else if(restart_frq ==0) then
        yr_next  = time_points(1, nd+1)
        mn_next  = time_points(2, nd+1)
        dy_next  = time_points(3, nd+1)
        tod_next = time_points(4, nd+1)
        if(is_end_curr_month) call restart_io("write", yr_next, mn_next, dy_next, tod_next)
     end if

  end do timeloop !! end time loop

  !--iLab::target vector should have been written completely (dimension consistency)!
  if( jcnt-1.ne.m ) then
     write(*,'(a,2(a,i4,1x))') ' FATAL::dimension inconsistency for target mapping',&
          'expected m=',m,'got jcnt=',jcnt
  endif

end subroutine evalf


!***********************************************************
!     x2beps
!
!> @brief maps one-dimensional (normalised) control vector
!>        to the respective physical BEPS parameter(s)
!>
!> @param[in]   n   length of control vector
!> @param[in]   x   control vector (in normalised units)
!
!> \authors Michael Vossbbeck, The Inversion Lab
!> \date    February 2020
!>
subroutine x2beps(n,x)
  use mo_prior
  use bepstype, only: assim
  use beps_par, only:PFT,texture
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(in) :: x(n)
  ! local declarations
  real(kind=8) :: xphys(n)
  integer :: i1,i2

  !-- convert to phyiscal units
  xphys = x*x_sigma

  !-- Vcmax
  !   (for the 'Vcmax' parameter)
  i1 = 1
  i2 = i1+PFT-1
  assim%p_Vcmax    = xphys(i1:i2)
  !-- VJ_slope
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  assim%p_VJ_slope = xphys(i1:i2)
  !-- Q10
  !   (parameter to adjust 'q10' in calculation plant respiration)
  i1 = i2+1
  i2 = i1+PFT-1
  assim%p_q10      = xphys(i1:i2)
  !-- SIF alpha
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  assim%p_sif_alpha = xphys(i1:i2)
  !-- SIF beta
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  assim%p_sif_beta = xphys(i1:i2)
  !-- D0
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  assim%p_D0     = xphys(i1:i2)
  !-- taueff
  !   (...)
  i1 = i2+1
  i2 = i1+PFT-1
  assim%p_taweff = xphys(i1:i2)
  !-- Ksat
  !   (parameter to adjust 'Ksat', the saturated hydraulic conductivity, in beps_soilMod.F90)
  i1 = i2+1
  i2 = i1+texture-1
  assim%p_Ksat_scalar = xphys(i1:i2)
  !-- b
  !   (parameter to adjust 'b', the parameters for calculating soil water characteristic, in beps_soilMode.F90)
  i1 = i2+1
  i2 = i1+texture-1
  assim%p_b_scalar    = xphys(i1:i2)
  !-- f_leaf
  !   (...)
  i1 = i2+1
  assim%p_f_leaf = xphys(i1)
  !-- kc25
  !   (...)
  i1 = i1+1
  assim%p_kc25 = xphys(i1)
  !-- ko25
  !   (...)
  i1 = i1+1
  assim%p_ko25 = xphys(i1)
  !-- tau25
  !   (...)
  i1 = i1+1
  assim%p_tau25 = xphys(i1)
  !-- agb2vod
  !   (...)
  i1 = i1+1
  assim%p_agb2vod = xphys(i1)
end subroutine x2beps
