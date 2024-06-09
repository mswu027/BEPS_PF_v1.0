!***************************************************************
! Monte Carlo
!! This module will calculate evaluation metrics for simulation and observation,
!! including Pearson correlation coefficient(R), R_Square, RMSE, ME, AIC, BIC
!! Edited by Lu Hu at Nanjing University
!! Nov. 2023
!***************************************************************
module MC_run
  use shr_kind_mod,only:r8=>shr_kind_r8
  use bepstype
  use bepstypeInit
  use controlInput_mod
  !--iLab::avoid beps_time_manager, all temporal information now passed as actual arguments
  ! use beps_time_manager
  use beps_par
  use beps_con
  implicit none



contains

  subroutine Performance_calc()

    real(r8) :: r2, rmse, me,r
    real(r8) :: y_mean, x_mean, y_std, x_std, cov_xy
    real(r8) :: mse, aic, bic
    real(r8) :: fill_value
    integer   :: i, j, k, effect_num,mc_len,ind_year
    integer  :: num_var, effect_flag
    character(len=*), parameter :: sub = 'MC performace calculation'
    real, dimension (:), allocatable :: y_sim, y_obs, y_sim1, y_obs1

    ! mc_length means the length of data for Monte Carlo analysis
    ! mc_length = sim_duration * 24 in beps_par.F90
    ! using three years data to force the model, the first two years for spinning up, the third year for comparison
    ! need to be check when model period was changed
    ! it takes too long for 2W MC, changed to 2 years,the first year for spinning up, the second year for comparison
    !ind_year = (365+366)*24
    ind_year = 365*24
    mc_len = mc_length-ind_year

    allocate ( y_sim(mc_len) )
    allocate ( y_obs(mc_len) )

    y_sim = mc_output%SM(ind_year+1:mc_length,1)
    y_obs = mc_output%obsSM(ind_year+1:mc_length,1)
    !write(*,*), 'y_sim=', y_sim(1:10)
    !write(*,*), 'y_obs=', y_obs(1:10)

    ! -99999.0 was set as missing values
    fill_value = -99999._r8
    !***** Count the number of valid paired elements between two variables***
    effect_num = 0
    i=1
    do
        if (y_sim(i) /= fill_value .and. y_obs(i) /= fill_value) then
            effect_num = effect_num + 1
        end if
        i=i+1
        if (i > mc_len) exit
    end do
!avoid the case: all the elements are filling value or only one valid values!!!
    if (effect_num <= 1) then
        effect_num = 2
        effect_flag = 1
        print *, 'Number of valid values less than one!'
    end if
    !write(*,*), 'effect_num=', effect_num
    allocate ( y_sim1(effect_num) )
    allocate ( y_obs1(effect_num) )
  !***** extract valid elements between two variables***
    i=1
    j=1
    do
      if ( y_sim(i) /= fill_value .and. y_obs(i) /= fill_value) then
         y_sim1(j) = y_sim(i)
         y_obs1(j) = y_obs(i)
         j=j+1

      end if
      i=i+1
      if (i > mc_len) exit

    end do

    if (effect_flag == 1) then
        me = fill_value
        rmse = fill_value
        r2 = fill_value
        r = fill_value
        aic  = fill_value
        bic = fill_value
    else
        ! calculate the mean bias
        me = sum(y_sim1-y_obs1)/real(effect_num)

        ! calculate the root mean square error
        rmse = 0.
        do k=1, effect_num
           rmse=rmse + (y_sim1(k)-y_obs1(k))**2

        end do
        rmse=rmse/real(effect_num)
        ! avoid the case: floating-point exception
        rmse=sqrt(max(rmse,1.e-12))

    !****** calculate the R2 *******************************************
        y_mean = sum(y_obs1) / real(effect_num)
        x_mean = sum(y_sim1) / real(effect_num)
        y_std = sqrt(max(sum((y_obs1-y_mean)**2) / real(effect_num-1),1.e-12))
        x_std = sqrt(max(sum((y_sim1-x_mean)**2) / real(effect_num-1),1.e-12))
        cov_xy = sum((y_obs1-y_mean) * (y_sim1-x_mean)) / real(effect_num-1)
        r2 = (cov_xy/max(y_std*x_std,1.e-12))**2
        r = cov_xy/max(y_std*x_std,1.e-12)
    ! *****************calculate the AIC and BIC*******************
    ! *Definition: AIC = -2*log(L) + 2*k, BIC = -2*log(L) + k*log(n), L is the likelihood
    ! assuming that errors of the model follow an independent normal distribution
    ! using mean squared errors(MSE) instead of likelihood(L)
    ! AIC = n*LL+2*k, LL is the log-likelihood for the model('A new look at the statistical identification',1974)
    ! BIC = n*LL+k*log(n)
    ! n is the number of examples,k is the number of parameters
        mse = sum((y_obs1-y_sim1)**2)/real(effect_num)
        num_var = 29 ! 29 parameters used for SM and LWP analysis
        aic = real(effect_num)*log(max(mse,1.e-12)) + 2.0*real(num_var)
        bic = real(effect_num)*log(max(mse,1.e-12)) + real(num_var)*log(real(effect_num))
    end if

    !********************************************************************
    MC%SM_R2(1,1) = r2
    MC%SM_RMSE(1,1) = rmse
    MC%SM_ME(1,1) = me
    MC%SM_R(1,1) = r

    MC%SM_AIC(1,1) = aic
    MC%SM_BIC(1,1) = bic

    deallocate (y_sim)
    deallocate (y_obs)
    deallocate (y_sim1)
    deallocate (y_obs1)


  end subroutine Performance_calc


   subroutine LWP_Performance_calc()

    real(r8) :: r2, rmse, me,r
    real(r8) :: y_mean, x_mean, y_std, x_std, cov_xy
    real(r8) :: mse, aic, bic
    real(r8) :: fill_value
    integer   :: i, j, k, effect_num,mc_len,ind_year
    integer  :: num_var, effect_flag
    character(len=*), parameter :: sub = 'MC performace calculation'
    real, dimension (:), allocatable :: y_sim, y_obs, y_sim1, y_obs1

    ! mc_length means the length of data for Monte Carlo analysis
    ! mc_length = sim_duration * 24 in beps_par.F90
    ! using three years data to force the model, the first two years for spinning up, the third year for comparison
    ! need to be check when model period was changed
    ! it takes too long for 2W MC, changed to 2 years,the first two years for spinning up, the second year for comparison
    !ind_year = (365+366)*24
    ind_year = 365*24
    mc_len = mc_length-ind_year

    allocate ( y_sim(mc_len) )
    allocate ( y_obs(mc_len) )

    y_sim = (mc_output%LWP(ind_year+1:mc_length,1)) !  m water
    y_obs = mc_output%obsLWP(ind_year+1:mc_length,1)
    !write(*,*), 'y_sim=', y_sim(1:10)
    !write(*,*), 'y_obs=', y_obs(1:10)

    ! zero was set as missing values
    fill_value = -99999._r8
    !***** Count the number of valid paired elements between two variables***
    effect_num = 0
    i=1
    do
        if (y_sim(i) /= fill_value .and. y_obs(i) /= fill_value) then
            effect_num = effect_num + 1
        end if
        i=i+1
        if (i > mc_len) exit
    end do
!avoid the case: all the elements are filling value or only one valid values!!!
    if (effect_num <= 1) then
        effect_num = 2
        effect_flag = 1
        print *, 'Number of valid values less than one!'
    end if
    write(*,*), 'valid LWP observations=', effect_num
    allocate ( y_sim1(effect_num) )
    allocate ( y_obs1(effect_num) )
  !***** extract valid elements between two variables***
    i=1
    j=1
    do
      if ( y_sim(i) /= fill_value .and. y_obs(i) /= fill_value) then
         y_sim1(j) = y_sim(i)
         y_obs1(j) = y_obs(i)
         j=j+1

      end if
      i=i+1
      if (i > mc_len) exit

    end do

    if (effect_flag == 1) then
        me = fill_value
        rmse = fill_value
        r2 = fill_value
        r = fill_value
        aic  = fill_value
        bic = fill_value
    else
        ! calculate the mean bias
        ! for leaf water potential, unit should be changed!
        y_sim1 = y_sim1/101.*(-1.) ! m water to -MPa
        !y_obs1 = y_obs1 *101.*(-1.) ! -MPa to m water
        me = sum(y_sim1-y_obs1)/real(effect_num)

        ! calculate the root mean square error
        rmse = 0.
        do k=1, effect_num
           rmse=rmse + (y_sim1(k)-y_obs1(k))**2

        end do
        rmse=rmse/real(effect_num)
        ! avoid the case: floating-point exception
        rmse=sqrt(max(rmse,1.e-12))

    !****** calculate the R2 *******************************************
        y_mean = sum(y_obs1) / real(effect_num)
        x_mean = sum(y_sim1) / real(effect_num)
        y_std = sqrt(max(sum((y_obs1-y_mean)**2) / real(effect_num-1),1.e-12))
        x_std = sqrt(max(sum((y_sim1-x_mean)**2) / real(effect_num-1),1.e-12))
        cov_xy = sum((y_obs1-y_mean) * (y_sim1-x_mean)) / real(effect_num-1)
        r2 = (cov_xy/max(y_std*x_std,1.e-12))**2
        r = cov_xy/max(y_std*x_std,1.e-12)
    ! *****************calculate the AIC and BIC*******************
    ! *Definition: AIC = -2*log(L) + 2*k, BIC = -2*log(L) + k*log(n), L is the likelihood
    ! assuming that errors of the model follow an independent normal distribution
    ! using mean squared errors(MSE) instead of likelihood(L)
    ! AIC = n*LL+2*k, LL is the log-likelihood for the model('A new look at the statistical identification',1974)
    ! BIC = n*LL+k*log(n)
    ! n is the number of examples,k is the number of parameters
        mse = sum((y_obs1-y_sim1)**2)/real(effect_num)
        num_var = 29 ! 29 parameters used for SM and LWP analysis
        aic = real(effect_num)*log(max(mse,1.e-12)) + 2.0*real(num_var)
        bic = real(effect_num)*log(max(mse,1.e-12)) + real(num_var)*log(real(effect_num))
    end if

    !********************************************************************
    MC%LWP_R2(1,1) = r2
    MC%LWP_RMSE(1,1) = rmse
    MC%LWP_ME(1,1) = me
    MC%LWP_R(1,1) = r

    MC%LWP_AIC(1,1) = aic
    MC%LWP_BIC(1,1) = bic

    deallocate (y_sim)
    deallocate (y_obs)
    deallocate (y_sim1)
    deallocate (y_obs1)


  end subroutine LWP_Performance_calc

  subroutine write_MC_para(p)
    implicit none
    integer, intent(in) :: p ! corresponding to the mc number
    integer   :: fid
    integer   :: nt,status
    integer   :: ncid,dimid_site,varid !,dimid_mc

    character(len=*), parameter :: sub = 'write_MC_para'

    character(len=255)  :: fln1
    character(len=8)    :: datestr,ppp
    integer :: i
    logical :: ldebug = .False.
    !
    real(r8),dimension(nlp)  :: Vcmax,VJ_slope,VN_slope, b_h2o, m_h2o, f_leaf, kc25, ko25,tau25
    real(r8),dimension(nlp)  :: sif_alpha,sif_beta,q10, f_resp, vod_a, vod_b, vod_c
    real(r8),dimension(nlp)  :: f_decay,Ksat_scalar,b_scalar,porosity_scalar,vfc_scalar,vwp_scalar,&
                                  psisat_scalar,drainage_scalar
    real(r8),dimension(nlp)  :: theta_Amin,pox, fei_c,spac_p1,spac_p2,tWA,tWB,Ttrig,r_xylem,r_r,Lr,&
                                deltal_min,deltal_max,p_delta,ppslh,fei_min,fei_th,p_excess,Tleaf_H,&
                                Tleaf_L,Tleaf_O
    real(r8),dimension(nlp)  ::  SM_R2,SM_RMSE,SM_ME, SM_AIC,SM_BIC,SM_R
    real(r8),dimension(nlp)  ::  LWP_R2,LWP_RMSE,LWP_ME, LWP_AIC,LWP_BIC,LWP_R
    !real(r8),dimension(mc_length) :: simVOD, obsVOD
    Vcmax = MC%Vcmax(1,1)
    VJ_slope = MC%VJ_slope(1,1)
    VN_slope = MC%VN_slope(1,1)
    b_h2o = MC%b_h2o(1,1)
    m_h2o = MC%m_h2o(1,1)
    f_leaf = MC%f_leaf(1,1)
    kc25 = MC%kc25(1,1)
    ko25 = MC%ko25(1,1)
    tau25 = MC%tau25(1,1)
    sif_alpha = MC%sif_alpha(1,1)
    sif_beta = MC%sif_beta(1,1)
    q10 = MC%q10(1,1)
    f_resp = MC%f_resp(1,1)
    vod_a = MC%vod_a(1,1)
    vod_b = MC%vod_b(1,1)
    vod_c = MC%vod_c(1,1)

    f_decay = MC%f_decay(1,1)
    Ksat_scalar=MC%Ksat_scalar(1,1)
    b_scalar=MC%b_scalar(1,1)
    porosity_scalar=MC%porosity_scalar(1,1)
    vfc_scalar=MC%vfc_scalar(1,1)
    vwp_scalar=MC%vwp_scalar(1,1)
    psisat_scalar=MC%psisat_scalar(1,1)
    drainage_scalar=MC%drainage_scalar(1,1)

    theta_Amin=MC%theta_Amin(1,1)
    pox=MC%pox(1,1)
    fei_c=MC%fei_c(1,1)
    spac_p1=MC%spac_p1(1,1)
    spac_p2=MC%spac_p2(1,1)
    tWA=MC%tWA(1,1)
    tWB=MC%tWB(1,1)
    Ttrig=MC%Ttrig(1,1)
    r_xylem=MC%r_xylem(1,1)
    r_r=MC%r_r(1,1)
    Lr=MC%Lr(1,1)
    deltal_min=MC%deltal_min(1,1)
    deltal_max=MC%deltal_max(1,1)
    p_delta=MC%p_delta(1,1)
    ppslh=MC%ppslh(1,1)
    fei_min=MC%fei_min(1,1)
    fei_th=MC%fei_th(1,1)
    p_excess=MC%p_excess(1,1)
    Tleaf_H=MC%Tleaf_H(1,1)
    Tleaf_L=MC%Tleaf_L(1,1)
    Tleaf_O=MC%Tleaf_O(1,1)
    !
    SM_R2=MC%SM_R2(1,1)
    SM_RMSE=MC%SM_RMSE(1,1)
    SM_ME=MC%SM_ME(1,1)
    SM_AIC=MC%SM_AIC(1,1)
    SM_BIC=MC%SM_BIC(1,1)
    SM_R=MC%SM_R(1,1)

    LWP_R2=MC%LWP_R2(1,1)
    LWP_RMSE=MC%LWP_RMSE(1,1)
    LWP_ME=MC%LWP_ME(1,1)
    LWP_AIC=MC%LWP_AIC(1,1)
    LWP_BIC=MC%LWP_BIC(1,1)
    LWP_R=MC%LWP_R(1,1)
    ! 2023/08/03 check the modeled vod and the point mc_out
    !simVOD=mc_output%VOD(1:mc_length,1)
    !obsVOD=mc_output%obsVOD(1:mc_length,1)

    ! 输出文件
    !if(nhtfrq <0) then
       !write(datestr,"(i8)") yy*10000+mm*100+dd
       !nt   = (tod/3600+1)/(-nhtfrq)
       ! !-- iLab::seconds elapsed since reference time (added for time-variable output)
       ! call timemgr_diff_secs(yy_ref*10000+mm_ref*100+dd_ref, tod_ref, yy*10000+mm*100+dd, tod,&
       !      secs_since_ref(1))
    !else if(nhtfrq ==0) then
       !write(datestr,"(i6)") yy*100+mm
       !nt   = 1
    !end if

    !-- iLab::make logging output depend on flag
    if(ldebug) then
       write(*,*) "Writing out MC para file now!"
    endif
    !write(ppp,"(i8)") nt
    !fln1  = trim(beps_out_dir)//"beps_site_MC"//trim(datestr)//"_"//trim(adjustl(ppp))//".nc"
    fln1  = trim(beps_out_dir)//"beps_site_"//"MC_performance1"//".nc"
    !write ()
    status =  nf90_open(fln1,nf90_write,ncid)
    if(status .ne. nf90_noerr) then

       call check(nf90_create(fln1,nf90_share,ncid))
       call check(nf90_def_dim(ncid,"nparameters",nparameters,dimid_site))
       !2023/08/03
       !call check(nf90_def_dim(ncid,"mc_length",mc_length,dimid_mc))
       !-- iLab::added NetCDF output initialisation with prescribed _FillValue
       !         *UPDATE*: disabled for now, since maybe not compatible on
       !                   platform which Mousong is using.
       call check(nf90_def_var(ncid,"Vcmax",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"VJ_slope",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"VN_slope",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"b_h2o",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"m_h2o",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"f_leaf",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"kc25",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"ko25",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"tau25",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"sif_alpha",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"sif_beta",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"q10",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"f_resp",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"vod_a",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"vod_b",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"vod_c",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"f_decay",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"Ksat_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"b_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"porosity_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"vfc_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"vwp_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"psisat_scalar",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"drainage_scalar",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"theta_Amin",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"pox",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"fei_c",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"spac_p1",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"spac_p2",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"tWA",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"tWB",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"Ttrig",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"r_xylem",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"r_r",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"Lr",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"deltal_min",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"deltal_max",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"p_delta",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"ppslh",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"fei_min",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"fei_th",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"p_excess",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"Tleaf_H",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"Tleaf_L",nf90_double,(/dimid_site/),varid))
       call check(nf90_def_var(ncid,"Tleaf_O",nf90_double,(/dimid_site/),varid))

       call check(nf90_def_var(ncid,"SM_R",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","Pearson correlation between observed and modelled SM"))
       call check(nf90_def_var(ncid,"SM_R2",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","R square between observed and modelled SM"))
       call check(nf90_def_var(ncid,"SM_RMSE",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","Root mean square error between observed and modelled SM"))
       call check(nf90_def_var(ncid,"SM_ME",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","Mean bias between observed and modelled SM"))

       call check(nf90_def_var(ncid,"SM_AIC",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","SM: Akaike Information Criterion"))
       call check(nf90_def_var(ncid,"SM_BIC",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","SM: Bayesian Information Criterion"))

       call check(nf90_def_var(ncid,"LWP_R",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","Pearson correlation between observed and modelled LWP"))
       call check(nf90_def_var(ncid,"LWP_R2",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","R square between observed and modelled LWP"))
       call check(nf90_def_var(ncid,"LWP_RMSE",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","Root mean square error between observed and modelled LWP"))
       call check(nf90_def_var(ncid,"LWP_ME",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","Mean bias between observed and modelled LWP"))

       call check(nf90_def_var(ncid,"LWP_AIC",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","LWP: Akaike Information Criterion"))
       call check(nf90_def_var(ncid,"LWP_BIC",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"name","LWP: Bayesian Information Criterion"))


       call check(nf90_enddef(ncid))
    end if

    call check(nf90_inq_varid(ncid,"Vcmax",varid))
    call check(nf90_put_var(ncid,varid,Vcmax,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"VJ_slope",varid))
    call check(nf90_put_var(ncid,varid,VJ_slope,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"VN_slope",varid))
    call check(nf90_put_var(ncid,varid,VN_slope,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"b_h2o",varid))
    call check(nf90_put_var(ncid,varid,b_h2o,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"m_h2o",varid))
    call check(nf90_put_var(ncid,varid,m_h2o,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"f_leaf",varid))
    call check(nf90_put_var(ncid,varid,f_leaf,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"kc25",varid))
    call check(nf90_put_var(ncid,varid,kc25,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"ko25",varid))
    call check(nf90_put_var(ncid,varid,ko25,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"tau25",varid))
    call check(nf90_put_var(ncid,varid,tau25,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"sif_alpha",varid))
    call check(nf90_put_var(ncid,varid,sif_alpha,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"sif_beta",varid))
    call check(nf90_put_var(ncid,varid,sif_beta,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"q10",varid))
    call check(nf90_put_var(ncid,varid,q10,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"f_resp",varid))
    call check(nf90_put_var(ncid,varid,f_resp,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"vod_a",varid))
    call check(nf90_put_var(ncid,varid,vod_a,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"vod_b",varid))
    call check(nf90_put_var(ncid,varid,vod_b,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"vod_c",varid))
    call check(nf90_put_var(ncid,varid,vod_c,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"f_decay",varid))
    call check(nf90_put_var(ncid,varid,f_decay,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Ksat_scalar",varid))
    call check(nf90_put_var(ncid,varid,Ksat_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"b_scalar",varid))
    call check(nf90_put_var(ncid,varid,b_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"porosity_scalar",varid))
    call check(nf90_put_var(ncid,varid,porosity_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"vfc_scalar",varid))
    call check(nf90_put_var(ncid,varid,vfc_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"vwp_scalar",varid))
    call check(nf90_put_var(ncid,varid,vwp_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"psisat_scalar",varid))
    call check(nf90_put_var(ncid,varid,psisat_scalar,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"drainage_scalar",varid))
    call check(nf90_put_var(ncid,varid,drainage_scalar,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"theta_Amin",varid))
    call check(nf90_put_var(ncid,varid,theta_Amin,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"pox",varid))
    call check(nf90_put_var(ncid,varid,pox,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"fei_c",varid))
    call check(nf90_put_var(ncid,varid,fei_c,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"spac_p1",varid))
    call check(nf90_put_var(ncid,varid,spac_p1,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"spac_p2",varid))
    call check(nf90_put_var(ncid,varid,spac_p2,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"tWA",varid))
    call check(nf90_put_var(ncid,varid,tWA,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"tWB",varid))
    call check(nf90_put_var(ncid,varid,tWB,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Ttrig",varid))
    call check(nf90_put_var(ncid,varid,Ttrig,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"r_xylem",varid))
    call check(nf90_put_var(ncid,varid,r_xylem,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"r_r",varid))
    call check(nf90_put_var(ncid,varid,r_r,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Lr",varid))
    call check(nf90_put_var(ncid,varid,Lr,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"deltal_min",varid))
    call check(nf90_put_var(ncid,varid,deltal_min,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"deltal_max",varid))
    call check(nf90_put_var(ncid,varid,deltal_max,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"p_delta",varid))
    call check(nf90_put_var(ncid,varid,p_delta,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"ppslh",varid))
    call check(nf90_put_var(ncid,varid,ppslh,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"fei_min",varid))
    call check(nf90_put_var(ncid,varid,fei_min,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"fei_th",varid))
    call check(nf90_put_var(ncid,varid,fei_th,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"p_excess",varid))
    call check(nf90_put_var(ncid,varid,p_excess,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Tleaf_H",varid))
    call check(nf90_put_var(ncid,varid,Tleaf_H,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Tleaf_L",varid))
    call check(nf90_put_var(ncid,varid,Tleaf_L,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"Tleaf_O",varid))
    call check(nf90_put_var(ncid,varid,Tleaf_O,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"SM_R2",varid))
    call check(nf90_put_var(ncid,varid,SM_R2,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"SM_R",varid))
    call check(nf90_put_var(ncid,varid,SM_R,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"SM_RMSE",varid))
    call check(nf90_put_var(ncid,varid,SM_RMSE,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"SM_ME",varid))
    call check(nf90_put_var(ncid,varid,SM_ME,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"SM_AIC",varid))
    call check(nf90_put_var(ncid,varid,SM_AIC,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"SM_BIC",varid))
    call check(nf90_put_var(ncid,varid,SM_BIC,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"LWP_R2",varid))
    call check(nf90_put_var(ncid,varid,LWP_R2,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"LWP_R",varid))
    call check(nf90_put_var(ncid,varid,LWP_R,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"LWP_RMSE",varid))
    call check(nf90_put_var(ncid,varid,LWP_RMSE,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"LWP_ME",varid))
    call check(nf90_put_var(ncid,varid,LWP_ME,start=(/p/),count=(/1/)))

    call check(nf90_inq_varid(ncid,"LWP_AIC",varid))
    call check(nf90_put_var(ncid,varid,LWP_AIC,start=(/p/),count=(/1/)))
    call check(nf90_inq_varid(ncid,"LWP_BIC",varid))
    call check(nf90_put_var(ncid,varid,LWP_BIC,start=(/p/),count=(/1/)))

    ! 2023/08/03
    !call check(nf90_inq_varid(ncid,"simVOD",varid))
    !call check(nf90_put_var(ncid,varid,simVOD,start=(/p,1/),count=(/1,mc_length/)))
    !call check(nf90_inq_varid(ncid,"obsVOD",varid))
    !call check(nf90_put_var(ncid,varid,obsVOD,start=(/p,1/),count=(/1,mc_length/)))

    call check(nf90_close(ncid))

  end subroutine write_MC_para

  ! for checking the modeling SM
  subroutine write_MCoutput_site(yy,mm,dd,tod,ref_date,secs_since_ref,p)
 ! Note: output the modeling SM for each pair of parameters
    use netcdf
    implicit none
    !--iLab::yy,mm,dd,tod turned into arguments
    integer, intent(in) :: yy,mm,dd,tod
    !integer, intent(in) :: yy,mm,dd,tod,p
    character(len=*), intent(in) :: ref_date
    real(r8), intent(in) :: secs_since_ref
    integer, intent(in) :: p ! corresponding to MC times
    character(len=*), parameter :: sub = 'write_MCoutput_site'
    real(r8),dimension(nlp)        :: Thetam1

    integer   :: ierr
    integer   :: ncid,dimid_site,dimid_time,dimid_mc,varid
    integer   :: nsite(nlp)
    !character(len=8)    :: datestr,ppp
    character(len=8)    :: datestr
    character(len=255)  :: fln1,fln2,name,unit
    integer   :: nt,status
    integer   :: i
    !-- iLab::reduce amount of terminal output (to be reactivated on purpose)
    logical :: ldebug = .False.
    !-- iLab::added for consistent initialisation of NetCDF variables
    real(r8), parameter :: fill_value = -99999._r8
    ! initialization
    Thetam1=0.
    Thetam1 = mc_output%tempSM

    !--iLab::yy,mm,dd,tod are arguments now
    ! call get_prev_date(yy,mm,dd,tod)
    if(nhtfrq <0) then
       !datestr="US_Ne1"
       datestr='US-MOz' ! site name as output NC file
       !write(datestr,"(i8)") yy
       !write(datestr,"(i8)") yy*10000+mm*100+dd
       nt   = (secs_since_ref/3600+1)/(-nhtfrq)
       !nt   = (tod/3600+1)/(-nhtfrq)
       ! !-- iLab::seconds elapsed since reference time (added for time-variable output)
       ! call timemgr_diff_secs(yy_ref*10000+mm_ref*100+dd_ref, tod_ref, yy*10000+mm*100+dd, tod,&
       !      secs_since_ref(1))
    else if(nhtfrq ==0) then
       write(datestr,"(i6)") yy*100+mm
       nt   = 1
    end if

    !-- iLab::make logging output depend on flag
    if(ldebug) then
       write(*,*) "Writing out simulation file now!"
    endif
    !write(ppp,"(i8)") p
    !fln1  = trim(beps_out_dir)//"beps_site_"//trim(adjustl(ppp))//"_"//trim(datestr)//".nc"
    fln1  = trim(beps_out_dir)//"beps_site_"//trim(datestr)//".nc"
    status =  nf90_open(fln1,nf90_write,ncid)
    if(status .ne. nf90_noerr) then

       call check(nf90_create(fln1,nf90_share,ncid))
       call check(nf90_def_dim(ncid,"nsite",nlp,dimid_site))
       call check(nf90_def_dim(ncid,"MCtimes",nparameters,dimid_mc))
       call check(nf90_def_dim(ncid,"time",nf90_unlimited,dimid_time))

       call check(nf90_def_var(ncid,"time",nf90_double,(/dimid_time/),varid))
       call check(nf90_put_att(ncid, varid, "long_name", "time"))
       !-- iLab::added for proper time-variable (hourly output only)
       if( nhtfrq<0 ) then
            call check(nf90_put_att(ncid, varid, "units", "seconds since "//ref_date))
            !call check(nf90_put_att(ncid, varid, "units", "hours since "//ref_date))
       endif
       call check(nf90_put_att(ncid,varid,"calendar","Gregorian"))

       !-- iLab::added NetCDF output initialisation with prescribed _FillValue
       !         *UPDATE*: disabled for now, since maybe not compatible on
       !                   platform which Mousong is using.
       call check(nf90_def_var(ncid,"nsite",nf90_double,(/dimid_site/),varid))
       call check(nf90_put_att(ncid,varid,"units","-"))
       call check(nf90_put_att(ncid,varid,"long_name","site_number"))
       call check(nf90_put_att(ncid,varid,"axis","X"))

       call check(nf90_def_var(ncid,"Thetam",nf90_double,(/dimid_site,dimid_mc,dimid_time/),varid))
       ! call check(nf90_def_var_fill(ncid, varid,  NF90_FILL, fill_value))
       call check(nf90_put_att(ncid,varid,"units","m3/m3"))
       call check(nf90_put_att(ncid,varid,"missing_value",fill_value))
       call check(nf90_put_att(ncid,varid,"name","Soil moisture in layer 1 (5 cm)"))

       call check(nf90_enddef(ncid))
    end if
    !! For temporary output  , Should be improved later @J.Wang
    do i = 1,nlp
       nsite(i) = i
    end do

    call check(nf90_inq_varid(ncid,"nsite",varid))
    call check(nf90_put_var(ncid,varid,nsite))

    !-- iLab::added writing of time-values (hourly output only)
    if( nhtfrq <0 ) then
       call check(nf90_inq_varid(ncid,"time",varid))
       call check(nf90_put_var(ncid,varid,(/secs_since_ref/),start=(/nt/),count=(/1/)))
    endif

    call check(nf90_inq_varid(ncid,"Thetam",varid))
    call check(nf90_put_var(ncid,varid,Thetam1,start=(/1,p,nt/),count=(/nlp,1,1/)))

    call check(nf90_close(ncid))
    !end if
    !call mpi_barrier(mpi_comm_world,ierr)


  end subroutine write_MCoutput_site

end module MC_run
