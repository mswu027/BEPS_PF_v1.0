!***************************************************
! Function : calculate the photosyntheis,leaf temp,
!            and the evapotraspiration
! created  : Jun Wang
! Date     : 2016/12/5
!***************************************************
module AnGsMod
use shr_kind_mod,only:r8=>shr_kind_r8
use beps_con
use beps_par
use mid_results
implicit none

type(meteorology)  :: met
type(factors)      :: fact
type(boundary_layer_resistances)   :: bound_layer_res

public :: lai_cal,       &
          Vcmax_Jmax,    &
          photosynthesis,&
          transpiration, &
          evaporation_canopy

private :: SFC_VPD, &
           TEMP_FUNC, &
           LAMBDA,  &
           ES, &
           TBOLTZ, &
           TBOLTZC4, &
           fluorescence, &
           VT_FUNC, &
           RT_FUNC, &
           KT_FUNC

contains

subroutine lai_cal ( stem_o,stem_u,LC,CosZs,lai_o,clumping,lai_u, &
             lai_o_sunlit,lai_o_shaded,lai_u_sunlit,lai_u_shaded, &
             PAI_o_sunlit,PAI_o_shaded,PAI_u_sunlit,PAI_u_shaded )
implicit none

    real(r8)  :: stem_o       !overstory woody area
    real(r8)  :: stem_u       !understory woody area
    integer   :: LC           !landcover type
    real(r8)  :: CosZs        !cosine solar zenith angle
    real(r8)  :: lai_o        !overstory lai
    real(r8)  :: clumping     !clumping index
    real(r8)  :: lai_u        !understory lai
    real(r8)  :: lai_o_sunlit !overstory sunlit lai
    real(r8)  :: lai_o_shaded !overstory shaded lai
    real(r8)  :: lai_u_sunlit ! understory sunlit lai
    real(r8)  :: lai_u_shaded !understory shaded lai

    real(r8)  :: PAI_o_sunlit !overstory sunlit lai
    real(r8)  :: PAI_o_shaded !overstory shaded lai
    real(r8)  :: PAI_u_sunlit !understory sunlit lai
    real(r8)  :: PAI_u_shaded !understory shaded lai
    real(r8)  :: temp

    if(CosZs>0)then
        temp = -0.5*clumping*(lai_o+stem_o)/CosZs
        if(temp < -10) then
           PAI_o_sunlit = 2*CosZs
        else
           PAI_o_sunlit=2*CosZs*(1-exp(-0.5*clumping*(lai_o+stem_o)/CosZs))
        end if
    else
        PAI_o_sunlit =0
    end if

    PAI_o_shaded =(lai_o+stem_o)-PAI_o_sunlit

    if(CosZs>0)then
      temp = -0.5*clumping*(lai_o+stem_o+lai_u+stem_u)/CosZs
      if(temp < -10) then
         PAI_u_sunlit=2*CosZs-PAI_o_sunlit
       else
         PAI_u_sunlit=2*CosZs*(1-exp(-0.5*clumping*(lai_o+stem_o+lai_u+stem_u)/CosZs))-PAI_o_sunlit
       end if
    else
        PAI_u_sunlit =0
    end if

    PAI_u_shaded =(lai_u+stem_u)-PAI_u_sunlit

    if(CosZs>0)then
        temp = -0.5*clumping*lai_o/CosZs
        if(temp < -10) then
           lai_o_sunlit=2*CosZs
        else
           lai_o_sunlit=2*CosZs*(1-exp(-0.5*clumping*lai_o/CosZs))
        end if
    else
        lai_o_sunlit=0
    end if

    lai_o_shaded=max(0.,lai_o-PAI_o_sunlit )

    if(CosZs>0)then
        temp = -0.5*clumping*(lai_o+lai_u)/CosZs
        if(temp < -10) then
          lai_u_sunlit=2*CosZs-lai_o_sunlit
        else
          lai_u_sunlit=2*CosZs*(1-exp(-0.5*clumping*(lai_o+lai_u)/CosZs))-lai_o_sunlit
         end if
    else
        lai_u_sunlit =0
    end if

    lai_u_shaded=max(0.,lai_u-PAI_u_sunlit )

end subroutine

!*************************************************************************************************
! Subroutine to calculate the Vcmax and Jmax for sunlit and shaded big-leaf

! Reference:
!       (1) Chen, J. M., G. Mo, J. Pisek, F. Deng, M. Ishozawa, D. Chan, 2012.
!       Effects of foliage clumping on global terrestrial gross primary
!       productivity. Global Biogeochemical Cycles, VOL. 26, GB1019, 18,
!       doi:10.1029/2010GB003996
!       (2) Medlyn, B.E. et al., 1999. Effects of elevated [CO2] on photosynthesis
!   in European forest species: a meta-analysis of model parameters.
!   Plant, Cell & Environment, 22(12): 1475-1495.

!*************************************************************************************************



subroutine Vcmax_Jmax(lai_o, clumping, Vcmax0, VJ_slope,slope_Vcmax_N, leaf_N, CosZs,&
                      Vcmax_sunlit, Vcmax_shaded, Jmax_sunlit, Jmax_shaded)
implicit none

    real(r8)  :: lai_o,clumping,Vcmax0,slope_Vcmax_N
    real(r8)  :: leaf_N,CosZs
    real(r8)  :: Vcmax_sunlit,Vcmax_shaded
    real(r8)  :: Jmax_sunlit,Jmax_shaded
    real(r8)  :: VJ_slope
    real(r8)  :: K,expr1,expr2,expr3
    real(r8),parameter:: Kn = 0.3
    real(r8),parameter:: G_theta = 0.5

    if(lai_o < 0.001)then
        lai_o = 0.001
    end if

    if (CosZs>0)then

        K = G_theta*clumping/CosZs

        if(K > 10) then
          expr1 = 1.
          expr2 = 1.
        else
          expr1 = 1 - exp(-K*lai_o)
          expr2 = 1 - exp(-(Kn + K)*lai_o)
        end if
        expr3 = 1 - exp(-Kn*lai_o)

        if(expr1>0)then
            Vcmax_sunlit = Vcmax0 * slope_Vcmax_N * leaf_N * K*expr2 /(Kn + K) / expr1
        else
            Vcmax_sunlit = Vcmax0
        end if


        if ((K > 0).and.(lai_o > expr1/K))then
            Vcmax_shaded = Vcmax0 * slope_Vcmax_N * leaf_N * (expr3 / Kn - expr2 / (Kn + K)) / (lai_o - expr1/K)
        else
            Vcmax_shaded = Vcmax0
        end if

    else
        Vcmax_sunlit = Vcmax0
        Vcmax_shaded = Vcmax0
    end if
    ! VJ_slope = 2.39, currently was not used!
    !Jmax_sunlit = Vcmax_sunlit * 2.39 - 14.2
    !Jmax_shaded = Vcmax_shaded * 2.39 - 14.2
    ! 2024/03/17, when water stress works on Vcmax, using the additive Eqs could avoid the negative values of jmax
    Jmax_sunlit = Vcmax_sunlit * 1.64 + 29.1
    Jmax_shaded = Vcmax_shaded * 1.64 + 29.1


end subroutine


subroutine leaf_temperatures(Tair, slope, psychrometer, VPD_air, Cp_ca,   &
               Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded,        &
               Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded,    &
               Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded,        &
               Xcs_o, Xcl_o, Xcs_u, Xcl_u, radiation_o_sun,               &
               radiation_o_shaded, radiation_u_sun, radiation_u_shaded,   &
               Tc_o_sunlit, Tc_o_shaded, Tc_u_sunlit, Tc_u_shaded )
implicit none
!*************************************************************************
! Subroutine to calculate the sunlit and shaded leaf temperatures for
! overstory and understory leave.
!*************************************************************************

    real(r8)  :: Tair, slope, psychrometer, VPD_air, Cp_ca
    real(r8)  :: Gw_o_sunlit, Gw_o_shaded, Gw_u_sunlit, Gw_u_shaded
    real(r8)  :: Gww_o_sunlit, Gww_o_shaded, Gww_u_sunlit, Gww_u_shaded
    real(r8)  :: Gh_o_sunlit, Gh_o_shaded, Gh_u_sunlit, Gh_u_shaded
    real(r8)  :: Xcs_o, Xcl_o, Xcs_u, Xcl_u
    real(r8)  :: radiation_o_sun, radiation_o_shaded, radiation_u_sun, radiation_u_shaded
    real(r8)  :: Tc_o_sunlit, Tc_o_shaded, Tc_u_sunlit, Tc_u_shaded

    call Leaf_Temperature(Tair, slope, psychrometer, VPD_air,Cp_ca, Gw_o_sunlit, Gww_o_sunlit,  &
                          Gh_o_sunlit, Xcs_o, Xcl_o, radiation_o_sun, Tc_o_sunlit )

    call Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw_o_shaded, Gww_o_shaded, &
                          Gh_o_shaded, Xcs_o, Xcl_o, radiation_o_shaded, Tc_o_shaded )

    call Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw_u_sunlit, Gww_u_sunlit, &
                          Gh_u_sunlit, Xcs_u, Xcl_u, radiation_u_sun, Tc_u_sunlit )

    call Leaf_Temperature(Tair, slope, psychrometer, VPD_air, Cp_ca, Gw_u_shaded, Gww_u_shaded, &
                          Gh_u_shaded, Xcs_u, Xcl_u, radiation_u_shaded, Tc_u_shaded  )


end subroutine


subroutine leaf_temperature( Tair, slope, psychrometer, VPD_air, Cp_ca, &
                             Gw, Gww, Gh, Xcs, Xcl, radiation, Tc )
implicit none
    real(r8)  :: Tair, slope, psychrometer, VPD_air, Cp_ca
    real(r8)  :: Gw, Gww, Gh, Xcs, Xcl, radiation
    real(r8)  :: p_star, Tc, R

!    R = 1. / Gw + 1. / (Gww* (Xcs + Xcl))
    p_star = (Gw + Gww* (Xcs + Xcl) ) / psychrometer

    Tc = Tair + (radiation - VPD_air * rho_a * Cp_ca * p_star) / (rho_a *Cp_ca * (Gh + slope * p_star) )
!   temperature of sunlit leaves, overstorey

    Tc = max(Tair - 3.0, Tc)
    Tc = min(Tair + 5.0, Tc)

end subroutine


subroutine transpiration(tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit,    &
               tempL_u_shaded, temp_air, rh_air, Gtrans_o_sunlit,           &
               Gtrans_o_shaded, Gtrans_u_sunlit, Gtrans_u_shaded,           &
               lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded,      &
               trans_o, trans_u,trans_o_sunlit, trans_o_shaded,&
               trans_u_sunlit,trans_u_shaded)
use meteoMod
implicit none

    real(r8)  ::  tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded
    real(r8)  ::  temp_air, rh_air
    real(r8)  ::  Gtrans_o_sunlit, Gtrans_o_shaded, Gtrans_u_sunlit, Gtrans_u_shaded
    real(r8)  ::  lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded
    real(r8)  ::  trans_o, trans_u
    ! 2024/03/30
    real(r8)  ::  trans_o_sunlit, trans_o_shaded,trans_u_sunlit,trans_u_shaded

    real(r8)  ::  LHt_o_sunlit, LHt_o_shaded, LHt_u_sunlit, LHt_u_shaded  !latent heat from leaves W/m2

    call meteo_pack(temp_air,rh_air)

    LHt_o_sunlit =(vpd+slope_vapor *(tempL_o_sunlit -temp_air))*density_air*cp_air*Gtrans_o_sunlit /psy
    LHt_o_shaded =(vpd+slope_vapor *(tempL_o_shaded -temp_air))*density_air*cp_air*Gtrans_o_shaded /psy

    LHt_u_sunlit =(vpd+slope_vapor *(tempL_u_sunlit -temp_air))*density_air*cp_air*Gtrans_u_sunlit /psy
    LHt_u_shaded =(vpd+slope_vapor *(tempL_u_shaded -temp_air))*density_air*cp_air*Gtrans_u_shaded /psy

    trans_o =1/latent_water *(LHt_o_sunlit *lai_o_sunlit +LHt_o_shaded*lai_o_shaded )
    trans_u =1/latent_water *(LHt_u_sunlit *lai_u_sunlit +LHt_u_shaded*lai_u_shaded )
    ! 2024/03/30
    trans_o_sunlit = 1/latent_water * LHt_o_sunlit * lai_o_sunlit
    trans_o_shaded = 1/latent_water * LHt_o_shaded * lai_o_shaded
    trans_u_sunlit = 1/latent_water * LHt_u_sunlit * lai_u_sunlit
    trans_u_shaded = 1/latent_water * LHt_u_shaded * lai_u_shaded

end subroutine


subroutine evaporation_canopy(tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit,        &
                tempL_u_shaded, temp_air, rh_air, Gwater_o_sunlit, Gwater_o_shaded,  &
                Gwater_u_sunlit, Gwater_u_shaded, lai_o_sunlit, lai_o_shaded,        &
                lai_u_sunlit, lai_u_shaded, percent_water_o, percent_water_u,        &
                percent_snow_o, percent_snow_u, evapo_water_o, evapo_water_u,        &
                evapo_snow_o, evapo_snow_u )
use meteoMod
implicit none
!*************************************************************************
! this module calculates evaporation and sublimation from canopy, from
!overstorey understorey sunlit and shaded

! input includes:
!temperature of sunlit and shaded leaves from other storey (leaf temperature module).
!temperature of air, relative humidity,
!aerodynamic conductance of water (snow) for sunlit shaded leaves from overstorey
!and understorey
!percentage of overstorey or understorey covered by water or snow
!leaf area index, sunlit and shaded, overstorey and understorey
!(from leaf area index module)

! output:
!evaporation of water and snow from overstorey and understorey
!*************************************************************************

    real(r8)  ::  tempL_o_sunlit, tempL_o_shaded, tempL_u_sunlit, tempL_u_shaded
    real(r8)  ::  temp_air, rh_air
    real(r8)  ::  Gwater_o_sunlit, Gwater_o_shaded, Gwater_u_sunlit, Gwater_u_shaded
    real(r8)  ::  lai_o_sunlit, lai_o_shaded, lai_u_sunlit, lai_u_shaded
    real(r8)  ::  percent_water_o, percent_water_u, percent_snow_o, percent_snow_u
    real(r8)  ::  evapo_water_o, evapo_water_u, evapo_snow_o, evapo_snow_u

    real(r8)  ::  LHw_o_sunlit, LHw_o_shaded, LHw_u_sunlit, LHw_u_shaded
    real(r8)  ::  LHs_o_sunlit, LHs_o_shaded, LHs_u_sunlit, LHs_u_shaded

    call meteo_pack(temp_air, rh_air)

    LHw_o_sunlit =percent_water_o*(vpd+slope_vapor *(tempL_o_sunlit  &
                  -temp_air ))*density_air*cp_air*Gwater_o_sunlit /psy
    LHw_o_shaded =percent_water_o*(vpd+slope_vapor *(tempL_o_shaded  &
                  -temp_air ))*density_air*cp_air*Gwater_o_shaded /psy

    LHw_u_sunlit =percent_water_u*(vpd+slope_vapor *(tempL_u_sunlit  &
                  -temp_air ))*density_air*cp_air*Gwater_u_sunlit /psy
    LHw_u_shaded =percent_water_u*(vpd+slope_vapor *(tempL_u_shaded  &
                  -temp_air ))*density_air*cp_air*Gwater_u_shaded /psy

    LHs_o_sunlit =percent_snow_o*(vpd+slope_vapor *(tempL_o_sunlit   &
                  -temp_air ))*density_air*cp_air*Gwater_o_sunlit /psy
    LHs_o_shaded =percent_snow_o*(vpd+slope_vapor *(tempL_o_shaded   &
                  -temp_air ))*density_air*cp_air*Gwater_o_shaded /psy

    LHs_u_sunlit =percent_snow_u*(vpd+slope_vapor *(tempL_u_sunlit   &
                  -temp_air ))*density_air*cp_air*Gwater_u_sunlit /psy
    LHs_u_shaded =percent_snow_u*(vpd+slope_vapor *(tempL_u_shaded   &
                  -temp_air ))*density_air*cp_air*Gwater_u_shaded /psy

    evapo_water_o =1./(latent_water )*(max(LHw_o_sunlit,1.e-6) *lai_o_sunlit + &
                   max(LHw_o_shaded,1.e-6) * lai_o_shaded )
    evapo_water_u =1./(latent_water )*(max(LHw_u_sunlit,1.e-6) *lai_u_sunlit + &
                   max(LHw_u_shaded,1.e-6) * lai_u_shaded )
    evapo_snow_o =1/(latent_snow )*(max(LHs_o_sunlit,1.e-6) *lai_o_sunlit + &
                  max(LHs_o_shaded,1.e-6) * lai_o_shaded )
    evapo_snow_u =1/(latent_snow )*(max(LHs_u_sunlit,1.e-6) *lai_u_sunlit + &
                   max(LHs_u_shaded,1.e-6) * lai_u_shaded )


    evapo_water_o=max(0.,evapo_water_o)
    evapo_water_u=max(0.,evapo_water_u)
    evapo_snow_o=max(0.,evapo_snow_o)
    evapo_snow_u=max(0.,evapo_snow_u)

end subroutine

subroutine photosynthesis(LC,temp_leaf_p,f_leaf,p_kc25,p_ko25,p_tau25,rad_leaf,e_air,&
                         g_lb_w,vc_opt,vj_slope,f_soilwater,b_h2o,m_h2o,cii,temp_leaf_c,&
                         LH_leaf,Gs_w,gs_h2o_mole,aphoto,ci,ffpa,sif_alpha,sif_beta,xSIF,&
                         cosii,cosi,cos_assim,f_feileaf,fw_flag)


!*************************************************************************
! This program solves a cubic equation to calculate
! leaf photosynthesis.  This cubic expression is derived from solving
! five simultaneous equations for A, PG, cs, CI and GS.
! Stomatal conductance is computed with the Ball-Berry model.
! The cubic derivation assumes that b', the intercept of the Ball-Berry
! stomatal conductance model, is non-zero.

! Gs = k A rh/cs + b'


! We also found that the solution for A can be obtained by a quadratic equation
! when Gs is constant or b' is zero.


! The derivation is published in:

! Baldocchi, D.D. 1994. An analytical solution for coupled leaf photosynthesis
! and stomatal conductance models. Tree Physiology 14: 1069-1079.


!-----------------------------------------------------------------------

! A Biochemical Model of C3 Photosynthesis

! After Farquhar, von Caemmerer and Berry (1980) Planta.
! 149: 78-90.

! The original program was modified to incorporate functions and parameters
! derived from gas exchange experiments of Harley, who paramertized Vc and J in
! terms of optimal temperature, rather than some reference temperature, eg 25C.

! Program calculates leaf photosynthesis from biochemical parameters

! rd25 - Dark respiration at 25 degrees C (umol m-2 s-1)
! tlk - leaf temperature, Kelvin
! jmax - optimal rate of electron transport
! vcopt - maximum rate of RuBP Carboxylase/oxygenase
! iphoton - incident photosynthetically active photon flux (umols m-2 s-1)

! note: Harley parameterized the model on the basis of incident PAR

! gs - stomatal conductance (mols m-2 s-1), typically 0.01-0.20
! pstat-station pressure, bars
! aphoto - net photosynthesis  (umol m-2 s-1)
! ps - gross photosynthesis (umol m-2 s-1)
! aps - net photosynthesis (mg m-2 s-1)
! aphoto (umol m-2 s-1)

!-----------------------------------------------------------------------

! iphoton is radiation incident on leaves

! The temperature dependency of the kinetic properties of
! RUBISCO are compensated for using the Arrhenius and
! Boltzmann equations.  From biochemistry, one observes that
! at moderate temperatures enzyme kinetic rates increase
! with temperature.  At extreme temperatures enzyme
! denaturization occurs and rates must decrease.

! Arrhenius Eq.

! f(T)=f(tk_25) exp(tk -298)eact/(298 R tk)), where eact is the
! activation energy.

! Boltzmann distribution

! F(T)=tboltz)

! Define terms for calculation of gross photosynthesis, PG

! PG is a function of the minimum of RuBP saturated rate of
! carboxylation, Wc, and the RuBP limited rate of carboxylation, Wj.
! Wj is limiting when light is low and electron transport, which
! re-generates RuBP, is limiting.  Wc is limiting when plenty of RuBP is
! available compared to the CO2 that is needed for carboxylation.

! Both equations take the form:

! PG-photorespiration= (a CI-a d)/(e CI + b)

! PG-photorespiration=min[Wj,Wc] (1-gamma_ps/Ci)

! Wc=Vcmax Ci/(Ci + Kc(1+O2/Ko))

! Wj=J Ci/(4 Ci + 8 gamma_ps)

! Ps kinetic coefficients from Harley at WBW.

! gamma_ps is the CO2 compensation point

! Jan 14, 1999 Updated the cubic solutions for photosynthesis.  There are
! times when the restriction that R^2 < Q^3 is violated.  I therefore need
! alternative algorithms to solve for the correct root.

!*************************************************************************
    implicit none

    real(r8)  :: temp_leaf_p     ! It is temporay, will be removed later
    real(r8)  :: f_leaf          ! Leaf respiration ratio, default value of 0.5
    real(r8)  :: rad_leaf        ! net shortwave radiation (W/m2)
    real(r8)  :: e_air           ! water vapor pressure above canopy (kPa)
    real(r8)  :: g_lb_w          ! leaf laminar boundary layer condunctance to H2O (m/s)
    real(r8)  :: vc_opt          ! the maximum velocities of carboxylation of Rubisco at 25 deg C (umol m-2 s-1)
    real(r8)  :: vj_slope        ! the ratio of Vmax/Jmax, default value of 2.39
    real(r8)  :: f_soilwater     ! an empirical scalar of soil water stress on stomatal conductance, dimensionless
    real(r8)  :: b_h2o           ! the intercept term in BWB model (mol H2O m-2 s-1)
    real(r8)  :: m_h2o           ! the slope in BWB model
    integer   :: LC              !landcover type
    real(r8)  :: p_kc25,p_ko25,p_tau25       ! replace, kc25,ko25,tau25 with three pramters, default values are 274.6,419.8,2904.12, respectively
    real(r8)  :: cii             ! initial intercellular co2 concentration (ppm)
    real(r8)  :: temp_leaf_c     ! leaf temperature (deg C)

    real(r8)  :: LH_leaf         ! leaf latent heat flux (W m-2)

    real(r8)  :: Gs_w           ! stomatal conductance to water vapor (m s-1)
    real(r8)  :: aphoto         ! net photosynthesis rate (umol CO2 m-2 s-1)
    real(r8)  :: ci             ! intercellular co2 concentration (ppm)
!    real(r8)  :: xSIF           ! SIF
    real(r8)  :: sif_alpha,sif_beta,xSIF           ! SIF
    real(r8),parameter  :: air_pres=101.325     !  air pressure (kPa)
    real(r8)  :: ca              ! atmospheric co2 concentration (ppm)
    real(r8)  :: iphoton         ! incident photosynthetic photon flux density (PPFD) umol m-2 s-1
!    real(r8)  :: g_lb_c          ! leaf laminar boundary layer condunctance to CO2 (mol m-2 s-1)
    real(r8)  :: g_lb_c,g_lb_cos          ! leaf laminar boundary layer condunctance to CO2,COS (mol m-2 s-1)
    real(r8)  :: rh_leaf         ! relative humidity at leaf surface (0-1)
    real(r8)  :: temp_leaf_K     ! leaf temperature (K)
    real(r8)  :: gs_co2_mole     ! stomatal conductance to CO2 (mol m-2 s-1)
    real(r8)  :: gs_h2o_mole     ! stomatal conductance to h2o (mol m-2 s-1)
    real(r8)  :: bc              ! temporary variable
    real(r8)  :: cs              ! CO2 concentration at leaf surface (ppm)

    real(r8)  :: b_co2           ! the intercept term in BWB model (mol CO2 m-2 s-1): b_h2o/1.6
    real(r8)  :: m_co2           ! the slope in BWB model: m_h2o/1.6

    real(r8)  :: gammac          ! CO2 compensation point (ppm)
    real(r8)  :: jmopt           ! the maximum potential electron transport rate at 25 deg C (umol m-2 s-1)
    real(r8)  :: jmax            ! the maximum potential electron transport rate (umol m-2 s-1)
    real(r8)  :: vcmax           ! the maximum velocities of carboxylation of Rubisco (umol m-2 s-1)
    real(r8)  :: km_co2          ! Michaelis-Menten constant for CO2 (?ol mol-1)
    real(r8)  :: km_o2           ! Michaelis-Menten constant for O2 (mmol mol-1)
    real(r8)  :: tau             ! the specifity of Rubisco for CO2 compared with O2
    real(r8)  :: resp_ld         ! leaf dark respiration (umol m-2 s-1)
    real(r8)  :: resp_ld25       ! leaf dark respiration at 25 deg C (umol m-2 s-1)


    real(r8)  :: j_photon        ! the flux of electrons through the thylakoid membrane (umol m-2 s-1)
    real(r8)  :: alpha_ps
    real(r8)  :: beta_ps
    real(r8)  :: gamma_ps
    real(r8)  :: theta_ps

    real(r8)  :: denom
    real(r8)  :: p_cubic
    real(r8)  :: q_cubic
    real(r8)  :: r_cubic

    real(r8)  :: Qroot
    real(r8)  :: Rroot
    real(r8)  :: root1,root2,root3
    real(r8)  :: ang_L

    real(r8)  :: j_sucrose        ! net photosynthesis rate limited by sucrose synthesis (umol m-2 s-1)
    real(r8)  :: wc,wj,psguess    ! gross photosynthesis rate limited by light (umol m-2 s-1)
    real(r8)  :: cosa,convf,cosi,cosii,coss  ! parameters used for COS calculations
    real(r8)  :: cos_assim
    real(r8)  :: kn,kf,kd,kp,ps_sif,fm
    real(r8)  :: je,xxn,fs,ffpa

    real(r8)  :: Aquad,Bquad,Cquad
    real(r8)  :: b_ps, a_ps, e_ps, d_ps
    real(r8)  :: product1,product
    real(r8)  :: ps_1
    real(r8)  :: delta_1
    real(r8)  :: r3q
    real(r8)  :: minroot, maxroot, midroot
    real(r8)  :: tprime25
    real(r8)  :: k_T_opt,k_T,M1,M2,M
    real(r8)  :: VT,RT,KT,kt25,Alpha,Beta,Theta,Qp,c4_a,c4_b,c4_c,c4_d
! 2024/01/07 flag for choosing the form of water stress on stomatal conductance
    real(r8)  :: f_feileaf       ! an empirical scalar of leaf water potential stress on stomatal conductance, dimensionless
    integer   :: fw_flag         !
    real(r8)  :: fws  ! factor of water stress
    !integer   :: quad2,OUTDAT2
    !quad2=500
    !OUTDAT2=600
    minroot=0.0
    maxroot=0.0
    midroot=0.0


    ca=CO2_air ! CO2_air was updated in get_CO2.F90
    cosa = COS_air
    !write(*,*), 'ca = ', ca
!    iphoton = 4.55*0.5*rad_leaf
    iphoton = 4.55*f_leaf*rad_leaf   ! replace 0.5 to a coefficient, f_leaf, leaf respiration rate, for optimization purpose, @MOUSONG WU, 2020-09-14

    if(2*iphoton < 1)then
        iphoton = 0
    end if

    temp_leaf_K = temp_leaf_c + 273.13


    call LAMBDA(temp_leaf_p, fact%latent)
    bound_layer_res%vapor = 1.0/g_lb_w

    met%press_bars  = 1.013
    met%pstat273    = 0.022624 / (273.16 * met%press_bars)

    met%T_Kelvin    = temp_leaf_c+273.13
    met%rhova_g     = e_air * 2165/met%T_Kelvin   ! absolute humidity, g m-3
    met%rhova_kg    = met%rhova_g / 1000.0        ! absolute humidity, kg m-3


    g_lb_c  = 1. / (1.0/g_lb_w*1.6 * temp_leaf_K * (met%pstat273))
    !write(*,*), 'g_lb_c = ', g_lb_c

    m_co2 = m_h2o/1.6
    b_co2 = b_h2o/1.6
    g_lb_cos =  1. / (1.0/g_lb_w*1.56 * temp_leaf_K * (met%pstat273))
    convf = temp_leaf_K * (met%pstat273)

    call SFC_VPD(temp_leaf_K, LH_leaf, rh_leaf)

    tprime25 = temp_leaf_K - tk_25                ! temperature difference

!   /** use Arrhenius Eq. to compute KC and km_o2 */

    call TEMP_FUNC(p_kc25, ekc, tprime25, tk_25, temp_leaf_K, km_co2)

    call TEMP_FUNC(p_ko25, eko, tprime25, tk_25,temp_leaf_K, km_o2)

    call TEMP_FUNC(p_tau25, ektau, tprime25, tk_25,temp_leaf_K, tau)

    bc = km_co2 * (1.0 + o2 / km_o2)

    gammac = 0.5 * o2/tau*1000.0                  ! umol mol-1

    resp_ld25=vc_opt * 0.004657

    if(2.0*iphoton > 10)then                      ! Bin Chen: check this later.
        resp_ld25 = resp_ld25 * 0.4               ! reduce respiration by 40% in light according to Amthor
    end if


    call TEMP_FUNC(resp_ld25, erd, tprime25, tk_25, temp_leaf_K, resp_ld)

!    jmopt   = 2.39*vc_opt - 14.2

! vj_slope default value is 2.39 when using subtractive Eqs
     !jmopt   = vj_slope*vc_opt - 14.2

! 2024/03/17, when water stress works on Vcmax, using the additive Eqs could avoid the negative values of jmax
! 2024/03/13 note: using jmax = 1.64*vcmax +29.1 not the TBOLTZC4 function for jmopt, so jmopt was not used behind
     !jmopt   = 1.64*vc_opt + 29.1

    ! 2024/03/13 using  j = 29.1+1.64*vcmax
    !call TBOLTZ(jmopt, ejm, toptjm, temp_leaf_K, jmax)  ! Apply temperature correction to JMAX

    if (LC == 40 .or. LC == 41) then
       call TBOLTZC4(vc_opt, toptvc, temp_leaf_K, vcmax) ! Apply temperature correction to vcmax
    else
       call TBOLTZ(vc_opt, evc, toptvc, temp_leaf_K, vcmax) ! Apply temperature correction to vcmax
    end if

!***************************************

!        APHOTO = PG - resp_ld, net photosynthesis is the difference
!        between gross photosynthesis and dark respiration. Note
!        photorespiration is already factored into PG.


!        Gs from Ball-Berry is for water vapor.  It must be divided
!        by the ratio of the molecular diffusivities to be valid
!        for A

!****************************************
! choose the form of water stress
    select case (fw_flag)
    case (0)
        fws = 1.0  ! non water stress on Gs
    case (1)
        fws = f_soilwater ! f_soilwater works on BWB slope
    case (2)
        vcmax = vcmax*max(f_feileaf,0.1)  !f_feileaf works on vcmax, and finally influnence the Gs
        fws = 1.0
    case (3)
        fws = f_feileaf ! f_feileaf works on BWB slope
    case (4)
        fws = 1.0 ! f_feileaf only works on Etp
    case (5)
        fws = 1.0
        vcmax = vcmax*max(f_soilwater,0.1) !f_soilwater works on vcmax, and finally influnence the Gs
    case (6)
        fws = 1.0
        vcmax = vcmax*max(min(f_soilwater,f_feileaf),0.1)
    case default
        fws = f_soilwater
    end select
! 2024/03/17, when water stress works on Vcmax, using the additive Eqs could avoid the negative values of jmax
    jmax = 29.1 + 1.64 * vcmax ! 2024/03/13, 1.64 could be set as vj_slope parameter when using additive Eqs




!*******************************************
    if (LC == 40 .or. LC==41) then   ! analytical solution for C4 photosynthesis, added by MOUSONG.WU@201907, refer to Chen et al.2019 AFM

        k_T_opt = 2.*1.e4*vc_opt

        k_T     = k_T_opt * 2.**((temp_leaf_K-298.15)/10.)
        j_photon = iphoton*0.05

        wj = j_photon
        wc = vcmax

        Theta = 0.80
        Beta = 0.95

        M1 = ((wj+wc)+sqrt((wj+wc)*(wj+wc)-4*Theta*wj*wc))/(2*Theta)
        M2 = ((wj+wc)+sqrt((wj+wc)*(wj+wc)+4*Theta*wj*wc))/(2*Theta)

        if (M1<M2) then
           M = M1
        else
           M = M2
        end if

        !alpha_ps = Beta*m_co2*f_soilwater*rh_leaf - Beta*b_co2/g_lb_c - 1.6*k_T/g_lb_c &
               !+ k_T*m_co2*f_soilwater*rh_leaf/g_lb_c - k_T*b_co2/g_lb_c/g_lb_c
        alpha_ps = Beta*m_co2*fws*rh_leaf - Beta*b_co2/g_lb_c - 1.6*k_T/g_lb_c &
               + k_T*m_co2*fws*rh_leaf/g_lb_c - k_T*b_co2/g_lb_c/g_lb_c

        !beta_ps  = Beta*b_co2*ca - Beta*m_co2*f_soilwater*rh_leaf*resp_ld + 1.6*k_T*ca &
               !- M*m_co2*f_soilwater*rh_leaf - k_T*m_co2*f_soilwater*rh_leaf*ca &
               !+ M*b_co2/g_lb_c - resp_ld*k_T*m_co2*f_soilwater*rh_leaf/g_lb_c &
               !+2*k_T*b_co2*ca/g_lb_c + (1.6*k_T*resp_ld+1.6*k_T*M- &
               !M*k_T*m_co2*f_soilwater*rh_leaf)/g_lb_c + M*k_T*b_co2/g_lb_c/g_lb_c

        beta_ps  = Beta*b_co2*ca - Beta*m_co2*fws*rh_leaf*resp_ld + 1.6*k_T*ca &
               - M*m_co2*fws*rh_leaf - k_T*m_co2*fws*rh_leaf*ca &
               + M*b_co2/g_lb_c - resp_ld*k_T*m_co2*fws*rh_leaf/g_lb_c &
               +2*k_T*b_co2*ca/g_lb_c + (1.6*k_T*resp_ld+1.6*k_T*M- &
               M*k_T*m_co2*fws*rh_leaf)/g_lb_c + M*k_T*b_co2/g_lb_c/g_lb_c

        !gamma_ps = resp_ld*M*m_co2*f_soilwater*rh_leaf - M*b_co2*ca + resp_ld*k_T &
               !*m_co2*f_soilwater*rh_leaf*ca - k_T*b_co2*ca*ca + (M*k_T*m_co2*f_soilwater &
               !*rh_leaf-1.6*resp_ld*k_T-1.6*k_T*M)*ca - 2*M*k_T*b_co2*ca/g_lb_c - (1.6*M &
               !*resp_ld*k_T-M*k_T*m_co2*f_soilwater*rh_leaf*resp_ld)/g_lb_c

        gamma_ps = resp_ld*M*m_co2*fws*rh_leaf - M*b_co2*ca + resp_ld*k_T &
               *m_co2*fws*rh_leaf*ca - k_T*b_co2*ca*ca + (M*k_T*m_co2*fws*rh_leaf &
               -1.6*resp_ld*k_T-1.6*k_T*M)*ca - 2*M*k_T*b_co2*ca/g_lb_c - (1.6*M &
               *resp_ld*k_T-M*k_T*m_co2*fws*rh_leaf*resp_ld)/g_lb_c

        !theta_ps = M*k_T*b_co2*ca*ca + 1.6*M*resp_ld*k_T*ca - M*resp_ld*k_T &
               !*m_co2*f_soilwater*rh_leaf*ca

        theta_ps = M*k_T*b_co2*ca*ca + 1.6*M*resp_ld*k_T*ca - M*resp_ld*k_T &
               *m_co2*fws*rh_leaf*ca

        if ((wj <= resp_ld).or.(wc <= resp_ld)) goto 300
        !        //cubic solution:
        !        // A^3 + p A^2 + q A + r = 0

        denom = alpha_ps

        p_cubic = beta_ps / denom

        q_cubic = gamma_ps / denom

        r_cubic = theta_ps / denom


!        // Use solution from Numerical Recipes from Press

        Qroot = (p_cubic*p_cubic - 3.0 * q_cubic) / 9.0
        Rroot = (2.0 * p_cubic*p_cubic*p_cubic  - 9.0 * p_cubic * q_cubic + 27.0 * r_cubic) / 54.0
!       write(*,*) Qroot
        Qroot = max(Qroot,1.e-6)
        r3q = Rroot / sqrt(Qroot*Qroot*Qroot)
        if (r3q>1) r3q=1
        if (r3q<-1) r3q=-1

        ang_L = acos(r3q)

        root1 = -2.0 * sqrt(Qroot) * cos(ang_L / 3.0) - p_cubic / 3.0
        root2 = -2.0 * sqrt(Qroot) * cos((ang_L + PI2) / 3.0) - p_cubic / 3.0
        root3 = -2.0 * sqrt(Qroot) * cos((ang_L - PI2) / 3.0) - p_cubic / 3.0

!       Here A = x - p / 3, allowing the cubic expression to be expressed
!       as: x^3 + ax + b = 0
!       rank roots #1,#2 and #3 according to the minimum, intermediate and maximum value

        if(( root1 <= root2 ).and.( root1 <= root3 ))then
            minroot=root1
            if (root2 <= root3)then
                midroot=root2
                maxroot=root3
            else
                midroot=root3
                maxroot=root2
            end if
        end if


        if(( root2 <= root1 ).and.( root2 <= root3 ))then
            minroot=root2
            if (root1 <= root3)then
                midroot=root1
                maxroot=root3
            else
                midroot=root3
                maxroot=root1
            end if
        end if


        if(( root3 <= root1 ).and.( root3 <= root2 ))then
            minroot=root3
            if (root1 < root2) then
                midroot=root1
                maxroot=root2
            else
                midroot=root2
                maxroot=root1
            end if
        end if

        aphoto=0.0
!     find out where roots plop down relative to the x-y axis

        if (minroot > 0 .and. midroot > 0 .and. maxroot > 0) aphoto=minroot


        if (minroot < 0 .and. midroot < 0 .and. maxroot > 0) aphoto=maxroot


        if (minroot < 0 .and. midroot > 0 .and. maxroot > 0) aphoto=midroot

!        //also test for sucrose limitation of photosynthesis, as suggested by
!        //Collatz.  Js=Vmax/2

        cs = ca - aphoto / g_lb_c

        !gs_h2o_mole = (f_soilwater *m_h2o* rh_leaf * aphoto / cs) + b_h2o   !mol m-2 s-1
        gs_h2o_mole = (fws * m_h2o*rh_leaf * aphoto / cs) + b_h2o   !mol m-2 s-1
        gs_co2_mole = gs_h2o_mole /1.6

        ci = cs - aphoto / gs_co2_mole

        j_sucrose = k_T*ci - resp_ld


        if(j_sucrose < aphoto) aphoto = j_sucrose

!         //Stomatal conductance for water vapor

!        //forest are hypostomatous.
!        //Hence we don't divide the total resistance
!        //by 2 since transfer is going on only one side of a leaf.


!        // if A < 0 then gs should go to cuticular value and recalculate A
!         //using quadratic solution

        if(aphoto <= 0.0)then
            goto 300
        else
            goto 400
        end if

!        // if aphoto < 0  set stomatal conductance to cuticle value


!        //a quadratic solution of A is derived if gs=b, but a cubic form occurs
!        //if gs =ax + b.  Use quadratic case when A <=0


!  // Bin Chen:
!         r_tot = 1.0/b_co2 + 1.0/g_lb_c  // total resistance to CO2 (m2 s mol-1)
!         denom = g_lb_c * b_co2

!         Aquad = r_tot * e_ps
!         Bquad = (e_ps*resp_ld + a_ps)*r_tot - b_ps - e_ps*ca
!         Cquad = a_ps*(ca-d_ps) - resp_ld*(e_ps*ca+b_ps)


300 Aquad = Beta + k_T/g_lb_c + 1.6*k_T/b_co2

    Bquad = 2*Beta*resp_ld - M - ca*k_T - (k_T/g_lb_c + 1.6*k_T/b_co2)*(M-resp_ld)

    Cquad = Beta*resp_ld*resp_ld - M*resp_ld + ca*k_T*(M-resp_ld)

    product1=Bquad * Bquad - 4.0 * Aquad * Cquad

    if (product1 >= 0)then
        aphoto = (-Bquad - sqrt(product1)) / (2.0 * Aquad)
    end if

400 aphoto = max(0., aphoto)

    cs = ca - aphoto / g_lb_c

    !gs_h2o_mole = (f_soilwater *m_h2o* rh_leaf * aphoto / cs) + b_h2o   !mol m-2 s-1
    gs_h2o_mole = (fws * m_h2o*rh_leaf * aphoto / cs) + b_h2o   !mol m-2 s-1
    gs_co2_mole = gs_h2o_mole /1.6

    ci = cs - aphoto / gs_co2_mole

    Gs_w = gs_h2o_mole * temp_leaf_K * (met%pstat273)     !m s-1

    ! 2024/03/31 Constraint the ci >=0.1*ca, cs >=0.8*ca, avoid the negative value in Ci when using water stress on BWB slope @Lu Hu
    ! ci>a*ca, cs>b*ca, a,b should be explored
    ! ci>0.1*ca, (ci/ca)min =[0.10-0.36],ref 'dynamics of changing intercellular co2 concentration during drought and determination of minimum functional ci'
    cs = max(cs,0.8*ca)
    ci = max(ci,0.1*ca)


    else

        !alpha_ps = 1.0 + (b_co2 / g_lb_c) - m_co2*rh_leaf*f_soilwater
        !beta_ps = ca * (g_lb_c*m_co2*rh_leaf*f_soilwater - 2.0 * b_co2 - g_lb_c)
        alpha_ps = 1.0 + (b_co2 / g_lb_c) - m_co2*rh_leaf*fws
        beta_ps = ca * (g_lb_c*m_co2*rh_leaf*fws - 2.0 * b_co2 - g_lb_c)
        gamma_ps = ca * ca * g_lb_c * b_co2
        !theta_ps = g_lb_c*m_co2*rh_leaf*f_soilwater - b_co2
        theta_ps = g_lb_c*m_co2*rh_leaf*fws - b_co2

    !***************************************

!        Test for the minimum of Wc and Wj.  Both have the form:

!        W = (a ci - ad)/(e ci + b)

!        after the minimum is chosen set a, b, e and d for the cubic solution.

!        estimate of J according to Farquhar and von Cammerer (1981)


!***************************************

        j_photon = jmax * iphoton/ (iphoton+ 2.1*jmax)

!    if(LC == 40) then
!        j_photon = jmax
!    else
!      j_photon =jmax * iphoton/ (iphoton+ 2.1*jmax)
!    end if

!  initial guess of intercellular CO2 concentration to estimate Wc and Wj:

        wj = j_photon * (cii - gammac) / (4. * cii + 8.0*gammac)
        !write(*,*), 'cii = ', cii
        !write(*,*), 'wj = ', wj

        wc = vcmax * (cii - gammac) / (cii + bc)
        !write(*,*), 'wc = ', wc

!    if(LC == 40) then
!        wc = vcmax
!    else
!      wc = vcmax * (cii - gammac) / (cii + bc)
!    end if

        if(wj < wc)then
!        // for Harley and Farquhar type model for Wj

            psguess=wj
            a_ps = j_photon
            b_ps = 8.0*gammac
            e_ps = 4.0
            d_ps = gammac

        else

            psguess=wc
            a_ps = vcmax
            b_ps = bc
            e_ps = 1.0
            d_ps = gammac

        end if

!***************************************

! if wj or wc are less than resp_ld then A would probably be less than
! zero.  This would yield a
! negative stomatal conductance.  In this case, assume gs equals the
! cuticular value. This
! assumptions yields a quadratic rather than cubic solution for A

!***************************************

        if ((wj <= resp_ld).or.(wc <= resp_ld)) goto 100

!        //cubic solution:
!        // A^3 + p A^2 + q A + r = 0

        denom = e_ps * alpha_ps

        p_cubic = (e_ps * beta_ps + b_ps * theta_ps - a_ps * alpha_ps + e_ps * resp_ld * alpha_ps)
        p_cubic = p_cubic / denom

        q_cubic = (e_ps * gamma_ps + (b_ps * gamma_ps / ca) - a_ps * beta_ps +   &
               a_ps * d_ps * theta_ps + e_ps * resp_ld * beta_ps + resp_ld * b_ps * theta_ps)
        q_cubic = q_cubic / denom

        r_cubic = -a_ps * gamma_ps + a_ps * d_ps * gamma_ps/ca + e_ps * resp_ld  &
              * gamma_ps + resp_ld * b_ps * gamma_ps/ca
        r_cubic = r_cubic / denom


!        // Use solution from Numerical Recipes from Press

        Qroot = (p_cubic*p_cubic - 3.0 * q_cubic) / 9.0
        Rroot = (2.0 * p_cubic*p_cubic*p_cubic  - 9.0 * p_cubic * q_cubic + 27.0 * r_cubic) / 54.0
!    write(*,*) Qroot
        Qroot = max(Qroot,1.e-6)
        r3q = Rroot / sqrt(Qroot*Qroot*Qroot)
        if (r3q>1) r3q=1
        if (r3q<-1) r3q=-1

        ang_L = acos(r3q)

        root1 = -2.0 * sqrt(Qroot) * cos(ang_L / 3.0) - p_cubic / 3.0
        root2 = -2.0 * sqrt(Qroot) * cos((ang_L + PI2) / 3.0) - p_cubic / 3.0
        root3 = -2.0 * sqrt(Qroot) * cos((ang_L - PI2) / 3.0) - p_cubic / 3.0

!   Here A = x - p / 3, allowing the cubic expression to be expressed
!   as: x^3 + ax + b = 0
!   rank roots #1,#2 and #3 according to the minimum, intermediate and maximum value

        if(( root1 <= root2 ).and.( root1 <= root3 ))then
            minroot=root1
            if (root2 <= root3)then
                midroot=root2
                maxroot=root3
            else
                midroot=root3
                maxroot=root2
            end if
        end if


        if(( root2 <= root1 ).and.( root2 <= root3 ))then
            minroot=root2
            if (root1 <= root3)then
                midroot=root1
                maxroot=root3
            else
                midroot=root3
                maxroot=root1
            end if
        end if


        if(( root3 <= root1 ).and.( root3 <= root2 ))then
            minroot=root3
            if (root1 < root2) then
                midroot=root1
                maxroot=root2
            else
                midroot=root2
                maxroot=root1
            end if
        end if

        aphoto=0.0
!     find out where roots plop down relative to the x-y axis

        if (minroot > 0 .and. midroot > 0 .and. maxroot > 0) aphoto=minroot


        if (minroot < 0 .and. midroot < 0 .and. maxroot > 0) aphoto=maxroot


        if (minroot < 0 .and. midroot > 0 .and. maxroot > 0) aphoto=midroot

!        //also test for sucrose limitation of photosynthesis, as suggested by
!        //Collatz.  Js=Vmax/2

        j_sucrose = vcmax / 2. - resp_ld

!    if(LC == 40) then
!      j_sucrose = 2.*1.e4*vcmax*cii/air_pres/1000. - resp_ld
!    else
!      j_sucrose = vcmax / 2. - resp_ld
!    end if

        if(j_sucrose < aphoto) aphoto = j_sucrose

!         //Stomatal conductance for water vapor

!        //forest are hypostomatous.
!        //Hence we don't divide the total resistance
!        //by 2 since transfer is going on only one side of a leaf.


!        // if A < 0 then gs should go to cuticular value and recalculate A
!         //using quadratic solution

        if(aphoto <= 0.0)then
            goto 100
        else
            goto 200
        end if

!        // if aphoto < 0  set stomatal conductance to cuticle value


!        //a quadratic solution of A is derived if gs=b, but a cubic form occurs
!        //if gs =ax + b.  Use quadratic case when A <=0


!  // Bin Chen:
!         r_tot = 1.0/b_co2 + 1.0/g_lb_c  // total resistance to CO2 (m2 s mol-1)
!         denom = g_lb_c * b_co2

!         Aquad = r_tot * e_ps
!         Bquad = (e_ps*resp_ld + a_ps)*r_tot - b_ps - e_ps*ca
!         Cquad = a_ps*(ca-d_ps) - resp_ld*(e_ps*ca+b_ps)

100 ps_1    = ca * g_lb_c * b_co2
    delta_1 = b_co2 + g_lb_c
    denom   = g_lb_c * b_co2

    Aquad = delta_1 * e_ps
    Bquad = -ps_1 * e_ps - a_ps * delta_1 + e_ps * resp_ld * delta_1- b_ps * denom
    Cquad = a_ps * ps_1 - a_ps * d_ps * denom - e_ps * resp_ld * ps_1 - resp_ld * b_ps * denom

    product1=Bquad * Bquad - 4.0 * Aquad * Cquad

    if (product1 >= 0)then
        aphoto = (-Bquad - sqrt(product1)) / (2.0 * Aquad)
    end if

200 aphoto = max(0., aphoto)

    cs = ca - aphoto / g_lb_c


    !gs_h2o_mole = (f_soilwater *m_h2o* rh_leaf * aphoto / cs) + b_h2o   !mol m-2 s-1
    gs_h2o_mole = (fws * m_h2o*rh_leaf * aphoto / cs) + b_h2o   !mol m-2 s-1
    gs_co2_mole = gs_h2o_mole /1.6

    ci = cs - aphoto / gs_co2_mole

    Gs_w = gs_h2o_mole * temp_leaf_K * (met%pstat273)     !m s-1

    ! 2024/03/31 Constraint the ci >=0.1*ca, cs >=0.8*ca, avoid the negative value in Ci when using water stress on BWB slope @Lu Hu
    ! ci>a*ca, cs>b*ca, a,b should be explored
    ! ci>0.1*ca, (ci/ca)min =[0.10-0.36],ref 'dynamics of changing intercellular co2 concentration during drought and determination of minimum functional ci'
    cs = max(cs,0.8*ca)
    ci = max(ci,0.1*ca)
    ! for US-MOZ site test,
    ! = max(ci,0.5*ca)

    end if


    if (LC == 1001) then   ! C4 model refer to Xiaorong Wang
! 2024/01/07 some codes were repeated compared with those before 'if (LC == 40 .or. LC==41) then '
        ca=CO2_air
!       iphoton = 4.55*0.5*rad_leaf

        iphoton = 4.55*f_leaf*rad_leaf   ! replace 0.5 to a coefficient, f_leaf, leaf respiration rate, for optimization purpose, @MOUSONG WU, 2020-09-14

        if(2*iphoton < 1)then
            iphoton = 0
        end if

        temp_leaf_K = temp_leaf_c + 273.13  !repeated in the front 2024/03/06 @Lu Hu


        call LAMBDA(temp_leaf_p, fact%latent)
        bound_layer_res%vapor = 1.0/g_lb_w

        met%press_bars  = 1.013
        met%pstat273    = 0.022624 / (273.16 * met%press_bars)

        met%T_Kelvin    = temp_leaf_c+273.13
        met%rhova_g     = e_air * 2165/met%T_Kelvin   ! absolute humidity, g m-3
        met%rhova_kg    = met%rhova_g / 1000.0        ! absolute humidity, kg m-3


        g_lb_c  = 1. / (1.0/g_lb_w*1.6 * temp_leaf_K * (met%pstat273))

        m_co2 = m_h2o/1.6
        b_co2 = b_h2o/1.6

        call SFC_VPD(temp_leaf_K, LH_leaf, rh_leaf)

        tprime25 = temp_leaf_K - tk_25                ! temperature difference

        resp_ld25=vc_opt * 0.021

        if(2.0*iphoton > 10) then                      ! Bin Chen: check this later.
            resp_ld25 = resp_ld25 * 0.4               ! reduce respiration by 40% in light according to Amthor
        end if

        kt25=20000*vc_opt    !Collatz 1992

        call VT_FUNC(temp_leaf_K,vc_opt,VT)      !VT=VT_FUNC(temp_leaf_K,vc_opt) //VT   //Collatz 1992
        call RT_FUNC(temp_leaf_K,resp_ld25,RT)   !RT=RT_FUNC(temp_leaf_K,resp_ld25) //RT   //Collatz 1992
        call KT_FUNC(temp_leaf_K,kt25,KT)        !KT=KT_FUNC(temp_leaf_K,kt25) //KT   //Collatz 1992

        !M is the intermediate variable that gives the minimum of wc and wj
        !Theta M^2 - M(wj + wc) + wj wc = 0
        ! a x^2+b x+c=0
        ! x = [-b +/- sqrt(b^2 - 4 a c)]/2a

        Alpha=0.05   ! in relation to the quantum efficiency          //Collatz 1992
        Qp=iphoton   ! absorbed photosynthetically active radiation   //Collatz 1992

        wc = VT         !Collatz 1992
        wj = Alpha*Qp   !Collatz 1992


        Theta=0.83   !Collatz 1992
        Beta=0.93    !Collatz 1992


		M1=((wj+wc)+sqrt((wj+wc)*(wj+wc)-4*Theta*wj*wc))/(2*Theta)
        M2=((wj+wc)-sqrt((wj+wc)*(wj+wc)-4*Theta*wj*wc))/(2*Theta)

        if(M1<M2) then
           M=M1
        else
           M=M2
        end if

        !c4_a = Beta*m_h2o*rh_leaf*f_soilwater -Beta*b_h2o/g_lb_c &
               !- 1.6*KT/g_lb_c + KT*m_h2o*rh_leaf*f_soilwater &
               !/g_lb_c - KT*b_h2o/(g_lb_c*g_lb_c)

        c4_a = Beta*m_h2o*rh_leaf*fws -Beta*b_h2o/g_lb_c &
               - 1.6*KT/g_lb_c + KT*m_h2o*rh_leaf*fws &
               /g_lb_c - KT*b_h2o/(g_lb_c*g_lb_c)

        !c4_b = Beta*b_h2o*ca - Beta*m_h2o*rh_leaf*RT*f_soilwater &
               !+  1.6*KT*ca - M*m_h2o*rh_leaf*f_soilwater &
               !- KT*m_h2o*rh_leaf*ca*f_soilwater - M*b_h2o/g_lb_c &
               !- RT*KT*m_co2*rh_leaf*f_soilwater /g_lb_c &
               !+2*KT*b_h2o*ca /g_lb_c + (1.6*KT*RT+1.6*KT*M-M*KT*m_h2o*rh_leaf*f_soilwater) &
               !/g_lb_c +  M*KT*b_h2o/(g_lb_c*g_lb_c)

        c4_b = Beta*b_h2o*ca - Beta*m_h2o*rh_leaf*RT*fws &
               +  1.6*KT*ca - M*m_h2o*rh_leaf*fws &
               - KT*m_h2o*rh_leaf*ca*fws - M*b_h2o/g_lb_c &
               - RT*KT*m_co2*rh_leaf*fws /g_lb_c &
               +2*KT*b_h2o*ca /g_lb_c + (1.6*KT*RT+1.6*KT*M-M*KT*m_h2o*rh_leaf*fws) &
               /g_lb_c +  M*KT*b_h2o/(g_lb_c*g_lb_c)

        !c4_c = RT*M*m_h2o*rh_leaf*f_soilwater - M*b_h2o*ca + RT*KT*m_h2o*rh_leaf*ca*f_soilwater &
               !- KT*b_h2o*ca*ca+(M*KT*m_h2o*rh_leaf*f_soilwater -1.6*KT*RT-1.6*KT*M)*ca &
               !- 2*M*KT*b_h2o*ca/g_lb_c - (1.6*M*RT*KT - M*KT*m_h2o*rh_leaf*RT*f_soilwater )/g_lb_c

        c4_c = RT*M*m_h2o*rh_leaf*fws - M*b_h2o*ca + RT*KT*m_h2o*rh_leaf*ca*fws &
               - KT*b_h2o*ca*ca+(M*KT*m_h2o*rh_leaf*fws -1.6*KT*RT-1.6*KT*M)*ca &
               - 2*M*KT*b_h2o*ca/g_lb_c - (1.6*M*RT*KT - M*KT*m_h2o*rh_leaf*RT*fws )/g_lb_c

        !c4_d = M*KT*b_h2o*ca*ca + 1.6*M*RT*KT*ca - M*RT*KT*m_h2o*ca*rh_leaf*f_soilwater

        c4_d = M*KT*b_h2o*ca*ca + 1.6*M*RT*KT*ca - M*RT*KT*m_h2o*ca*rh_leaf*fws


        denom = c4_a
        p_cubic = c4_b/c4_a
        q_cubic = c4_c/c4_a
        r_cubic = c4_d/c4_a


		!Use solution from Numerical Recipes from Press
        Qroot = (p_cubic*p_cubic - 3.0 * q_cubic) / 9.0
        Rroot = (2.0 * p_cubic*p_cubic*p_cubic  - 9.0 * p_cubic * q_cubic + 27.0 * r_cubic) / 54.0
       	r3q = Rroot / sqrt(Qroot*Qroot*Qroot)

	    if (r3q>1) then
            r3q=1 	!  by G. Mo
        end if
	    if (r3q<-1) then
            r3q=-1 	!  by G. Mo
        end if

	    ang_L = acos(r3q)

        root1 = -2.0 * sqrt(Qroot) * cos(ang_L / 3.0) - p_cubic / 3.0   ! real roots
        root2 = -2.0 * sqrt(Qroot) * cos((ang_L + PI2) / 3.0) - p_cubic / 3.0
        root3 = -2.0 * sqrt(Qroot) * cos((ang_L - PI2) / 3.0) - p_cubic / 3.0


		! Here A = x - p / 3, allowing the cubic expression to be expressed
        !as: x^3 + ax + b = 0
		! rank roots #1,#2 and #3 according to the minimum, intermediate and maximum value

	   if(root1 <= root2 .and. root1 <= root3) then
		   minroot=root1
		   if (root2 <= root3) then
			   midroot=root2
               maxroot=root3
           else

			   midroot=root3
               maxroot=root2
            end if
        end if


       if(root2 <= root1 .and. root2 <= root3) then
		   minroot=root2
		   if (root1 <= root3) then
			   midroot=root1
               maxroot=root3
           else
			   midroot=root3
               maxroot=root1
            end if
        end if


       if(root3 <= root1 .and. root3 <= root2) then
		   minroot=root3
           if (root1 < root2) then
			   midroot=root1
               maxroot=root2
           else
			   midroot=root2
               maxroot=root1
            end if
        end if

		aphoto=0.0

        if (minroot > 0 .and. midroot > 0 .and. maxroot > 0) aphoto=minroot


        if (minroot < 0 .and. midroot < 0 .and. maxroot > 0) aphoto=maxroot


        if (minroot < 0 .and. midroot > 0 .and. maxroot > 0) aphoto=midroot

        if(aphoto <= 0.0) then
			goto 500
		else
			goto 600
        end if
!a quadratic solution of A is derived if gs=b, but a cubic form occurs
!if gs =ax + b.  Use quadratic case when A <=0
500 Aquad = Beta+KT/g_lb_c+1.6*KT/b_h2o
	Bquad = 2*Beta*RT-M-ca*KT-(KT/g_lb_c+1.6*KT/b_h2o)*(M-RT)
	Cquad = Beta*RT*RT-M*RT+ca*KT*(M-RT)

    product=Bquad * Bquad - 4.0 * Aquad * Cquad

    if (product >= 0) then
	!	*aphoto = (-Bquad + sqrt(product)) / (2.0 * Aquad)
		aphoto = (-Bquad - sqrt(product)) / (2.0 * Aquad)
    end if

    !Tests suggest that APHOTO2 is the correct photosynthetic root when
    !light is zero because root 2, not root 1 yields the dark respiration
    !value rd.

600 aphoto =max(0., aphoto)

    cs = ca - aphoto / g_lb_c

	!gs_h2o_mole = (f_soilwater *m_h2o* rh_leaf *aphoto / cs) + b_h2o  ! mol m-2 s-1
	gs_h2o_mole = (fws * m_h2o*rh_leaf * aphoto / cs) + b_h2o  ! mol m-2 s-1
	gs_co2_mole = gs_h2o_mole /1.6

	ci = cs - aphoto / gs_co2_mole

	Gs_w = gs_h2o_mole * temp_leaf_K * (met%pstat273)  ! m s-1

    ! 2024/03/31 Constraint the ci >=0.1*ca, cs >=0.8*ca, avoid the negative value in Ci when using water stress on BWB slope @Lu Hu
    ! ci>a*ca, cs>b*ca, a,b should be explored in the future
    ! ci>0.1*ca, (ci/ca)min =[0.10-0.36],ref 'dynamics of changing intercellular co2 concentration during drought and determination of minimum functional ci'
    cs = max(cs,0.8*ca)
    ci = max(ci,0.1*ca)

    end if


   !! QBO
!   ffpa = 0.8    ! 0.6~0.9

    if(iphoton <= 0.) then
        xSIF = 0.
    else    ! C4 je is added, @MOUSONG WU, 2020-09-10
        if (LC == 40 .or. LC==41) then
            je = aphoto
        else
            je = max(aphoto*(ci+2.0*gammac)/(ci-gammac),0.0)
        end if

        xxn=1.0-je/(iphoton*ffpa*0.05)
        if(xxn < 0) xxn = 0.

        kf  = 0.05
        kd  = 0.95
        kp  = 4.0

        ps_sif = kp/(kf+kp+kd)*(1.-xxn)
        kn     =  (sif_alpha*xxn + sif_beta)*xxn
        ! kn     = (6.2473*xxn - 0.5994)*xxn
        fm     = kf/(kf+kd+kn)
        fs     = fm*(1.0-ps_sif)

        if(xxn < 0.26) fs = -0.0075*xxn+0.0181

        xSIF   = fs*iphoton*ffpa
!       xSIF = iphoton
    end if

    ! add the module for calculating COS uptake by plants, @MOUSONG.WU 2020-09-17
    ! for COS modeling, whether replacing f_soilwater with f_feileaf should be explored in the future @Lu Hu 2024/01/07
    call cos_plant(LC,cosii,convf,gs_h2o_mole,g_lb_cos,vcmax,ffpa,f_soilwater,cos_assim)

    coss = cosa - cos_assim/g_lb_cos                           ! ppb
    cosi = coss - cos_assim/(gs_h2o_mole/1.56)                 ! ppb

end subroutine


subroutine SFC_VPD (temp_leaf_K, leleafpt, y)
implicit none
!****************************************************

! this function computes the relative humidity at the leaf surface for
! application in the Ball Berry Equation
! latent heat flux, LE, is passed through the function, mol m-2 s-1
! and it solves for the humidity at leaf surface

! es_leaf : saturation vapor pressure at leaf temperature.

!****************************************************
    real(r8)  :: temp_leaf_K, leleafpt
    real(r8)  :: y, rhov_sfc,e_sfc,vpd_sfc,rhum_leaf
    real(r8)  :: es_leaf

    call ES(temp_leaf_K,es_leaf)
    rhov_sfc = (leleafpt / (fact%latent)) * bound_layer_res%vapor + met%rhova_kg  !  kg m-3

    e_sfc = rhov_sfc * temp_leaf_K / 0.2165                                       ! mb
    vpd_sfc = es_leaf - e_sfc                                                     ! mb
    rhum_leaf = 1.0 - vpd_sfc / es_leaf                                           ! 0 to 1.0
    y = rhum_leaf

end subroutine


subroutine TEMP_FUNC(rate, eact, tprime, tref, t_lk, y)
implicit none
    real(r8)  :: rate, eact, tprime, tref, t_lk, y

!   Arhennius temperature function

    y = rate * exp(tprime * eact / (tref * rugc*t_lk))

end subroutine


subroutine LAMBDA (tak, y)
implicit none
!   Latent heat of Vaporiation, J kg-1

    real(r8)  :: tak, y

    y = 3149000.0 - 2370.0 * tak

!   add heat of fusion for melting ice
    if(tak < 273.0)then
        y = y + 333
    end if

end subroutine


subroutine ES(t, y)
implicit none
! saturation vapor pressure function (mb)
! T is temperature in Kelvin

    real(r8)  :: t, y, y1

    if(t > 0)then
        y1 = (54.8781919 - 6790.4985 / t - 5.02808 * log(t))
        y = exp(y1)
    else
        write(*,*)'bad es calc'
    end if

end subroutine


subroutine TBOLTZ(rate, eakin, topt, tl, y)
implicit none
!   Boltzmann temperature distribution for photosynthesis
    real(r8)  :: rate, eakin, topt, tl
    real(r8)  :: y, dtlopt,prodt,numm,denom

    dtlopt = tl   - topt
    prodt  = rugc * topt * tl
    numm   = hkin * exp(eakin * (dtlopt) / (prodt))
    denom  = hkin - eakin * (1.0 - exp(hkin * (dtlopt) / (prodt)))
    y      = rate * numm / denom
end subroutine

subroutine TBOLTZC4(rate, topt, tl, y)
implicit none
! temperature correction for C4 plants
    real(r8)  :: rate, topt, tl
    real(r8)  :: y, dtlopt,q
    real(r8)  :: s1,s2,s3,s4,fh,fl

    dtlopt = tl - topt
    q = 2.
    s1 = 0.3
    s2 = 313.15
    s3 = 0.2
    s4 = 288.15

    fh = 1 + exp(s1*(tl-s2))
    fl = 1 + exp(s3*(s4-tl))
    y = rate*q**dtlopt/fh/fl
end subroutine

subroutine fluorescence(x,fs)
     implicit none
     real, intent(in)    :: x
     real, intent(out)   :: fs
     real:: Kn
     real:: Kf
     real:: Kd
     real:: Kp
     real:: po0
     real:: ps
     real:: fo0
     real:: fo
     real:: fm
     real:: fm0
     real:: eta
     real:: qQ
     real:: qE

     Kf          = 0.05
     Kd          = 0.95
     Kp          = 4.0

     po0         = Kp/(Kf+Kd+Kp)
     ps          = po0*(1.0-x)
     Kn          = (6.2473 * x - 0.5944)*x


     fo0         = Kf/(Kf+Kp+Kd)
     fo          = Kf/(Kf+Kp+Kd+Kn)
     fm          = Kf/(Kf   +Kd+Kn)
     fm0         = Kf/(Kf   +Kd)
     fs          = fm*(1.0-ps)
     eta         = fs/fo0


     qQ          = 1.0-(fs-fo)/(fm-fo)
     qE          = 1.0-(fm-fo)/(fm0-fo0)

end subroutine

!*----------------------------------------------------------------*!
!xiaorong 21/07/2019
!VT_FUNC :The temperature dependence of the substrate saturated rubisco capacity(vmax)
!//Q10:(umol/m2/s)   Vcmax25:vcmax at 25 Degrees celsiusï¼Œ  Tl:leave temperature(K)
!/*----------------------------------------------------------------*/
subroutine VT_FUNC(Tl,Vcmax25,y)
   implicit none
   real(r8)  :: Tl, Vcmax25, y,Q10
   Q10=2.4
   y=(Vcmax25*Q10**((Tl-298.13)/10))/( (1+exp (0.3*(286.13-Tl)))*(1+exp(0.3*(Tl-309.13))))
end subroutine


!/*----------------------------------------------------------------*/
!//xiaorong 21/07/2019
!//RT_FUNC :the pseudo-first order rate constant with respect to co2 (k) and the leaf respiration (RT) are goven as Q10 function
!//Q10:(umol/m2/s)   Rd25:respiration at 25 Degrees celsiusï¼Œ  Tl:leave temperature(K)
!/*----------------------------------------------------------------*/
subroutine RT_FUNC(Tl,Rd25,y)
   implicit none
   real(r8)  :: Tl,Rd25, Q10,y
   Q10=2
   y=(Rd25*Q10**((Tl-298.13)/10))/(1+exp(1.3*(Tl-328.13)))
end subroutine

!/*----------------------------------------------------------------*/
!//xiaorong 21/07/2019
!//Q10:(mol/m2/s)  Tl:leave temperature(K)
!/*----------------------------------------------------------------*/
subroutine KT_FUNC(Tl,KT25,y)
   implicit none
   real(r8)  :: Tl,KT25, Q10,y
   Q10=2.4  !//mol/m2/s-->umol/m2/s
   y=(KT25*Q10**((Tl-298.13)/10))
end subroutine

end module

