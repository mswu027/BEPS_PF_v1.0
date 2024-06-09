module ecoRespMod
use shr_kind_mod, only: r8=>shr_kind_r8
use beps_par
use mid_results
implicit none
real(r8),parameter:: sec_per_day = 86400.

contains

subroutine  plant_resp(f_q10,lc,mid_res,lai_yr,lai,temp_air,temp_soil,CosZs)
implicit none
integer,intent(in)   :: lc
type(results),intent(inout)  :: mid_res
real(r8),intent(in)  :: f_q10,lai_yr
real(r8),intent(in)  :: lai,temp_air,temp_soil,CosZs

real(r8) :: temp_opt25  = 25.0
real(r8) :: biomass,biomass_leaf_o,biomass_stem_o,biomass_root_o,biomass_leaf_u,biomass_stem_u,biomass_root_u
real(r8) :: respir_croot_o,respir_root_o,respir_stem_o,respir_leaf_o
real(r8) :: respir_croot_u,respir_root_u,respir_stem_u,respir_leaf_u
real(r8) :: q10
real(r8) :: exponent1
real(r8) :: lai_u,lai_max_o,lai_max_u
real(r8) :: ra
real(r8) ::coef_leaf_respir,coef_stem_respir,coef_root_respir,coef_fineroot_respir
real(r8) :: gpp_o,gpp_u,gpp_r,rg,ratio_froot

if(lc == 25 .or. lc ==40 .or. lc==41 .or. lc==1001) then
  lai_u  = 0.01
else
  lai_u  = 1.18*exp(-0.99*lai)
end if

if(lai_u > lai) lai_u = 0.01

if(lc ==6) then
   ra  = 0.6
else
   ra  = 1.0
end if

q10  = 3.22 - f_q10*temp_air      ! f_q10 default value 0.046

if(lc >=1 .and. lc <= 5) then
           !/*  calculating aboveground biomass based on LAI  J. Liu 2002 */
        biomass=0.9097*lai_yr+0.125*lai_yr*lai_yr
        biomass_leaf_o=0.05*biomass    !/* leaf C of overstory */
        biomass_stem_o=0.95*biomass    !/* stem C of overstory */
        biomass_root_o=0.454*biomass
        !/*biomass_root_o=0.232*biomass; // root C of overstoryKurz 1996 */
        biomass_leaf_u=0.3*biomass_leaf_o  !/* leaf C of understory */
        biomass_stem_u=0.02*biomass_stem_o     !/* stem C of understory */
        biomass_root_u=0.05*biomass_root_o !/* root C of understory */

        coef_leaf_respir=0.0015/sec_per_day  !/*  leaf resp co.  kg C-1 d-1 kg-1    */
        coef_stem_respir=0.0020/sec_per_day  !/*  stem resp co.   kg C-1 d-1 kg-1      */
        coef_root_respir=0.0020/sec_per_day  !/*  root resp co.   kg C-1 d-1 kg-1   */
        coef_fineroot_respir=0.003/sec_per_day   !/*  fine root resp co.   kg C-1 d-1 kg-1  */

        lai_max_o=4.5     ! /* LAI_max_overstory                         3.3*/
        lai_max_u=2.4     ! /* LAI_max_understroy                        2.4*/
else if(lc ==6 .or. lc ==9) then
        !!/*  calculating aboveground biomass based on LAI  J. Liu 2002 */
        biomass=1.545*lai_yr+0.183*lai_yr*lai_yr
        biomass_leaf_o=0.04*biomass    !/* leaf C of overstory */
        biomass_stem_o=0.96*biomass    !/* stem C of overstory */
        biomass_root_o=1.432*biomass**0.639    !/* root C of overstory  Kurz 1996 */
        biomass_leaf_u=0.3*biomass_leaf_o  !/* leaf C of understory */
        biomass_stem_u=0.01*biomass_stem_o     !/* stem C of understory */
        biomass_root_u=0.01*biomass_root_o  !/* root C of understory */

        coef_leaf_respir=0.015/sec_per_day  !/*  leaf resp co.  kg C-1 d-1 kg-1    */
        coef_stem_respir=0.0035/sec_per_day  !/*  stem resp co.   kg C-1 d-1 kg-1      */
        coef_root_respir=0.0025/sec_per_day  !/*  root resp co.   kg C-1 d-1 kg-1   */
        coef_fineroot_respir=0.003/sec_per_day    !/*  fine root resp co.   kg C-1 d-1 kg-1  */

        lai_max_o=4.5      !/* LAI_max_overstory                         3.3*/
        lai_max_u=2.4      !/* LAI_max_understroy                        2.4*/
else if (lc == 10) then
        biomass = 1.227*lai_yr+0.154*lai_yr*lai_yr
        biomass_leaf_o  = 0.045*biomass
        biomass_stem_o  = 0.95*biomass
        biomass_root_o  = (0.454*biomass+1.432*biomass**0.639)/2.
        biomass_leaf_u  = 0.3*biomass_leaf_o
        biomass_stem_u  = 0.015*biomass_stem_o
        biomass_root_u  = 0.03*biomass_root_o

        coef_leaf_respir = 0.008/sec_per_day
        coef_stem_respir = 0.0028/sec_per_day
        coef_root_respir = 0.0023/sec_per_day
        coef_fineroot_respir = 0.003/sec_per_day

        lai_max_o  = 4.5
        lai_max_u  = 2.4

else if (lc ==13) then
        biomass=1.545*lai_yr+0.183*lai_yr*lai_yr
        biomass_leaf_o=0.1*biomass    !/* leaf C of overstory */
        biomass_stem_o=0.90*biomass    !/* stem C of overstory */
        biomass_root_o=1.432*biomass**0.639    !/* root C of overstory  Kurz 1996 */
        biomass_leaf_u=0.3*biomass_leaf_o     !/* leaf C of understory */
        biomass_stem_u=0.01*biomass_stem_o    ! /* stem C of understory */
        biomass_root_u=0.01*biomass_root_o    !/* root C of understory */

        coef_leaf_respir=0.001/sec_per_day  !/*  leaf resp co.  kg C-1 d-1 kg-1    */
        coef_stem_respir=0.002/sec_per_day  !/*  stem resp co.   kg C-1 d-1 kg-1      */
        coef_root_respir=0.0015/sec_per_day  !/*  root resp co.   kg C-1 d-1 kg-1   */
        coef_fineroot_respir=0.003/sec_per_day    !/*  fine root resp co.   kg C-1 d-1 kg-1  */

        lai_max_o=3.3    ! /*   LAI_max_overstory   */
        lai_max_u=0.01   !/*  LAI_max_understroy   */
else if(lc == 14 .or. lc == 15 .or. lc ==25 .or. lc ==40 .or. lc==41) then
        biomass_leaf_o=0.05*lai_yr  ! /* leaf C = lai/20  from W.Ju 05y11*/
        biomass_stem_o=0.0          !/* stem C */
        biomass_root_o=0.061*lai_yr    !/* root C = lai/20*0.55/0.45  from W.Ju 05y11*/
        biomass_leaf_u=0.0
        biomass_stem_u=0.0
        biomass_root_u=0.0

        coef_leaf_respir=0.001/sec_per_day
        coef_stem_respir=0.002/sec_per_day  !/*  stem resp co.   kg C-1 d-1 kg-1      */
        coef_root_respir=0.0015/sec_per_day !/* root resp co.   kg C-1 d-1 kg-1   */
        coef_fineroot_respir=0.003/sec_per_day  !/*  fine root resp co.   kg C-1 d-1 kg-1  */

        lai_max_o=3.3  !  /*   LAI_max_overstory     */
        lai_max_u=0.01 !/*  LAI_max_understroy     */  !/* leaf resp co.  kg C-1 d-1 kg-1    */
else if(lc ==1001 .or. lc ==1002 .or. lc ==1003 .or. lc ==1004) then
        biomass_leaf_o=0.05*lai_yr  ! /* leaf C = lai/20  from W.Ju 05y11*/
        biomass_stem_o=0.0          !/* stem C */
        biomass_root_o=0.061*lai_yr    !/* root C = lai/20*0.55/0.45  from W.Ju 05y11*/
        biomass_leaf_u=0.0
        biomass_stem_u=0.0
        biomass_root_u=0.0

        coef_leaf_respir=0.001/sec_per_day
        coef_stem_respir=0.002/sec_per_day  !/*  stem resp co.   kg C-1 d-1 kg-1      */
        coef_root_respir=0.0015/sec_per_day !/* root resp co.   kg C-1 d-1 kg-1   */
        coef_fineroot_respir=0.003/sec_per_day  !/*  fine root resp co.   kg C-1 d-1 kg-1  */

        lai_max_o=3.3  !  /*   LAI_max_overstory     */
        lai_max_u=0.01 !/*  LAI_max_understroy     */  !/* leaf resp co.  kg C-1 d-1 kg-1    */
end if

!! calculation for overstorey
! stem maintenance respiration
exponent1 = (temp_air-temp_opt25)/10.0
respir_stem_o = (biomass_stem_o*0.35/(biomass_stem_o+0.35))*coef_stem_respir*q10**exponent1*ra
respir_stem_o = max(respir_stem_o, 0.0)

! root maintenance
exponent1=(temp_soil-temp_opt25)/10.0
if(lc ==14 .or. lc == 15 .or. lc ==25 .or. lc ==40 .or. lc==41 .or. lc ==1001 .or. lc ==1002 .or. lc ==1003 .or. lc ==1004) then
    respir_root_o = biomass_root_o*coef_root_respir*q10**exponent1*ra
else
    ratio_froot=exp(1.007)*biomass_root_o**(-0.841)
    ratio_froot=min(0.9, ratio_froot)

    respir_croot_o=0.05*biomass_root_o*(1-ratio_froot)*coef_root_respir*q10**exponent1      !/*coarse root */
    respir_root_o=respir_croot_o+0.05*biomass_root_o*ratio_froot*coef_fineroot_respir*q10**exponent1   ! /* coarse + fine root */
end if
 respir_root_o = max  (respir_root_o, 0.0)

!/* leaf day/night time maintenance respiration */
if (CosZs>0.01) then
    respir_leaf_o=0
else
    exponent1=(temp_air-temp_opt25)/10.0
    respir_leaf_o =lai/lai_max_o*biomass_leaf_o*coef_leaf_respir*q10**exponent1*ra   !kgC/m2/s
end if
respir_leaf_o =max( respir_leaf_o, 0.0)

!   /*     changed in Nov.2005 by W. Ju */
gpp_o = (mid_res%gpp_o_sunlit + mid_res%gpp_o_shaded)    !kgC/m2/s
gpp_r = gpp_o - (respir_leaf_o+respir_stem_o+respir_root_o)
if(gpp_r <=0) then
   rg  = 0.
else
   rg  = 0.35*gpp_r
end if

!mid_res%npp_o  = gpp_r - rg    !kgC/m2/s
mid_res%npp_o   = gpp_o*0.45    !kgC/m2/s @J.Wang

!! calculation for understorey

! /* stem maintenance respiration */
exponent1=(temp_air-temp_opt25)/10.0
respir_stem_u =(biomass_stem_u*0.35/(biomass_stem_u+0.35))*coef_stem_respir*q10**exponent1*ra
respir_stem_u = max(respir_stem_u, 0.0)

!/* root maintenance respiration      changed in Nov.2005 by W. Ju */
exponent1=(temp_soil-temp_opt25)/10.0
if(lc == 14 .or. lc == 15 .or. lc ==25 .or. lc ==40 .or. lc==41 .or. lc ==1001 .or. lc ==1002 .or. lc ==1003 .or. lc ==1004) then
  respir_root_u = biomass_root_u*coef_root_respir*q10**exponent1*ra
else
  ratio_froot=exp(1.007)*biomass_root_u**(-(0.841))
  ratio_froot=min(0.9, ratio_froot)

  respir_croot_u=0.05*biomass_root_u*(1-ratio_froot)*coef_root_respir*q10**exponent1
  respir_root_u=respir_croot_u+0.05*biomass_root_u*ratio_froot*coef_fineroot_respir*q10**exponent1
end if
respir_root_u = max(respir_root_u, 0.0)

if (CosZs>0.01) then
   respir_leaf_u=0
else
    exponent1=(temp_air-temp_opt25)/10.0
    respir_leaf_u =lai_u/lai_max_u*biomass_leaf_u*coef_leaf_respir*q10**exponent1*ra*0.5
end if
 respir_leaf_u =max(respir_leaf_u, 0.0)

!!/*     changed in Nov.2005 by W. Ju */
gpp_u = (mid_res%gpp_u_sunlit + mid_res%gpp_u_shaded)
gpp_r = gpp_u - (respir_leaf_u+respir_stem_u+respir_root_u)

if(gpp_r <=0) then
    rg = 0
else
    rg = 0.35*gpp_r
end if

!mid_res%npp_u = gpp_r - rg  !kgC/m2/s
mid_res%npp_u  = gpp_u*0.45  !kgC/m2/s

mid_res%NPP   = mid_res%npp_u + mid_res%npp_o
end subroutine


subroutine  soil_resp(f_slow,Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,Cp,npp_yr,coef,soiltype,soilp,mid_res)
use beps_soilMod
use beps_par
implicit none
real(r8),intent(inout) :: Ccd(0:4),Cssd(0:4),Csmd(0:4),Cfsd(0:4),Cfmd(0:4),Csm(0:4),Cm(0:4),Cs(0:4),Cp(0:4)
real(r8),intent(in)    :: npp_yr,f_slow
real(r8),intent(in)    :: coef(0:49)
integer,intent(in)     :: soiltype
type(soil),intent(in)  :: soilp
type(results),intent(inout) :: mid_res

real(r8) :: fw, fcr, fl, ffr, kw_cd, kcr_cd, kl_sl, kfr_fl, km_p, ks_p
real(r8) :: kssd_a, kssd_sm, kssd_s, ksmd_a, ksmd_sm,kfsd_a, kfsd_m, kfsd_s, kfmd_a, kfmd_m
real(r8) :: kcd_a, kcd_m
real(r8) :: kcd_s,ksm_a,ksm_s, km_a, km_s, ks_a, ks_m,kp_a, kp_m
real(r8) :: Cw(0:9),Ccr(0:9),Cl(0:9),Cfr(0:9),dCw(0:9),dCcr(0:9),dCl(0:9),DCfr(0:9)
real(r8) :: dCcd(0:9),dCssd(0:9),dCsmd(0:9),dCfsd(0:9),dCfmd(0:9),dCsm(0:9),dCm(0:9),dCs(0:9),dCp(0:9)
real(r8) :: part1,part2
real(r8) :: Fm(0:9),npp
real(r8) :: lambda(0:layer),lambda_t(0:layer),lambda_w(0:layer)
real(r8) :: lam_u,lam_d
integer  :: ii

do ii= 1,layer
   if (308.56*(1/(35.0+46.032)-1/(46.032+soilp%temp_soil_c(ii - 1))) < -2.3) then
      lambda_t(ii) = 0.1                     !! to get rid of enormous value @MOUSONG.WU
   else if (308.56*(1/(35.0+46.032)-1/(46.032+soilp%temp_soil_c(ii - 1))) > 0.) then
      lambda_t(ii) = 1.
   else
      lambda_t(ii) = exp(308.56*(1/(35.0+46.032)-1/(46.032+soilp%temp_soil_c(ii - 1))))  ! Arrenius Equation
   end if
   lambda_t(ii) = min(1.0,lambda_t(ii))
   lambda_t(ii) = max(0.3,lambda_t(ii))
end do

do ii=1,layer
  if(soiltype >=6) then
    lambda_w(ii) = 5.44*soilp%thetam(ii-1)/soilp%fei(ii-1)-5.03*(soilp%thetam(ii-1)/soilp%fei(ii-1))**2-0.472
  else
    lambda_w(ii) = 5.63*soilp%thetam(ii-1)/soilp%fei(ii-1)-4.64*(soilp%thetam(ii-1)/soilp%fei(ii-1))**2-0.710
  end if
  lambda_w(ii)=max(0.3,lambda_w(ii))
end do

do ii=1,layer
    lambda(ii)=lambda_t(ii)*lambda_w(ii)
end do


lam_u  = lambda(1)  ! for surface pool
lam_d  = lambda(2)  ! for soil pool

fw     = coef(0)
fcr    = coef(1)
fl     = coef(2)
ffr    = coef(3)
kw_cd  = coef(4)/8760      ! units?? @J.Wang
kcr_cd = coef(5)/8760
kl_sl  = coef(6)/8760
kfr_fl = coef(7)/8760
kssd_a = coef(8)/8760
kssd_sm= coef(9)/8760
kssd_s = coef(10)/8760
ksmd_a = coef(11)/8760
ksmd_sm= coef(12)/8760
kfsd_a = coef(13)/8760
kfsd_m = coef(14)/8760
kfsd_s = coef(15)/8760
kfmd_a = coef(16)/8760
kfmd_m = coef(17)/8760
kcd_a  = coef(18)/8760
kcd_m  = coef(19)/8760
kcd_s  = coef(20)/8760
km_a   = coef(21)/8760
km_p   = coef(22)/8760
km_s   = coef(23)/8760
ksm_a  = coef(24)/8760
ksm_s  = coef(25)/8760
ks_a   = coef(26)/8760
ks_p   = coef(27)/8760
ks_m   = coef(28)/8760
kp_a   = coef(29)/8760
kp_m   = coef(30)/8760

Cw(0)  = coef(0)/coef(4)*npp_yr  !for stem gC.m2
Ccr(0) = coef(1)/coef(5)*npp_yr  ! for coast root
Cl(0)  = coef(2)/coef(6)*npp_yr  ! for leaf
Cfr(0) = coef(3)/coef(7)*npp_yr  ! for fine root

Fm(1)  = 0.2
npp    = mid_res%npp_o + mid_res%npp_u   !kg/m2/s
npp    = npp*1000*step         !gC/m2/step return to original units @J.Wang

dCw(1) = fw*npp  - kw_cd*Cw(0)
dCcr(1)= fcr*npp - kcr_cd*Ccr(0)
dCl(1) = fl*npp  - kl_sl*Cl(0)
dCfr(1)= ffr*npp - kfr_fl*Cfr(0)

Cw(1)  = Cw(0) + dCw(1)
Ccr(1) = Ccr(0)+ dCcr(1)
Cl(1)  = Cl(0) + dCl(1)
Cfr(1) = Cfr(0)+ dCfr(1)

part1  = (kw_cd*Cw(1) +kcr_cd * Ccr(1))/(1+lam_d*(kcd_a + kcd_m + kcd_s))
part2  = Ccd(0)*lam_d*(Kcd_a+ kcd_m + kcd_s)
dCcd(1)= part1 - part2
Ccd(1) = Ccd(0)+dCcd(1)
!Coarse detrius from woody and coarse root

part1  = (1-Fm(1))*kl_sl*Cl(1)/(1+lam_u*(kssd_a+kssd_sm + kssd_s))
part2  = Cssd(0)* lam_u * (kssd_a + kssd_sm + kssd_s)
dCssd(1)  = part1 - part2
Cssd(1)   = Cssd(0)+dCssd(1)
!surface structural litter

part1  = Fm(1)*kl_sl*Cl(1)/(1+lam_u*(ksmd_a+ksmd_sm))
part2  = Csmd(0)*lam_u*(ksmd_a+ksmd_sm)
dCsmd(1) = part1 - part2
Csmd(1)  = Csmd(0) + dCsmd(1)
!surface metobolic litter

part1  = (1-Fm(1))*kfr_fl*Cfr(1)/(1+lam_d*(kfsd_a + kfsd_m + kfsd_s))
part2  = Cfsd(0)*lam_d*(kfsd_a + kfsd_m+kfsd_s)
dCfsd(1) = part1 - part2
Cfsd(1)  = Cfsd(0) + dCfsd(1)
!for soil strutural litter pool

part1  = Fm(1)*kfr_fl*Cfr(1)/(1+lam_d*(kfmd_a + kfmd_m))
part2  = lam_d*(kfmd_a + kfmd_m)*Cfmd(0)
dCfmd(1) = part1 - part2
Cfmd(1)  = Cfmd(0) + dCfmd(1)
! soil metobolic pool

part1  = lam_u*(Cssd(1)*kssd_sm+Csmd(1)*ksmd_sm)
part2  = lam_u*Csm(0)*(ksm_a + ksm_s)
dCsm(1)  = part1 - part2
Csm(1)   = Csm(0) + dCsm(1)
! surface microbe pool

part1  = (lam_d*(kfsd_m*Cfsd(1)+kfmd_m*Cfmd(1)+Ccd(1)*kcd_m)+lam_d*(Cs(0)*ks_m+Cp(0)*kp_m))
part2  = Cm(0)*lam_d*(km_a+km_s+km_p)
dCm(1) = part1 - part2
Cm(1)  = Cm(0)+dCm(1)
!soil microbe pool

part1  = (lam_d*(Cm(1)*km_s+Ccd(1)*kcd_s+Cfsd(1)*kfsd_s)+lam_u*(Csm(1)*ksm_s+Cssd(1)*kssd_s))
part2  = Cs(0)*lam_d*(ks_a+ks_p+ks_m)
dCs(1) = f_slow*(part1 - part2)
Cs(1)  =Cs(0)+dCs(1)
!for slow carbon pool

dCp(1) = (lam_d*(km_p*Cm(1)+Ks_p*Cs(1))-lam_d*(kp_m * Cp(0) + kp_a * Cp(0)))
Cp(1)  = Cp(0)+dCp(1)
! passive carbon pool

!NEP
mid_res%NEP  = npp+(dCsmd(1)+dCssd(1)+dCfsd(1)+dCfmd(1)+dCcd(1)+dCm(1)+dCsm(1)+dCs(1)+dCp(1))
mid_res%NEP  = mid_res%NEP*1e-3/step    !kgC/m2/s
mid_res%NPP  = npp*1e-3/step
mid_res%soil_resp   = mid_res%NPP - mid_res%NEP
return

end subroutine
end module

