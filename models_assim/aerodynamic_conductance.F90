!! calculate aerodynamic resistance/conductance
!! Written by: J. Liu, and W. Ju, Modified by G. MO
!! Last update : May 2015

subroutine aerodynamic_conductance(canopy_height_o,canopy_height_u,zz,clumping,temp_air,wind_sp,SH_o_p,lai_o,&
                                   lai_u,rm,ra_u,ra_g,G_o_a,G_o_b,G_u_a,G_u_b)
use shr_kind_mod,only: r8=>shr_kind_r8
implicit none
real(r8),intent(in)  :: canopy_height_o,canopy_height_u,zz,clumping,temp_air,wind_sp
real(r8),intent(in)  :: SH_o_p,lai_o,lai_u
real(r8),intent(out) :: rm,ra_u,ra_g,G_o_a,G_o_b,G_u_a,G_u_b

real(r8)  :: kh_o
real(r8)  :: lw = 0.3   !leaf characteristic width  = 0.3 for BS
real(r8)  :: sigma = 5  !shelter factor = 5 for BS
real(r8)  :: rb_o,rb_u
real(r8)  :: k = 0.4    ! von Karman's constant
real(r8)  :: beta = 0.5 ! Bowen's ratio
real(r8)  :: cp   = 1010!specific heat of air  J/kg/K
real(r8)  :: density_air  = 1.225
real(r8)  :: gg   = 9.8
real(r8)  :: n    = 5.0
real(r8)  :: nu_lower   ! viscosity (cm2/s)
real(r8)  :: uf,psi
real(r8)  :: d    !displacement height (m)
real(r8)  :: z0     !roughness length    (m)
real(r8)  :: ustar  !friction velocity (m/s)
real(r8)  :: L,Le
real(r8)  :: uh     !wind speed at height h
real(r8)  :: ud     !wind speed at height d
real(r8)  :: u_50   ! wind speed at height 50 m, @mousong.wu, to correct the calculation of ustar
                    ! from 2 m wind speed but the canopy height is around 30 m,which leads to negative value in log function.
real(r8)  :: gamma1,Re,Nu,alfac,alfaw,ram,un_d,un_t,kh_u,z_50

nu_lower    = (13.3+temp_air*0.07)/1000000.
alfac       = 0.15    ! for CO2
alfaw       = (18.9+temp_air*0.07)/1000000.
z_50         = 50.     ! for calculation of windspeed at 50 m height

if(wind_sp ==0.) then
   uh   = 0.
   uf   = 0.
   psi  = 6.
   G_o_a   = 1./200.0
   G_o_b   = 1./200.0
   G_u_a   = 1./200.0
   G_u_b   = 1./200.0
   ra_g    = 300.
else
   d       = 0.8*canopy_height_o
   z0      = 0.08*canopy_height_o
   u_50    = wind_sp*log((z_50-zz)/z0)
   ustar   = u_50*k/log((z_50-d)/z0)
   L       = -(k*gg*SH_o_p)/(density_air*cp*(temp_air+273.3)*ustar**3)
   L       = max(-2.,L)

   ram     = 1./(k*ustar)*(log((z_50-d)/z0)+(n*(z_50-d)*L))
   ram     = max(2.,ram)
   ram     = min(100.,ram)

   if(L>0.) then
      psi = 1.+5.*(z_50-d)*L
   else
      psi = (1.-16.*(z_50-d)*L)**(-0.5)
   end if
   psi = min(10.0,psi)

   !! Leaf boundary layer resistance
   !! wind speed at tree top
    uh    = 1.1*ustar/k
    Le    = lai_o*clumping
    gamma1 = (0.167+0.179*uh)*Le**(1./3.)

   !! wind speed at d ,taking as the mean wind speed inside a stand
   ud     = uh*exp(-gamma1*(1.-d/canopy_height_o))

   !! Reynold's number
   Re     = (ud*0.1)/nu_lower

   !!Nusselt number
   Nu     = 1.0*Re**0.5

   !!leaf boudnary resistance
   rb_o   = min(40.,0.5*0.1/(alfaw*Nu))

   uf     = ustar
   rm     = ram
   G_o_a  = 1./ram
   G_o_b  = 1./rb_o

   kh_o   = 0.41*ustar*(canopy_height_o-canopy_height_o*0.8)/psi
   gamma1  = 0.1+lai_o**0.75

   !!wind speed at the zero displancement of canopy
   un_d=uh*exp(-gamma1*(1.-canopy_height_u*0.8/canopy_height_o))
   un_t=uh*exp(-gamma1*(1.-canopy_height_u/canopy_height_o))

   Re  = (un_d*0.1)/nu_lower
   Nu  = 1.0*Re**0.5

   rb_u  = 0.5*0.1/(alfaw*Nu)
   rb_u  = min(40.,rb_u)
   G_u_b = 1./rb_u
   ra_u  = canopy_height_o/(gamma1* kh_o)*(exp(gamma1*(1.-canopy_height_u/canopy_height_o))-1.)
   G_u_a = 1./(ram+ra_u)

   gamma1 = 4.0
   kh_u=kh_o*exp(-gamma1*(1.-canopy_height_u/canopy_height_o))
   ! 2024/05/14 revisal
   ra_g=canopy_height_o/(gamma1* kh_o)*(exp(gamma1*(1.-0.0/canopy_height_o))-exp(gamma1*(1.-canopy_height_u/canopy_height_o)))
   ra_g = ra_g+ra_u+ram
   ra_g = max(120.,ra_g)

end if

   return
end subroutine
