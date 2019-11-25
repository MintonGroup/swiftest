!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_viscosity
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : solves
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_viscosity(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_viscosity(Gsigma, Gm_pdisk, v2_pdisk, r_pdisk, r_hstar, Q, tau, w) result(nu)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_viscosity
      implicit none

! Arguments
      real(DP),intent(in) :: Gsigma, Gm_pdisk, v2_pdisk, r_pdisk, r_hstar, Q, tau, w
      real(DP) :: nu

! Internals
      real(DP)       :: kappa_Q,eta_Q,y
      real(DP)       :: nu_trans_stable,nu_grav_stable,nu_trans_unstable,nu_grav_unstable
      real(DP)       :: nu_trans,nu_grav,nu_coll

! Executable code


   nu_trans_unstable = 13 * r_hstar**5 * Gsigma**2 / w**3
   nu_trans_stable = (v2_pdisk / (2 * w)) * (0.46_DP * tau / (1._DP + tau**2))

   nu_grav_stable = 0.0_DP
   nu_grav_unstable = nu_trans_unstable

   y = 0.25_DP * Q 
   kappa_Q = ringmoons_transition_function(y)
   eta_Q = 1._DP - kappa_Q

   nu_trans = kappa_Q * nu_trans_stable + eta_Q * nu_trans_unstable
   nu_grav  = kappa_Q * nu_grav_stable  + eta_Q * nu_grav_unstable
   nu_coll = r_pdisk**2 * w * tau

   nu = nu_trans + nu_grav + nu_coll


end function ringmoons_viscosity
