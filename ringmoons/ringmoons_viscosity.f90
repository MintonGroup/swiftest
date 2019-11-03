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
SUBROUTINE ringmoons_viscosity(ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_viscosity
      IMPLICIT NONE

! Arguments
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      real(DP) :: Q,tau,sigma_r,kappa,eta,Gsigma
      integer(I4B) :: i
      real(DP) :: nu_trans_stable,nu_grav_stable,nu_trans_unstable,nu_grav_unstable,y,nu_trans,nu_grav,nu_coll

! Executable code
 
   do concurrent (i = 1:ring%N) 

      if (ring%Gsigma(i) <= sqrt(TINY)) then
         ring%nu(i) = 0.0_DP
         cycle
      end if

      kappa = ringmoons_transition_function(ring%r_hstar(i))
      eta = 1._DP - kappa
      sigma_r = kappa * sqrt(ring%Gm_pdisk / ring%r_pdisk) + eta * (2 * ring%r_pdisk * ring%w(i))

      Q = ring%w(i) * sigma_r / (3.36_DP * ring%Gsigma(i))
      tau = PI * ring%r_pdisk**2 * ring%Gsigma(i) / ring%Gm_pdisk

      nu_trans_unstable = 13 * ring%r_hstar(i)**5 * ring%Gsigma(i)**2 / ring%w(i)**3
      nu_trans_stable = sigma_r**2 / (2 * ring%w(i)) * (0.46_DP * tau / (1._DP + tau**2))

      nu_grav_stable = 0.0_DP

      nu_grav_unstable = nu_trans_unstable
      y = 0.25_DP * Q 
      kappa = ringmoons_transition_function(y)
      eta = 1._DP - kappa
      nu_trans = kappa * nu_trans_stable + eta * nu_trans_unstable
      nu_grav  = kappa * nu_grav_stable  + eta * nu_grav_unstable
      nu_coll = ring%r_pdisk**2 * ring%w(i) * tau
      ring%nu(i) = nu_trans + nu_grav + nu_coll
      ! Add in the "effective" viscosity arising from Torques
      ring%nu(i) = ring%nu(i)  + ring%Torque(i) / (3 * PI * ring%Gsigma(i) * ring%w(i) * ring%r(i)**2)
   end do


END SUBROUTINE ringmoons_viscosity
