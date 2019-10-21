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
SUBROUTINE ringmoons_viscosity(GM_Planet,R_Planet,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_viscosity
      IMPLICIT NONE

! Arguments
      real(DP),intent(in) :: GM_Planet,R_Planet
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      real(DP) :: Q,tau,sigma_r,kappa,eta
      integer(I4B) :: i
      real(DP) :: nu_trans_stable,nu_grav_stable,nu_trans_unstable,nu_grav_unstable,y,nu_trans,nu_grav,nu_coll
   

! Executable code

   !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) &
   !$OMP SHARED(ring)
   do i = 1, ring%N 
      tau = PI * ring%r_pdisk**2 * ring%Gsigma(i) / ring%Gm_pdisk

      kappa = ringmoons_transition_function(ring%r_hstar(i))
      eta = 1._DP - kappa
      sigma_r = kappa * sqrt(ring%Gm_pdisk / ring%r_pdisk) + eta * (2 * ring%r_pdisk * ring%w(i))

      Q = ring%w(i) * sigma_r / (3.36_DP * (ring%Gsigma(i) + TINY))
      nu_trans_stable = sigma_r**2 / (2 * ring%w(i)) * (0.46_DP * tau / (1._DP + tau**2))
      nu_grav_stable = 0.0_DP

      nu_trans_unstable = 13 * ring%r_hstar(i)**5 * ring%Gsigma(i)**2 / ring%w(i)**3
      nu_grav_unstable = nu_trans_unstable
      y = Q / 4._DP
      kappa = ringmoons_transition_function(y)
      eta = 1._DP - kappa

      nu_trans = kappa * nu_trans_stable + eta * nu_trans_unstable
      nu_grav  = kappa * nu_grav_stable  + eta * nu_grav_unstable
      nu_coll = ring%r_pdisk**2 * ring%w(i) * tau
      ring%nu(i) = nu_trans + nu_grav + nu_coll
   end do
   !$OMP END PARALLEL DO


END SUBROUTINE ringmoons_viscosity
