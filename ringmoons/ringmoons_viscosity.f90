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
      real(DP),dimension(0:ring%N+1)  :: tau,sigma_r,kappa_rhstar,eta_rhstar,kappa_Q,eta_Q
      integer(I4B) :: i
      real(DP),dimension(0:ring%N+1) :: nu_trans_stable,nu_grav_stable,nu_trans_unstable,nu_grav_unstable,y,nu_trans,nu_grav,nu_coll

! Executable code

   !where(ring%Gsigma(:) > VSMALL_SQRT)
   where(ring%Gm(:) > N_DISK_FACTOR * ring%Gm_pdisk)
      kappa_rhstar(:) = ringmoons_transition_function(ring%r_hstar(:))
      eta_rhstar(:) = 1._DP - kappa_rhstar(:)
      sigma_r(:) = kappa_rhstar(:) * sqrt(ring%Gm_pdisk / ring%r_pdisk) + eta_rhstar(:) * (2 * ring%r_pdisk * ring%w(:))
      ring%Q(:) = ring%w(:) * sigma_r(:) / (3.36_DP * ring%Gsigma(:))
      tau(:) = PI * ring%r_pdisk**2 * ring%Gsigma(:) / ring%Gm_pdisk

      nu_trans_unstable(:) = 13 * ring%r_hstar(:)**5 * ring%Gsigma(:)**2 / ring%w(:)**3
      nu_trans_stable(:) = (sigma_r(:))**2 / (2 * ring%w(:)) * (0.46_DP * tau(:) / (1._DP + (tau(:))**2))

      nu_grav_stable(:) = 0.0_DP

      nu_grav_unstable(:) = nu_trans_unstable(:)
      y(:) = 0.25_DP * ring%Q(:) 
      kappa_Q(:) = ringmoons_transition_function(y(:))
      eta_Q(:) = 1._DP - kappa_Q(:)
      nu_trans = kappa_Q(:) * nu_trans_stable(:) + eta_Q(:) * nu_trans_unstable(:)
      nu_grav  = kappa_Q(:) * nu_grav_stable(:)  + eta_Q(:) * nu_grav_unstable(:)
      nu_coll(:) = ring%r_pdisk**2 * ring%w(:) * tau(:)
      ring%nu(:) = nu_trans(:) + nu_grav(:) + nu_coll(:)

   elsewhere
      ring%nu(:) = 0.0_DP
      ring%Q(:) = 0.0_DP
   end where


END SUBROUTINE ringmoons_viscosity
