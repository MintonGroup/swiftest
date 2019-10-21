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
      real(DP) :: Q,r_hstar,tau,sigma_r,kappa,eta
      integer(I4B) :: i
      real(DP) :: sigsmall
      real(DP) :: nu_trans_stable,nu_grav_stable,nu_trans_unstable,nu_grav_unstable,y,nu_trans,nu_grav,nu_coll
   

! Executable code
   sigsmall = GU * 1e-6_DP / MU2GM

   r_hstar = R_Planet / (2 * ring%r_pdisk) * (2 *ring%m_pdisk /(3._DP * GM_Planet))**(1._DP/3._DP)
   !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(AUTO) &
   !$OMP SHARED(ring,GU,GM_Planet,r_hstar,sigsmall)
   do i = 1, ring%N 
      if (ring%sigma(i) <= sigsmall) then
         ring%nu(i) = 0.0_DP
      else
         tau = PI * ring%r_pdisk**2 * ring%sigma(i) / ring%m_pdisk
         if (r_hstar < 1._DP) then
             kappa = ringmoons_transition_function(r_hstar)
             eta = 1._DP - kappa
             sigma_r = kappa * sqrt(ring%m_pdisk / ring%r_pdisk) + eta * (2 * ring%r_pdisk * ring%w(i))
         else
             sigma_r = sqrt(ring%m_pdisk / ring%r_pdisk)
         end if

         Q = ring%w(i) * sigma_r / (3.36_DP * (ring%sigma(i) + TINY))
         nu_trans_stable = sigma_r**2 / (2 * ring%w(i)) * (0.46_DP * tau / (1._DP + tau**2))
         nu_grav_stable = 0.0_DP

         if (Q <= 4._DP) then
             nu_trans_unstable = 13 * r_hstar**5 * ring%sigma(i)**2 / ring%w(i)**3
             nu_grav_unstable = nu_trans_unstable
             y = Q / 4._DP
             kappa = ringmoons_transition_function(y)
             eta = 1._DP - kappa

             nu_trans = kappa * nu_trans_stable + eta * nu_trans_unstable
             nu_grav  = kappa * nu_grav_stable  + eta * nu_grav_unstable
         else
             nu_trans = nu_trans_stable
             nu_grav = nu_grav_stable
         end if
         nu_coll = ring%r_pdisk**2 * ring%w(i) * tau
         ring%nu(i) = nu_trans + nu_grav + nu_coll
      end if
   end do
   !$OMP END PARALLEL DO

END SUBROUTINE ringmoons_viscosity
