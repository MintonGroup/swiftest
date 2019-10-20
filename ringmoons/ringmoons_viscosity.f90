!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_viscocity
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
!  Invocation  : CALL ringmoons_viscocity(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
SUBROUTINE ringmoons_viscocity(GM_Planet,R_Planet,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_viscocity
      IMPLICIT NONE

! Arguments
      real(DP),intent(in) :: GM_Planet,R_Planet
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      real(DP) :: Q,r_hstar,tau,sigma_r
      integer(I4B) :: i
      real(DP) :: sigsmall
      real(DP) :: nu_trans_stable,nu_grav_stable,nu_trans_unstable,nu_grav_unstable,y,nu_trans,nu_grav,nu_coll
   

! Executable code
   sigsmall = 1e-6_DP * MU2GM

   r_hstar = R_Planet / (2 * ring%r_pdisk) * (2 *ring%m_pdisk /(3.0_DP * M_Planet))**(1.0/3.0)
   do i = 1, ring%N 
      if (ring%sigma(i) <= sigsmall) then
         ring%nu(i) = 0.0_DP
      else
         tau = PI * ring%r_pdisk**2 * ring%sigma(i) / ring%m_pdisk
         if (r_hstar < 1.0_DP) then
             sigma_r = 0.5_DP * (1 + tanh((2 * r_hstar - 1.0_DP) / (r_hstar * (r_hstar - 1.0_DP)))) &
                      * sqrt(ring%m_pdisk / ring%r_pdisk) + 0.5_DP * (1.0_DP - tanh((2 * r_hstar - 1.0_DP) &
                      / (r_hstar * (r_hstar - 1._DP)))) * (2 * ring%r_pdisk * ring%w(i))
         else
             sigma_r = sqrt(ring%m_pdisk / ring%r_pdisk)
         end if

         Q = ring%w(i) * sigma_r / (3.36_DP * (ring%sigma(i)))

         if (Q <= 4.0_DP) then
             nu_trans_stable = sigma_r**2 / (2 *ring%w(i)) * (0.46 * tau /(1.0 + tau**2))
             nu_grav_stable = 0.0_DP
             nu_trans_unstable = 13 * r_hstar**5 * ring%sigma(i)**2 / ring%w(i)**3
             nu_grav_unstable = nu_trans_unstable

             y = Q / 4.0_DP
             nu_trans = 0.5_DP * (1._DP + tanh((2 * y - 1._DP) / (y * (y - 1._DP)))) * nu_trans_stable &
                      + 0.5_DP * (1._DP - tanh((2 * y - 1._DP) / (y * (y - 1._DP)))) * nu_trans_unstable
             nu_grav  = 0.5_DP * (1._DP + tanh((2 * y - 1._DP) / (y * (y - 1._DP)))) * nu_grav_stable &
                      + 0.5_DP * (1._DP - tanh((2 * y - 1._DP) / (y * (y - 1._DP)))) * nu_grav_unstable
         else
             nu_trans = sigma_r**2 / (2 * ring%w(i)) * (0.46_DP * tau / (1._DP + tau**2))
             nu_grav = 0.0_DP
         end if

         nu_coll = ring%r_pdisk**2 * ring%w(i) * tau
         nu(i) = nu_trans + nu_grav + nu_coll
      end if
   end do

END SUBROUTINE ringmoons_viscocity
