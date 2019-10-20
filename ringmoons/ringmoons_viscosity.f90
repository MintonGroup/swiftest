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
SUBROUTINE ringmoons_viscocity(M_Planet,R_Planet,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_viscocity
      IMPLICIT NONE

! Arguments
      real(DP),intent(in) :: M_Planet,R_Planet
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      real(DP) :: Q,r_hstar,tau
      integer(I4B) :: i
      real(DP),parameter :: sigsmall = 1e-6_DP * MU2GM
   

! Executable code


    do i = 1, ring%N 
      if (ring%sigma(i) <= sigsmall) then
         ring%nu(i) = 0.0_DP
      else
            r_hstar = /(2.0 * ring%r_pdisk)*(2.0*ring%m_pdisk/(3.0*M_Planet))**(1.0/3.0)
            tau = PI*r_pdisk**2*sigma[a]/m_pdisk


            if r_hstar < 1.0:
                sigma_r = 0.5*(1+tanh((2.*r_hstar-1)/(r_hstar*(r_hstar-1))))*sqrt(G*m_pdisk/r_pdisk) + 0.5*(1-tanh((2.*r_hstar-1)/(r_hstar*(r_hstar-1))))*(2.0*r_pdisk*w[a])
            else:
                sigma_r = sqrt(G*m_pdisk/r_pdisk)


            Q = w[a]*sigma_r/(3.36*G*(sigma[a]+0.000001))

            if Q <= 4.0:

                nu_trans_stable = sigma_r**2.0/(2.0*w[a])*(0.46*tau/(1.0+tau**2.0))
                nu_grav_stable = 0.0
                nu_trans_unstable = 13.0*r_hstar**5.0*G**2.0*sigma[a]**2.0/w[a]**3.0
                nu_grav_unstable = nu_trans_unstable

                y = Q/4.
                nu_trans = 0.5*(1+tanh((2.*y-1)/(y*(y-1))))*(nu_trans_stable) +  0.5*(1-tanh((2.*y-1)/(y*(y-1))))*(nu_trans_unstable)
                nu_grav = 0.5*(1+tanh((2.*y-1)/(y*(y-1))))*(nu_grav_stable) +  0.5*(1-tanh((2.*y-1)/(y*(y-1))))*(nu_grav_unstable)

            else:
                nu_trans = sigma_r**2.0/(2.0*w[a])*(0.46*tau/(1.0+tau**2.0))
                nu_grav = 0.0

            nu_coll = r_pdisk**2.0*w[a]*tau
            nu[a] = nu_trans + nu_grav + nu_coll

END SUBROUTINE ringmoons_viscocity
