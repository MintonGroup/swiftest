!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_ring_dvdt
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Models the rate of change of velocity dispersion of the aggregates in unstable regions of the disk following the 
!                predator/prey model of Esposito et al. (2012) 
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
!  Invocation  : CALL ringmoons_ring_dvdt(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_ring_dvdt(Gm_pdisk,v2_pdisk,tau,nu,w) result(v2dot)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_dvdt
      implicit none

! Arguments
      real(DP), intent(in)                   :: Gm_pdisk,v2_pdisk,tau,nu,w 
      real(DP)                               :: v2dot

! Internals
      real(DP)                               :: Torb,Tacc,Tcoll,S2
      real(DP)                               :: eps2 ! coefficient of restitution
      real(DP), parameter                    :: Vc = 0.0077 ! See Schmidt et al. (2006) eqn. 14.14
      real(DP), parameter                    :: eps_exponent = -0.234_DP
      real(DP), parameter                    :: eps_constant = 0.89_DP ! See Brisset et al. (2019)

! Executable code
  
      Torb = 2 * PI / w
      Tcoll = Torb / (4 * tau) 
      !eps2 = (v2_pdisk * (DU2CM / TU2S)**2 / Vc**2)**(eps_exponent)
      !eps2 = 0.9_DP * exp(-0.22_DP * (sqrt(v2_pdisk) * DU2CM / TU2S)) + 0.01_DP * (sqrt(v2_pdisk) * DU2CM / TU2S)**(-0.6_DP)
      eps2 = eps_constant**2
      S2 = 9 * w**2 / 4._DP ! Shear rate squared
      
      v2dot = nu * S2 - v2_pdisk * (1._DP - eps2) / Tcoll

    
      return
end function ringmoons_ring_dvdt
