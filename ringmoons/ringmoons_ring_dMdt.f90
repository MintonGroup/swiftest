!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_ring_dMdt
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Models the growth rate of the aggregates in unstable regions of the disk following the predator/prey model of Esposito et al. (2012)
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
!  Invocation  : CALL ringmoons_ring_dMdt(dt,ring,ring)
!
!  Notes       : 
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_ring_dMdt(Gm_pdisk,v2_pdisk,r_pdisk,tau,w) result(Gmdot)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_dMdt
      implicit none

! Arguments
      real(DP), intent(in)                   :: Gm_pdisk,v2_pdisk,r_pdisk,tau,w
      real(DP)                               :: Gmdot

! Internals
      real(DP)                               :: Torb,Tacc,Tcoll,vth2
      real(DP),parameter                     :: MBOUNCE_G = 2.1e-10_DP !See Weidling et al. (2012) eq. 8

! Executable code

      !Threshold velocity for sticking
      !See Esposito et al. (2012) and Blum (2006)
     ! vth2 = 1.0_DP

      !See Weidling et al. (2012)
      vth2 = ((Gm_pdisk * MU2GM / GU) / MBOUNCE_G)**(-5._DP / 18._DP)
      ! Convert from m/s to system units
      vth2 = vth2 * 100._DP * TU2S / DU2CM
      vth2 = max(vth2**2,2 * Gm_pdisk / r_pdisk)

      Torb = 2 * PI / w
      Tcoll = Torb / (4 * tau)
      Tacc = Tcoll

      Gmdot = (Gm_pdisk/Tcoll) * (1._DP - (v2_pdisk / vth2)) 
      
     
    
      return
end function ringmoons_ring_dMdt
