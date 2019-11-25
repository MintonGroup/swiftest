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
elemental function ringmoons_ring_dMdt(Gm_pdisk,v2_pdisk,tau,w) result(Gmdot)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_dMdt
      implicit none

! Arguments
      real(DP), intent(in)                   :: Gm_pdisk,v2_pdisk,tau,w
      real(DP)                               :: Gmdot

! Internals
      real(DP)                               :: Torb,Tacc,Tcoll,vth2

! Executable code

      !Threshold velocity for sticking
      !See Esposito et al. (2012) and Blum (2006)
      vth2 = (1.0e2_DP * TU2S / DU2CM)**2

      Torb = 2 * PI / w
      Tcoll = Torb / (4 * tau)
      Tacc = Tcoll

      Gmdot = Gm_pdisk / Tacc - (v2_pdisk / vth2) * (Gm_pdisk / Tcoll)
      
     
    
      return
end function ringmoons_ring_dMdt
