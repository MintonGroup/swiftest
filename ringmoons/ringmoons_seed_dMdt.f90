!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_dMdt
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Models the growth rate of a seed embedded within the ring using a runaway growth model with a factor that corrects
!                for distance from the FRL
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
!  Invocation  : CALL ringmoons_seed_dMdt(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_seed_dMdt(ring,GMP,Gsigma,Gmseed,rhoseed,a) result(Gmdot)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_dMdt
      implicit none

! Arguments
      type(ringmoons_ring), intent(in)       :: ring
      real(DP), intent(in)                   :: GMP,Gsigma,Gmseed,rhoseed,a
      real(DP)                               :: Gmdot

! Internals
      integer(I4B)                           :: i
      real(DP)                               :: C
      real(DP),parameter                     :: eff2   = 1e-10_DP ! This term gets the growth rate to match up closely to Andy's
      real(DP),parameter                     :: growth_exponent = 4._DP / 3._DP
      real(DP), parameter                    :: SIGLIMIT = 1e-100_DP

! Executable code
     
      C = 12 * PI**(2._DP / 3._DP) * (3._DP / (4 * rhoseed))**(1._DP / 3._DP) / sqrt(GMP) 
      Gmdot = C * Gsigma / (eff2 * sqrt(a)) * Gmseed**(growth_exponent)
    
      return
end function ringmoons_seed_dMdt
