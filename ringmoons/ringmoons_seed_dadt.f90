!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_dadt
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Models the migration rate of a seed from its torques
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
!  Invocation  : CALL ringmoons_seed_dadt(GMP,Gmseed,a,Torque)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_seed_dadt(GMP,Gmseed,a,Torque,mdot) result(adot)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_dadt
      implicit none

! Arguments
      real(DP), intent(in)                   :: GMP,Gmseed,a,Torque,mdot
      real(DP)                               :: adot

! Internals

! Executable code
      adot = 2 * (Torque / (sqrt((GMP + Gmseed) * a))  &
               - mdot * (1._DP  + Gmseed / (2 *(GMP + Gmseed))))
      adot = adot * a / Gmseed
      return
end function ringmoons_seed_dadt
