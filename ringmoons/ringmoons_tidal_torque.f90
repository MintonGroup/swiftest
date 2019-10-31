!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_tidal_lindblad
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Calculates the lindblad torques between each ring element and a given satellite. Function returns total torque on the
!                satellite, and stores the torques acting on each ring element in the ring
!
!  Input
!    Arguments : 
!                
!    Terminal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : Torque = ringmoons_tidal_torque(swifter_pl1P,ring,Gm,a,e,inc)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
function ringmoons_tidal_torque(swifter_pl1P,Gm,n,a,e,inc) result(Torque)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_tidal_torque
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      real(DP),intent(in)                    :: Gm, a, n, e, inc
      real(DP)                               :: Torque
      

! Internals
      integer(I4B)                           :: i

! Executable code

     
      Torque = sign(1._DP,swifter_pl1P%rot(3) - n) * &
               1.5_DP * swifter_pl1P%k2 / swifter_pl1P%Q * &
               Gm**2 * swifter_pl1P%radius**5 / a**6

      return
end function ringmoons_tidal_torque
