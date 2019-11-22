!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_transition_function
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Transition function for viscosity regime parameters
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
!  Invocation  : kappa = ringmoons_transition_function(y)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_transition_function(yin) result(kappa)

! Modules
      use module_parameters
      IMPLICIT NONE

! Arguments
      real(DP),intent(in) ::yin
      real(DP) :: kappa

! Internals
      real(DP) :: y
   

! Executable code

   y = min(max(yin,epsilon(1._DP)),1.0_DP-epsilon(1._DP))
   kappa =  0.5_DP * (1._DP + tanh((2 * y - 1._DP) / (y * (1._DP - y)))) 
   return

end function ringmoons_transition_function
