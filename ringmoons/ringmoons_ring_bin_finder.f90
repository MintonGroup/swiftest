!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_ring_bin_finder
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Finds the ring bin corresponding to a distance r
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
!  Invocation  : CALL ringmoons_ring_bin_finder(ring,r) result(bin)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_ring_bin_finder(ring,r) result(bin)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_bin_finder
      implicit none

! Arguments
      type(ringmoons_ring), intent(in)       :: ring
      real(DP), intent(in)                   :: r
      integer(I4B)                           :: bin

! Internals

! Executable code

      bin = ceiling(2 * (sqrt(r) - sqrt(ring%r_I)) / ring%deltaX) !+ 1 
      bin = min(max(1,bin),ring%N)
    
      return
end function ringmoons_ring_bin_finder
