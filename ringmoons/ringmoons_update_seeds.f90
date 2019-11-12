!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_update_seeds
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Updates seed and ring values (Torques, feeding zones, etc) between RK steps
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
!  Invocation  : CALL ringmoons_update_seeds(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_update_seeds(swifter_pl1P,ring,seeds)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_update_seeds
   implicit none

! Arguments
   type(swifter_pl),pointer               :: swifter_pl1P
   type(ringmoons_ring), intent(in)       :: ring
   type(ringmoons_seeds), intent(inout)   :: seeds

! Internals
   real(DP),dimension(seeds%N)            :: fz_width
   integer(I4B)                           :: i
   

! Executable code

   ! Compute feeding zone info
   do i = 1,seeds%N
      seeds%Rhill(i) = seeds%a(i) * (seeds%Gm(i) / (3 * swifter_pl1P%mass))**(1.0_DP / 3.0_DP)
      fz_width(i) = FEEDING_ZONE_FACTOR * seeds%Rhill(i)
      seeds%rbin(i) = ringmoons_ring_bin_finder(ring,seeds%a(i))
      seeds%fz_bin_inner(i) = ringmoons_ring_bin_finder(ring,seeds%a(i) - fz_width(i))
      seeds%fz_bin_outer(i) = ringmoons_ring_bin_finder(ring,seeds%a(i) + fz_width(i))
   end do

   return
end subroutine ringmoons_update_seeds
