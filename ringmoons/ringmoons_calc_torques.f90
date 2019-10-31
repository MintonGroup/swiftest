!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_calc_torques
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Calculates all the torques on the seeds and the ring
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
!  Invocation  : CALL ringmoons_calc_torques(swifter_pl1P,ring,seeds)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_calc_torques(swifter_pl1P,ring,seeds)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_calc_torques
   implicit none

! Arguments
   type(swifter_pl),pointer               :: swifter_pl1P
   type(ringmoons_ring), intent(inout)    :: ring
   type(ringmoons_seeds), intent(inout)   :: seeds

! Internals
   integer(I4B)                           :: i
   real(DP)                               :: e, inc, n, Ttide

! Executable code

   e = 0.0_DP
   inc = 0.0_DP
   ring%Torque = 0.0_DP
   do i = 1, seeds%N
      if (seeds%active(i)) then 
         seeds%Torque(i) = ringmoons_lindblad_torque(swifter_pl1P,ring,seeds%Gm(i),seeds%a(i),e,inc)
         n = sqrt((swifter_pl1P%mass + seeds%Gm(i)) / seeds%a(i)**3)
         Ttide = ringmoons_tidal_torque(swifter_pl1P,seeds%Gm(i),n,seeds%a(i),e,inc) 
         seeds%Torque(i) = seeds%Torque(i) + Ttide
         swifter_pl1P%rot(3) = swifter_pl1P%rot(3) - Ttide / (swifter_pl1P%mass * swifter_pl1P%Ip(3) * swifter_pl1P%radius**2)
      end if
   end do
         

   return
end subroutine ringmoons_calc_torques
