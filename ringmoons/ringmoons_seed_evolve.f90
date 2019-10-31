!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_evolve
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Grows the seeds by accreting mass from within their local feeding zones from either the disk or other seeds
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
!  Invocation  : CALL ringmoons_seed_evolve(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_evolve(swifter_pl1P,ring,seeds,dt)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_evolve
   implicit none

! Arguments
   type(swifter_pl),pointer               :: swifter_pl1P
   type(ringmoons_ring), intent(inout)    :: ring
   type(ringmoons_seeds), intent(inout)   :: seeds
   real(DP), intent(in)                   :: dt

! Internals
   integer(I4B)                           :: i,iRRL
   real(DP)                               :: dadt, e, inc,Loriginal,deltaL

! Executable code

   !e = 0.0_DP
   !inc = 0.0_DP
   do i = 1, seeds%N
      if (.not.seeds%active(i)) cycle
      dadt = 2 * seeds%Torque(i) / seeds%Gm(i) * sqrt(seeds%a(i) / swifter_pl1P%mass)
      seeds%a(i) = seeds%a(i) + dadt
      if (seeds%a(i) <= ring%RRL) then   ! Destroy the satellite!
         seeds%active(i) = .false.
         Loriginal = seeds%Gm(i) * sqrt(swifter_pl1P%mass * seeds%a(i)) 
         iRRL = ring%iRRL
         ring%Gm(iRRL) = ring%Gm(iRRL) + seeds%Gm(i)
         ring%Gsigma(iRRL) = ring%Gm(iRRL) / ring%deltaA(iRRL)
         deltaL = Loriginal - seeds%Gm(i) * ring%Iz(iRRL) * ring%w(iRRL) ! TODO: conserve angular momentum properly
      end if
         
   end do
         

   return
end subroutine ringmoons_seed_evolve
