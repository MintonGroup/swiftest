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
   real(DP)                               :: dadt, e, inc,Loriginal,deltaL,P,Q,R,S
   real(DP),dimension(seeds%N)            :: daseeds

! Executable code

   !e = 0.0_DP
   !inc = 0.0_DP

   do concurrent(i=1:seeds%N,seeds%active(i)) 
      P = dt * ringmoons_seed_dadt(swifter_pl1P%mass,seeds%Gm(i),seeds%a(i),         seeds%Torque(i))
      Q = dt * ringmoons_seed_dadt(swifter_pl1P%mass,seeds%Gm(i),seeds%a(i)+0.5_DP*P,seeds%Torque(i))
      R = dt * ringmoons_seed_dadt(swifter_pl1P%mass,seeds%Gm(i),seeds%a(i)+0.5_DP*Q,seeds%Torque(i))
      S = dt * ringmoons_seed_dadt(swifter_pl1P%mass,seeds%Gm(i),seeds%a(i)+       R,seeds%Torque(i))
      daseeds(i) = (P + 2 * Q + 2 * R + S) / 6._DP
   end do

   do i = 1, seeds%N
      if (seeds%active(i)) then
         seeds%a(i) = seeds%a(i) + daseeds(i)
         if (seeds%a(i) <= ring%RRL) then   ! Destroy the satellite!
            write(*,*) 'We are on our way to destruction!'
            seeds%active(i) = .false.
            Loriginal = seeds%Gm(i) * sqrt(swifter_pl1P%mass * seeds%a(i)) 
            iRRL = ring%iRRL
            ring%Gm(iRRL) = ring%Gm(iRRL) + seeds%Gm(i)
            ring%Gsigma(iRRL) = ring%Gm(iRRL) / ring%deltaA(iRRL)
            deltaL = Loriginal - seeds%Gm(i) * ring%Iz(iRRL) * ring%w(iRRL) ! TODO: conserve angular momentum properly
         end if
      end if
   end do
         

   return
end subroutine ringmoons_seed_evolve
