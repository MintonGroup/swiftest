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
   integer(I4B)                              :: i,j,iRRL
   real(DP)                                  :: dadt, e, inc, P,Q,R,S,ai
   real(DP),dimension(seeds%N)               :: daseeds
   type(ringmoons_ring)                      :: iring
   type(ringmoons_seeds)                     :: iseeds

! Executable code

   !e = 0.0_DP
   !inc = 0.0_DP
   

   !do concurrent(i=1:seeds%N,seeds%active(i)) 
   !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE (AUTO) &
   !$OMP SHARED(dt,seeds,ring,swifter_pl1P,daseeds)
   do i = 1, seeds%N
      if (seeds%active(i)) then
         iring = ring 
         iseeds = seeds

         ai = seeds%a(i)
         call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)
         P = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))

         iseeds%a(i) = ai + 0.5_DP * P
         call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)
         Q = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))

         iseeds%a(i) = ai + 0.5_DP * Q
         call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)
         R = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))
         
         iseeds%a(i) = ai + R
         call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)
         S = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))

         daseeds(i) = (P + 2 * Q + 2 * R + S) / 6._DP
      else
         daseeds(i) = 0.0_DP
      end if
   end do
   !$OMP END PARALLEL DO


   do i = 1, seeds%N

      if (seeds%active(i)) then
         !if (daseeds(i) > (ring%deltaX / 2._DP)**2) then
         !   write(*,*) 'big migration. timestep too long?'
         !   write(*,*) i,daseeds(i) / (ring%deltaX / 2._DP)**2
         !end if
         seeds%a(i) = seeds%a(i) + daseeds(i)
         if (seeds%a(i) <= ring%RRL) then   ! Destroy the satellite!
            write(*,*) 'We are on our way to destruction!'
            DESTRUCTION_EVENT = .true.
            DESTRUCTION_COUNTER = 0
            seeds%active(i) = .false.
         end if
      end if
   end do

   ! Adjust seed parameters 
   do concurrent(i=1:seeds%N,seeds%active(i))
      seeds%Rhill(i) = seeds%a(i) * (seeds%Gm(i) / (3 * swifter_pl1P%mass))**(1.0_DP / 3.0_DP)
      seeds%rbin(i) = ringmoons_ring_bin_finder(ring,seeds%a(i))
      seeds%fz_bin_inner(i) = ringmoons_ring_bin_finder(ring,seeds%a(i) - FEEDING_ZONE_FACTOR * seeds%Rhill(i))
      seeds%fz_bin_outer(i) = ringmoons_ring_bin_finder(ring,seeds%a(i) + FEEDING_ZONE_FACTOR * seeds%Rhill(i))
   end do 
         

   return
end subroutine ringmoons_seed_evolve
