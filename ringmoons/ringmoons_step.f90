!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step the fluid ring in time
!
!  Input
!    Arguments   swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!                ring           : Ringmoons data structure
!                dtin           : input time step
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                ring           : Updated ring
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL ringmoons_step(swifter_pl1P,ring,dtin,lfirst)
!
!  Notes       : Adapted from Andy Hesselbrock's RING-MOONS Python scripts
!
!**********************************************************************************************************************************
subroutine ringmoons_step(t,swifter_pl1P,ring,seeds,dtin,lfirst,Merror,Lerror)

! Modules
     use module_parameters
     use module_swifter
     use module_ringmoons
     use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_step
     implicit none

! Arguments
      type(swifter_pl), pointer                       :: swifter_pl1P
      real(DP), intent(in)                            :: t,dtin
      type(ringmoons_ring),intent(inout)              :: ring
      type(ringmoons_seeds),intent(inout)             :: seeds
      logical(LGT), intent(inout)                     :: lfirst
      real(DP),intent(out)                            :: Merror,Lerror

! Internals
      integer(I4B) :: i,loop,seedloop,subcount,loopcount
      integer(I4B), parameter :: submax = 2
      real(DP) :: dt,dtleft,dtnewring
      real(DP),save :: Mtot_orig,Mtot_now,Ltot_orig,Ltot_now,Lerr_now,dLerror
      CHARACTER(*),parameter :: ring_outfile = "ring.dat"   ! Name of ringmoons output binary file
      type(ringmoons_ring)                            :: old_ring
      type(ringmoons_seeds)                           :: old_seeds
      logical(LGT)                                    :: stepfail
      type(swifter_pl)                                :: old_swifter_pl1P
      real(DP),dimension(0:ring%N+1)                  :: dTorque
      real(DP),save                                   :: dLerror_max = 0._DP
      real(DP),save                                   :: Lerr_old
      real(DP),save :: dL_planet_accrete = 0.0_DP, dL_seed_evolve = 0.0_DP, dL_sigma_solver = 0.0_DP, dL_seed_construct = 0.0_DP

! Executable code
      !if (lfirst) then
      dtleft = dtin

      !TESTING
      !write(*,*) 'Ring mass: ',sum(ring%Gm)
      if (lfirst) then
         Mtot_orig = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
         Ltot_orig = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active)
         Ltot_orig = Ltot_orig + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
         Ltot_orig = Ltot_orig + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
         lfirst = .false.
         Lerr_old = 0.0_DP
         !write(*,*) 'calc_torques'
         !call ringmoons_calc_torques(swifter_pl1P,ring,seeds)
      end if
      !^^^^^^^^  
!write(*,*)
!write(*,*) (t + (dtin - dtleft)) * 1e-6_DP
!write(*,*) 'viscosity'
      seeds%Torque(:) = 0.0_DP

      old_ring%N = ring%N
      old_seeds%N = seeds%N
      call ringmoons_allocate(old_ring,old_seeds)
      old_ring = ring
      old_seeds = seeds
      old_swifter_pl1P%mass = swifter_pl1P%mass
      old_swifter_pl1P%radius = swifter_pl1P%radius
      old_swifter_pl1P%rot = swifter_pl1P%rot

      subcount = 0
      dt = dtleft
      do loop = 1, LOOPMAX
         if (loop == LOOPMAX) then
            write(*,*) 'LOOPMAX reached in seed evolution. Ringmoons_step failed'
            call util_exit(FAILURE)
         end if

         ring%Torque(:) = 0.0_DP
         seeds%Torque(:) = 0.0_DP
!write(*,*) 'timestep'
         call ringmoons_viscosity(ring)
         dt = ringmoons_ring_timestep(swifter_pl1P,ring,dt)
         dt = ringmoons_seed_timestep(swifter_pl1P,ring,seeds,dt) 

!write(*,*) 'planet_accrete'
         call ringmoons_planet_accrete(swifter_pl1P,ring,seeds,dt)
!Ltot_now = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active)
!Ltot_now = Ltot_now + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
!Ltot_now = Ltot_now + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
!Lerr_now = (Ltot_now - Ltot_orig) / Ltot_orig
!dLerror = Lerr_now - Lerr_old
!dL_planet_accrete = dL_planet_accrete + dLerror
!!write(*,*) 'dL/L0 = ',Lerr_now
!!write(*,*) 'change= ',dLerror
!Lerr_old = Lerr_now


!write(*,*) 'seed_evolve'
         call ringmoons_seed_evolve(swifter_pl1P,ring,seeds,dt,stepfail)

         if (stepfail) then
            dt = dt / submax
            subcount = 0
            ring = old_ring
            seeds = old_seeds
            swifter_pl1P%mass = old_swifter_pl1P%mass
            swifter_pl1P%radius = old_swifter_pl1P%radius
            swifter_pl1P%rot = old_swifter_pl1P%rot
            cycle
         end if

!Ltot_now = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active)
!Ltot_now = Ltot_now + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
!Ltot_now = Ltot_now + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
!Ltot_now = Ltot_now + sum(ring%Torque(:)) * dt
!Lerr_now = (Ltot_now - Ltot_orig) / Ltot_orig
!dLerror = Lerr_now - Lerr_old
!dL_seed_evolve = dL_seed_evolve + dLerror
!write(*,*) 'dL/L0 = ',Lerr_now
!!write(*,*) 'change= ',dLerror
!write(*,*) 'dtT/L0= ',sum(ring%Torque(:)) * dt / Ltot_orig
!Lerr_old = Lerr_now


!write(*,*) 'sigma_solver'
         call ringmoons_viscosity(ring)
         dtnewring = ringmoons_ring_timestep(swifter_pl1P,ring,dt)
         if (dtnewring < dt) then
            !Step size is now too small
            dt = dtnewring
            ring = old_ring
            seeds = old_seeds
            swifter_pl1P%mass = old_swifter_pl1P%mass
            swifter_pl1P%radius = old_swifter_pl1P%radius
            swifter_pl1P%rot = old_swifter_pl1P%rot
            cycle
         end if
         call ringmoons_sigma_solver(ring,swifter_pl1P%mass,dt)
!Ltot_now = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active)
!Ltot_now = Ltot_now + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
!Ltot_now = Ltot_now + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
!Lerr_now = (Ltot_now - Ltot_orig) / Ltot_orig
!dLerror = Lerr_now - Lerr_old
!dL_sigma_solver = dL_sigma_solver + dLerror
!!write(*,*) 'dL/L0 = ',Lerr_now
!!write(*,*) 'change= ',dLerror
!Lerr_old = Lerr_now


!write(*,*) 'seed_construct'
         call ringmoons_seed_construct(swifter_pl1P,ring,seeds) ! Spawn new seeds in any available bins outside the FRL where there is ring material
!Ltot_now = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active)
!Ltot_now = Ltot_now + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
!Ltot_now = Ltot_now + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
!Lerr_now = (Ltot_now - Ltot_orig) / Ltot_orig
!dLerror = Lerr_now - Lerr_old
!dL_seed_construct = dL_seed_construct + dLerror
!!   write(*,*) 'dL/L0 = ',Lerr_now
!!   write(*,*) 'change= ',dLerror
!Lerr_old = Lerr_now

         subcount = subcount + 1
         if (DESTRUCTION_EVENT) then
            call ringmoons_io_write_frame(t + (dtin - dtleft - dtleft), ring, seeds, ring_outfile, out_stat = "APPEND")
            DESTRUCTION_COUNTER = DESTRUCTION_COUNTER + 1
            if (DESTRUCTION_COUNTER > 1000) then
               DESTRUCTION_EVENT = .false.
               DESTRUCTION_COUNTER = 0
            end if
         end if

         dtleft = dtleft - dt
         ! Scale the change in the ring torques by the step size reduction in order to get the time-averaged Torque
         if (dtleft <= 0.0_DP) exit
         loopcount = loop
         if (subcount == 2 * submax) then
            dt = max(min(dtleft, submax * dt),dtin)
            subcount = 0
         end if
         old_ring = ring
         old_seeds = seeds
         old_swifter_pl1P%mass = swifter_pl1P%mass
         old_swifter_pl1P%radius = swifter_pl1P%radius
         old_swifter_pl1P%rot = swifter_pl1P%rot
      end do
      call ringmoons_deallocate(old_ring,old_seeds)

      Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
      Ltot_now = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active)
      Ltot_now = Ltot_now + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
      Ltot_now = Ltot_now + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
      Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
      Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
      Lerror = (Ltot_now - Ltot_orig) / Ltot_orig
      !if ((Mtot_now /= Mtot_now).or.(Ltot_now /= Ltot_now)) then
      !   write(*,*) 'ERROR at time: ',t + (dtin - dtleft)
      !   write(*,*) 'Seeds:'
      !   do i = 1,seeds%N
      !      if (.not.seeds%active(i)) cycle
      !      write(*,*)  i,seeds%a(i),seeds%Gm(i)
      !   end do
      !   write(*,*) 'Ring:'
      !   do i = 0,ring%N + 1
      !      write(*,*) i,ring%r(i),ring%Gsigma(i),ring%Gm(i)
      !   end do
      !   call util_exit(FAILURE)
      !end if 
!!      write(*,*) 'Loop count: ',loopcount
!      write(*,*) 'Lerror: ',Lerror
!      write(*,*) 'dL_planet_accrete = ',dL_planet_accrete
!      write(*,*) 'dL_seed_evolve = ',dL_seed_evolve
!      write(*,*) 'dL_sigma_solver = ',dL_sigma_solver
!      write(*,*) 'dL_seed_construct = ',dL_seed_construct



      RETURN

END SUBROUTINE ringmoons_step
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
