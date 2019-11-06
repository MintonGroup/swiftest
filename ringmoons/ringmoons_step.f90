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
      integer(I4B) :: i,loop
      real(DP) :: dtstab,dtleft,dt,seedmass
      real(DP),save :: Mtot_orig,Mtot_now,Ltot_orig,Ltot_now
      CHARACTER(*),parameter :: ring_outfile = "ring.dat"   ! Name of ringmoons output binary file

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
         end if
         !call ringmoons_viscosity(ring)
         Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
         Ltot_now = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active)
         Ltot_now = Ltot_now + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
         Ltot_now = Ltot_now + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
         Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
         Lerror = (Ltot_now - Ltot_orig) / Ltot_orig
      !^^^^^^^^  
      do loop = 1, LOOPMAX
         !write(*,*)
         !write(*,*) (t + (dtin - dtleft)) * 1e-6_DP
         !write(*,*) 'calc_torques'
         call ringmoons_calc_torques(swifter_pl1P,ring,seeds)
         !write(*,*) 'timestep'
         dt = ringmoons_timestep(swifter_pl1P,ring,seeds,dtleft)
         !write(*,*) 'dt = ',dt
         !write(*,*) 'seed_grow'
         call ringmoons_seed_grow(swifter_pl1P,ring,seeds,dt)
         !write(*,*) 'seed_evolve'
         call ringmoons_seed_evolve(swifter_pl1P,ring,seeds,dt)
         !write(*,*) 'viscosity'
         call ringmoons_viscosity(ring)
         !write(*,*) 'sigma_solver'
         call ringmoons_sigma_solver(ring,swifter_pl1P%mass,dt)
         !write(*,*) 'planet_accrete'
         call ringmoons_planet_accrete(swifter_pl1P,ring,seeds)
         !write(*,*) 'seed_construct'
         call ringmoons_seed_construct(swifter_pl1P,ring,seeds) ! Spawn new seeds in any available bins outside the FRL where there is ring material
         dtleft = dtleft - dt
         if (DESTRUCTION_EVENT) then
            call ringmoons_io_write_frame(t + (dtin - dtleft), ring, seeds, ring_outfile, out_stat = "APPEND")
            DESTRUCTION_COUNTER = DESTRUCTION_COUNTER + 1
            if (DESTRUCTION_COUNTER > 1000) then
               DESTRUCTION_EVENT = .false.
               DESTRUCTION_COUNTER = 0
            end if
         end if
         if (dtleft <= 0.0_DP) exit
         Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
         if (Mtot_now /= Mtot_now) then
            write(*,*) 'ERROR at time: ',t + (dtin - dtleft)
            write(*,*) 'Seeds:'
            do i = 1,seeds%N
               if (.not.seeds%active(i)) cycle
               write(*,*)  i,seeds%a(i),seeds%Gm(i)
            end do
            write(*,*) 'Ring:'
            do i = 0,ring%N + 1
               write(*,*) i,ring%r(i),ring%Gsigma(i),ring%Gm(i)
            end do
            call util_exit(FAILURE)
         end if 
         
      end do 


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
