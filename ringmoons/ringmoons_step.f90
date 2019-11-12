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
      integer(I4B) :: i,loop,seedloop,subcount
      integer(I4B), parameter :: submax = 2
      real(DP) :: dtstab,dtleft,dtring,seedmass,dtseed,dtseedleft
      real(DP),save :: Mtot_orig,Mtot_now,Ltot_orig,Ltot_now
      CHARACTER(*),parameter :: ring_outfile = "ring.dat"   ! Name of ringmoons output binary file
      type(ringmoons_ring)                            :: old_ring
      type(ringmoons_seeds)                           :: old_seeds
      logical(LGT)                                    :: stepfail
      type(swifter_pl), pointer                       :: tmp_swifter_pl1P
      real(DP),dimension(0:ring%N+1)                  :: dTorque
      real(DP)                                        :: dMerror,Merror_old
      real(DP),save                                   :: dMerror_max = 0._DP

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
         !write(*,*) 'calc_torques'
         call ringmoons_calc_torques(swifter_pl1P,ring,seeds)
      end if

      !^^^^^^^^  

      do loop = 1, LOOPMAX
         if (loop == LOOPMAX) then
            write(*,*) 'LOOPMAX reached. Ringmoons_step failed'
            call util_exit(FAILURE)
         end if
         !write(*,*)
         !write(*,*) (t + (dtin - dtleft)) * 1e-6_DP

         !write(*,*) 'viscosity'
         call ringmoons_viscosity(ring)

         !write(*,*) 'timestep'
         dtring = ringmoons_ring_timestep(swifter_pl1P,ring,dtleft)
         !write(*,*) 'dtring = ',dtring

         !write(*,*) 'planet_accrete'
         !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
         !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
         !Merror_old = Merror
         !write(*,*) 'Merror = ',Merror
         call ringmoons_planet_accrete(swifter_pl1P,ring,seeds,dtring)
         !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
         !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
         !write(*,*) 'Merror = ',Merror
         !dMerror = (Merror - Merror_old) 
         !if (abs(dMerror) > dMerror_max) then
         !   write(*,*) 'A max in the mass error occurred!'
         !   write(*,*) 'after planet_accrete'
         !   write(*,*) 'change = ',dMerror
         !   dMerror_max = abs(dMerror)
         !end if
         !if (abs(dMerror) > 10 * epsilon(1._DP)) call util_exit(FAILURE)
         !Merror_old = Merror

         !write(*,*) 'before sigma_solver'
         !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
         !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
         !write(*,*) 'Merror = ',Merror
         call ringmoons_sigma_solver(ring,swifter_pl1P%mass,dtring)
         !write(*,*) 'sigma_solver complete'
         !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
         !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
         !write(*,*) 'Merror = ',Merror
         !dMerror = (Merror - Merror_old) 
         !if (abs(dMerror) > dMerror_max) then
         !   write(*,*) 'A max in the mass error occurred!'
         !   write(*,*) 'after sigma_solver'
         !   write(*,*) 'change = ',dMerror
         !   dMerror_max = abs(dMerror)
         !end if
         !if (abs(dMerror) > 10 * epsilon(1._DP)) call util_exit(FAILURE)
         !Merror_old = Merror
         ring%Torque(:) = 0.0_DP


         old_ring%N = ring%N
         old_seeds%N = seeds%N
         call ringmoons_allocate(old_ring,old_seeds)
         old_ring = ring
         old_seeds = seeds

         dtseedleft = dtring
         dtseed = dtring 
         subcount = 0
         do seedloop = 1, LOOPMAX 
            if (seedloop == LOOPMAX) then
               write(*,*) 'LOOPMAX reached in seed evolution. Ringmoons_step failed'
               call util_exit(FAILURE)
            end if
            !write(*,*) 'before seed_construct'
            !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
            !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
            !Merror_old = Merror
            !write(*,*) 'Merror = ',Merror
            call ringmoons_seed_construct(swifter_pl1P,ring,seeds) ! Spawn new seeds in any available bins outside the FRL where there is ring material
            !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
            !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
            !write(*,*) 'Merror = ',Merror
            !dMerror = (Merror - Merror_old) 
            !if (abs(dMerror) > dMerror_max) then
            !   write(*,*) 'A max in the mass error occurred!'
            !   write(*,*) 'after seed_construct'
            !   write(*,*) 'change = ',dMerror
            !   dMerror_max = abs(dMerror)
               !read(*,*)
            !end if
            !if (abs(dMerror) > 10 * epsilon(1._DP)) call util_exit(FAILURE)
            !Merror_old = Merror
            call ringmoons_update_seeds(swifter_pl1P,ring,seeds)
            dtseed = ringmoons_seed_timestep(swifter_pl1P,ring,seeds,dtseed) 
            !ring%Torque(:) = 0.0_DP  
            !write(*,*) 'evolve'
            !write(*,*) seedloop,'evolve',dtring/dtseed
            !write(*,*) 'before seed_evolve'
            !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
            !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
            !write(*,*) 'Merror = ',Merror
            call ringmoons_seed_evolve(swifter_pl1P,ring,seeds,dtseed,stepfail)

            if (stepfail) then
               dtseed = dtseed / submax
               subcount = 0
               ring = old_ring
               seeds = old_seeds
               !write(*,*) 'Failed the step',dtring/dtseed
            else
               !Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm,seeds%active)
               !Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
               !write(*,*) 'Merror = ',Merror
               !dMerror = (Merror - Merror_old) 
               !if (abs(dMerror) > dMerror_max) then
               !   write(*,*) 'A max in the mass error occurred!'
               !   write(*,*) 'after seed_evolve'
               !   write(*,*) 'change = ',dMerror
               !   dMerror_max = abs(dMerror)
               !end if
               !if (abs(dMerror) > 10 * epsilon(1._DP)) call util_exit(FAILURE)
               !Merror_old = Merror
               subcount = subcount + 1
               if (DESTRUCTION_EVENT) then
                  call ringmoons_io_write_frame(t + (dtin - dtleft - dtseedleft), ring, seeds, ring_outfile, out_stat = "APPEND")
                  DESTRUCTION_COUNTER = DESTRUCTION_COUNTER + 1
                  if (DESTRUCTION_COUNTER > 1000) then
                     DESTRUCTION_EVENT = .false.
                     DESTRUCTION_COUNTER = 0
                  end if
               end if


               dtseedleft = dtseedleft - dtseed
               ! Scale the change in the ring torques by the step size reduction in order to get the time-averaged Torque
               dTorque(:) = ring%Torque(:) - old_ring%Torque(:) 
               ring%Torque(:) = old_ring%Torque(:) + dTorque(:) * (dtseed / dtring)
               if (dtseedleft <= 0.0_DP) exit
               if (subcount == 2* submax) then
                  dtseed = min(dtseedleft, submax * dtseed)
                  subcount = 0
               end if
               old_ring = ring
               old_seeds = seeds
            end if
         end do
         call ringmoons_deallocate(old_ring,old_seeds)

         dtleft = dtleft - dtring

         if (dtleft <= 0.0_DP) exit
      end do 
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
