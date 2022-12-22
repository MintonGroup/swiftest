!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle) s_fraggle_resolve
   use swiftest
   use symba
contains

   module function fraggle_resolve_disruption(collision_system, nbody_system, param, t)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic disruption collision
      !! 
      implicit none
      ! Arguments
      class(fraggle_system),        intent(inout) :: collision_system
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
      real(DP),                     intent(in)    :: t            !! Time of collision
      ! Result
      integer(I4B)                                :: status    !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)          :: i, ibiggest, nfrag
      logical               :: lfailure
      character(len=STRMAX) :: message 
      real(DP) :: dpe

      associate(impactors => collision_system%impactors, fragments => collision_system%fragments,  pl => nbody_system%pl)
         select case(impactors%regime)
         case(COLLRESOLVE_REGIME_DISRUPTION)
            message = "Disruption between"
         case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
            message = "Supercatastrophic disruption between"
         end select
         call collision_io_collider_message(nbody_system%pl, impactors%id, message)
         call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)

         ! Collisional fragments will be uniformly distributed around the pre-impact barycenter
         call collision_system%set_mass_dist(param)

         ! Generate the position and velocity distributions of the fragments
         call fraggle_generate_fragments(collision_system, nbody_system, param, lfailure)

         dpe = collision_system%pe(2) - collision_system%pe(1) 
         nbody_system%Ecollisions = nbody_system%Ecollisions - dpe 
         nbody_system%Euntracked  = nbody_system%Euntracked + dpe 

         if (lfailure) then
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, "No fragment solution found, so treat as a pure hit-and-run")
            status = ACTIVE 
            nfrag = 0
            pl%status(impactors%id(:)) = status
            pl%ldiscard(impactors%id(:)) = .false.
            pl%lcollision(impactors%id(:)) = .false.
            select type(before => collision_system%before)
            class is (swiftest_nbody_system)
            select type(after => collision_system%after)
            class is (swiftest_nbody_system)
               allocate(after%pl, source=before%pl) ! Be sure to save the pl so that snapshots still work 
            end select
            end select
         else
            ! Populate the list of new bodies
            nfrag = fragments%nbody
            write(message, *) nfrag
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")
            select case(impactors%regime)
            case(COLLRESOLVE_REGIME_DISRUPTION)
               status = DISRUPTED
               ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
               fragments%id(1) = pl%id(ibiggest)
               fragments%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
               param%maxid = fragments%id(nfrag)
            case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
               status = SUPERCATASTROPHIC
               fragments%id(1:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag)]
               param%maxid = fragments%id(nfrag)
            end select

            call collision_resolve_mergeaddsub(nbody_system, param, t, status)
         end if
      end associate

      return
   end function fraggle_resolve_disruption


   module function fraggle_resolve_hitandrun(collision_system, nbody_system, param, t)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(fraggle_system),        intent(inout) :: collision_system !! Fraggle collision system object
      class(swiftest_nbody_system), intent(inout) :: nbody_system     !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
      real(DP),                     intent(in)    :: t                !! Time of collision
      ! Result
      integer(I4B)                                :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, ibiggest, nfrag, jtarg, jproj
      logical                                 :: lpure 
      character(len=STRMAX) :: message
      real(DP) :: dpe

      select type(before => collision_system%before)
      class is (swiftest_nbody_system)
      select type(after => collision_system%after)
      class is (swiftest_nbody_system)
         associate(impactors => collision_system%impactors, fragments => collision_system%fragments,  pl => nbody_system%pl)
            message = "Hit and run between"
            call collision_io_collider_message(nbody_system%pl, impactors%id, message)
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, trim(adjustl(message)))

            if (impactors%mass(1) > impactors%mass(2)) then
               jtarg = 1
               jproj = 2
            else
               jtarg = 2
               jproj = 1
            end if

            if (impactors%mass_dist(2) > 0.9_DP * impactors%mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Pure hit and run. No new fragments generated.")
               nfrag = 0
               lpure = .true.
            else ! Imperfect hit and run, so we'll keep the largest body and destroy the other
               lpure = .false.
               call collision_system%set_mass_dist(param)

               ! Generate the position and velocity distributions of the fragments
               call fraggle_generate_fragments(collision_system, nbody_system, param, lpure)

               dpe = collision_system%pe(2) - collision_system%pe(1) 
               nbody_system%Ecollisions = nbody_system%Ecollisions - dpe 
               nbody_system%Euntracked  = nbody_system%Euntracked + dpe 

               if (lpure) then
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Should have been a pure hit and run instead")
                  nfrag = 0
               else
                  nfrag = fragments%nbody
                  write(message, *) nfrag
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")
               end if
            end if
            if (lpure) then ! Reset these bodies back to being active so that nothing further is done to them
               status = HIT_AND_RUN_PURE
               pl%status(impactors%id(:)) = ACTIVE
               pl%ldiscard(impactors%id(:)) = .false.
               pl%lcollision(impactors%id(:)) = .false.
               allocate(after%pl, source=before%pl) ! Be sure to save the pl so that snapshots still work 
            else
               ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
               fragments%id(1) = pl%id(ibiggest)
               fragments%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
               param%maxid = fragments%id(nfrag)
               status = HIT_AND_RUN_DISRUPT
               call collision_resolve_mergeaddsub(nbody_system, param, t, status)
            end if
         end associate
      end select
      end select


      return
   end function fraggle_resolve_hitandrun

end submodule s_fraggle_resolve