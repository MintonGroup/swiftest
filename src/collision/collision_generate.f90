
!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(collision) s_collision_generate
   use swiftest
contains

   module subroutine collision_generate_bounce(self, nbody_system, param, t)
      !! author: David A. Minton
      !!
      !! In this collision model, if the collision would result in a disruption, the bodies are instead "bounced" off 
      !! of the center of mass. This is done as a reflection in the 2-body equivalent distance vector direction.
      implicit none
      ! Arguments
      class(collision_bounce),  intent(inout) :: self         !! Bounce fragment system object 
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! The time of the collision
      ! Internals
      integer(I4B) :: i,j,nfrag
      real(DP), dimension(NDIM) :: vcom, rnorm

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type (pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(impactors => nbody_system%collider%impactors, fragments => nbody_system%collider%fragments)
            select case (impactors%regime) 
            case (COLLRESOLVE_REGIME_DISRUPTION, COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
               nfrag = size(impactors%id(:))
               do i = 1, nfrag
                  j = impactors%id(i)
                  vcom(:) = pl%vb(:,j) - impactors%vbcom(:)
                  rnorm(:) = .unit. (impactors%rb(:,2) - impactors%rb(:,1))
                  ! Do the reflection
                  vcom(:) = vcom(:) - 2 * dot_product(vcom(:),rnorm(:)) * rnorm(:)
                  pl%vb(:,j) = impactors%vbcom(:) + vcom(:)
                  self%status = DISRUPTED
                  pl%status(j) = ACTIVE
                  pl%ldiscard(j) = .false.
                  pl%lcollision(j) = .false.
               end do
               select type(before => self%before)
               class is (swiftest_nbody_system)
               select type(after => self%after)
               class is (swiftest_nbody_system)
                  allocate(after%pl, source=before%pl) ! Be sure to save the pl so that snapshots still work   
               end select
               end select
            case (COLLRESOLVE_REGIME_HIT_AND_RUN)
               call collision_generate_hitandrun(self, nbody_system, param, t)
            case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
               call self%collision_merge%generate(nbody_system, param, t) ! Use the default collision model, which is merge
            case default 
               write(*,*) "Error in swiftest_collision, unrecognized collision regime"
               call util_exit(FAILURE)
            end select
            end associate
      end select
      end select

      return
   end subroutine collision_generate_bounce


   module subroutine collision_generate_hitandrun(collider, nbody_system, param, t) 
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(collision_merge),   intent(inout) :: collider !! Fraggle collision system object
      class(base_nbody_system), intent(inout) :: nbody_system     !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
      real(DP),                 intent(in)    :: t                !! Time of collision
      ! Internals
      integer(I4B)                            :: i, ibiggest, nfrag, jtarg, jproj
      logical                                 :: lpure 
      character(len=STRMAX) :: message
      real(DP) :: dpe

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type (pl => nbody_system%pl)
      class is (swiftest_pl)

         associate(impactors => collider%impactors, fragments => collider%fragments,  pl => nbody_system%pl)
            message = "Hit and run between"
            call collision_io_collider_message(nbody_system%pl, impactors%id, message)
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, trim(adjustl(message)))

            collider%status = HIT_AND_RUN_PURE
            pl%status(impactors%id(:)) = ACTIVE
            pl%ldiscard(impactors%id(:)) = .false.
            pl%lcollision(impactors%id(:)) = .false.

            ! Be sure to save the pl so that snapshots still work 
            select type(before => collider%before)
            class is (swiftest_nbody_system)
            select type(after => collider%after)
            class is (swiftest_nbody_system)
               allocate(after%pl, source=before%pl) 
            end select
            end select
         end associate
      end select
      end select



      return
   end subroutine collision_generate_hitandrun


   module subroutine collision_generate_merge(self, nbody_system, param, t)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Merge massive bodies in any collisionals ystem.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routines symba_merge_pl.f90 and symba_discard_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routines symba5_merge.f and discard_mass_merge.f
      implicit none
      ! Arguments
      class(collision_merge),  intent(inout) :: self         !! Merge fragment system object 
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! The time of the collision
      ! Internals
      integer(I4B)                              :: i, j, k, ibiggest
      real(DP), dimension(NDIM)                 :: Lspin_new
      real(DP)                                  :: volume, dpe
      character(len=STRMAX) :: message

      select type(nbody_system)
      class is (swiftest_nbody_system)
         associate(impactors => nbody_system%collider%impactors, fragments => nbody_system%collider%fragments)
            message = "Merging"
            call collision_io_collider_message(nbody_system%pl, impactors%id, message)
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)

            select type(pl => nbody_system%pl)
            class is (swiftest_pl)
               ! Generate the merged body as a single fragment
               call self%setup_fragments(1)

               ! Calculate the initial energy of the nbody_system without the collisional family
               call self%get_energy_and_momentum(nbody_system, param, lbefore=.true.)

               ! The new body's metadata will be taken from the largest of the two impactor bodies, so we need 
               ! its index in the main pl structure
               ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
               fragments%id(1) = pl%id(ibiggest)
               allocate(fragments%info, source=pl%info(ibiggest:ibiggest))

               ! Compute the physical properties of the new body after the merge.
               volume = 4._DP / 3._DP * PI * sum(impactors%radius(:)**3)
               fragments%mass(1) = impactors%mass_dist(1)
               fragments%density(1) = fragments%mass(1) / volume
               fragments%radius(1) = (3._DP * volume / (4._DP * PI))**(THIRD)
               if (param%lrotation) then
                  do concurrent(i = 1:NDIM)
                     fragments%Ip(i,1) = sum(impactors%mass(:) * impactors%Ip(i,:)) 
                     Lspin_new(i) = sum(impactors%Lorbit(i,:) + impactors%Lorbit(i,:))
                  end do
                  fragments%Ip(:,1) = fragments%Ip(:,1) / fragments%mass(1)
                  fragments%rot(:,1) = Lspin_new(:) / (fragments%Ip(3,1) * fragments%mass(1) * fragments%radius(1)**2)
               else ! If spin is not enabled, we will consider the lost pre-collision angular momentum as "escaped" and add it to our bookkeeping variable
                  nbody_system%Lescape(:) = nbody_system%Lescape(:) + impactors%Lorbit(:,1) + impactors%Lorbit(:,2) 
               end if

               ! The fragment trajectory will be the barycentric trajectory
               fragments%rb(:,1) = impactors%rbcom(:)
               fragments%vb(:,1) = impactors%vbcom(:)

               ! Get the energy of the system after the collision
               call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)

               ! Keep track of the component of potential energy that is now not considered because two bodies became one
               dpe = self%pe(2) - self%pe(1)  
               nbody_system%Ecollisions = nbody_system%Ecollisions - dpe 
               nbody_system%Euntracked  = nbody_system%Euntracked + dpe 

               ! Update any encounter lists that have the removed bodies in them so that they instead point to the new body
               do k = 1, nbody_system%plpl_encounter%nenc
                  do j = 1, impactors%ncoll
                     i = impactors%id(j)
                     if (i == ibiggest) cycle
                     if (nbody_system%plpl_encounter%id1(k) == pl%id(i)) then
                        nbody_system%plpl_encounter%id1(k) = pl%id(ibiggest)
                        nbody_system%plpl_encounter%index1(k) = i
                     end if
                     if (nbody_system%plpl_encounter%id2(k) == pl%id(i)) then
                        nbody_system%plpl_encounter%id2(k) = pl%id(ibiggest)
                        nbody_system%plpl_encounter%index2(k) = i
                     end if
                     if (nbody_system%plpl_encounter%id1(k) == nbody_system%plpl_encounter%id2(k)) nbody_system%plpl_encounter%status(k) = INACTIVE
                  end do
               end do

               self%status = MERGED
               
               call collision_resolve_mergeaddsub(nbody_system, param, t, self%status) 

            end select
         end associate
      end select
      return 
   end subroutine collision_generate_merge


   module subroutine collision_generate_simple_disruption(self, nbody_system, param, t)
      implicit none
      class(collision_simple_disruption),  intent(inout) :: self         !! Simple fragment system object 
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! The time of the collision
   end subroutine collision_generate_simple_disruption

end submodule s_collision_generate