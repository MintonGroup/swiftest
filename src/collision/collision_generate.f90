
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

   module subroutine collision_generate_basic(self, nbody_system, param, t)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Merge massive bodies no matter the regime
      !! 
      !! Adapted from David E. Kaufmann's Swifter routines symba_merge_pl.f90 and symba_discard_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routines symba5_merge.f and discard_mass_merge.f
      implicit none
      ! Arguments
      class(collision_basic),   intent(inout) :: self         !! Merge fragment system object 
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! The time of the collision

      call self%merge(nbody_system, param, t)

      return
   end subroutine collision_generate_basic


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
               call self%hitandrun(nbody_system, param, t)
            case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
               call self%merge(nbody_system, param, t) ! Use the default collision model, which is merge
            case default 
               write(*,*) "Error in swiftest_collision, unrecognized collision regime"
               call util_exit(FAILURE)
            end select
            end associate
      end select
      end select

      return
   end subroutine collision_generate_bounce


   module subroutine collision_generate_hitandrun(self, nbody_system, param, t) 
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(collision_basic),   intent(inout) :: self         !! Collision system object
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
      real(DP),                 intent(in)    :: t            !! Time of collision
      ! Result
      integer(I4B)                            :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, ibiggest, nfrag, jtarg, jproj
      logical                                 :: lpure 
      character(len=STRMAX) :: message
      real(DP) :: dpe


      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(impactors => self%impactors)
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

            ! The simple disruption model (and its extended types allow for non-pure hit and run. 
            !For the basic merge model, all hit and runs are pure
            select type(self) 
            class is (collision_disruption)
               if (impactors%mass_dist(2) > 0.9_DP * impactors%mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Pure hit and run. No new fragments generated.")
                  nfrag = 0
                  lpure = .true.
               else ! Imperfect hit and run, so we'll keep the largest body and destroy the other
                  lpure = .false.
                  call self%set_mass_dist(param)

                  ! Generate the position and velocity distributions of the fragments
                  call self%get_energy_and_momentum(nbody_system, param, lbefore=.true.)
                  call self%disrupt(nbody_system, param, t, lpure)
                  call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)

                  dpe = self%pe(2) - self%pe(1) 
                  nbody_system%Ecollisions = nbody_system%Ecollisions - dpe 
                  nbody_system%Euntracked  = nbody_system%Euntracked + dpe 

                  if (lpure) then
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Should have been a pure hit and run instead")
                     nfrag = 0
                  else
                     nfrag = self%fragments%nbody
                     write(message, *) nfrag
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")
                  end if
               end if
            class default
               lpure = .true.
            end select 

            if (lpure) then ! Reset these bodies back to being active so that nothing further is done to them
               status = HIT_AND_RUN_PURE
               pl%status(impactors%id(:)) = ACTIVE
               pl%ldiscard(impactors%id(:)) = .false.
               pl%lcollision(impactors%id(:)) = .false.
               ! Be sure to save the pl so that snapshots still work 
               select type(before => self%before)
               class is (swiftest_nbody_system)
               select type(after => self%after)
               class is (swiftest_nbody_system)
                  allocate(after%pl, source=before%pl) 
               end select
               end select  
            else
               ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
               self%fragments%id(1) = pl%id(ibiggest)
               self%fragments%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
               param%maxid = self%fragments%id(nfrag)
               status = HIT_AND_RUN_DISRUPT
               call collision_resolve_mergeaddsub(nbody_system, param, t, status)
            end if
         end associate
      end select
      end select

      return
   end subroutine collision_generate_hitandrun


   module subroutine collision_generate_simple(self, nbody_system, param, t)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic disruption collision
      !! 
      implicit none
      ! Arguments
      class(collision_disruption), intent(inout) :: self
      class(base_nbody_system),           intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),             intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
      real(DP),                           intent(in)    :: t            !! Time of collision
      ! Internals
      integer(I4B)          :: i, ibiggest, nfrag
      character(len=STRMAX) :: message 
      real(DP) :: dpe

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(impactors => self%impactors, status => self%status)

            select case (impactors%regime) 
            case (COLLRESOLVE_REGIME_HIT_AND_RUN)
               call self%hitandrun(nbody_system, param, t)
               return
            case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
               call self%merge(nbody_system, param, t) ! Use the default collision model, which is merge
               return
            case(COLLRESOLVE_REGIME_DISRUPTION)
               message = "Disruption between"
            case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
               message = "Supercatastrophic disruption between"
            case default 
               write(*,*) "Error in swiftest_collision, unrecognized collision regime"
               call util_exit(FAILURE)
            end select
            call self%set_mass_dist(param) 
            call self%get_energy_and_momentum(nbody_system, param, lbefore=.true.)
            call self%set_budgets()
            call self%set_coordinate_system()
            call self%disrupt(nbody_system, param, t)
            call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)

            dpe = self%pe(2) - self%pe(1) 
            nbody_system%Ecollisions = nbody_system%Ecollisions - dpe 
            nbody_system%Euntracked  = nbody_system%Euntracked + dpe 

            associate (fragments => self%fragments)
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
            end associate
         end associate
      end select
      end select
      return
   end subroutine collision_generate_simple


   module subroutine collision_generate_merge(self, nbody_system, param, t)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Merge massive bodies in any collisional system.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routines symba_merge_pl.f90 and symba_discard_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routines symba5_merge.f and discard_mass_merge.f
      implicit none
      ! Arguments
      class(collision_basic),   intent(inout) :: self         !! Merge fragment system object 
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
         associate(impactors => nbody_system%collider%impactors)
            message = "Merging"
            call collision_io_collider_message(nbody_system%pl, impactors%id, message)
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)

            select type(pl => nbody_system%pl)
            class is (swiftest_pl)
               ! Get coordinate system
               call impactors%set_coordinate_system()

               ! Generate the merged body as a single fragment
               call self%setup_fragments(1)
               associate(fragments => nbody_system%collider%fragments)

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
                  fragments%rc(:,1) = 0.0_DP
                  fragments%vc(:,1) = 0.0_DP

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

               end associate
            end select
         end associate
      end select
      return 
   end subroutine collision_generate_merge


   module subroutine collision_generate_disrupt(self, nbody_system, param, t, lfailure)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Generates a simple fragment position and velocity distribution based on the collision 
      !! regime. 
      implicit none
      ! Arguments
      class(collision_disruption),  intent(inout) :: self         !! Simple fragment system object 
      class(base_nbody_system),            intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),              intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                            intent(in)    :: t            !! The time of the collision
      logical, optional,                   intent(out)   :: lfailure
      ! Internals

      call collision_generate_simple_pos_vec(self)
      call collision_generate_simple_rot_vec(self)
      call collision_generate_simple_vel_vec(self)

      return
   end subroutine collision_generate_disrupt


   module subroutine collision_generate_simple_pos_vec(collider)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the position vectors of the fragments around the center of mass based on the collision style.
      !! For hit and run with disruption, the fragments are generated in a random cloud around the smallest of the two colliders (body 2)
      !! For disruptive collisions, the fragments are generated in a random cloud around the impact point. Bodies are checked for overlap and
      !! regenerated if they overlap.
      implicit none
      ! Arguments
      class(collision_disruption), intent(inout) :: collider !! Fraggle collision system object
      ! Internals
      real(DP)  :: dis
      real(DP), dimension(NDIM,2) :: fragment_cloud_center
      real(DP), dimension(2) :: fragment_cloud_radius
      logical, dimension(collider%fragments%nbody) :: loverlap
      integer(I4B) :: i, j, loop
      logical :: lcat, lhitandrun
      integer(I4B), parameter :: MAXLOOP = 10000
      real(DP) :: rdistance
      real(DP), parameter :: fail_scale = 1.1_DP ! Scale factor to apply to cloud radius and distance if cloud generation fails


      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)
         lcat = (impactors%regime == COLLRESOLVE_REGIME_SUPERCATASTROPHIC) 
         lhitandrun = (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) 

         ! We will treat the first two fragments of the list as special cases. 
         ! Place the first two bodies at the centers of the two fragment clouds, but be sure they are sufficiently far apart to avoid overlap
         rdistance = .mag. (impactors%rc(:,2) - impactors%rc(:,1)) - sum(fragments%radius(1:2))
         rdistance = min(0.5_DP*rdistance, 1e-6_DP*impactors%radius(2))

         fragment_cloud_radius(:) = impactors%radius(:)

         loverlap(:) = .true.
         do loop = 1, MAXLOOP
            if (.not.any(loverlap(:))) exit
            fragment_cloud_center(:,1) = impactors%rc(:,1) + rdistance * impactors%bounce_unit(:)
            fragment_cloud_center(:,2) = impactors%rc(:,2) - rdistance * impactors%bounce_unit(:)
            do concurrent(i = 1:nfrag, loverlap(i))
               if (i < 3) then
                  fragments%rc(:,i) = fragment_cloud_center(:,i)
               else
                  ! Make a random cloud
                  call random_number(fragments%rc(:,i))
   
                  ! Make the fragment cloud symmertic about 0
                  fragments%rc(:,i) = 2 *(fragments%rc(:,i) - 0.5_DP)
   
                  j = fragments%origin_body(i)
   
                  ! Scale the cloud size
                  fragments%rc(:,i) = fragment_cloud_radius(j) * fragments%rc(:,i)
   
                  ! Shift to the cloud center coordinates
                  fragments%rc(:,i) = fragments%rc(:,i) + fragment_cloud_center(:,j)
               end if
            end do

            ! Check for any overlapping bodies.
            loverlap(:) = .false.
            do j = 1, nfrag
               do i = j + 1, nfrag
                  dis = .mag.(fragments%rc(:,j) - fragments%rc(:,i))
                  loverlap(i) = loverlap(i) .or. (dis <= (fragments%radius(i) + fragments%radius(j))) 
               end do
            end do
            rdistance = rdistance * fail_scale
            fragment_cloud_radius(:) = fragment_cloud_radius(:) * fail_scale
         end do

         call collision_util_shift_vector_to_origin(fragments%mass, fragments%rc)
         call collider%set_coordinate_system()

         do concurrent(i = 1:nfrag)
            fragments%rb(:,i) = fragments%rc(:,i) + impactors%rbcom(:)
         end do

         impactors%rbcom(:) = 0.0_DP
         do concurrent(i = 1:nfrag)
            impactors%rbcom(:) = impactors%rbcom(:) + fragments%mass(i) * fragments%rb(:,i) 
         end do
         impactors%rbcom(:) = impactors%rbcom(:) / fragments%mtot
      end associate

      return
   end subroutine collision_generate_simple_pos_vec


   module subroutine collision_generate_simple_rot_vec(collider)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Calculates the spins of a collection of fragments such that they conserve angular momentum 
      implicit none
      ! Arguments
      class(collision_basic), intent(inout) :: collider !! Fraggle collision system object
      ! Internals
      real(DP), dimension(NDIM) :: Lbefore, Lafter, Lspin, rotdir
      real(DP) :: v_init, v_final, mass_init, mass_final, rotmag
      integer(I4B) :: i

      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)

         ! Torque the first body based on the change in angular momentum betwen the pre- and post-impact system assuming a simple bounce
         mass_init = impactors%mass(2)
         mass_final = sum(fragments%mass(2:nfrag))
         v_init = .mag.(impactors%vb(:,2) - impactors%vb(:,1))
         v_final = sqrt(mass_init / mass_final * v_init**2 - 2 * impactors%Qloss / mass_final)

         Lbefore(:) = mass_init * (impactors%rb(:,2) - impactors%rb(:,1)) .cross. (impactors%vb(:,2) - impactors%vb(:,1))
          
         Lafter(:) = mass_final * (impactors%rb(:,2) - impactors%rb(:,1)) .cross. (v_final * impactors%bounce_unit(:))
         Lspin(:) = impactors%Lspin(:,1) + (Lbefore(:) - Lafter(:))
         fragments%rot(:,1) = Lspin(:) / (fragments%mass(1) * fragments%radius(1)**2 * fragments%Ip(3,1))

         ! Add in some random spin noise. The magnitude will be scaled by the before-after amount and the direction will be random
         do concurrent(i = 2:nfrag)
            call random_number(rotdir)
            rotdir = rotdir - 0.5_DP
            rotdir = .unit. rotdir
            call random_number(rotmag)
            rotmag = rotmag * .mag. (Lbefore(:) - Lafter(:)) / ((nfrag - 1) * fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i))
            fragments%rot(:,i) = rotmag * rotdir
         end do

      end associate

      return
   end subroutine collision_generate_simple_rot_vec


   module subroutine collision_generate_simple_vel_vec(collider)
      !! Author:  David A. Minton
      !!
      !! Generates an initial velocity distribution. For disruptions, the velocity magnitude is set to be
      !! 2x the escape velocity of the colliding pair. For hit and runs the velocity magnitude is set to be
      !! 2x the escape velocity of the smallest of the two bodies.
      implicit none
      ! Arguments
      class(collision_disruption), intent(inout) :: collider !! Fraggle collision system object
      ! Internals
      integer(I4B) :: i,j
      logical :: lhitandrun, lnoncat
      real(DP), dimension(NDIM) :: vimp_unit, rimp, vrot, Lresidual, vshear
      real(DP), dimension(2) :: vimp
      real(DP) :: vmag, vdisp, Lmag
      integer(I4B), dimension(collider%fragments%nbody) :: vsign
      real(DP), dimension(collider%fragments%nbody) :: vscale, mass_vscale

      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)
         lhitandrun = (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) 
         lnoncat = (impactors%regime /= COLLRESOLVE_REGIME_SUPERCATASTROPHIC) ! For non-catastrophic impacts, make the fragments act like ejecta and point away from the impact point
         ! "Bounce" the first two bodies
         where(fragments%origin_body(:) == 1)
            vsign(:) = -1
         elsewhere
            vsign(:) = 1
         end where

         ! Compute the velocity dispersion based on the escape velocity
         if (lhitandrun) then
            vdisp = 2 * sqrt(2 * impactors%Gmass(2) / impactors%radius(2))
         else
            vdisp = 2 * sqrt(2 * sum(impactors%Gmass(:)) / sum(impactors%radius(:)))
         end if

         vimp(:) = .mag.impactors%vc(:,:)

         ! Scale the magnitude of the velocity by the distance from the impact point and add a bit of shear
         do concurrent(i = 1:nfrag)
            rimp(:) = fragments%rc(:,i) - impactors%rbimp(:) 
            vscale(i) = .mag. rimp(:) / (.mag. (impactors%rb(:,2) - impactors%rb(:,1)))
         end do
         vscale(:) = vscale(:)/maxval(vscale(:))

         ! Give the fragment velocities a random value that is somewhat scaled with fragment mass
         call random_number(mass_vscale)
         mass_vscale(:) = (mass_vscale(:) + 1.0_DP) / 2
         mass_vscale(:) = mass_vscale(:) * (fragments%mtot / fragments%mass(:))**(0.125_DP)
         mass_vscale(:) = 2*mass_vscale(:) / maxval(mass_vscale(:))
         fragments%vc(:,1) = .mag.impactors%vc(:,1) * impactors%bounce_unit(:) 
         do concurrent(i = 2:nfrag)
            j = fragments%origin_body(i)
            vrot(:) = impactors%rot(:,j) .cross. (fragments%rc(:,i) - impactors%rc(:,j))
            vmag = .mag.impactors%vc(:,j) * vscale(i) * mass_vscale(i)
            if (lhitandrun) then
               fragments%vc(:,i) = vmag * 0.5_DP * impactors%bounce_unit(:) * vsign(i) + vrot(:)
            else
               ! Add more velocity dispersion to disruptions vs hit and runs.
               rimp(:) = fragments%rc(:,i) - impactors%rbimp(:)
               vimp_unit(:) = .unit. rimp(:)
               fragments%vc(:,i) = vmag * 0.5_DP * (impactors%bounce_unit(:) + vimp_unit(:)) * vsign(i) + vrot(:)
            end if
         end do
         call collider%set_coordinate_system()
         ! Check for any residual angular momentum, and if there is any, put it into velocity shear of the fragments
         call fragments%get_angular_momentum()
         if (all(fragments%L_budget(:) / (fragments%Lorbit(:) + fragments%Lspin(:)) > epsilon(1.0_DP))) then
            Lresidual(:) = fragments%L_budget(:) - (fragments%Lorbit(:) + fragments%Lspin(:))
            do concurrent(i = 3:nfrag)
               vshear(:) = Lresidual(:) / (nfrag - 2) / (fragments%mass(i) * fragments%rc(:,i) .cross. fragments%v_t_unit(:,i))
               vshear(:) = .mag.vshear(:) * fragments%v_t_unit(:,i)
               fragments%vc(:,i) = fragments%vc(:,i) + vshear(:)
            end do
         end if

         do concurrent(i = 1:nfrag)
            fragments%vb(:,i) = fragments%vc(:,i) + impactors%vbcom(:)
         end do

         impactors%vbcom(:) = 0.0_DP
         do concurrent(i = 1:nfrag)
            impactors%vbcom(:) = impactors%vbcom(:) + fragments%mass(i) * fragments%vb(:,i) 
         end do
         impactors%vbcom(:) = impactors%vbcom(:) / fragments%mtot

         ! Distribute any remaining angular momentum into fragments pin
         call fragments%set_spins()

      end associate
      return
   end subroutine collision_generate_simple_vel_vec


end submodule s_collision_generate