!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (collision) s_collision_util
   use swiftest
contains

   module subroutine collision_util_add_fragments_to_system(self, nbody_system, param)
      !! Author: David A. Minton
      !!
      !! Adds fragments to the temporary system pl object
      implicit none
      ! Arguments
      class(collision_merge),      intent(in)    :: self      !! Collision system system object
      class(base_nbody_system), intent(inout) :: nbody_system    !! Swiftest nbody system object
      class(base_parameters),   intent(in)    :: param     !! Current swiftest run configuration parameters
      ! Internals
      integer(I4B) :: i, npl_before, npl_after
      logical, dimension(:), allocatable :: lexclude

      select type(nbody_system)
      class is (swiftest_nbody_system)
         associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody, pl => nbody_system%pl, cb => nbody_system%cb)
            npl_after = pl%nbody
            npl_before = npl_after - nfrag
            allocate(lexclude(npl_after))

            pl%status(npl_before+1:npl_after) = ACTIVE
            pl%mass(npl_before+1:npl_after) = fragments%mass(1:nfrag)
            pl%Gmass(npl_before+1:npl_after) = fragments%mass(1:nfrag) * param%GU
            pl%radius(npl_before+1:npl_after) = fragments%radius(1:nfrag)
            do concurrent (i = 1:nfrag)
               pl%rb(:,npl_before+i) =  fragments%rb(:,i) 
               pl%vb(:,npl_before+i) =  fragments%vb(:,i) 
               pl%rh(:,npl_before+i) =  fragments%rb(:,i) - cb%rb(:)
               pl%vh(:,npl_before+i) =  fragments%vb(:,i) - cb%vb(:)
            end do
            if (param%lrotation) then
               pl%Ip(:,npl_before+1:npl_after) = fragments%Ip(:,1:nfrag)
               pl%rot(:,npl_before+1:npl_after) = fragments%rot(:,1:nfrag)
            end if
            ! This will remove the impactors from the system since we've replaced them with fragments
            lexclude(1:npl_after) = .false.
            lexclude(impactors%id(1:impactors%ncoll)) = .true.
            where(lexclude(1:npl_after)) 
               pl%status(1:npl_after) = INACTIVE
            elsewhere
               pl%status(1:npl_after) = ACTIVE
            endwhere

         end associate
      end select

      return
   end subroutine collision_util_add_fragments_to_system


   module subroutine collision_util_construct_temporary_system(self, nbody_system, param, tmpsys, tmpparam)
      !! Author: David A. Minton
      !!
      !! Constructs a temporary internal system consisting of active bodies and additional fragments. This internal temporary system is used to calculate system energy with and without fragments
      implicit none
      ! Arguments
      class(collision_merge),                intent(inout) :: self         !! Fraggle collision system object
      class(base_nbody_system),               intent(in)    :: nbody_system !! Original swiftest nbody system object
      class(base_parameters),                 intent(in)    :: param        !! Current swiftest run configuration parameters
      class(base_nbody_system), allocatable,  intent(out)   :: tmpsys       !! Output temporary swiftest nbody system object
      class(base_parameters),   allocatable,  intent(out)   :: tmpparam     !! Output temporary configuration run parameters
      ! Internals
      logical, dimension(:), allocatable :: linclude
      integer(I4B) :: npl_tot
      ! The following are needed in order to deal with typing requirements
      class(swiftest_nbody_system), allocatable :: tmpsys_local
      class(swiftest_parameters), allocatable :: tmpparam_local

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         associate(fragments => self%fragments, nfrag => self%fragments%nbody, pl => nbody_system%pl, npl => nbody_system%pl%nbody, cb => nbody_system%cb)
            ! Set up a new system based on the original
            if (allocated(tmpparam)) deallocate(tmpparam)
            if (allocated(tmpsys)) deallocate(tmpsys)
            allocate(tmpparam, source=param)
            call swiftest_setup_construct_system(tmpsys_local, tmpparam_local)

            ! No test particles necessary for energy/momentum calcs
            call tmpsys_local%tp%setup(0, tmpparam_local)

            ! Replace the empty central body object with a copy of the original
            deallocate(tmpsys_local%cb)
            allocate(tmpsys_local%cb, source=cb)

            ! Make space for the fragments
            npl_tot = npl + nfrag
            call tmpsys_local%pl%setup(npl_tot, tmpparam_local)
            allocate(linclude(npl_tot))

            ! Fill up the temporary system with all of the original bodies, leaving the spaces for fragments empty until we add them in later
            linclude(1:npl) = .true.
            linclude(npl+1:npl_tot) = .false.
            call tmpsys_local%pl%fill(pl, linclude)

            call move_alloc(tmpsys_local, tmpsys) 
            call move_alloc(tmpparam_local, tmpparam) 

         end associate
      end select
      end select

      return
   end subroutine collision_util_construct_temporary_system


   module subroutine collision_util_get_idvalues_snapshot(self, idvals)
      !! author: David A. Minton
      !!
      !! Returns an array of all id values saved in this snapshot
      implicit none
      ! Arguments
      class(collision_snapshot),               intent(in)  :: self   !! Fraggle snapshot object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      ! Internals
      integer(I4B) :: npl_before, ntp_before, npl_after, ntp_after, ntot, nlo, nhi

      select type(before => self%collider%before)
      class is (swiftest_nbody_system)
      select type(after => self%collider%after)
      class is (swiftest_nbody_system)
         npl_before = 0; ntp_before = 0; npl_after = 0; ntp_after = 0
         if (allocated(before%pl)) then
            npl_before = before%pl%nbody
         endif

         if (allocated(before%tp)) then
            ntp_before = before%tp%nbody
         end if 

         if (allocated(after%pl)) then
            npl_after = after%pl%nbody
         end if

         if (allocated(after%tp)) then
            ntp_after = after%tp%nbody
         end if 

         ntot = npl_before + ntp_before + npl_after + ntp_after
         if (ntot == 0) return
         allocate(idvals(ntot))

         nlo = 1; nhi = npl_before
         if (npl_before > 0) idvals(nlo:nhi) = before%pl%id(1:npl_before)
         nlo = nhi + 1; nhi = nhi + ntp_before
         if (ntp_before > 0) idvals(nlo:nhi) = before%tp%id(1:ntp_before)

         nlo = nhi + 1; nhi = nhi + npl_after
         if (npl_after > 0) idvals(nlo:nhi) = after%pl%id(1:npl_after)
         nlo = nhi + 1; nhi = nhi + ntp_after
         if (ntp_after > 0) idvals(nlo:nhi) = after%tp%id(1:ntp_after)
      end select
      end select

      return

   end subroutine collision_util_get_idvalues_snapshot


   module subroutine collision_util_get_energy_momentum(self,  nbody_system, param, lbefore)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision outcome state (lbefore = .false.)
      !! This subrourtine works by building a temporary internal massive body object out of the non-excluded bodies and optionally with fragments appended. 
      !! This will get passed to the energy calculation subroutine so that energy is computed exactly the same way is it is in the main program. 
      !! This will temporarily expand the massive body object in a temporary system object called tmpsys to feed it into symba_energy
      implicit none
      ! Arguments
      class(collision_merge),  intent(inout) :: self    !! Encounter collision system object
      class(base_nbody_system), intent(inout) :: nbody_system  !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param   !! Current swiftest run configuration parameters
      logical,                  intent(in)    :: lbefore !! Flag indicating that this the "before" state of the nbody_system, with impactors included and fragments excluded or vice versa
      ! Internals
      class(base_nbody_system), allocatable, save :: tmpsys
      class(base_parameters), allocatable, save   :: tmpparam
      integer(I4B)  :: npl_before, npl_after, stage

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody, pl => nbody_system%pl, cb => nbody_system%cb)

            ! Because we're making a copy of the massive body object with the excludes/fragments appended, we need to deallocate the
            ! big k_plpl array and recreate it when we're done, otherwise we run the risk of blowing up the memory by
            ! allocating two of these ginormous arrays simulteouously. This is not particularly efficient, but as this
            ! subroutine should be called relatively infrequently, it shouldn't matter too much.

            npl_before = pl%nbody
            npl_after = npl_before + nfrag

            if (lbefore) then
               call self%construct_temporary_system(nbody_system, param, tmpsys, tmpparam)
               select type(tmpsys)
               class is (swiftest_nbody_system)
                  ! Build the exluded body logical mask for the *before* case: Only the original bodies are used to compute energy and momentum
                  tmpsys%pl%status(impactors%id(1:impactors%ncoll)) = ACTIVE
                  tmpsys%pl%status(npl_before+1:npl_after) = INACTIVE
               end select
            else
               if (.not.allocated(tmpsys)) then
                  write(*,*) "Error in collision_util_get_energy_momentum. " // &
                           " This must be called with lbefore=.true. at least once before calling it with lbefore=.false."
                  call util_exit(FAILURE)
               end if
               select type(tmpsys)
               class is (swiftest_nbody_system)
                  ! Build the exluded body logical mask for the *after* case: Only the new bodies are used to compute energy and momentum
                  call self%add_fragments(tmpsys, tmpparam)
                  tmpsys%pl%status(impactors%id(1:impactors%ncoll)) = INACTIVE
                  tmpsys%pl%status(npl_before+1:npl_after) = ACTIVE
               end select
            end if 
            select type(tmpsys)
            class is (swiftest_nbody_system)

               if (param%lflatten_interactions) call tmpsys%pl%flatten(param)

               call tmpsys%get_energy_and_momentum(param) 

               ! Calculate the current fragment energy and momentum balances
               if (lbefore) then
                  stage = 1
               else
                  stage = 2
               end if
               self%Lorbit(:,stage) = tmpsys%Lorbit(:)
               self%Lspin(:,stage) = tmpsys%Lspin(:)
               self%Ltot(:,stage) = tmpsys%Ltot(:)
               self%ke_orbit(stage) = tmpsys%ke_orbit
               self%ke_spin(stage) = tmpsys%ke_spin
               self%pe(stage) = tmpsys%pe
               self%Etot(stage) = tmpsys%te 
               if (stage == 2) self%Etot(stage) = self%Etot(stage) - (self%pe(2) - self%pe(1)) ! Gotta be careful with PE when number of bodies changes.
            end select
         end associate
      end select
      end select

      return
   end subroutine collision_util_get_energy_momentum


   module subroutine collision_util_index_map(self)
      !! author: David A. Minton
      !!
      !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      implicit none
      ! Arguments
      class(collision_storage(*)), intent(inout) :: self  !! Swiftest storage object
      ! Internals
      integer(I4B), dimension(:), allocatable :: idvals
      real(DP), dimension(:), allocatable :: tvals

      call self%get_index_values(idvals, tvals)

      ! Consolidate ids to only unique values
      call swiftest_util_unique(idvals,self%idvals,self%idmap)
      self%nid = size(self%idvals)

      ! Don't consolidate time values (multiple collisions can happen in a single time step)
      self%nt = size(self%tvals)

      return
   end subroutine collision_util_index_map


   module subroutine collision_util_reset_impactors(self)
      !! author: David A. Minton
      !!
      !! Resets the collider object variables to 0 and deallocates the index and mass distributions
      implicit none
      ! Arguments
      class(collision_impactors),  intent(inout) :: self

      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%mass_dist)) deallocate(self%mass_dist)
      self%ncoll = 0
      self%rb(:,:) = 0.0_DP
      self%vb(:,:) = 0.0_DP
      self%rot(:,:) = 0.0_DP
      self%Lspin(:,:) = 0.0_DP
      self%Lorbit(:,:) = 0.0_DP
      self%Ip(:,:) = 0.0_DP
      self%mass(:) = 0.0_DP
      self%radius(:) = 0.0_DP
      self%Qloss = 0.0_DP
      self%regime = 0

      self%x_unit(:) = 0.0_DP
      self%y_unit(:) = 0.0_DP
      self%z_unit(:) = 0.0_DP
      self%v_unit(:) = 0.0_DP
      self%rbcom(:) = 0.0_DP
      self%vbcom(:) = 0.0_DP
      self%rbimp(:) = 0.0_DP

      return
   end subroutine collision_util_reset_impactors


   module subroutine collision_util_reset_fragments(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_fragments(*)),  intent(inout) :: self

      if (allocated(self%info)) deallocate(self%info) 
      self%mtot = 0.0_DP
      self%status = 0
      self%rh(:,:) = 0.0_DP
      self%vh(:,:) = 0.0_DP
      self%rb(:,:) = 0.0_DP
      self%vb(:,:) = 0.0_DP
      self%rot(:,:) = 0.0_DP
      self%Ip(:,:) = 0.0_DP
      self%mass(:) = 0.0_DP
      self%radius(:) = 0.0_DP
      self%density(:) = 0.0_DP
      self%rc(:,:) = 0.0_DP
      self%vc(:,:) = 0.0_DP
      self%v_r_unit(:,:) = 0.0_DP
      self%v_t_unit(:,:) = 0.0_DP
      self%v_n_unit(:,:) = 0.0_DP

      return
   end subroutine collision_util_reset_fragments


   module subroutine collision_util_reset_system(self)
      !! author: David A. Minton
      !!
      !! Resets the collider nbody_system and deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_merge),    intent(inout) :: self  !! Collision system object

      if (allocated(self%before)) deallocate(self%before)
      if (allocated(self%after)) deallocate(self%after)

      self%Lorbit(:,:) = 0.0_DP
      self%Lspin(:,:) = 0.0_DP
      self%Ltot(:,:) = 0.0_DP
      self%ke_orbit(:) = 0.0_DP
      self%ke_spin(:) = 0.0_DP
      self%pe(:) = 0.0_DP
      self%Etot(:) = 0.0_DP

      if (allocated(self%impactors)) call self%impactors%reset()
      if (allocated(self%fragments)) deallocate(self%fragments)

      return
   end subroutine collision_util_reset_system


   module subroutine collision_util_set_coordinate_system(self)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual fragments.
      implicit none
      ! Arguments
      class(collision_merge),    intent(inout) :: self      !! Collisional nbody_system
      ! Internals
      integer(I4B) :: i
      real(DP), dimension(NDIM) ::  delta_r, delta_v, Ltot
      real(DP)   ::  L_mag
      real(DP), dimension(NDIM, self%fragments%nbody) :: L_sigma

      associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody)
         delta_v(:) = impactors%vb(:, 2) - impactors%vb(:, 1)
         delta_r(:) = impactors%rb(:, 2) - impactors%rb(:, 1)
   
         ! We will initialize fragments on a plane defined by the pre-impact nbody_system, with the z-axis aligned with the angular momentum vector
         ! and the y-axis aligned with the pre-impact distance vector.

         ! y-axis is the separation distance
         impactors%y_unit(:) = .unit.delta_r(:) 
         Ltot = impactors%Lorbit(:,1) + impactors%Lorbit(:,2) + impactors%Lspin(:,1) + impactors%Lspin(:,2)

         L_mag = .mag.Ltot(:)
         if (L_mag > sqrt(tiny(L_mag))) then
            impactors%z_unit(:) = .unit.Ltot(:) 
         else ! Not enough angular momentum to determine a z-axis direction. We'll just pick a random direction
            call random_number(impactors%z_unit(:))
            impactors%z_unit(:) = .unit.impactors%z_unit(:) 
         end if

         ! The cross product of the y- by z-axis will give us the x-axis
         impactors%x_unit(:) = impactors%y_unit(:) .cross. impactors%z_unit(:)

         impactors%v_unit(:) = .unit.delta_v(:)

         if (.not.any(fragments%rc(:,:) > 0.0_DP)) return
         fragments%rmag(:) = .mag. fragments%rc(:,:)
  
         ! Randomize the tangential velocity direction. 
         ! This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector, otherwise we can get an ill-conditioned nbody_system
         call random_number(L_sigma(:,:)) 
         do concurrent(i = 1:nfrag, fragments%rmag(i) > 0.0_DP)
            fragments%v_n_unit(:, i) = impactors%z_unit(:) + 2e-1_DP * (L_sigma(:,i) - 0.5_DP)
         end do

         ! Define the radial, normal, and tangential unit vectors for each individual fragment
         fragments%v_r_unit(:,:) = .unit. fragments%rc(:,:) 
         fragments%v_n_unit(:,:) = .unit. fragments%v_n_unit(:,:) 
         fragments%v_t_unit(:,:) = .unit. (fragments%v_n_unit(:,:) .cross. fragments%v_r_unit(:,:))

      end associate

      return
   end subroutine collision_util_set_coordinate_system


   subroutine collision_util_save_snapshot(collision_history, snapshot)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter storage against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
      !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
      !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
      implicit none
      ! Arguments
      class(collision_storage(*)), allocatable, intent(inout) :: collision_history  !! Collision history object
      class(encounter_snapshot),               intent(in)    :: snapshot           !! Encounter snapshot object
      ! Internals
      type(collision_storage(nframes=:)), allocatable :: tmp
      integer(I4B) :: i, nnew, nold, nbig

      ! Advance the snapshot frame counter
      collision_history%iframe = collision_history%iframe + 1

      ! Check to make sure the current encounter_history object is big enough. If not, grow it by a factor of 2
      nnew = collision_history%iframe
      nold = collision_history%nframes

      if (nnew > nold) then
         nbig = nold
         do while (nbig < nnew)
            nbig = nbig * 2
         end do
         allocate(collision_storage(nbig) :: tmp) 
         tmp%iframe = collision_history%iframe
         call move_alloc(collision_history%nc, tmp%nc)

         do i = 1, nold
            if (allocated(collision_history%frame(i)%item)) call move_alloc(collision_history%frame(i)%item, tmp%frame(i)%item)
         end do
         deallocate(collision_history)
         call move_alloc(tmp,collision_history)
         nnew = nbig
      end if

      collision_history%frame(nnew) = snapshot

      return
   end subroutine collision_util_save_snapshot


   module subroutine collision_util_snapshot(self, param, nbody_system, t, arg)
      !! author: David A. Minton
      !!
      !! Takes a minimal snapshot of the state of the nbody_system during an encounter so that the trajectories
      !! can be played back through the encounter
      implicit none
      ! Internals
      class(collision_storage(*)), intent(inout)          :: self   !! Swiftest storage object
      class(base_parameters),      intent(inout)          :: param  !! Current run configuration parameters
      class(base_nbody_system),    intent(inout)          :: nbody_system !! Swiftest nbody system object to store
      real(DP),                    intent(in),   optional :: t      !! Time of snapshot if different from nbody_system time
      character(*),                intent(in),   optional :: arg    !! "before": takes a snapshot just before the collision. "after" takes the snapshot just after the collision.
      ! Arguments
      class(collision_snapshot), allocatable :: snapshot
      class(swiftest_pl), allocatable                    :: pl
      character(len=:), allocatable :: stage

      if (present(arg)) then
         stage = arg
      else
         stage = ""
      end if 

      select type (nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         select case(stage)
         case("before")
            ! Saves the states of the bodies involved in the collision before the collision is resolved
            associate (idx => nbody_system%collider%impactors%id, ncoll => nbody_system%collider%impactors%ncoll)
               allocate(pl, mold=nbody_system%pl)
               call pl%setup(ncoll, param)
               pl%id(:) = nbody_system%pl%id(idx(:))
               pl%Gmass(:) = nbody_system%pl%Gmass(idx(:))
               pl%radius(:) = nbody_system%pl%radius(idx(:))
               pl%rot(:,:) = nbody_system%pl%rot(:,idx(:))
               pl%Ip(:,:) = nbody_system%pl%Ip(:,idx(:))
               pl%rh(:,:) = nbody_system%pl%rh(:,idx(:))
               pl%vh(:,:) = nbody_system%pl%vh(:,idx(:))
               pl%info(:) = nbody_system%pl%info(idx(:))
               select type (before => nbody_system%collider%before)
               class is (swiftest_nbody_system)
                  allocate(before%pl, source=pl)
               end select
            end associate
         case("after")
            allocate(collision_snapshot :: snapshot)
            allocate(snapshot%collider, source=nbody_system%collider) 
            snapshot%t = t
            call collision_util_save_snapshot(nbody_system%collision_history,snapshot)
         case default
            write(*,*) "collision_util_snapshot requies either 'before' or 'after' passed to 'arg'"
         end select
      end select
      end select

      return
   end subroutine collision_util_snapshot


end submodule s_collision_util