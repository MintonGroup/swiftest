!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (collision_classes) s_collision_util
   use swiftest
contains

   module subroutine collision_util_dealloc_fragments(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_fragments),  intent(inout) :: self

      call util_dealloc_pl(self)

      if (allocated(self%rc)) deallocate(self%rc) 
      if (allocated(self%vc)) deallocate(self%vc)
      if (allocated(self%rmag)) deallocate(self%rmag)
      if (allocated(self%rotmag)) deallocate(self%rotmag)
      if (allocated(self%v_r_unit)) deallocate(self%v_r_unit)
      if (allocated(self%v_t_unit)) deallocate(self%v_t_unit)
      if (allocated(self%v_n_unit)) deallocate(self%v_n_unit)

      return
   end subroutine collision_util_dealloc_fragments

   module subroutine collision_util_final_impactors(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(collision_impactors),  intent(inout) :: self !! Collision impactors storage object

      call self%reset()

      return
   end subroutine collision_util_final_impactors


   module subroutine collision_util_final_snapshot(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(collision_snapshot),  intent(inout) :: self !! Fraggle encountar storage object

      call encounter_util_final_snapshot(self%encounter_snapshot)

      return
   end subroutine collision_util_final_snapshot


   module subroutine collision_util_final_storage(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(collision_storage(*)),  intent(inout) :: self !! Collision storage object

      call util_final_storage(self%swiftest_storage)

      return
   end subroutine collision_util_final_storage


   module subroutine collision_util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(collision_system),  intent(inout) :: self !!  Collision system object

      call self%reset()

      return
   end subroutine collision_util_final_system


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

      npl_before = 0; ntp_before = 0; npl_after = 0; ntp_after = 0
      if (allocated(self%collision_system%before%pl)) then
         npl_before = self%collision_system%before%pl%nbody
      endif

      if (allocated(self%collision_system%before%tp)) then
         ntp_before = self%collision_system%before%tp%nbody
      end if 

      if (allocated(self%collision_system%after%pl)) then
         npl_after = self%collision_system%after%pl%nbody
      end if

      if (allocated(self%collision_system%after%tp)) then
         ntp_after = self%collision_system%after%tp%nbody
      end if 

      ntot = npl_before + ntp_before + npl_after + ntp_after
      if (ntot == 0) return
      allocate(idvals(ntot))

      nlo = 1; nhi = npl_before
      if (npl_before > 0) idvals(nlo:nhi) = self%collision_system%before%pl%id(1:npl_before)
      nlo = nhi + 1; nhi = nhi + ntp_before
      if (ntp_before > 0) idvals(nlo:nhi) = self%collision_system%before%tp%id(1:ntp_before)

      nlo = nhi + 1; nhi = nhi + npl_after
      if (npl_after > 0) idvals(nlo:nhi) = self%collision_system%after%pl%id(1:npl_after)
      nlo = nhi + 1; nhi = nhi + ntp_after
      if (ntp_after > 0) idvals(nlo:nhi) = self%collision_system%after%tp%id(1:ntp_after)

      return

   end subroutine collision_util_get_idvalues_snapshot


   module subroutine collision_util_get_energy_momentum(self,  system, param, lbefore)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision outcome state (lbefore = .false.)
      !! This subrourtine works by building a temporary internal massive body object out of the non-excluded bodies and optionally with fragments appended. 
      !! This will get passed to the energy calculation subroutine so that energy is computed exactly the same way is it is in the main program. 
      !! This will temporarily expand the massive body object in a temporary system object called tmpsys to feed it into symba_energy
      implicit none
      ! Arguments
      class(collision_system),      intent(inout) :: self    !! Encounter collision system object
      class(swiftest_nbody_system), intent(inout) :: system  !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param   !! Current swiftest run configuration parameters
      logical,                      intent(in)    :: lbefore !! Flag indicating that this the "before" state of the system, with impactors included and fragments excluded or vice versa
      ! Internals
      class(swiftest_nbody_system), allocatable, save :: tmpsys
      class(swiftest_parameters), allocatable, save   :: tmpparam
      integer(I4B)  :: npl_before, npl_after, stage

      associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody, pl => system%pl, cb => system%cb)

         ! Because we're making a copy of the massive body object with the excludes/fragments appended, we need to deallocate the
         ! big k_plpl array and recreate it when we're done, otherwise we run the risk of blowing up the memory by
         ! allocating two of these ginormous arrays simulteouously. This is not particularly efficient, but as this
         ! subroutine should be called relatively infrequently, it shouldn't matter too much.

         npl_before = pl%nbody
         npl_after = npl_before + nfrag

         if (lbefore) then
            call encounter_util_construct_temporary_system(fragments, system, param, tmpsys, tmpparam)
            ! Build the exluded body logical mask for the *before* case: Only the original bodies are used to compute energy and momentum
            tmpsys%pl%status(impactors%idx(1:impactors%ncoll)) = ACTIVE
            tmpsys%pl%status(npl_before+1:npl_after) = INACTIVE
         else
            if (.not.allocated(tmpsys)) then
               write(*,*) "Error in collision_util_get_energy_momentum. " // &
                         " This must be called with lbefore=.true. at least once before calling it with lbefore=.false."
               call util_exit(FAILURE)
            end if
            ! Build the exluded body logical mask for the *after* case: Only the new bodies are used to compute energy and momentum
            call encounter_util_add_fragments_to_system(fragments, impactors, tmpsys, tmpparam)
            tmpsys%pl%status(impactors%idx(1:impactors%ncoll)) = INACTIVE
            tmpsys%pl%status(npl_before+1:npl_after) = ACTIVE
         end if 

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
      end associate

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

      call encounter_util_get_vals_storage(self, idvals, tvals)

      ! Consolidate ids to only unique values
      call util_unique(idvals,self%idvals,self%idmap)
      self%nid = size(self%idvals)

      ! Don't consolidate time values (multiple collisions can happen in a single time step)
      self%nt = size(self%tvals)

      return
   end subroutine collision_util_index_map

   !> The following interfaces are placeholders intended to satisfy the required abstract methods given by the parent class
   module subroutine collision_util_placeholder_accel(self, system, param, t, lbeg)
      implicit none
      class(collision_fragments),     intent(inout) :: self   !! Fraggle fragment system object 
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      write(*,*) "The type-bound procedure 'accel' is not defined for the collision_fragments class"
      return
   end subroutine collision_util_placeholder_accel

   module subroutine collision_util_placeholder_kick(self, system, param, t, dt, lbeg)
      implicit none
      class(collision_fragments),     intent(inout) :: self   !! Fraggle fragment system object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system objec
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 

      write(*,*) "The type-bound procedure 'kick' is not defined for the collision_fragments class"
      return
   end subroutine collision_util_placeholder_kick

   module subroutine collision_util_placeholder_step(self, system, param, t, dt)
      implicit none
      class(collision_fragments),     intent(inout) :: self   !! Swiftest body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Simulation time
      real(DP),                     intent(in)    :: dt     !! Current stepsize

      write(*,*) "The type-bound procedure 'step' is not defined for the collision_fragments class"
      return
   end subroutine collision_util_placeholder_step


   module subroutine collision_util_reset_impactors(self)
      !! author: David A. Minton
      !!
      !! Resets the collider object variables to 0 and deallocates the index and mass distributions
      implicit none
      ! Arguments
      class(collision_impactors),  intent(inout) :: self

      if (allocated(self%idx)) deallocate(self%idx)
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


   module subroutine collision_util_reset_system(self)
      !! author: David A. Minton
      !!
      !! Resets the collider system and deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_system),  intent(inout) :: self

      if (allocated(self%impactors)) deallocate(self%impactors)
      if (allocated(self%fragments)) deallocate(self%fragments)
      if (allocated(self%before)) deallocate(self%before)
      if (allocated(self%after)) deallocate(self%after)

      self%Lorbit(:,:) = 0.0_DP
      self%Lspin(:,:) = 0.0_DP
      self%Ltot(:,:) = 0.0_DP
      self%ke_orbit(:) = 0.0_DP
      self%ke_spin(:) = 0.0_DP
      self%pe(:) = 0.0_DP
      self%Etot(:) = 0.0_DP

      return
   end subroutine collision_util_reset_system




   module subroutine collision_util_set_coordinate_system(self)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate system, including the unit vectors of both the system and individual fragments.
      implicit none
      ! Arguments
      class(collision_system),    intent(inout) :: self      !! Collisional system
      ! Internals
      integer(I4B) :: i
      real(DP), dimension(NDIM) ::  delta_r, delta_v, Ltot
      real(DP)   ::  L_mag
      real(DP), dimension(NDIM, self%fragments%nbody) :: L_sigma

      associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody)
         delta_v(:) = impactors%vb(:, 2) - impactors%vb(:, 1)
         delta_r(:) = impactors%rb(:, 2) - impactors%rb(:, 1)
   
         ! We will initialize fragments on a plane defined by the pre-impact system, with the z-axis aligned with the angular momentum vector
         ! and the y-axis aligned with the pre-impact distance vector.

         ! y-axis is the separation distance
         fragments%y_coll_unit(:) = .unit.delta_r(:) 
         Ltot = impactors%Lorbit(:,1) + impactors%Lorbit(:,2) + impactors%Lspin(:,1) + impactors%Lspin(:,2)

         L_mag = .mag.Ltot(:)
         if (L_mag > sqrt(tiny(L_mag))) then
            fragments%z_coll_unit(:) = .unit.Ltot(:) 
         else ! Not enough angular momentum to determine a z-axis direction. We'll just pick a random direction
            call random_number(fragments%z_coll_unit(:))
            fragments%z_coll_unit(:) = .unit.fragments%z_coll_unit(:) 
         end if

         ! The cross product of the y- by z-axis will give us the x-axis
         fragments%x_coll_unit(:) = fragments%y_coll_unit(:) .cross. fragments%z_coll_unit(:)

         fragments%v_coll_unit(:) = .unit.delta_v(:)

         if (.not.any(fragments%r_coll(:,:) > 0.0_DP)) return
         fragments%rmag(:) = .mag. fragments%r_coll(:,:)
  
         ! Randomize the tangential velocity direction. 
         ! This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector, otherwise we can get an ill-conditioned system
         call random_number(L_sigma(:,:)) 
         do concurrent(i = 1:nfrag, fragments%rmag(i) > 0.0_DP)
            fragments%v_n_unit(:, i) = fragments%z_coll_unit(:) + 2e-1_DP * (L_sigma(:,i) - 0.5_DP)
         end do

         ! Define the radial, normal, and tangential unit vectors for each individual fragment
         fragments%v_r_unit(:,:) = .unit. fragments%r_coll(:,:) 
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
      type(collision_storage(*)), allocatable, intent(inout) :: collision_history  !! Collision history object
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


end submodule s_collision_util