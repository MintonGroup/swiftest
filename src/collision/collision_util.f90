!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (encounter_classes) s_encounter_util
   use swiftest
contains

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


   module subroutine collision_util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(collision_system),  intent(inout) :: self !! Collision impactors storage object

      call self%reset()

      return
   end subroutine collision_util_final_system


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
   end subroutine collision_util_final_storage


   module subroutine collision_util_get_idvalues_snapshot(self, idvals)
      !! author: David A. Minton
      !!
      !! Returns an array of all id values saved in this snapshot
      implicit none
      ! Arguments
      class(collision_snapshot),                 intent(in)  :: self   !! Fraggle snapshot object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      ! Internals
      integer(I4B) :: ncoll, nfrag

      if (allocated(self%impactors)) then
         ncoll = self%impactors%pl%nbody
      else
         ncoll = 0
      end if 

      if (allocated(self%fragments)) then
         nfrag = self%fragments%pl%nbody
      else
         nfrag = 0
      end if

      if (ncoll + nfrag == 0) return
      allocate(idvals(ncoll+nfrag))

      if (ncoll > 0) idvals(1:ncoll) = self%impactors%pl%id(:)
      if (nfrag > 0) idvals(ncoll+1:ncoll+nfrag) = self%fragments%pl%id(:)

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
      class(collision_system), intent(inout) :: self    !! Encounter collision system object
      class(swiftest_nbody_system),      intent(inout) :: system  !! Swiftest nbody system object
      class(swiftest_parameters),        intent(inout) :: param   !! Current swiftest run configuration parameters
      logical,                           intent(in)    :: lbefore !! Flag indicating that this the "before" state of the system, with impactors included and fragments excluded or vice versa
      ! Internals
      class(swiftest_nbody_system), allocatable, save :: tmpsys
      class(swiftest_parameters), allocatable, save   :: tmpparam
      integer(I4B)  :: npl_before, npl_after

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
            fragments%Lorbit_before(:) = tmpsys%Lorbit(:)
            fragments%Lspin_before(:) = tmpsys%Lspin(:)
            fragments%Ltot_before(:) = tmpsys%Ltot(:)
            fragments%ke_orbit_before = tmpsys%ke_orbit
            fragments%ke_spin_before = tmpsys%ke_spin
            fragments%pe_before = tmpsys%pe
            fragments%Etot_before = tmpsys%te 
         else
            fragments%Lorbit_after(:) = tmpsys%Lorbit(:)
            fragments%Lspin_after(:) = tmpsys%Lspin(:)
            fragments%Ltot_after(:) = tmpsys%Ltot(:)
            fragments%ke_orbit_after = tmpsys%ke_orbit
            fragments%ke_spin_after = tmpsys%ke_spin
            fragments%pe_after = tmpsys%pe
            fragments%Etot_after = tmpsys%te - (fragments%pe_after - fragments%pe_before) ! Gotta be careful with PE when number of bodies changes.
         end if
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


   module subroutine collision_util_reset_fragments(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_fragments),  intent(inout) :: self

      call self%swiftest_pl%dealloc()

      if (allocated(self%rc)) deallocate(self%rc) 
      if (allocated(self%vc)) deallocate(self%vc)
      if (allocated(self%rmag)) deallocate(self%rmag)
      if (allocated(self%rotmag)) deallocate(self%rotmag)
      if (allocated(self%v_r_unit)) deallocate(self%v_r_unit)
      if (allocated(self%v_t_unit)) deallocate(self%v_t_unit)
      if (allocated(self%v_n_unit)) deallocate(self%v_n_unit)

      return
   end subroutine collision_util_reset_fragments


   module subroutine collision_util_reset_impactors(self)
      !! author: David A. Minton
      !!
      !! Resets the collider object variables to 0 and deallocates the index and mass distributions
      implicit none
      ! Arguments
      class(collision_impactors),  intent(inout) :: self

      if (allocated(self%idx)) deallocate(self%idx)
      if (allocated(self%mass_dist)) deallocate(self%mass_dist)
      ncoll = 0
      rb(:,:) = 0.0_DP
      vb(:,:) = 0.0_DP
      rot(:,:) = 0.0_DP
      L_spin(:,:) = 0.0_DP
      L_orbit(:,:) = 0.0_DP
      Ip(:,:) = 0.0_DP
      mass(:) = 0.0_DP
      radius(:) = 0.0_DP
      Qloss = 0.0_DP
      regime = 0

      x_unit(:) = 0.0_DP
      y_unit(:) = 0.0_DP
      z_unit(:) = 0.0_DP
      v_unit(:) = 0.0_DP
      rbcom(:) = 0.0_DP
      vbcom(:) = 0.0_DP
      rbimp(:) = 0.0_DP

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

      Lorbit(:,:) = 0.0_DP
      Lspin(:,:) = 0.0_DP
      Ltot(:,:) = 0.0_DP
      ke_orbit(:) = 0.0_DP
      ke_spin(:) = 0.0_DP
      pe(:) = 0.0_DP
      Etot(:) = 0.0_DP

      return
   end subroutine collision_util_reset_impactors


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


end submodule s_encounter_util