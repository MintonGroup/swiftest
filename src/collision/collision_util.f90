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

   module subroutine collision_util_add_fragments_to_collider(self, nbody_system, param)
      !! Author: David A. Minton
      !!
      !! Adds fragments to the temporary system pl object
      implicit none
      ! Arguments
      class(collision_basic),   intent(in)    :: self      !! Collision system system object
      class(base_nbody_system), intent(inout) :: nbody_system    !! Swiftest nbody system object
      class(base_parameters),   intent(in)    :: param     !! Current Swiftest run configuration parameters
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
   end subroutine collision_util_add_fragments_to_collider


   module subroutine collision_util_construct_constraint_system(collider, nbody_system, param, constraint_system, phase)
      !! Author: David A. Minton
      !!
      !! Constructs a temporary system that is used to evaluate the convergence on energy and angular momentum constraints. 
      !! The rotations of all bodies other than those involved in the collision are set to 0 so that the changes in spin kinetic
      !! energy and momentum are isolated to the collision system.
      implicit none
      ! Arguments
      class(collision_basic),                 intent(inout) :: collider          !! Collision system object
      class(base_nbody_system),               intent(in)    :: nbody_system      !! Original Swiftest nbody system object
      class(base_parameters),                 intent(inout) :: param             !! Current Swiftest run configuration parameters
      class(base_nbody_system), allocatable,  intent(out)   :: constraint_system !! Output temporary Swiftest nbody system object
      character(len=*),                       intent(in)    :: phase             !! One of "before" or "after", indicating which phase of the calculation this needs to be done
      ! Internals
      integer(I4B) :: i, status
      class(swiftest_nbody_system), allocatable :: tmpsys

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         associate(pl => nbody_system%pl, npl => nbody_system%pl%nbody, cb => nbody_system%cb, nfrag => collider%fragments%nbody)

            ! Set up a new temporary system based on the original
            allocate(tmpsys, mold=nbody_system)
            allocate(tmpsys%cb, source=cb)
            allocate(tmpsys%pl, source=pl)
            allocate(tmpsys%collider, source=collider)
            call tmpsys%collider%set_original_scale()

            ! Remove spins and velocities from all bodies other than the new fragments so that we can isolate the kinetic energy and momentum of the collision system, but still be able to compute
            ! the potential energy correctly
            tmpsys%cb%rot(:) = 0.0_DP
            tmpsys%pl%rot(:,:) = 0.0_DP
            tmpsys%pl%vb(:,:) = 0.0_DP

            if (phase == "before") then
               ! Put back the spins and velocities of the colliding bodies to compute pre-impact KE and L
               tmpsys%pl%rot(:,collider%impactors%id(:)) = pl%rot(:,collider%impactors%id(:))
               tmpsys%pl%vb(:,collider%impactors%id(:)) = pl%vb(:,collider%impactors%id(:))
            else if (phase == "after") then
               allocate(tmpsys%pl_adds, mold=pl)
               allocate(tmpsys%pl_discards, mold=pl)
               associate(impactors => tmpsys%collider%impactors, fragments => tmpsys%collider%fragments) ! Be sure to select the temporary version because its unit system has been updated
                  ! Update barycentric vector values
                  do concurrent(i = 1:nfrag)
                     fragments%rb(:,i) = fragments%rc(:,i) + impactors%rbcom(:)
                     fragments%vb(:,i) = fragments%vc(:,i) + impactors%vbcom(:)
                  end do

                  select case(impactors%regime)
                  case(COLLRESOLVE_REGIME_DISRUPTION)
                     status = DISRUPTED
                  case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
                     status = SUPERCATASTROPHIC
                  case(COLLRESOLVE_REGIME_HIT_AND_RUN)
                     status = HIT_AND_RUN_DISRUPT
                  end select
               end associate

               call collision_resolve_mergeaddsub(tmpsys, param, nbody_system%t, status)
               call tmpsys%pl%rearray(tmpsys, param)
            end if
            call move_alloc(tmpsys, constraint_system)

         end associate
      end select
      end select

      return
   end subroutine collision_util_construct_constraint_system


   module subroutine collision_util_get_idvalues_snapshot(self, idvals)
      !! author: David A. Minton
      !!
      !! Returns an array of all id values saved in this snapshot
      implicit none
      ! Arguments
      class(collision_snapshot),               intent(in)  :: self   !! Collision snapshot object
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


   module subroutine collision_util_get_energy_and_momentum(self, nbody_system, param, phase)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy in either the pre-collision outcome state (phase = "before") or the post-collision outcome state (lbefore = .false.)
      !! This subrourtine works by building a temporary internal massive body object out of the non-excluded bodies and optionally with fragments appended. 
      !! This will get passed to the energy calculation subroutine so that energy is computed exactly the same way is it is in the main program. 
      !! This will temporarily expand the massive body object in a temporary system object called constraint_system to feed it into symba_energy
      implicit none
      ! Arguments
      class(collision_basic),   intent(inout) :: self         !! Encounter collision system object
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current Swiftest run configuration parameters
      character(len=*),         intent(in)    :: phase        !! One of "before" or "after", indicating which phase of the calculation this needs to be done
      ! Internals
      class(base_nbody_system), allocatable :: constraint_system
      integer(I4B) :: i, phase_val

      select case(phase)
      case("before")
         phase_val = 1
      case("after")
         phase_val = 2
      case default
         write(*,*) "Unknown value of phase argument passed to collision_util_get_energy_and_momentum: ",trim(adjustl(phase))
         return
      end select

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody, pl => nbody_system%pl, cb => nbody_system%cb)
            call collision_util_construct_constraint_system(self, nbody_system, param, constraint_system, phase)
            select type(constraint_system)
            class is (swiftest_nbody_system)
               call constraint_system%get_energy_and_momentum(param)
               self%L_orbit(:,phase_val) = constraint_system%L_orbit(:) / self%Lscale
               self%L_spin(:,phase_val) = constraint_system%L_spin(:) / self%Lscale
               self%L_total(:,phase_val) = constraint_system%L_total(:) / self%Lscale
               self%ke_orbit(phase_val) = constraint_system%ke_orbit / self%Escale
               self%ke_spin(phase_val) = constraint_system%ke_spin / self%Escale
               self%pe(phase_val) = constraint_system%pe / self%Escale
               self%be(phase_val) = constraint_system%be / self%Escale
               self%te(phase_val) = constraint_system%te / self%Escale

               if (phase_val == 2) then
                  do concurrent(i = 1:nfrag)
                     fragments%ke_orbit(i) = 0.5_DP * fragments%mass(i) * dot_product(fragments%vc(:,i), fragments%vc(:,i))
                     fragments%ke_spin(i) = 0.5_DP * fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i) * dot_product(fragments%rot(:,i), fragments%rot(:,i))
                     fragments%L_orbit(:,i) = fragments%mass(i) * fragments%rc(:,i) .cross. fragments%vc(:,i)
                     fragments%L_spin(:,i) = fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i) * fragments%rot(:,i)
                  end do
                  call swiftest_util_get_potential_energy(nfrag, [(.true., i = 1, nfrag)], constraint_system%cb%Gmass, fragments%Gmass, fragments%mass, fragments%rb, fragments%pe)
                  fragments%be = sum(-3*fragments%Gmass(:)*fragments%mass(:)/(5*fragments%radius(:)))
                  fragments%L_orbit_tot(:) = sum(fragments%L_orbit(:,:),dim=2)
                  fragments%L_spin_tot(:) = sum(fragments%L_spin(:,:),dim=2)
                  fragments%ke_orbit_tot = sum(fragments%ke_orbit(:))
                  fragments%ke_spin_tot = sum(fragments%ke_spin(:))
               end if
            end select

         end associate
      end select
      end select

      return
   end subroutine collision_util_get_energy_and_momentum


   module subroutine collision_util_index_map(self)
      !! author: David A. Minton
      !!
      !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      implicit none
      ! Arguments
      class(collision_storage), intent(inout) :: self  !! Swiftest storage object
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


   module subroutine collision_util_dealloc_snapshot(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      class(collision_snapshot),  intent(inout) :: self !! Collsion snapshot object

      if (allocated(self%collider)) then
         call self%collider%dealloc()
         deallocate(self%collider)
      end if

      call self%encounter_snapshot%dealloc()

      return
   end subroutine collision_util_dealloc_snapshot

   module subroutine collision_util_dealloc_impactors(self)
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
      self%L_spin(:,:) = 0.0_DP
      self%L_orbit(:,:) = 0.0_DP
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
   end subroutine collision_util_dealloc_impactors


   module subroutine collision_util_dealloc_fragments(self)
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
      self%r_unit(:,:) = 0.0_DP
      self%t_unit(:,:) = 0.0_DP
      self%n_unit(:,:) = 0.0_DP

      return
   end subroutine collision_util_dealloc_fragments


   module subroutine collision_util_dealloc_basic(self)
      !! author: David A. Minton
      !!
      !! Resets the collider nbody_system and deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_basic), intent(inout) :: self  !! Collision system object

      if (allocated(self%impactors)) then 
         call self%impactors%dealloc()
         deallocate(self%impactors)
      end if

      if (allocated(self%fragments)) then
         call self%fragments%dealloc()
         deallocate(self%fragments)
      end if

      if (allocated(self%before)) then
         select type(before => self%before)
         class is (swiftest_nbody_system)
            if (allocated(before%pl)) deallocate(before%pl)
            if (allocated(before%tp)) deallocate(before%tp)
         end select
         deallocate(self%before)
      end if

      if (allocated(self%after)) then
         select type(after => self%after)
         class is (swiftest_nbody_system)
            if (allocated(after%pl)) deallocate(after%pl)
            if (allocated(after%tp)) deallocate(after%tp)
         end select
         deallocate(self%after)
      end if

      self%L_orbit(:,:) = 0.0_DP
      self%L_spin(:,:) = 0.0_DP
      self%L_total(:,:) = 0.0_DP
      self%ke_orbit(:) = 0.0_DP
      self%ke_spin(:) = 0.0_DP
      self%pe(:) = 0.0_DP
      self%te(:) = 0.0_DP

      self%dscale = 1.0_DP 
      self%mscale = 1.0_DP 
      self%tscale = 1.0_DP 
      self%vscale = 1.0_DP 
      self%Escale = 1.0_DP 
      self%Lscale = 1.0_DP


      return
   end subroutine collision_util_dealloc_basic


   module subroutine collision_util_reset_fragments(self)
      !! author: David A. Minton
      !!
      !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, radius, or other values that get set prior to the call to fraggle_generate)
      implicit none
      ! Arguments
      class(collision_fragments(*)), intent(inout) :: self

      self%rc(:,:) = 0.0_DP
      self%vc(:,:) = 0.0_DP
      self%rh(:,:) = 0.0_DP
      self%vh(:,:) = 0.0_DP
      self%rb(:,:) = 0.0_DP
      self%vb(:,:) = 0.0_DP
      self%rot(:,:) = 0.0_DP
      self%r_unit(:,:) = 0.0_DP
      self%t_unit(:,:) = 0.0_DP
      self%n_unit(:,:) = 0.0_DP

      self%rmag(:) = 0.0_DP
      self%vmag(:) = 0.0_DP
      self%rotmag(:) = 0.0_DP

      self%L_orbit_tot(:) = 0.0_DP 
      self%L_spin_tot(:)  = 0.0_DP 
      self%L_orbit(:,:)   = 0.0_DP 
      self%L_spin(:,:)    = 0.0_DP 
      self%ke_orbit_tot   = 0.0_DP 
      self%ke_spin_tot    = 0.0_DP 
      self%pe             = 0.0_DP 
      self%be             = 0.0_DP 
      self%ke_orbit(:)    = 0.0_DP 
      self%ke_spin(:)     = 0.0_DP 

      return
   end subroutine collision_util_reset_fragments

   module subroutine collision_util_set_coordinate_collider(self)
      
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual fragments.
      implicit none
      ! Arguments
      class(collision_basic),    intent(inout) :: self      !! Collisional nbody_system

      associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody)
         call impactors%set_coordinate_system() 

         if (.not.allocated(self%fragments)) return
         call fragments%set_coordinate_system()


      end associate

      return
   end subroutine collision_util_set_coordinate_collider


   module subroutine collision_util_set_coordinate_fragments(self)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual fragments.
      implicit none
      ! Arguments
      class(collision_fragments(*)), intent(inout) :: self      !! Collisional nbody_system

      associate(fragments => self, nfrag => self%nbody)
         if ((nfrag == 0) .or. (.not.any(fragments%rc(:,:) > 0.0_DP))) return

         fragments%rmag(:) = .mag. fragments%rc(:,:)
         fragments%vmag(:) = .mag. fragments%vc(:,:)
         fragments%rotmag(:) = .mag. fragments%rot(:,:)
   
         ! Define the radial, normal, and tangential unit vectors for each individual fragment
         fragments%r_unit(:,:) = .unit. fragments%rc(:,:) 
         fragments%v_unit(:,:) = .unit. fragments%vc(:,:) 
         fragments%n_unit(:,:) = .unit. (fragments%rc(:,:) .cross. fragments%vc(:,:))
         fragments%t_unit(:,:) = -.unit. (fragments%r_unit(:,:) .cross. fragments%n_unit(:,:))
      end associate

      return
   end subroutine collision_util_set_coordinate_fragments


   module subroutine collision_util_set_coordinate_impactors(self)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual fragments.
      implicit none
      ! Arguments
      class(collision_impactors), intent(inout) :: self      !! Collisional nbody_system
      ! Internals
      real(DP), dimension(NDIM) ::  delta_r, delta_v, L_total
      real(DP)   ::  L_mag, mtot

      associate(impactors => self)
         delta_v(:) = impactors%vb(:, 2) - impactors%vb(:, 1)
         delta_r(:) = impactors%rb(:, 2) - impactors%rb(:, 1)
   
         ! We will initialize fragments on a plane defined by the pre-impact nbody_system, with the z-axis aligned with the angular momentum vector
         ! and the y-axis aligned with the pre-impact distance vector.

         ! y-axis is the separation distance
         impactors%y_unit(:) = .unit.delta_r(:) 
         L_total = impactors%L_orbit(:,1) + impactors%L_orbit(:,2) + impactors%L_spin(:,1) + impactors%L_spin(:,2)

         L_mag = .mag.L_total(:)
         if (L_mag > sqrt(tiny(L_mag))) then
            impactors%z_unit(:) = .unit.L_total(:) 
         else ! Not enough angular momentum to determine a z-axis direction. We'll just pick a random direction
            call random_number(impactors%z_unit(:))
            impactors%z_unit(:) = .unit.impactors%z_unit(:) 
         end if

         ! The cross product of the y- by z-axis will give us the x-axis
         impactors%x_unit(:) = impactors%y_unit(:) .cross. impactors%z_unit(:)
         impactors%v_unit(:) = .unit.delta_v(:)

         ! Find the center of mass of the collisional system	
         mtot = sum(impactors%mass(:))
         impactors%rbcom(:) = (impactors%mass(1) * impactors%rb(:,1) + impactors%mass(2) * impactors%rb(:,2)) / mtot 
         impactors%vbcom(:) = (impactors%mass(1) * impactors%vb(:,1) + impactors%mass(2) * impactors%vb(:,2)) / mtot

         ! The center of mass coordinate position and velocities
         impactors%rc(:,1) = impactors%rb(:,1) - impactors%rbcom(:)
         impactors%rc(:,2) = impactors%rb(:,2) - impactors%rbcom(:)
         impactors%vc(:,1) = impactors%vb(:,1) - impactors%vbcom(:)
         impactors%vc(:,2) = impactors%vb(:,2) - impactors%vbcom(:)
   
         ! Find the point of impact between the two bodies
         impactors%rbimp(:) = impactors%rb(:,1) + impactors%radius(1) * impactors%y_unit(:) - impactors%rbcom(:)

         ! Set the velocity direction as the "bounce" direction" for disruptions, and body 2's direction for hit and runs
         if (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) then
            impactors%bounce_unit(:) = .unit. impactors%vc(:,2)
         else
            impactors%bounce_unit(:) = .unit. (impactors%vc(:,2) - 2 * dot_product(impactors%vc(:,2),impactors%y_unit(:)) * impactors%y_unit(:))
         end if

      end associate

      return
   end subroutine collision_util_set_coordinate_impactors


   module subroutine collision_util_setup_collider(self, nbody_system)
      !! author: David A. Minton
      !!
      !! Initializer for the encounter collision system. Sets up impactors and the before/after snapshots,
      !! but not fragments. Those are setup later when the number of fragments is known.
      implicit none
      ! Arguments
      class(collision_basic),   intent(inout) :: self         !! Encounter collision system object
      class(base_nbody_system), intent(in)    :: nbody_system !! Current nbody system. Used as a mold for the before/after snapshots

      call self%setup_impactors()
      if (allocated(self%before)) deallocate(self%before)
      if (allocated(self%after)) deallocate(self%after)

      allocate(self%before, mold=nbody_system)
      allocate(self%after,  mold=nbody_system)

      return
   end subroutine collision_util_setup_collider


   module subroutine collision_util_setup_impactors_collider(self)
      !! author: David A. Minton
      !!
      !! Initializer for the impactors for the encounter collision system. Deallocates old impactors before creating new ones
      implicit none
      ! Arguments
      class(collision_basic), intent(inout) :: self   !! Encounter collision system object

      if (allocated(self%impactors)) deallocate(self%impactors)
      allocate(collision_impactors :: self%impactors)

      return
   end subroutine collision_util_setup_impactors_collider


   module subroutine collision_util_setup_fragments_collider(self, nfrag)
      !! author: David A. Minton
      !!
      !! Initializer for the fragments of the collision system. 
      implicit none
      ! Arguments
      class(collision_basic), intent(inout) :: self  !! Encounter collision system object
      integer(I4B),            intent(in)    :: nfrag !! Number of fragments to create

      if (allocated(self%fragments)) deallocate(self%fragments)
      allocate(collision_fragments(nfrag) :: self%fragments)
      self%fragments%nbody = nfrag
      self%fragments%nbody = nfrag
      self%fragments%status(:) = ACTIVE
      self%fragments%rh(:,:) = 0.0_DP
      self%fragments%vh(:,:) = 0.0_DP
      self%fragments%rb(:,:) = 0.0_DP
      self%fragments%vb(:,:) = 0.0_DP
      self%fragments%rc(:,:) = 0.0_DP
      self%fragments%vc(:,:) = 0.0_DP
      self%fragments%rot(:,:) = 0.0_DP
      self%fragments%Ip(:,:) = 0.0_DP
      self%fragments%r_unit(:,:) = 0.0_DP
      self%fragments%t_unit(:,:) = 0.0_DP
      self%fragments%n_unit(:,:) = 0.0_DP
      self%fragments%mass(:) = 0.0_DP
      self%fragments%radius(:) = 0.0_DP
      self%fragments%density(:) = 0.0_DP
      self%fragments%rmag(:) = 0.0_DP
      self%fragments%vmag(:) = 0.0_DP
      self%fragments%L_orbit_tot(:) = 0.0_DP
      self%fragments%L_spin_tot(:) = 0.0_DP
      self%fragments%ke_orbit_tot = 0.0_DP
      self%fragments%ke_spin_tot = 0.0_DP

      return
   end subroutine collision_util_setup_fragments_collider


   module subroutine collision_util_shift_vector_to_origin(m_frag, vec_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the position or velocity of the fragments as needed to align them with the center of mass origin
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: m_frag    !! Fragment masses
      real(DP), dimension(:,:), intent(inout) :: vec_frag  !! Fragment positions or velocities in the center of mass frame

      ! Internals
      real(DP), dimension(NDIM) :: mvec_frag, COM_offset
      integer(I4B) :: i, nfrag
      real(DP) :: mtot

      mvec_frag(:) = 0.0_DP
      mtot = sum(m_frag)
      nfrag = size(m_frag)

      do i = 1, nfrag
         mvec_frag = mvec_frag(:) + vec_frag(:,i) * m_frag(i)
      end do
      COM_offset(:) = -mvec_frag(:) / mtot
      do i = 1, nfrag 
         vec_frag(:, i) = vec_frag(:, i) + COM_offset(:)
      end do

      return
   end subroutine collision_util_shift_vector_to_origin


   module subroutine collision_util_snapshot(self, param, nbody_system, t, arg)
      !! author: David A. Minton
      !!
      !! Takes a minimal snapshot of the state of the nbody_system during an encounter so that the trajectories
      !! can be played back through the encounter
      implicit none
      ! Internals
      class(collision_storage), intent(inout)          :: self   !! Swiftest storage object
      class(base_parameters),   intent(inout)          :: param  !! Current run configuration parameters
      class(base_nbody_system), intent(inout)          :: nbody_system !! Swiftest nbody system object to store
      real(DP),                 intent(in),   optional :: t      !! Time of snapshot if different from nbody_system time
      character(*),             intent(in),   optional :: arg    !! "before": takes a snapshot just before the collision. "after" takes the snapshot just after the collision.
      ! Arguments
      class(collision_snapshot), allocatable, save :: snapshot
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
            allocate(collision_snapshot :: snapshot)
            allocate(snapshot%collider, source=nbody_system%collider) 
            snapshot%t = t

            ! Get and record the energy of the system before the collision
            call nbody_system%get_energy_and_momentum(param)
            snapshot%collider%L_orbit(:,1) = nbody_system%L_orbit(:)
            snapshot%collider%L_spin(:,1) = nbody_system%L_spin(:)
            snapshot%collider%L_total(:,1) = nbody_system%L_total(:)
            snapshot%collider%ke_orbit(1) = nbody_system%ke_orbit
            snapshot%collider%ke_spin(1) = nbody_system%ke_spin
            snapshot%collider%pe(1) = nbody_system%pe
            snapshot%collider%be(1) = nbody_system%be
            snapshot%collider%te(1) = nbody_system%te

         case("after")
            ! Get record the energy of the sytem after the collision
            call nbody_system%get_energy_and_momentum(param)
            snapshot%collider%L_orbit(:,2) = nbody_system%L_orbit(:)
            snapshot%collider%L_spin(:,2) = nbody_system%L_spin(:)
            snapshot%collider%L_total(:,2) = nbody_system%L_total(:)
            snapshot%collider%ke_orbit(2) = nbody_system%ke_orbit
            snapshot%collider%ke_spin(2) = nbody_system%ke_spin
            snapshot%collider%pe(2) = nbody_system%pe
            snapshot%collider%be(2) = nbody_system%be
            snapshot%collider%te(2) = nbody_system%te

            select type(before_snap => snapshot%collider%before )
            class is (swiftest_nbody_system)
            select type(before_orig => nbody_system%collider%before)
            class is (swiftest_nbody_system)
               call move_alloc(before_orig%pl, before_snap%pl)
            end select
            end select

            select type(after_snap => snapshot%collider%after )
            class is (swiftest_nbody_system)
            select type(after_orig => nbody_system%collider%after)
            class is (swiftest_nbody_system)
               call move_alloc(after_orig%pl, after_snap%pl)
            end select
            end select

            ! Save the snapshot for posterity
            call self%save(snapshot)
            deallocate(snapshot)
         case default
            write(*,*) "collision_util_snapshot requies either 'before' or 'after' passed to 'arg'"
         end select
      end select
      end select

      return
   end subroutine collision_util_snapshot


   module subroutine collision_util_set_natural_scale_factors(self)
      !! author: David A. Minton
      !!
      !! Scales dimenional quantities to ~O(1) with respect to the collisional system. 
      !! This scaling makes it it easier to converge on a solution without having floating point issues
      implicit none
      ! Arguments
      class(collision_basic), intent(inout) :: self  !! Collision system object
      ! Internals
      integer(I4B) :: i
      real(DP) :: vesc

      associate(collider => self, fragments => self%fragments, impactors => self%impactors)
         ! Set primary scale factors (mass, length, and time) based on the impactor properties at the time of collision
         collider%mscale = minval(fragments%mass(:))
         collider%dscale = minval(fragments%radius(:))

         vesc = sqrt(2 * sum(impactors%Gmass(:)) / sum(impactors%radius(:)))
         collider%tscale = collider%dscale / vesc

         ! Set secondary scale factors for convenience
         collider%vscale = collider%dscale / collider%tscale
         collider%Escale = collider%mscale * collider%vscale**2
         collider%Lscale = collider%mscale * collider%dscale * collider%vscale

         ! Scale all dimensioned quantities of impactors and fragments
         impactors%rbcom(:)     = impactors%rbcom(:)      / collider%dscale
         impactors%vbcom(:)     = impactors%vbcom(:)      / collider%vscale
         impactors%rbimp(:)     = impactors%rbimp(:)      / collider%dscale
         impactors%rb(:,:)      = impactors%rb(:,:)       / collider%dscale
         impactors%vb(:,:)      = impactors%vb(:,:)       / collider%vscale
         impactors%rc(:,:)      = impactors%rc(:,:)       / collider%dscale
         impactors%vc(:,:)      = impactors%vc(:,:)       / collider%vscale
         impactors%mass(:)      = impactors%mass(:)       / collider%mscale
         impactors%Gmass(:)     = impactors%Gmass(:)      / (collider%dscale**3/collider%tscale**2)
         impactors%Mcb          = impactors%Mcb           / collider%mscale
         impactors%radius(:)    = impactors%radius(:)     / collider%dscale
         impactors%L_spin(:,:)  = impactors%L_spin(:,:)   / collider%Lscale
         impactors%L_orbit(:,:) = impactors%L_orbit(:,:)  / collider%Lscale

         do concurrent(i = 1:2)
            impactors%rot(:,i) = impactors%L_spin(:,i) / (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i))
         end do

         fragments%mtot      = fragments%mtot      / collider%mscale
         fragments%mass(:)   = fragments%mass(:)   / collider%mscale
         fragments%Gmass(:)  = fragments%Gmass(:)  / (collider%dscale**3/collider%tscale**2)
         fragments%radius(:) = fragments%radius(:) / collider%dscale
         impactors%Qloss     = impactors%Qloss     / collider%Escale
      end associate

      return
   end subroutine collision_util_set_natural_scale_factors


   module subroutine collision_util_set_original_scale_factors(self)
      !! author: David A. Minton
      !!
      !! Restores dimenional quantities back to the system units
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(collision_basic),      intent(inout) :: self      !! Fragment system object
      ! Internals
      integer(I4B) :: i
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      call ieee_set_halting_mode(IEEE_ALL,.false.)

      associate(collider => self, fragments => self%fragments, impactors => self%impactors)

         ! Restore scale factors
         impactors%rbcom(:)  = impactors%rbcom(:) * collider%dscale
         impactors%vbcom(:)  = impactors%vbcom(:) * collider%vscale
         impactors%rbimp(:)  = impactors%rbimp(:) * collider%dscale

         impactors%mass      = impactors%mass      * collider%mscale
         impactors%Gmass(:)  = impactors%Gmass(:)  * (collider%dscale**3/collider%tscale**2)
         impactors%Mcb       = impactors%Mcb       * collider%mscale
         impactors%mass_dist = impactors%mass_dist * collider%mscale
         impactors%radius    = impactors%radius    * collider%dscale
         impactors%rb        = impactors%rb        * collider%dscale
         impactors%vb        = impactors%vb        * collider%vscale
         impactors%rc        = impactors%rc        * collider%dscale
         impactors%vc        = impactors%vc        * collider%vscale
         impactors%L_spin    = impactors%L_spin    * collider%Lscale
         impactors%L_orbit   = impactors%L_orbit   * collider%Lscale
         do concurrent(i = 1:2)
            impactors%rot(:,i) = impactors%L_spin(:,i) * (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i))
         end do
   
         fragments%mtot      = fragments%mtot      * collider%mscale
         fragments%mass(:)   = fragments%mass(:)   * collider%mscale
         fragments%Gmass(:)  = fragments%Gmass(:)  * (collider%dscale**3/collider%tscale**2)
         fragments%radius(:) = fragments%radius(:) * collider%dscale
         fragments%rot(:,:)  = fragments%rot(:,:)  / collider%tscale
         fragments%rc(:,:)   = fragments%rc(:,:)   * collider%dscale
         fragments%vc(:,:)   = fragments%vc(:,:)   * collider%vscale
         fragments%rb(:,:)   = fragments%rb(:,:)   * collider%dscale
         fragments%vb(:,:)   = fragments%vb(:,:)   * collider%vscale

         impactors%Qloss = impactors%Qloss * collider%Escale

         collider%L_orbit(:,:) = collider%L_orbit(:,:) * collider%Lscale
         collider%L_spin(:,:)  = collider%L_spin(:,:)  * collider%Lscale
         collider%L_total(:,:) = collider%L_total(:,:) * collider%Lscale
         collider%ke_orbit(:)  = collider%ke_orbit(:)  * collider%Escale
         collider%ke_spin(:)   = collider%ke_spin(:)   * collider%Escale
         collider%pe(:)        = collider%pe(:)        * collider%Escale
         collider%be(:)        = collider%be(:)        * collider%Escale
         collider%te(:)        = collider%te(:)        * collider%Escale
   
         collider%mscale = 1.0_DP
         collider%dscale = 1.0_DP
         collider%vscale = 1.0_DP
         collider%tscale = 1.0_DP
         collider%Lscale = 1.0_DP
         collider%Escale = 1.0_DP
      end associate
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)
   
      return
   end subroutine collision_util_set_original_scale_factors


end submodule s_collision_util