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
      class(collision_basic),      intent(in)    :: self      !! Collision system system object
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
   end subroutine collision_util_add_fragments_to_collider


   module subroutine collision_util_construct_temporary_system(self, nbody_system, param, tmpsys, tmpparam)
      !! Author: David A. Minton
      !!
      !! Constructs a temporary internal system consisting of active bodies and additional fragments. This internal temporary system is used to calculate system energy with and without fragments
      implicit none
      ! Arguments
      class(collision_basic),                intent(inout) :: self         !! Fraggle collision system object
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
            allocate(tmpparam_local, source=param)
            select type(tmpparam_local)
            class is (swiftest_parameters)
               tmpparam_local%system_history%nc%lfile_is_open = .false.
            end select
            call swiftest_util_setup_construct_system(tmpsys_local, tmpparam_local)

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


   module subroutine collision_util_get_angular_momentum(self) 
      !! Author: David A. Minton
      !!
      !! Calculates the current angular momentum of the fragments
      implicit none
      ! Arguments
      class(collision_fragments(*)), intent(inout)  :: self !! Fraggle fragment system object
      ! Internals
      integer(I4B) :: i

      associate(fragments => self, nfrag => self%nbody)
   
         do i = 1, nfrag
            fragments%Lorbit(:,i) = fragments%mass(i) * (fragments%rc(:,i) .cross. fragments%vc(:, i))
            fragments%Lspin(:,i)  = fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(:,i) * fragments%rot(:,i)
         end do

         fragments%Lorbit_tot(:) = sum(fragments%Lorbit, dim=2)
         fragments%Lspin_tot(:) = sum(fragments%Lspin, dim=2)
      end associate

      return
   end subroutine collision_util_get_angular_momentum


   module subroutine collision_util_get_kinetic_energy(self) 
      !! Author: David A. Minton
      !!
      !! Calculates the current kinetic energy of the fragments
      implicit none
      ! Argument
      class(collision_fragments(*)), intent(inout)  :: self !! Fragment system object
      ! Internals
      integer(I4B) :: i

      associate(fragments => self, nfrag => self%nbody)
   
         do concurrent(i = 1:nfrag)
            fragments%ke_orbit(i) = fragments%mass(i) * dot_product(fragments%vc(:,i), fragments%vc(:,i))
            fragments%ke_spin(i) =  fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i) * dot_product(fragments%rot(:,i),fragments%rot(:,i) )
         end do

         fragments%ke_orbit(:) = fragments%ke_orbit(:) / 2
         fragments%ke_spin(:) = fragments%ke_spin(:) / 2
         fragments%ke_orbit_tot = sum(fragments%ke_orbit(:))
         fragments%ke_spin_tot = sum(fragments%ke_spin(:))

      end associate

      return
   end subroutine collision_util_get_kinetic_energy


   module subroutine collision_util_get_energy_momentum(self,  nbody_system, param, lbefore)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision outcome state (lbefore = .false.)
      !! This subrourtine works by building a temporary internal massive body object out of the non-excluded bodies and optionally with fragments appended. 
      !! This will get passed to the energy calculation subroutine so that energy is computed exactly the same way is it is in the main program. 
      !! This will temporarily expand the massive body object in a temporary system object called tmpsys to feed it into symba_energy
      implicit none
      ! Arguments
      class(collision_basic),  intent(inout) :: self    !! Encounter collision system object
      class(base_nbody_system), intent(inout) :: nbody_system  !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param   !! Current swiftest run configuration parameters
      logical,                  intent(in)    :: lbefore !! Flag indicating that this the "before" state of the nbody_system, with impactors included and fragments excluded or vice versa
      ! Internals
      integer(I4B)  :: stage,i
      real(DP), dimension(NDIM) :: Lorbit, Lspin
      real(DP) :: ke_orbit, ke_spin

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody, pl => nbody_system%pl, cb => nbody_system%cb)


            if (lbefore) then

               Lorbit(:) = sum(impactors%Lorbit(:,:), dim=2)
               Lspin(:) = sum(impactors%Lspin(:,:), dim=2)
               ke_orbit = 0.0_DP
               ke_spin = 0.0_DP
               do concurrent(i = 1:2)
                  ke_orbit = ke_orbit + impactors%mass(i) * dot_product(impactors%vc(:,i), impactors%vc(:,i))
                  ke_spin = ke_spin + impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i) * dot_product(impactors%rot(:,i), impactors%rot(:,i))
               end do
               ke_orbit = ke_orbit / 2
               ke_spin = ke_spin / 2
               
            else
               call fragments%get_angular_momentum()
               Lorbit(:) = fragments%Lorbit_tot(:) 
               Lspin(:) = fragments%Lspin_tot(:) 

               call fragments%get_kinetic_energy()
               ke_orbit = fragments%ke_orbit_tot
               ke_spin = fragments%ke_spin_tot

            end if 
            ! Calculate the current fragment energy and momentum balances
            if (lbefore) then
               stage = 1
            else
               stage = 2
            end if
            self%Lorbit(:,stage) = Lorbit(:)
            self%Lspin(:,stage) = Lspin(:)
            self%Ltot(:,stage) = Lorbit(:) + Lspin(:)
            self%ke_orbit(stage) = ke_orbit
            self%ke_spin(stage) = ke_spin
            self%Etot(stage) = ke_orbit + ke_spin
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
      self%r_unit(:,:) = 0.0_DP
      self%t_unit(:,:) = 0.0_DP
      self%n_unit(:,:) = 0.0_DP

      return
   end subroutine collision_util_reset_fragments


   module subroutine collision_util_reset_system(self)
      !! author: David A. Minton
      !!
      !! Resets the collider nbody_system and deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_basic),    intent(inout) :: self  !! Collision system object

      select type(before => self%before)
      class is (swiftest_nbody_system)
         if (allocated(before%pl)) deallocate(before%pl)
         if (allocated(before%tp)) deallocate(before%tp)
      end select
      select type(after => self%after)
      class is (swiftest_nbody_system)
         if (allocated(after%pl)) deallocate(after%pl)
         if (allocated(after%tp)) deallocate(after%tp)
      end select

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


   module subroutine collision_util_set_budgets(self)
      !! author: David A. Minton
      !!
      !! Sets the energy and momentum budgets of the fragments based on the collider values and the before/after values of energy and momentum
      implicit none
      ! Arguments
      class(collision_basic), intent(inout) :: self !! Fraggle collision system object

      associate(impactors => self%impactors, fragments => self%fragments)

         fragments%L_budget(:) = self%Ltot(:,1)
         fragments%ke_budget = self%Etot(1) - impactors%Qloss

      end associate
      
      return
   end subroutine collision_util_set_budgets


   module subroutine collision_util_set_coordinate_collider(self)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual fragments.
      implicit none
      ! Arguments
      class(collision_basic),    intent(inout) :: self      !! Collisional nbody_system

      associate(fragments => self%fragments, impactors => self%impactors, nfrag => self%fragments%nbody)
         call impactors%set_coordinate_system() 

         if (.not.allocated(self%fragments)) return
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
   end subroutine collision_util_set_coordinate_collider


   module subroutine collision_util_set_coordinate_impactors(self)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual fragments.
      implicit none
      ! Arguments
      class(collision_impactors),    intent(inout) :: self      !! Collisional nbody_system
      ! Internals
      real(DP), dimension(NDIM) ::  delta_r, delta_v, Ltot
      real(DP)   ::  L_mag, mtot

      associate(impactors => self)
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
      class(collision_basic),  intent(inout) :: self         !! Encounter collision system object
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
      self%fragments%Lorbit_tot(:) = 0.0_DP
      self%fragments%Lspin_tot(:) = 0.0_DP
      self%fragments%L_budget(:) = 0.0_DP
      self%fragments%ke_orbit_tot = 0.0_DP
      self%fragments%ke_spin_tot = 0.0_DP
      self%fragments%ke_budget = 0.0_DP

      return
   end subroutine collision_util_setup_fragments_collider


   module subroutine collision_util_shift_vector_to_origin(m_frag, vec_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the position or velocity of the fragments as needed to align them with the origin
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
      class(swiftest_pl), allocatable :: pl
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
                  call move_alloc(pl, before%pl)
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