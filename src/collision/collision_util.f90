! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

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
      integer(I4B) :: i, npl_before, npl_after, nfrag
      logical, dimension(:), allocatable :: lexclude

      select type(nbody_system)
      class is (swiftest_nbody_system)
         associate(fragments => self%fragments, impactors => self%impactors, pl => nbody_system%pl, cb => nbody_system%cb)
            npl_after = pl%nbody
            npl_before = npl_after - nfrag
            nfrag = self%fragments%nbody
            allocate(lexclude(npl_after))

            pl%status(npl_before+1:npl_after) = ACTIVE
            pl%mass(npl_before+1:npl_after) = fragments%mass(1:nfrag)
            pl%Gmass(npl_before+1:npl_after) = fragments%mass(1:nfrag) * param%GU
            pl%radius(npl_before+1:npl_after) = fragments%radius(1:nfrag)
#ifdef DOCONLOC
            do concurrent (i = 1:nfrag) shared(cb,pl,fragments)
#else
            do concurrent (i = 1:nfrag)
#endif
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

   module subroutine collision_util_bounce_one(r,v,rcom,vcom,radius)
      !! Author: David A. Minton
      !!
      !! Performs a "bounce" operation on a single body by reversing its velocity in a center of mass frame.
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(inout) :: r,v
      real(DP), dimension(:), intent(in)    :: rcom,vcom
      real(DP),               intent(in)    :: radius
      ! Internals
      real(DP), dimension(NDIM) :: rrel, vrel, rnorm

      rrel(:) = r(:) - rcom(:)
      vrel(:) = v(:) - vcom(:)
      rnorm(:) = .unit. rrel(:)
      ! Do the reflection
      vrel(:) = vrel(:) - 2 * dot_product(vrel(:),rnorm(:)) * rnorm(:)
      ! Shift the positions so that the collision doesn't immediately occur again
      r(:) = r(:) + 0.5_DP * radius * rnorm(:)
      v(:) = vcom(:) + vrel(:)

      return
   end subroutine collision_util_bounce_one


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
      !! Calculates total system energy in either the pre-collision outcome state (phase = "before") or the post-collision outcome 
      !! state (lbefore = .false.)
      implicit none
      ! Arguments
      class(collision_basic),   intent(inout) :: self         !! Encounter collision system object
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current Swiftest run configuration parameters
      character(len=*),         intent(in)    :: phase        !! One of "before" or "after", indicating which phase of the 
                                                              !!    calculation this needs to be done
      ! Internals
      integer(I4B) :: i, phase_val, nfrag

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
         associate(fragments => self%fragments, impactors => self%impactors, pl => nbody_system%pl, cb => nbody_system%cb)
            nfrag = self%fragments%nbody
            if (phase_val == 1) then
#ifdef DOCONLOC
               do concurrent(i = 1:2) shared(impactors)
#else
               do concurrent(i = 1:2)
#endif
                  impactors%ke_orbit(i) = 0.5_DP * impactors%mass(i) * dot_product(impactors%vc(:,i), impactors%vc(:,i))
                  impactors%ke_spin(i) = 0.5_DP * impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i)  & 
                                                * dot_product(impactors%rot(:,i), impactors%rot(:,i))
                  impactors%be(i) = -3 * impactors%Gmass(i) * impactors%mass(i) / (5 * impactors%radius(i))
                  impactors%L_orbit(1,i) = impactors%mass(i) * (impactors%rc(2,i) * impactors%vc(3,i) &
                                                              - impactors%rc(3,i) * impactors%vc(2,i))
                  impactors%L_orbit(2,i) = impactors%mass(i) * (impactors%rc(3,i) * impactors%vc(1,i) &
                                                              - impactors%rc(1,i) * impactors%vc(3,i))
                  impactors%L_orbit(3,i) = impactors%mass(i) * (impactors%rc(1,i) * impactors%vc(2,i) &
                                                              - impactors%rc(2,i) * impactors%vc(1,i))
                  impactors%L_spin(:,i) = impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i) * impactors%rot(:,i)
               end do
               self%L_orbit(:,phase_val) = sum(impactors%L_orbit(:,1:2),dim=2)
               self%L_spin(:,phase_val) = sum(impactors%L_spin(:,1:2),dim=2)
               self%L_total(:,phase_val) = self%L_orbit(:,phase_val) + self%L_spin(:,phase_val) 
               self%ke_orbit(phase_val) = sum(impactors%ke_orbit(1:2))
               self%ke_spin(phase_val) = sum(impactors%ke_spin(1:2))
               self%be(phase_val) = sum(impactors%be(1:2))
               call swiftest_util_get_potential_energy(2, [(.true., i = 1, 2)], 0.0_DP, impactors%Gmass, impactors%mass, & 
                                                            impactors%rb, self%pe(phase_val))
               self%te(phase_val) = self%ke_orbit(phase_val) + self%ke_spin(phase_val) + self%be(phase_val) + self%pe(phase_val)
            else if (phase_val == 2) then
#ifdef DOCONLOC
               do concurrent(i = 1:nfrag) shared(fragments)
#else
               do concurrent(i = 1:nfrag)
#endif
                  fragments%ke_orbit(i) = 0.5_DP * fragments%mass(i) * dot_product(fragments%vc(:,i), fragments%vc(:,i))
                  fragments%ke_spin(i) = 0.5_DP * fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i) & 
                                                * dot_product(fragments%rot(:,i), fragments%rot(:,i))
                  fragments%L_orbit(1,i) = fragments%mass(i) * (fragments%rc(2,i) * fragments%vc(3,i) - &
                                                                fragments%rc(3,i) * fragments%vc(2,i))
                  fragments%L_orbit(2,i) = fragments%mass(i) * (fragments%rc(3,i) * fragments%vc(1,i) - &
                                                                fragments%rc(1,i) * fragments%vc(3,i))
                  fragments%L_orbit(3,i) = fragments%mass(i) * (fragments%rc(1,i) * fragments%vc(2,i) - &
                                                                fragments%rc(2,i) * fragments%vc(1,i))
                  fragments%L_spin(:,i) = fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i) * fragments%rot(:,i)
               end do
               call swiftest_util_get_potential_energy(nfrag, [(.true., i = 1, nfrag)], 0.0_DP, fragments%Gmass, fragments%mass, &
                                                      fragments%rb, fragments%pe)
               fragments%be = sum(-3*fragments%Gmass(1:nfrag)*fragments%mass(1:nfrag)/(5*fragments%radius(1:nfrag)))
               fragments%L_orbit_tot(:) = sum(fragments%L_orbit(:,1:nfrag),dim=2)
               fragments%L_spin_tot(:) = sum(fragments%L_spin(:,1:nfrag),dim=2)
               fragments%ke_orbit_tot = sum(fragments%ke_orbit(1:nfrag))
               fragments%ke_spin_tot = sum(fragments%ke_spin(1:nfrag))
               self%L_orbit(:,phase_val) = fragments%L_orbit_tot(:)
               self%L_spin(:,phase_val) = fragments%L_spin_tot(:)
               self%L_total(:,phase_val) = self%L_orbit(:,phase_val) + self%L_spin(:,phase_val) 
               self%ke_orbit(phase_val) = fragments%ke_orbit_tot
               self%ke_spin(phase_val) = fragments%ke_spin_tot
               self%be(phase_val) = fragments%be
               self%pe(phase_val) = fragments%pe
               self%te(phase_val) = self%ke_orbit(phase_val) + self%ke_spin(phase_val) + self%be(phase_val) + self%pe(phase_val)
            end if

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
      call util_unique(idvals,self%idvals,self%idmap)
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
      self%ke_spin(:) = 0.0_DP
      self%ke_orbit(:) = 0.0_DP
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
      self%rcimp(:) = 0.0_DP

      return
   end subroutine collision_util_dealloc_impactors


   module subroutine collision_util_dealloc_fragments(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_fragments),  intent(inout) :: self

      self%nbody = 0
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%info)) deallocate(self%info) 
      if (allocated(self%status)) deallocate(self%status) 
      if (allocated(self%rh)) deallocate(self%rh)
      if (allocated(self%vh)) deallocate(self%vh)
      if (allocated(self%rb))  deallocate(self%rb)
      if (allocated(self%vb)) deallocate(self%vb)
      if (allocated(self%rc)) deallocate(self%rc)
      if (allocated(self%vc)) deallocate(self%vc)
      if (allocated(self%r_unit)) deallocate(self%r_unit)
      if (allocated(self%t_unit)) deallocate(self%t_unit)
      if (allocated(self%n_unit)) deallocate(self%n_unit)
      if (allocated(self%rot)) deallocate(self%rot)
      if (allocated(self%Ip)) deallocate(self%Ip)
      if (allocated(self%mass)) deallocate(self%mass)
      if (allocated(self%radius)) deallocate(self%radius)
      if (allocated(self%density)) deallocate(self%density)
      if (allocated(self%rmag)) deallocate(self%rmag)
      if (allocated(self%vmag)) deallocate(self%vmag)
      if (allocated(self%rotmag)) deallocate(self%rotmag)
      if (allocated(self%origin_body)) deallocate(self%origin_body)
      if (allocated(self%L_orbit)) deallocate(self%L_orbit)
      if (allocated(self%L_spin)) deallocate(self%L_spin)
      if (allocated(self%ke_orbit)) deallocate(self%ke_orbit)
      if (allocated(self%ke_spin)) deallocate(self%ke_spin)

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
      !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, 
      !! radius, or other values that get set prior to the call to fraggle_generate)
      implicit none
      ! Arguments
      class(collision_fragments), intent(inout) :: self

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


   module subroutine collision_util_setup_fragments(self, n)
      !! author: David A. Minton
      !!
      !! Constructor for fragment class. Allocates space for all particles and
      implicit none
      ! Arguments
      class(collision_fragments), intent(inout) :: self  !! Swiftest generic body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      ! Internals
      integer(I4B) :: i

      if (n < 0) return

      self%nbody = n
      if (n == 0) return

      if(allocated(self%info)) deallocate(self%info); allocate(swiftest_particle_info :: self%info(n))
      if (allocated(self%id)) deallocate(self%id); allocate(self%id(n))
      if (allocated(self%status)) deallocate(self%status); allocate(self%status(n))
      if (allocated(self%rh)) deallocate(self%rh); allocate(self%rh(NDIM, n))
      if (allocated(self%vh)) deallocate(self%vh); allocate(self%vh(NDIM, n))
      if (allocated(self%rb)) deallocate(self%rb); allocate(self%rb(NDIM, n))
      if (allocated(self%vb)) deallocate(self%vb); allocate(self%vb(NDIM, n))
      if (allocated(self%rc)) deallocate(self%rc); allocate(self%rc(NDIM, n))
      if (allocated(self%vc)) deallocate(self%vc); allocate(self%vc(NDIM, n))
      if (allocated(self%r_unit)) deallocate(self%r_unit); allocate(self%r_unit(NDIM, n))
      if (allocated(self%v_unit)) deallocate(self%v_unit); allocate(self%v_unit(NDIM, n))
      if (allocated(self%t_unit)) deallocate(self%t_unit); allocate(self%t_unit(NDIM, n))
      if (allocated(self%n_unit)) deallocate(self%n_unit); allocate(self%n_unit(NDIM, n))
      if (allocated(self%rot)) deallocate(self%rot); allocate(self%rot(NDIM, n))
      if (allocated(self%Ip)) deallocate(self%Ip); allocate(self%Ip(NDIM, n))
      if (allocated(self%Gmass)) deallocate(self%Gmass); allocate(self%Gmass(n))
      if (allocated(self%mass)) deallocate(self%mass); allocate(self%mass(n))
      if (allocated(self%radius)) deallocate(self%radius); allocate(self%radius(n))
      if (allocated(self%density)) deallocate(self%density); allocate(self%density(n))
      if (allocated(self%rmag)) deallocate(self%rmag); allocate(self%rmag(n))
      if (allocated(self%vmag)) deallocate(self%vmag); allocate(self%vmag(n))
      if (allocated(self%rotmag)) deallocate(self%rotmag); allocate(self%rotmag(n))
      if (allocated(self%origin_body)) deallocate(self%origin_body); allocate(self%origin_body(n))
      if (allocated(self%L_orbit)) deallocate(self%L_orbit); allocate(self%L_orbit(NDIM, n))
      if (allocated(self%L_spin)) deallocate(self%L_spin); allocate(self%L_spin(NDIM, n))
      if (allocated(self%ke_orbit)) deallocate(self%ke_orbit); allocate(self%ke_orbit(n))
      if (allocated(self%ke_spin)) deallocate(self%ke_spin); allocate(self%ke_spin(n))

      self%id(:) = 0
      select type(info => self%info)
      class is (swiftest_particle_info)
         do i = 1, n
            call info(i)%set_value(&
               name = "UNNAMED", &
               particle_type = "UNKNOWN", &
               status = "INACTIVE", & 
               origin_type = "UNKNOWN", &
               collision_id = 0, &
               origin_time = -huge(1.0_DP), & 
               origin_rh = [0.0_DP, 0.0_DP, 0.0_DP], &
               origin_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
               discard_time = huge(1.0_DP), & 
               discard_rh = [0.0_DP, 0.0_DP, 0.0_DP], &
               discard_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
               discard_body_id = -1  &
            )
         end do
      end select

      self%mtot = 0.0_DP
      self%status(:) = ACTIVE
      self%rh(:,:)   = 0.0_DP
      self%vh(:,:)   = 0.0_DP
      self%rb(:,:)   = 0.0_DP
      self%vb(:,:)   = 0.0_DP
      self%rc(:,:)   = 0.0_DP
      self%vc(:,:)   = 0.0_DP
      self%r_unit(:,:)   = 0.0_DP
      self%v_unit(:,:)   = 0.0_DP
      self%t_unit(:,:)   = 0.0_DP
      self%n_unit(:,:)   = 0.0_DP
      self%rot(:,:)   = 0.0_DP
      self%Ip(:,:)   = 0.0_DP
      self%Gmass(:)   = 0.0_DP
      self%mass(:)   = 0.0_DP
      self%radius(:)   = 0.0_DP
      self%density(:)   = 0.0_DP
      self%rmag(:)   = 0.0_DP
      self%vmag(:)   = 0.0_DP
      self%rotmag(:)   = 0.0_DP
      self%origin_body(:)   = 0
      self%L_orbit_tot(:)   = 0.0_DP
      self%L_spin_tot(:)   = 0.0_DP
      self%L_orbit(:,:)   = 0.0_DP
      self%L_spin(:,:)   = 0.0_DP
      self%ke_orbit_tot = 0.0_DP
      self%ke_spin_tot = 0.0_DP
      self%pe = 0.0_DP
      self%be = 0.0_DP
      self%ke_orbit(:)   = 0.0_DP
      self%ke_spin(:)   = 0.0_DP

      return
   end subroutine collision_util_setup_fragments


   module subroutine collision_util_set_coordinate_collider(self)
      !! author: David A. Minton
      !! 
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual 
      !! fragments.
      implicit none
      ! Arguments
      class(collision_basic),    intent(inout) :: self      !! Collisional nbody_system

      associate(fragments => self%fragments, impactors => self%impactors)
         call impactors%set_coordinate_system() 

         if (.not.allocated(self%fragments)) return
         call fragments%set_coordinate_system()


      end associate

      return
   end subroutine collision_util_set_coordinate_collider


   module subroutine collision_util_set_coordinate_fragments(self)
      !! author: David A. Minton
      !!
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual 
      !! fragments.
      implicit none
      ! Arguments
      class(collision_fragments), intent(inout) :: self      !! Collisional nbody_system

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
      !! Defines the collisional coordinate nbody_system, including the unit vectors of both the nbody_system and individual 
      !! fragments.
      implicit none
      ! Arguments
      class(collision_impactors), intent(inout) :: self      !! Collisional nbody_system
      ! Internals
      real(DP), dimension(NDIM) ::  delta_r, delta_v, L_total
      real(DP)   ::  L_mag, mtot

      associate(impactors => self)
         delta_v(:) = impactors%vb(:, 2) - impactors%vb(:, 1)
         delta_r(:) = impactors%rb(:, 2) - impactors%rb(:, 1)
   
         ! We will initialize fragments on a plane defined by the pre-impact nbody_system, with the z-axis aligned with the angular 
         ! momentum vector and the y-axis aligned with the pre-impact distance vector.

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
   
         ! Find the point of impact between the two bodies, defined as the location (in the collisional coordinate system) at the 
         ! surface of body 1 along the line connecting the two bodies.
         impactors%rcimp(:) = impactors%rb(:,1) + impactors%radius(1) * impactors%y_unit(:) - impactors%rbcom(:)

         ! Set the velocity direction as the "bounce" direction" for disruptions, and body 2's direction for hit and runs
         if (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) then
            impactors%bounce_unit(:) = .unit. impactors%vc(:,2)
         else
            impactors%bounce_unit(:) = .unit. (impactors%vc(:,2) - 2 * dot_product(impactors%vc(:,2),impactors%y_unit(:)) &
                                                                     * impactors%y_unit(:))
         end if

      end associate

      return
   end subroutine collision_util_set_coordinate_impactors


   module subroutine collision_util_setup_collider(self, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Initializer for the encounter collision system. Sets up impactors and the before/after snapshots,
      !! but not fragments. Those are setup later when the number of fragments is known.
      implicit none
      ! Arguments
      class(collision_basic),   intent(inout) :: self         !! Encounter collision system object
      class(base_nbody_system), intent(in)    :: nbody_system !! Current nbody system. Used as a mold for the before/after snapshots
      class(base_parameters),   intent(inout) :: param        !! Current Swiftest run configuration parameters

      call self%setup_impactors()
      if (allocated(self%before)) deallocate(self%before)
      if (allocated(self%after)) deallocate(self%after)

      allocate(self%before, mold=nbody_system)
      allocate(self%after,  mold=nbody_system)

      self%max_rot = MAX_ROT_SI * param%TU2S          

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
      allocate(collision_fragments :: self%fragments)
      call self%fragments%setup(nfrag)

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
      nfrag = count(m_frag > tiny(0.0_DP))

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
      !! Takes a minimal snapshot of the state of the nbody_system during a collision to record the before and after states of the
      !! system through the collision.
      implicit none
      ! Internals
      class(collision_storage), intent(inout)          :: self   !! Swiftest storage object
      class(base_parameters),   intent(inout)          :: param  !! Current run configuration parameters
      class(base_nbody_system), intent(inout)          :: nbody_system !! Swiftest nbody system object to store
      real(DP),                 intent(in),   optional :: t      !! Time of snapshot if different from nbody_system time
      character(*),             intent(in),   optional :: arg    !! "before": takes a snapshot just before the collision. "after" 
                                                                 !!    takes the snapshot just after the collision.
      ! Arguments
      class(collision_snapshot), allocatable, save :: snapshot
      character(len=:), allocatable :: stage
      integer(I4B) :: i,phase_val
      character(len=STRMAX) :: message

      if (present(arg)) then
         stage = arg
      else
         stage = ""
      end if 

      select type (nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)

         select case (stage)
         case ("before")
            phase_val = 1
            allocate(collision_snapshot :: snapshot)
            allocate(snapshot%collider, source=nbody_system%collider) 
            snapshot%t = t
         case ("after")
            phase_val = 2
         case ("particle")
            phase_val = -1
            allocate(collision_snapshot :: snapshot)
            allocate(snapshot%collider, source=nbody_system%collider) 
         case default
            write(*,*) "collision_util_snapshot requies either 'before', 'after', or 'particle' passed to 'arg'"
            return
         end select

         if (stage /= "particle" ) then
            ! Get and record the energy of the system before the collision
            call nbody_system%get_energy_and_momentum(param)
            snapshot%collider%L_orbit(:,phase_val) = nbody_system%L_orbit(:)
            snapshot%collider%L_spin(:,phase_val) = nbody_system%L_spin(:)
            snapshot%collider%L_total(:,phase_val) = nbody_system%L_total(:)
            snapshot%collider%ke_orbit(phase_val) = nbody_system%ke_orbit
            snapshot%collider%ke_spin(phase_val) = nbody_system%ke_spin
            snapshot%collider%pe(phase_val) = nbody_system%pe
            snapshot%collider%be(phase_val) = nbody_system%be
            snapshot%collider%te(phase_val) = nbody_system%te

            if (stage == "after") then
               select type(before_snap => snapshot%collider%before )
               class is (swiftest_nbody_system)
               select type(before_orig => nbody_system%collider%before)
                  class is (swiftest_nbody_system)
                  select type(plsub => before_orig%pl)
                  class is (swiftest_pl)
                     ! Log the properties of the old and new bodies
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Removing bodies:")
                     do i = 1, plsub%nbody
                        write(message,*) trim(adjustl(plsub%info(i)%name)), " (", trim(adjustl(plsub%info(i)%particle_type)),")"
                        call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                     end do

                     allocate(before_snap%pl, source=plsub)
                  end select
                  deallocate(before_orig%pl)
               end select
               end select

               select type(after_snap => snapshot%collider%after )
               class is (swiftest_nbody_system)
               select type(after_orig => nbody_system%collider%after)
               class is (swiftest_nbody_system)
                  select type(plnew => after_orig%pl)
                  class is (swiftest_pl)
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Adding bodies:")
                     do i = 1, plnew%nbody
                        write(message,*) trim(adjustl(plnew%info(i)%name)), " (", trim(adjustl(plnew%info(i)%particle_type)),")"
                        call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                     end do
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, & 
                        "***********************************************************" // &
                        "***********************************************************")
                     allocate(after_snap%pl, source=plnew)
                  end select
                  deallocate(after_orig%pl)
               end select
               end select
            end if
         end if

         if ((stage == "after") .or. (stage == "particle")) then
            ! Save the snapshot for posterity
            call self%save(snapshot)
            deallocate(snapshot)
         end if
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
         collider%mscale = minval(impactors%mass(:))
         collider%dscale = minval(impactors%radius(:))

         vesc = sqrt(2 * sum(impactors%Gmass(:)) / sum(impactors%radius(:)))
         collider%tscale = collider%dscale / vesc

         ! Set secondary scale factors for convenience
         collider%vscale = collider%dscale / collider%tscale
         collider%Escale = collider%mscale * collider%vscale**2
         collider%Lscale = collider%mscale * collider%dscale * collider%vscale

         ! Scale all dimensioned quantities of impactors and fragments
         impactors%rbcom(:)     = impactors%rbcom(:)      / collider%dscale
         impactors%vbcom(:)     = impactors%vbcom(:)      / collider%vscale
         impactors%rcimp(:)     = impactors%rcimp(:)      / collider%dscale
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
         impactors%ke_orbit(:)  = impactors%ke_orbit(:)   / collider%Escale
         impactors%ke_spin(:)   = impactors%ke_spin(:)    / collider%Escale

         do concurrent(i = 1:2)
            impactors%rot(:,i) = impactors%L_spin(:,i) / (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i))
         end do

         fragments%mtot      = fragments%mtot      / collider%mscale
         fragments%mass(:)   = fragments%mass(:)   / collider%mscale
         fragments%Gmass(:)  = fragments%Gmass(:)  / (collider%dscale**3/collider%tscale**2)
         fragments%radius(:) = fragments%radius(:) / collider%dscale
         impactors%Qloss     = impactors%Qloss     / collider%Escale

         collider%min_mfrag = collider%min_mfrag / collider%mscale
         collider%max_rot   = collider%max_rot   * collider%tscale
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
         impactors%rcimp(:)  = impactors%rcimp(:) * collider%dscale

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
         impactors%ke_orbit  = impactors%ke_orbit  * collider%Escale
         impactors%ke_spin   = impactors%ke_spin   * collider%Escale
#ifdef DOCONLOC
         do concurrent(i = 1:2) shared(impactors)
#else
         do concurrent(i = 1:2)
#endif
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
         collider%min_mfrag    = collider%min_mfrag    * collider%mscale
         collider%max_rot      = collider%max_rot      / collider%tscale
   
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


   module subroutine collision_util_velocity_torque(dL, mass, r, v)
      !! author: David A. Minton
      !!
      !! Applies a torque to a body's center of mass velocity given a change in angular momentum
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(in)    :: dL   !! Change in angular momentum to apply
      real(DP),               intent(in)    :: mass !! Mass of body
      real(DP), dimension(:), intent(in)    :: r    !! Position of body wrt system center of mass
      real(DP), dimension(:), intent(inout) :: v !! Velocity of body wrt system center of mass
      ! Internals
      real(DP), dimension(NDIM) :: dL_unit, r_unit, vapply
      real(DP) :: rmag, vmag, vapply_mag

      dL_unit(:) = .unit. dL
      r_unit(:) = .unit.r(:)
      rmag = .mag.r(:)
      vmag = .mag.v(:)
      vapply_mag = .mag.(dL(:)/(rmag * mass))
      if ((vapply_mag / vmag) > epsilon(1.0_DP)) then
         vapply(:) = vapply_mag * (dL_unit(:) .cross. r_unit(:))
         v(:) = v(:) + vapply(:)
      end if

      return
   end subroutine collision_util_velocity_torque

end submodule s_collision_util