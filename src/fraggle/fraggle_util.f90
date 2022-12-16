!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle_classes) s_fraggle_util
   use swiftest
contains

   module subroutine fraggle_util_add_fragments_to_system(fragments, impactors, system, param)
      !! Author: David A. Minton
      !!
      !! Adds fragments to the temporary system pl object
      implicit none
      ! Arguments
      class(fraggle_fragments),     intent(in)    :: fragments      !! Fraggle fragment system object
      class(collision_impactors),     intent(in)    :: impactors !! Fraggle collider system object
      class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param     !! Current swiftest run configuration parameters
      ! Internals
      integer(I4B) :: i, npl_before, npl_after
      logical, dimension(:), allocatable :: lexclude

      associate(nfrag => fragments%nbody, pl => system%pl, cb => system%cb)
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
         lexclude(impactors%idx(1:impactors%ncoll)) = .true.
         where(lexclude(1:npl_after)) 
            pl%status(1:npl_after) = INACTIVE
         elsewhere
            pl%status(1:npl_after) = ACTIVE
         endwhere

      end associate

      return
   end subroutine fraggle_util_add_fragments_to_system
   

   module subroutine fraggle_util_get_angular_momentum(self) 
      !! Author: David A. Minton
      !!
      !! Calcualtes the current angular momentum of the fragments
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout)  :: self !! Fraggle fragment system object
      ! Internals
      integer(I4B) :: i

      associate(fragments => self, nfrag => self%nbody)
         fragments%Lorbit(:) = 0.0_DP
         fragments%Lspin(:) = 0.0_DP
   
         do i = 1, nfrag
            fragments%Lorbit(:) = fragments%Lorbit(:) + fragments%mass(i) * (fragments%rc(:, i) .cross. fragments%vc(:, i))
            fragments%Lspin(:) = fragments%Lspin(:) + fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(:, i) * fragments%rot(:, i)
         end do
      end associate

      return
   end subroutine fraggle_util_get_angular_momentum


   module subroutine fraggle_util_construct_temporary_system(fragments, system, param, tmpsys, tmpparam)
      !! Author: David A. Minton
      !!
      !! Constructs a temporary internal system consisting of active bodies and additional fragments. This internal temporary system is used to calculate system energy with and without fragments
      implicit none
      ! Arguments
      class(fraggle_fragments),                   intent(in)  :: fragments     !! Fraggle fragment system object
      class(swiftest_nbody_system),               intent(in)  :: system   !! Original swiftest nbody system object
      class(swiftest_parameters),                 intent(in)  :: param    !! Current swiftest run configuration parameters
      class(swiftest_nbody_system), allocatable,  intent(out) :: tmpsys   !! Output temporary swiftest nbody system object
      class(swiftest_parameters),   allocatable,  intent(out) :: tmpparam !! Output temporary configuration run parameters
      ! Internals
      logical, dimension(:), allocatable :: linclude
      integer(I4B) :: npl_tot

      associate(nfrag => fragments%nbody, pl => system%pl, npl => system%pl%nbody, cb => system%cb)
         ! Set up a new system based on the original
         if (allocated(tmpparam)) deallocate(tmpparam)
         if (allocated(tmpsys)) deallocate(tmpsys)
         allocate(tmpparam, source=param)
         call setup_construct_system(tmpsys, tmpparam)

         ! No test particles necessary for energy/momentum calcs
         call tmpsys%tp%setup(0, param)

         ! Replace the empty central body object with a copy of the original
         deallocate(tmpsys%cb)
         allocate(tmpsys%cb, source=cb)

         ! Make space for the fragments
         npl_tot = npl + nfrag
         call tmpsys%pl%setup(npl_tot, tmpparam)
         allocate(linclude(npl_tot))

         ! Fill up the temporary system with all of the original bodies, leaving the spaces for fragments empty until we add them in later
         linclude(1:npl) = .true.
         linclude(npl+1:npl_tot) = .false.
         call tmpsys%pl%fill(pl, linclude)

         ! Scale the temporary system to the natural units of the current Fraggle calculation
         call tmpsys%rescale(tmpparam, fragments%mscale, fragments%dscale, fragments%tscale)

      end associate

      return
   end subroutine fraggle_util_construct_temporary_system


   module subroutine fraggle_util_dealloc_fragments(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables
      implicit none
      ! Arguments
      class(fraggle_fragments),  intent(inout) :: self

      call collision_util_deallocate_fragments(self)

      if (allocated(self%v_r_mag)) deallocate(self%v_r_mag)
      if (allocated(self%v_t_mag)) deallocate(self%v_t_mag)
      if (allocated(self%v_n_mag)) deallocate(self%v_n_mag)

      return
   end subroutine fraggle_util_dealloc_fragments



   module subroutine fraggle_util_final_impactors(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(collision_impactors),  intent(inout) :: self !! Fraggle encountar storage object

      if (allocated(self%idx)) deallocate(self%idx)

      return
   end subroutine fraggle_util_final_impactors


   module subroutine fraggle_util_final_fragments(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(fraggle_fragments),  intent(inout) :: self !! Fraggle encountar storage object

      call self%dealloc()

      return
   end subroutine fraggle_util_final_fragments


   module subroutine fraggle_util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalizer will deallocate all allocatables
      implicit none
      ! Arguments
      type(fraggle_system),  intent(inout) :: self !! Collision impactors storage object

      call self%reset()

      return
   end subroutine fraggle_util_final_system


   module subroutine fraggle_util_restructure(self, impactors, try, f_spin, r_max_start)
      !! Author: David A. Minton
      !!
      !! Restructure the inputs after a failed attempt failed to find a set of positions and velocities that satisfy the energy and momentum constraints
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: self        !! Fraggle fragment system object
      class(collision_impactors), intent(in)    :: impactors   !! Fraggle collider system object
      integer(I4B),             intent(in)    :: try         !! The current number of times Fraggle has tried to find a solution
      real(DP),                 intent(inout) :: f_spin      !! Fraction of energy/momentum that goes into spin. This decreases ater a failed attempt
      real(DP),                 intent(inout) :: r_max_start !! The maximum radial distance that the position calculation starts with. This increases after a failed attempt
      ! Internals
      real(DP), save :: ke_tot_deficit, r_max_start_old, ke_avg_deficit_old
      real(DP) :: delta_r, delta_r_max, ke_avg_deficit
      real(DP), parameter :: ke_avg_deficit_target = 0.0_DP 

      ! Introduce a bit of noise in the radius determination so we don't just flip flop between similar failed positions
      associate(fragments => self)
         call random_number(delta_r_max)
         delta_r_max = sum(impactors%radius(:)) * (1.0_DP + 2e-1_DP * (delta_r_max - 0.5_DP))
         if (try == 1) then
            ke_tot_deficit = - (fragments%ke_budget - fragments%ke_orbit - fragments%ke_spin)
            ke_avg_deficit = ke_tot_deficit
            delta_r = delta_r_max
         else
            ! Linearly interpolate the last two failed solution ke deficits to find a new distance value to try
            ke_tot_deficit = ke_tot_deficit - (fragments%ke_budget - fragments%ke_orbit - fragments%ke_spin)
            ke_avg_deficit = ke_tot_deficit / try
            delta_r = (r_max_start - r_max_start_old) * (ke_avg_deficit_target - ke_avg_deficit_old) &
                                                      / (ke_avg_deficit - ke_avg_deficit_old)
            if (abs(delta_r) > delta_r_max) delta_r = sign(delta_r_max, delta_r)
         end if
         r_max_start_old = r_max_start
         r_max_start = r_max_start + delta_r ! The larger lever arm can help if the problem is in the angular momentum step
         ke_avg_deficit_old = ke_avg_deficit
   
         if (f_spin > epsilon(1.0_DP)) then ! Try reducing the fraction in spin
            f_spin = f_spin / 2
         else
            f_spin = 0.0_DP
         end if
      end associate 

      return
   end subroutine fraggle_util_restructure


   module subroutine fraggle_util_shift_vector_to_origin(m_frag, vec_frag)
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
   end subroutine fraggle_util_shift_vector_to_origin


   module function fraggle_util_vmag_to_vb(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) result(vb) 
      !! Author: David A. Minton
      !!
      !! Converts radial and tangential velocity magnitudes into barycentric velocity
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)  :: v_r_mag   !! Unknown radial component of fragment velocity vector
      real(DP), dimension(:),   intent(in)  :: v_t_mag   !! Tangential component of velocity vector set previously by angular momentum constraint
      real(DP), dimension(:,:), intent(in)  :: v_r_unit, v_t_unit !! Radial and tangential unit vectors for each fragment
      real(DP), dimension(:),   intent(in)  :: m_frag    !! Fragment masses
      real(DP), dimension(:),   intent(in)  :: vcom      !! Barycentric velocity of collisional system center of mass
      ! Result
      real(DP), dimension(:,:), allocatable   :: vb
      ! Internals
      integer(I4B) :: i, nfrag

      allocate(vb, mold=v_r_unit)
      ! Make sure the velocity magnitude stays positive
      nfrag = size(m_frag)
      do i = 1, nfrag
         vb(:,i) = abs(v_r_mag(i)) * v_r_unit(:, i)
      end do
      ! In order to keep satisfying the kinetic energy constraint, we must shift the origin of the radial component of the velocities to the center of mass
      call fraggle_util_shift_vector_to_origin(m_frag, vb)
      
      do i = 1, nfrag
         vb(:, i) = vb(:, i) + v_t_mag(i) * v_t_unit(:, i) + vcom(:)
      end do

      return
   end function fraggle_util_vmag_to_vb


end submodule s_fraggle_util
