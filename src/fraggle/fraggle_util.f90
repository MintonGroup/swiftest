!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle) s_fraggle_util
   use swiftest
contains

   module subroutine fraggle_util_construct_temporary_system(self, nbody_system, param, tmpsys, tmpparam)
      !! Author: David A. Minton
      !!
      !! Constructs a temporary internal system consisting of active bodies and additional fragments. This internal temporary system is used to calculate system energy with and without fragments
      implicit none
      ! Arguments
      class(collision_fraggle),                  intent(inout) :: self         !! Fraggle collision system object
      class(base_nbody_system),               intent(in)    :: nbody_system !! Original swiftest nbody system object
      class(base_parameters),                 intent(in)    :: param        !! Current swiftest run configuration parameters
      class(base_nbody_system), allocatable,  intent(out)   :: tmpsys       !! Output temporary swiftest nbody system object
      class(base_parameters),   allocatable,  intent(out)   :: tmpparam     !! Output temporary configuration run parameters

      call collision_util_construct_temporary_system(self, nbody_system, param, tmpsys, tmpparam)

      select type(tmpsys)
      class is (swiftest_nbody_system)
      select type(tmpparam)
      class is (swiftest_parameters)
         call tmpsys%rescale(tmpparam, self%mscale, self%dscale, self%tscale)
      end select
      end select

      return
   end subroutine fraggle_util_construct_temporary_system


   module subroutine fraggle_util_reset_fragments(self)
      !! author: David A. Minton
      !!
      !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, radius, or other values that get set prior to the call to fraggle_generate)
      implicit none
      ! Arguments
      class(fraggle_fragments(*)), intent(inout) :: self

      self%rc(:,:) = 0.0_DP
      self%vc(:,:) = 0.0_DP
      self%rh(:,:) = 0.0_DP
      self%vh(:,:) = 0.0_DP
      self%rb(:,:) = 0.0_DP
      self%vb(:,:) = 0.0_DP
      self%rot(:,:) = 0.0_DP
      self%v_r_unit(:,:) = 0.0_DP
      self%v_t_unit(:,:) = 0.0_DP
      self%v_n_unit(:,:) = 0.0_DP

      self%rmag(:) = 0.0_DP
      self%rotmag(:) = 0.0_DP
      self%v_r_mag(:) = 0.0_DP
      self%v_t_mag(:) = 0.0_DP
      self%v_n_mag(:) = 0.0_DP

      return
   end subroutine fraggle_util_reset_fragments


   module subroutine fraggle_util_reset_system(self)
      !! author: David A. Minton
      !!
      !! Resets the collider system and deallocates all allocatables
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self  !! Collision system object

      self%dscale = 1.0_DP
      self%mscale = 1.0_DP
      self%tscale = 1.0_DP
      self%vscale = 1.0_DP
      self%Escale = 1.0_DP
      self%Lscale = 1.0_DP

      call collision_util_reset_system(self)

      return
   end subroutine fraggle_util_reset_system


   module subroutine fraggle_util_set_natural_scale_factors(self)
      !! author: David A. Minton
      !!
      !! Scales dimenional quantities to ~O(1) with respect to the collisional system. 
      !! This scaling makes it easier for the non-linear minimization to converge on a solution
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self  !! Fraggle collision system object
      ! Internals
      integer(I4B) :: i

      associate(collider => self, fragments => self%fragments, impactors => self%impactors)
         ! Set scale factors
         collider%Escale = 0.5_DP * ( impactors%mass(1) * dot_product(impactors%vb(:,1), impactors%vb(:,1)) &
                                            + impactors%mass(2) * dot_product(impactors%vb(:,2), impactors%vb(:,2)))
         collider%dscale = sum(impactors%radius(:))
         collider%mscale = fragments%mtot 
         collider%vscale = sqrt(collider%Escale / collider%mscale) 
         collider%tscale = collider%dscale / collider%vscale 
         collider%Lscale = collider%mscale * collider%dscale * collider%vscale

         ! Scale all dimensioned quantities of impactors and fragments
         impactors%rbcom(:)    = impactors%rbcom(:)    / collider%dscale
         impactors%vbcom(:)    = impactors%vbcom(:)    / collider%vscale
         impactors%rbimp(:)    = impactors%rbimp(:)    / collider%dscale
         impactors%bounce_unit(:)    = impactors%bounce_unit(:)    / collider%vscale
         impactors%rb(:,:)     = impactors%rb(:,:)     / collider%dscale
         impactors%vb(:,:)     = impactors%vb(:,:)     / collider%vscale
         impactors%rc(:,:)     = impactors%rc(:,:)     / collider%dscale
         impactors%vc(:,:)     = impactors%vc(:,:)     / collider%vscale
         impactors%mass(:)     = impactors%mass(:)     / collider%mscale
         impactors%Mcb         = impactors%Mcb         / collider%mscale
         impactors%radius(:)   = impactors%radius(:)   / collider%dscale
         impactors%Lspin(:,:)  = impactors%Lspin(:,:)  / collider%Lscale
         impactors%Lorbit(:,:) = impactors%Lorbit(:,:) / collider%Lscale

         do i = 1, 2
            impactors%rot(:,i) = impactors%Lspin(:,i) / (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3, i))
         end do

         fragments%mtot    = fragments%mtot   / collider%mscale
         fragments%mass    = fragments%mass   / collider%mscale
         fragments%radius  = fragments%radius / collider%dscale
         impactors%Qloss   = impactors%Qloss  / collider%Escale
      end associate

      return
   end subroutine fraggle_util_set_natural_scale_factors


   module subroutine fraggle_util_set_original_scale_factors(self)
      !! author: David A. Minton
      !!
      !! Restores dimenional quantities back to the system units
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(collision_fraggle),      intent(inout) :: self      !! Fraggle fragment system object
      ! Internals
      integer(I4B) :: i
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      call ieee_set_halting_mode(IEEE_ALL,.false.)

      associate(collider => self, fragments => self%fragments, impactors => self%impactors)

         ! Restore scale factors
         impactors%rbcom(:) = impactors%rbcom(:) * collider%dscale
         impactors%vbcom(:) = impactors%vbcom(:) * collider%vscale
         impactors%rbimp(:) = impactors%rbimp(:) * collider%dscale
         impactors%bounce_unit(:) = impactors%bounce_unit(:) * collider%vscale
   
         impactors%mass      = impactors%mass      * collider%mscale
         impactors%Mcb       = impactors%Mcb       * collider%mscale
         impactors%mass_dist = impactors%mass_dist * collider%mscale
         impactors%radius    = impactors%radius    * collider%dscale
         impactors%rb        = impactors%rb        * collider%dscale
         impactors%vb        = impactors%vb        * collider%vscale
         impactors%rc        = impactors%rc        * collider%dscale
         impactors%vc        = impactors%vc        * collider%vscale
         impactors%Lspin     = impactors%Lspin     * collider%Lscale
         impactors%Lorbit    = impactors%Lorbit    * collider%Lscale
         do i = 1, 2
            impactors%rot(:,i) = impactors%Lspin(:,i) * (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3, i))
         end do
   
         fragments%mtot   = fragments%mtot   * collider%mscale
         fragments%mass   = fragments%mass   * collider%mscale
         fragments%radius = fragments%radius * collider%dscale
         fragments%rot    = fragments%rot    / collider%tscale
         fragments%rc     = fragments%rc     * collider%dscale
         fragments%vc     = fragments%vc     * collider%vscale
   
         do i = 1, fragments%nbody
            fragments%rb(:, i) = fragments%rc(:, i) + impactors%rbcom(:)
            fragments%vb(:, i) = fragments%vc(:, i) + impactors%vbcom(:)
         end do

         impactors%Qloss = impactors%Qloss * collider%Escale

         collider%Lorbit(:,:) = collider%Lorbit(:,:) * collider%Lscale
         collider%Lspin(:,:)  = collider%Lspin(:,:)  * collider%Lscale
         collider%Ltot(:,:)   = collider%Ltot(:,:)   * collider%Lscale
         collider%ke_orbit(:) = collider%ke_orbit(:) * collider%Escale
         collider%ke_spin(:)  = collider%ke_spin(:)  * collider%Escale
         collider%pe(:)       = collider%pe(:)       * collider%Escale
         collider%Etot(:)     = collider%Etot(:)     * collider%Escale
   
         collider%mscale = 1.0_DP
         collider%dscale = 1.0_DP
         collider%vscale = 1.0_DP
         collider%tscale = 1.0_DP
         collider%Lscale = 1.0_DP
         collider%Escale = 1.0_DP
      end associate
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)
   
      return
   end subroutine fraggle_util_set_original_scale_factors


   module subroutine fraggle_util_setup_fragments_system(self, nfrag)
      !! author: David A. Minton
      !!
      !! Initializer for the fragments of the collision system. 
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self  !! Encounter collision system object
      integer(I4B),          intent(in)    :: nfrag !! Number of fragments to create

      if (allocated(self%fragments)) deallocate(self%fragments)
      allocate(fraggle_fragments(nbody=nfrag) :: self%fragments)
      self%fragments%nbody = nfrag

      return
   end subroutine fraggle_util_setup_fragments_system


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
      call collision_util_shift_vector_to_origin(m_frag, vb)
      
      do i = 1, nfrag
         vb(:, i) = vb(:, i) + v_t_mag(i) * v_t_unit(:, i) + vcom(:)
      end do

      return
   end function fraggle_util_vmag_to_vb


end submodule s_fraggle_util
