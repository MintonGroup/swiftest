!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle) s_fraggle_set
   use swiftest
   use symba
contains

   module subroutine fraggle_set_budgets(self)
      !! author: David A. Minton
      !!
      !! Sets the energy and momentum budgets of the fragments based on the collider values and the before/after values of energy and momentum
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self !! Fraggle collision system object
      ! Internals
      real(DP) :: dEtot
      real(DP), dimension(NDIM) :: dL

      associate(impactors => self%impactors)
         select type(fragments => self%fragments)
         class is (fraggle_fragments(*))

            dEtot = self%Etot(2) - self%Etot(1)
            dL(:) = self%Ltot(:,2) - self%Ltot(:,1)

            fragments%L_budget(:) = -dL(:)
            fragments%ke_budget = -(dEtot - 0.5_DP * fragments%mtot * dot_product(impactors%vbcom(:), impactors%vbcom(:))) - impactors%Qloss 

         end select
      end associate
      
      return
   end subroutine fraggle_set_budgets


   module subroutine fraggle_set_natural_scale_factors(self)
      !! author: David A. Minton
      !!
      !! Scales dimenional quantities to ~O(1) with respect to the collisional system. 
      !! This scaling makes it easier for the non-linear minimization to converge on a solution
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self  !! Fraggle collision system object
      ! Internals
      integer(I4B) :: i

      associate(collision_merge => self, fragments => self%fragments, impactors => self%impactors)
         ! Set scale factors
         collision_merge%Escale = 0.5_DP * ( impactors%mass(1) * dot_product(impactors%vb(:,1), impactors%vb(:,1)) &
                                            + impactors%mass(2) * dot_product(impactors%vb(:,2), impactors%vb(:,2)))
         collision_merge%dscale = sum(impactors%radius(:))
         collision_merge%mscale = fragments%mtot 
         collision_merge%vscale = sqrt(collision_merge%Escale / collision_merge%mscale) 
         collision_merge%tscale = collision_merge%dscale / collision_merge%vscale 
         collision_merge%Lscale = collision_merge%mscale * collision_merge%dscale * collision_merge%vscale

         ! Scale all dimensioned quantities of impactors and fragments
         impactors%rbcom(:)    = impactors%rbcom(:)    / collision_merge%dscale
         impactors%vbcom(:)    = impactors%vbcom(:)    / collision_merge%vscale
         impactors%rbimp(:)    = impactors%rbimp(:)    / collision_merge%dscale
         impactors%rb(:,:)     = impactors%rb(:,:)     / collision_merge%dscale
         impactors%vb(:,:)     = impactors%vb(:,:)     / collision_merge%vscale
         impactors%mass(:)     = impactors%mass(:)     / collision_merge%mscale
         impactors%radius(:)   = impactors%radius(:)   / collision_merge%dscale
         impactors%Lspin(:,:)  = impactors%Lspin(:,:)  / collision_merge%Lscale
         impactors%Lorbit(:,:) = impactors%Lorbit(:,:) / collision_merge%Lscale

         do i = 1, 2
            impactors%rot(:,i) = impactors%Lspin(:,i) / (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3, i))
         end do

         fragments%mtot    = fragments%mtot   / collision_merge%mscale
         fragments%mass    = fragments%mass   / collision_merge%mscale
         fragments%radius  = fragments%radius / collision_merge%dscale
         impactors%Qloss   = impactors%Qloss  / collision_merge%Escale
      end associate

      return
   end subroutine fraggle_set_natural_scale_factors


   module subroutine fraggle_set_original_scale_factors(self)
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

      associate(collision_merge => self, fragments => self%fragments, impactors => self%impactors)

         ! Restore scale factors
         impactors%rbcom(:) = impactors%rbcom(:) * collision_merge%dscale
         impactors%vbcom(:) = impactors%vbcom(:) * collision_merge%vscale
         impactors%rbimp(:) = impactors%rbimp(:) * collision_merge%dscale
   
         impactors%mass   = impactors%mass   * collision_merge%mscale
         impactors%radius = impactors%radius * collision_merge%dscale
         impactors%rb     = impactors%rb     * collision_merge%dscale
         impactors%vb     = impactors%vb     * collision_merge%vscale
         impactors%Lspin  = impactors%Lspin  * collision_merge%Lscale
         do i = 1, 2
            impactors%rot(:,i) = impactors%Lspin(:,i) * (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3, i))
         end do
   
         fragments%mtot   = fragments%mtot   * collision_merge%mscale
         fragments%mass   = fragments%mass   * collision_merge%mscale
         fragments%radius = fragments%radius * collision_merge%dscale
         fragments%rot    = fragments%rot    / collision_merge%tscale
         fragments%rc     = fragments%rc     * collision_merge%dscale
         fragments%vc     = fragments%vc     * collision_merge%vscale
   
         do i = 1, fragments%nbody
            fragments%rb(:, i) = fragments%rc(:, i) + impactors%rbcom(:)
            fragments%vb(:, i) = fragments%vc(:, i) + impactors%vbcom(:)
         end do

         impactors%Qloss = impactors%Qloss * collision_merge%Escale

         collision_merge%Lorbit(:,:) = collision_merge%Lorbit(:,:) * collision_merge%Lscale
         collision_merge%Lspin(:,:)  = collision_merge%Lspin(:,:)  * collision_merge%Lscale
         collision_merge%Ltot(:,:)   = collision_merge%Ltot(:,:)   * collision_merge%Lscale
         collision_merge%ke_orbit(:) = collision_merge%ke_orbit(:) * collision_merge%Escale
         collision_merge%ke_spin(:)  = collision_merge%ke_spin(:)  * collision_merge%Escale
         collision_merge%pe(:)       = collision_merge%pe(:)       * collision_merge%Escale
         collision_merge%Etot(:)     = collision_merge%Etot(:)     * collision_merge%Escale
   
         collision_merge%mscale = 1.0_DP
         collision_merge%dscale = 1.0_DP
         collision_merge%vscale = 1.0_DP
         collision_merge%tscale = 1.0_DP
         collision_merge%Lscale = 1.0_DP
         collision_merge%Escale = 1.0_DP
      end associate
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)
   
      return
   end subroutine fraggle_set_original_scale_factors


end submodule s_fraggle_set