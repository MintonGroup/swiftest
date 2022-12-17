!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle_classes) s_fraggle_set
   use swiftest
contains

   module subroutine fraggle_set_budgets(self)
      !! author: David A. Minton
      !!
      !! Sets the energy and momentum budgets of the fragments based on the collider values and the before/after values of energy and momentum
      implicit none
      ! Arguments
      class(fraggle_system), intent(inout) :: self !! Fraggle collision system object
      ! Internals
      real(DP) :: dEtot
      real(DP), dimension(NDIM) :: dL

      associate(impactors => self%impactors)
      select type(fragments => self%fragments)
      class is (fraggle_fragments)

         dEtot = self%Etot(2) - self%Etot(1)
         dL(:) = self%Ltot(:,2) - self%Ltot(:,1)

         fragments%L_budget(:) = -dL(:)
         fragments%ke_budget = -(dEtot - 0.5_DP * fragments%mtot * dot_product(impactors%vbcom(:), impactors%vbcom(:))) - impactors%Qloss 

      end select
      end associate
      return
   end subroutine fraggle_set_budgets


   module subroutine fraggle_set_mass_dist(self, param)
      !! author: David A. Minton
      !!
      !! Sets the mass of fragments based on the mass distribution returned by the regime calculation.
      !! This subroutine must be run after the the setup routine has been run on the fragments
      !!
      implicit none
      ! Arguments
      class(fraggle_system),        intent(inout) :: self  !! Fraggle collision system object
      class(swiftest_parameters),   intent(in)    :: param !! Current Swiftest run configuration parameters
      ! Internals
      integer(I4B)              :: i, jproj, jtarg, nfrag, istart
      real(DP), dimension(2)    :: volume
      real(DP), dimension(NDIM) :: Ip_avg
      real(DP) :: mfrag, mremaining, min_mfrag
      real(DP), parameter :: BETA = 2.85_DP
      integer(I4B), parameter :: NFRAGMAX = 100  !! Maximum number of fragments that can be generated
      integer(I4B), parameter :: NFRAGMIN = 7 !! Minimum number of fragments that can be generated (set by the fraggle_generate algorithm for constraining momentum and energy)
      integer(I4B), parameter :: NFRAG_SIZE_MULTIPLIER = 3 !! Log-space scale factor that scales the number of fragments by the collisional system mass
      integer(I4B), parameter :: iMlr = 1
      integer(I4B), parameter :: iMslr = 2
      integer(I4B), parameter :: iMrem = 3
     
      associate(impactors => self%impactors)
      select type(fragments => self%fragments)
      class is (fraggle_fragments)
         ! Get mass weighted mean of Ip and density
         volume(1:2) = 4._DP / 3._DP * PI * impactors%radius(1:2)**3
         Ip_avg(:) = (impactors%mass(1) * impactors%Ip(:,1) + impactors%mass(2) * impactors%Ip(:,2)) / fragments%mtot
         if (impactors%mass(1) > impactors%mass(2)) then
            jtarg = 1
            jproj = 2
         else
            jtarg = 2
            jproj = 1
         end if
  
         select case(impactors%regime)
         case(COLLRESOLVE_REGIME_DISRUPTION, COLLRESOLVE_REGIME_SUPERCATASTROPHIC, COLLRESOLVE_REGIME_HIT_AND_RUN)
            ! The first two bins of the mass_dist are the largest and second-largest fragments that came out of collision_regime.
            ! The remainder from the third bin will be distributed among nfrag-2 bodies. The following code will determine nfrag based on
            ! the limits bracketed above and the model size distribution of fragments.
            ! Check to see if our size distribution would give us a smaller number of fragments than the maximum number

            select type(param)
            class is (symba_parameters)
               min_mfrag = (param%min_GMfrag / param%GU) 
               ! The number of fragments we generate is bracked by the minimum required by fraggle_generate (7) and the 
               ! maximum set by the NFRAG_SIZE_MULTIPLIER which limits the total number of fragments to prevent the nbody
               ! code from getting an overwhelmingly large number of fragments
               nfrag = ceiling(NFRAG_SIZE_MULTIPLIER  * log(fragments%mtot / min_mfrag))
               nfrag = max(min(nfrag, NFRAGMAX), NFRAGMIN)
            class default
               min_mfrag = 0.0_DP
               nfrag = NFRAGMAX
            end select

            i = iMrem
            mremaining = impactors%mass_dist(iMrem)
            do while (i <= nfrag)
               mfrag = (1 + i - iMslr)**(-3._DP / BETA) * impactors%mass_dist(iMslr)
               if (mremaining - mfrag < 0.0_DP) exit
               mremaining = mremaining - mfrag
               i = i + 1
            end do
            if (i < nfrag) nfrag = max(i, NFRAGMIN)  ! The sfd would actually give us fewer fragments than our maximum
    
            call fragments%setup(nfrag, param)
         case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE) 
            call fragments%setup(1, param)
            fragments%mass(1) = impactors%mass_dist(1)
            fragments%radius(1) = impactors%radius(jtarg)
            fragments%density(1) = impactors%mass_dist(1) / volume(jtarg)
            if (param%lrotation) fragments%Ip(:, 1) = impactors%Ip(:,1)
            return
         case default
            write(*,*) "fraggle_set_mass_dist_fragments error: Unrecognized regime code",impactors%regime
         end select

         ! Make the first two bins the same as the Mlr and Mslr values that came from collision_regime
         fragments%mass(1) = impactors%mass_dist(iMlr) 
         fragments%mass(2) = impactors%mass_dist(iMslr) 

         ! Distribute the remaining mass the 3:nfrag bodies following the model SFD given by slope BETA 
         mremaining = impactors%mass_dist(iMrem)
         do i = iMrem, nfrag
            mfrag = (1 + i - iMslr)**(-3._DP / BETA) * impactors%mass_dist(iMslr)
            fragments%mass(i) = mfrag
            mremaining = mremaining - mfrag
         end do

         ! If there is any residual mass (either positive or negative) we will distribute remaining mass proportionally among the the fragments
         if (mremaining < 0.0_DP) then ! If the remainder is negative, this means that that the number of fragments required by the SFD is smaller than our lower limit set by fraggle_generate. 
            istart = iMrem ! We will reduce the mass of the 3:nfrag bodies to prevent the second-largest fragment from going smaller
         else ! If the remainder is postiive, this means that the number of fragments required by the SFD is larger than our upper limit set by computational expediency. 
            istart = iMslr ! We will increase the mass of the 2:nfrag bodies to compensate, which ensures that the second largest fragment remains the second largest
         end if
         mfrag = 1._DP + mremaining / sum(fragments%mass(istart:nfrag))
         fragments%mass(istart:nfrag) = fragments%mass(istart:nfrag) * mfrag

         ! There may still be some small residual due to round-off error. If so, simply add it to the last bin of the mass distribution.
         mremaining = fragments%mtot - sum(fragments%mass(1:nfrag))
         fragments%mass(nfrag) = fragments%mass(nfrag) + mremaining

         ! Compute physical properties of the new fragments
         select case(impactors%regime)
         case(COLLRESOLVE_REGIME_HIT_AND_RUN)  ! The hit and run case always preserves the largest body intact, so there is no need to recompute the physical properties of the first fragment
            fragments%radius(1) = impactors%radius(jtarg)
            fragments%density(1) = impactors%mass_dist(iMlr) / volume(jtarg)
            fragments%Ip(:, 1) = impactors%Ip(:,1)
            istart = 2
         case default
            istart = 1
         end select
         fragments%density(istart:nfrag) = fragments%mtot / sum(volume(:))
         fragments%radius(istart:nfrag) = (3 * fragments%mass(istart:nfrag) / (4 * PI * fragments%density(istart:nfrag)))**(1.0_DP / 3.0_DP)
         do i = istart, nfrag
            fragments%Ip(:, i) = Ip_avg(:)
         end do
      end select
      end associate

      return
   end subroutine fraggle_set_mass_dist


   module subroutine fraggle_set_natural_scale_factors(self)
      !! author: David A. Minton
      !!
      !! Scales dimenional quantities to ~O(1) with respect to the collisional system. 
      !! This scaling makes it easier for the non-linear minimization to converge on a solution
      implicit none
      ! Arguments
      class(fraggle_system), intent(inout) :: self  !! Fraggle collision system object
      ! Internals
      integer(I4B) :: i

      associate(collision_system => self, fragments => self%fragments, impactors => self%impactors)
         ! Set scale factors
         collision_system%Escale = 0.5_DP * ( impactors%mass(1) * dot_product(impactors%vb(:,1), impactors%vb(:,1)) &
                                            + impactors%mass(2) * dot_product(impactors%vb(:,2), impactors%vb(:,2)))
         collision_system%dscale = sum(impactors%radius(:))
         collision_system%mscale = fragments%mtot 
         collision_system%vscale = sqrt(collision_system%Escale / collision_system%mscale) 
         collision_system%tscale = collision_system%dscale / collision_system%vscale 
         collision_system%Lscale = collision_system%mscale * collision_system%dscale * collision_system%vscale

         ! Scale all dimensioned quantities of impactors and fragments
         impactors%rbcom(:)    = impactors%rbcom(:)    / collision_system%dscale
         impactors%vbcom(:)    = impactors%vbcom(:)    / collision_system%vscale
         impactors%rbimp(:)    = impactors%rbimp(:)    / collision_system%dscale
         impactors%rb(:,:)     = impactors%rb(:,:)     / collision_system%dscale
         impactors%vb(:,:)     = impactors%vb(:,:)     / collision_system%vscale
         impactors%mass(:)     = impactors%mass(:)     / collision_system%mscale
         impactors%radius(:)   = impactors%radius(:)   / collision_system%dscale
         impactors%Lspin(:,:)  = impactors%Lspin(:,:)  / collision_system%Lscale
         impactors%Lorbit(:,:) = impactors%Lorbit(:,:) / collision_system%Lscale

         do i = 1, 2
            impactors%rot(:,i) = impactors%Lspin(:,i) / (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3, i))
         end do

         fragments%mtot    = fragments%mtot   / collision_system%mscale
         fragments%mass    = fragments%mass   / collision_system%mscale
         fragments%radius  = fragments%radius / collision_system%dscale
         impactors%Qloss   = impactors%Qloss  / collision_system%Escale
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
      class(fraggle_system),      intent(inout) :: self      !! Fraggle fragment system object
      ! Internals
      integer(I4B) :: i
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      call ieee_set_halting_mode(IEEE_ALL,.false.)

      associate(collision_system => self, fragments => self%fragments, impactors => self%impactors)

         ! Restore scale factors
         impactors%rbcom(:) = impactors%rbcom(:) * collision_system%dscale
         impactors%vbcom(:) = impactors%vbcom(:) * collision_system%vscale
         impactors%rbimp(:) = impactors%rbimp(:) * collision_system%dscale
   
         impactors%mass   = impactors%mass   * collision_system%mscale
         impactors%radius = impactors%radius * collision_system%dscale
         impactors%rb     = impactors%rb     * collision_system%dscale
         impactors%vb     = impactors%vb     * collision_system%vscale
         impactors%Lspin  = impactors%Lspin  * collision_system%Lscale
         do i = 1, 2
            impactors%rot(:,i) = impactors%Lspin(:,i) * (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3, i))
         end do
   
         fragments%mtot   = fragments%mtot   * collision_system%mscale
         fragments%mass   = fragments%mass   * collision_system%mscale
         fragments%radius = fragments%radius * collision_system%dscale
         fragments%rot    = fragments%rot    / collision_system%tscale
         fragments%rc     = fragments%rc     * collision_system%dscale
         fragments%vc     = fragments%vc     * collision_system%vscale
   
         do i = 1, fragments%nbody
            fragments%rb(:, i) = fragments%rc(:, i) + impactors%rbcom(:)
            fragments%vb(:, i) = fragments%vc(:, i) + impactors%vbcom(:)
         end do

         impactors%Qloss = impactors%Qloss * collision_system%Escale

         collision_system%Lorbit(:,:) = collision_system%Lorbit(:,:) * collision_system%Lscale
         collision_system%Lspin(:,:)  = collision_system%Lspin(:,:)  * collision_system%Lscale
         collision_system%Ltot(:,:)   = collision_system%Ltot(:,:)   * collision_system%Lscale
         collision_system%ke_orbit(:) = collision_system%ke_orbit(:) * collision_system%Escale
         collision_system%ke_spin(:)  = collision_system%ke_spin(:)  * collision_system%Escale
         collision_system%pe(:)       = collision_system%pe(:)       * collision_system%Escale
         collision_system%Etot(:)     = collision_system%Etot(:)     * collision_system%Escale
   
         collision_system%mscale = 1.0_DP
         collision_system%dscale = 1.0_DP
         collision_system%vscale = 1.0_DP
         collision_system%tscale = 1.0_DP
         collision_system%Lscale = 1.0_DP
         collision_system%Escale = 1.0_DP
      end associate
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)
   
      return
   end subroutine fraggle_set_original_scale_factors


end submodule s_fraggle_set