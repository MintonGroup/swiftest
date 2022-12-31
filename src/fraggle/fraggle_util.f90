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
      self%r_unit(:,:) = 0.0_DP
      self%t_unit(:,:) = 0.0_DP
      self%n_unit(:,:) = 0.0_DP

      self%rmag(:) = 0.0_DP
      self%rotmag(:) = 0.0_DP

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


   module subroutine fraggle_util_set_mass_dist(self, param)
      !! author: David A. Minton
      !!
      !! Sets the mass of fragments based on the mass distribution returned by the regime calculation.
      !! This subroutine must be run after the the setup routine has been run on the fragments
      !!
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self  !! Fraggle collision system object
      class(base_parameters),             intent(in)    :: param !! Current Swiftest run configuration parameters
      ! Internals
      integer(I4B)              :: i, j, jproj, jtarg, nfrag, istart
      real(DP), dimension(2)    :: volume
      real(DP), dimension(NDIM) :: Ip_avg
      real(DP) :: mfrag, mremaining, min_mfrag, mtot, mcumul
      real(DP), parameter :: BETA = 2.85_DP
      integer(I4B), parameter :: NFRAGMAX = 100  !! Maximum number of fragments that can be generated
      integer(I4B), parameter :: NFRAGMIN = 7 !! Minimum number of fragments that can be generated (set by the fraggle_generate algorithm for constraining momentum and energy)
      integer(I4B), parameter :: NFRAG_SIZE_MULTIPLIER = 3 !! Log-space scale factor that scales the number of fragments by the collisional system mass
      integer(I4B), parameter :: iMlr = 1
      integer(I4B), parameter :: iMslr = 2
      integer(I4B), parameter :: iMrem = 3
      logical :: flipper
     
      associate(impactors => self%impactors)
         ! Get mass weighted mean of Ip and density
         volume(1:2) = 4._DP / 3._DP * PI * impactors%radius(1:2)**3
         mtot = sum(impactors%mass(:))
         Ip_avg(:) = (impactors%mass(1) * impactors%Ip(:,1) + impactors%mass(2) * impactors%Ip(:,2)) / mtot

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
            class is (swiftest_parameters)
               min_mfrag = (param%min_GMfrag / param%GU) 
               ! The number of fragments we generate is bracked by the minimum required by fraggle_generate (7) and the 
               ! maximum set by the NFRAG_SIZE_MULTIPLIER which limits the total number of fragments to prevent the nbody
               ! code from getting an overwhelmingly large number of fragments
               nfrag = ceiling(NFRAG_SIZE_MULTIPLIER  * log(mtot / min_mfrag))
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
            call self%setup_fragments(nfrag)

         case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE) 

            call self%setup_fragments(1)
            select type(fragments => self%fragments)
            class is (collision_fragments(*))
               fragments%mass(1) = impactors%mass_dist(1)
               fragments%radius(1) = impactors%radius(jtarg)
               fragments%density(1) = impactors%mass_dist(1) / volume(jtarg)
               if (param%lrotation) fragments%Ip(:, 1) = impactors%Ip(:,1)
            end select
            return
         case default
            write(*,*) "collision_util_set_mass_dist_fragments error: Unrecognized regime code",impactors%regime
         end select

         select type(fragments => self%fragments)
         class is (collision_fragments(*))
            fragments%mtot = mtot

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

            ! For catastrophic impacts, we will assign each of the n>2 fragments to one of the two original bodies so that the fragment cloud occupies 
            ! roughly the same space as both original bodies. For all other disruption cases, we use body 2 as the center of the cloud.
               fragments%origin_body(1) = 1
               fragments%origin_body(2) = 2
               if (impactors%regime == COLLRESOLVE_REGIME_SUPERCATASTROPHIC) then
                  mcumul = fragments%mass(1)
                  flipper = .true.
                  j = 2
                  do i = 1, nfrag
                     if (flipper .and. (mcumul < impactors%mass(1))) then
                        flipper = .false.
                        j = 1
                     else
                        j = 2
                        flipper = .true.
                     end if
                     fragments%origin_body(i) = j
                  end do
               else
                  fragments%origin_body(3:nfrag) = 2
               end if

         end select


      end associate

      return
   end subroutine fraggle_util_set_mass_dist


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
         impactors%rbcom(:)       = impactors%rbcom(:)       / collider%dscale
         impactors%vbcom(:)       = impactors%vbcom(:)       / collider%vscale
         impactors%rbimp(:)       = impactors%rbimp(:)       / collider%dscale
         impactors%rb(:,:)        = impactors%rb(:,:)        / collider%dscale
         impactors%vb(:,:)        = impactors%vb(:,:)        / collider%vscale
         impactors%rc(:,:)        = impactors%rc(:,:)        / collider%dscale
         impactors%vc(:,:)        = impactors%vc(:,:)        / collider%vscale
         impactors%mass(:)        = impactors%mass(:)        / collider%mscale
         impactors%Gmass(:)       = impactors%Gmass(:)       / (collider%dscale**3/collider%tscale**2)
         impactors%Mcb            = impactors%Mcb            / collider%mscale
         impactors%radius(:)      = impactors%radius(:)      / collider%dscale
         impactors%Lspin(:,:)     = impactors%Lspin(:,:)     / collider%Lscale
         impactors%Lorbit(:,:)    = impactors%Lorbit(:,:)    / collider%Lscale
         impactors%bounce_unit(:) = impactors%bounce_unit(:) / collider%vscale

         do i = 1, 2
            impactors%rot(:,i) = impactors%Lspin(:,i) / (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i))
         end do

         fragments%mtot   = fragments%mtot   / collider%mscale
         fragments%mass   = fragments%mass   / collider%mscale
         fragments%radius = fragments%radius / collider%dscale
         impactors%Qloss  = impactors%Qloss  / collider%Escale
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
         impactors%rbcom(:)       = impactors%rbcom(:) * collider%dscale
         impactors%vbcom(:)       = impactors%vbcom(:) * collider%vscale
         impactors%rbimp(:)       = impactors%rbimp(:) * collider%dscale
         impactors%bounce_unit(:) = impactors%bounce_unit(:) * collider%vscale

         impactors%mass      = impactors%mass      * collider%mscale
         impactors%Gmass(:)  = impactors%Gmass(:)  * (collider%dscale**3/collider%tscale**2)
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
            impactors%rot(:,i) = impactors%Lspin(:,i) * (impactors%mass(i) * impactors%radius(i)**2 * impactors%Ip(3,i))
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


end submodule s_fraggle_util
