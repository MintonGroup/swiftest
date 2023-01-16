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
      real(DP) :: mfrag, mremaining, mtot, mcumul, G
      real(DP), dimension(:), allocatable :: mass
      real(DP), parameter :: BETA = 2.85_DP
      integer(I4B), parameter :: NFRAGMAX = 100  !! Maximum number of fragments that can be generated
      integer(I4B), parameter :: NFRAGMIN = 7 !! Minimum number of fragments that can be generated (set by the fraggle_generate algorithm for constraining momentum and energy)
      integer(I4B), parameter :: NFRAG_SIZE_MULTIPLIER = 3 !! Log-space scale factor that scales the number of fragments by the collisional system mass
      integer(I4B), parameter :: iMlr = 1
      integer(I4B), parameter :: iMslr = 2
      integer(I4B), parameter :: iMrem = 3
      integer(I4B), dimension(:), allocatable :: ind
      logical :: flipper
     
      associate(impactors => self%impactors, min_mfrag => self%min_mfrag)
         ! Get mass weighted mean of Ip and density
         volume(1:2) = 4._DP / 3._DP * PI * impactors%radius(1:2)**3
         mtot = sum(impactors%mass(:))
         G = impactors%Gmass(1) / impactors%mass(1)
         Ip_avg(:) = (impactors%mass(1) * impactors%Ip(:,1) + impactors%mass(2) * impactors%Ip(:,2)) / mtot

         if (impactors%mass(1) >= impactors%mass(2)) then
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
            associate(fragments => self%fragments)
               fragments%mass(1) = impactors%mass_dist(1)
               fragments%Gmass(1) = G * impactors%mass_dist(1)
               fragments%radius(1) = impactors%radius(jtarg)
               fragments%density(1) = impactors%mass_dist(1) / volume(jtarg)
               if (param%lrotation) fragments%Ip(:, 1) = impactors%Ip(:,1)
            end associate
            return
         case default
            write(*,*) "collision_util_set_mass_dist_fragments error: Unrecognized regime code",impactors%regime
         end select

         associate(fragments => self%fragments)
            fragments%mtot = mtot
            allocate(mass, mold=fragments%mass)

            ! Make the first two bins the same as the Mlr and Mslr values that came from collision_regime
            mass(1) = impactors%mass_dist(iMlr) 
            mass(2) = impactors%mass_dist(iMslr) 

            ! Distribute the remaining mass the 3:nfrag bodies following the model SFD given by slope BETA 
            mremaining = impactors%mass_dist(iMrem)
            do i = iMrem, nfrag
               mfrag = (1 + i - iMslr)**(-3._DP / BETA) * impactors%mass_dist(iMslr)
               mass(i) = mfrag
               mremaining = mremaining - mfrag
            end do

            ! If there is any residual mass (either positive or negative) we will distribute remaining mass proportionally among the the fragments
            if (mremaining < 0.0_DP) then ! If the remainder is negative, this means that that the number of fragments required by the SFD is smaller than our lower limit set by fraggle_generate. 
               istart = iMrem ! We will reduce the mass of the 3:nfrag bodies to prevent the second-largest fragment from going smaller
            else ! If the remainder is postiive, this means that the number of fragments required by the SFD is larger than our upper limit set by computational expediency. 
               istart = iMslr ! We will increase the mass of the 2:nfrag bodies to compensate, which ensures that the second largest fragment remains the second largest
            end if
            mfrag = 1._DP + mremaining / sum(mass(istart:nfrag))
            mass(istart:nfrag) = mass(istart:nfrag) * mfrag

            ! There may still be some small residual due to round-off error. If so, simply add it to the last bin of the mass distribution.
            mremaining = fragments%mtot - sum(mass(1:nfrag))
            mass(nfrag) = mass(nfrag) + mremaining

            ! Sort the distribution in descending order by mass so that the largest fragment is always the first
            call swiftest_util_sort(-mass, ind)
            call swiftest_util_sort_rearrange(mass, ind, nfrag)
            fragments%mass(:) = mass(:)
            deallocate(mass)

            fragments%Gmass(:) = G * fragments%mass(:)

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
            do concurrent(i = istart:nfrag)
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

         end associate


      end associate

      return
   end subroutine fraggle_util_set_mass_dist


end submodule s_fraggle_util
