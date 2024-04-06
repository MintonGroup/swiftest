! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle) s_fraggle_util
   use swiftest
contains

   module subroutine fraggle_util_restructure(self, nbody_system, param, lfailure)
      !! author: David A. Minton
      !!
      !! Restructures the fragment distribution after a failure to converge on a solution.
      implicit none
      ! Arguments
      class(collision_fraggle),     intent(inout) :: self         !! Fraggle collision system object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      logical,                      intent(out)   :: lfailure     !! Did the computation fail?
      ! Internals
      class(collision_fragments), allocatable :: new_fragments
      integer(I4B) :: i,nnew, nold
      real(DP) :: volume

      associate(old_fragments => self%fragments)
         lfailure = .false.
         ! Merge the second-largest fragment into the first and shuffle the masses up
         nold = old_fragments%nbody
         if (nold <= 3) then
            lfailure = .true.
            return
         end if
         nnew = nold - 1
         allocate(new_fragments, mold=old_fragments)
         call new_fragments%setup(nnew)

         ! Merge the first and second bodies
         volume = 4._DP / 3._DP * PI * sum(old_fragments%radius(1:2)**3)
         new_fragments%mass(1) = sum(old_fragments%mass(1:2))
         new_fragments%Gmass(1) =sum(old_fragments%Gmass(1:2))
         new_fragments%density(1) = new_fragments%mass(1) / volume
         new_fragments%radius(1) = (3._DP * volume / (4._DP * PI))**(THIRD)
#ifdef DOCONLOC
         do concurrent(i = 1:NDIM) shared(old_fragments, new_fragments)
#else
         do concurrent(i = 1:NDIM)
#endif
            new_fragments%Ip(i,1) = sum(old_fragments%mass(1:2) * old_fragments%Ip(i,1:2)) 
         end do
         new_fragments%Ip(:,1) = new_fragments%Ip(:,1) / new_fragments%mass(1)
         new_fragments%origin_body(1) = old_fragments%origin_body(1)
         
         ! Copy the n>2 old fragments to the n>1 new fragments
         new_fragments%mass(2:nnew) = old_fragments%mass(3:nold)
         new_fragments%Gmass(2:nnew) = old_fragments%Gmass(3:nold)
         new_fragments%density(2:nnew) = old_fragments%density(3:nold)
         new_fragments%radius(2:nnew) = old_fragments%radius(3:nold)
#ifdef DOCONLOC
         do concurrent(i = 1:NDIM) shared(old_fragments,new_fragments)
#else
         do concurrent(i = 1:NDIM)
#endif
            new_fragments%Ip(i,2:nnew) = old_fragments%Ip(i,3:nold) 
         end do
         new_fragments%origin_body(2:nnew) = old_fragments%origin_body(3:nold)

         new_fragments%mtot = old_fragments%mtot
      end associate

      deallocate(self%fragments)
      call move_alloc(new_fragments, self%fragments)

      call fraggle_generate_pos_vec(self, nbody_system, param, lfailure)
      if (lfailure) return
      call fraggle_generate_rot_vec(self, nbody_system, param)

      ! Increase the spatial size factor to get a less dense cloud
      self%fail_scale = self%fail_scale * 1.001_DP

   end subroutine fraggle_util_restructure

   module subroutine fraggle_util_set_mass_dist(self, param)
      !! author: David A. Minton
      !!
      !! Sets the mass of fragments based on the mass distribution returned by the regime calculation.
      !! This subroutine must be run after the the setup routine has been run on the fragments
      !!
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(collision_fraggle),   intent(inout) :: self  !! Fraggle collision system object
      class(swiftest_parameters), intent(in)    :: param !! Current Swiftest run configuration parameters
      ! Internals
      integer(I4B)              :: i, j, jproj, jtarg, istart, nfragmax
      real(DP), dimension(2)    :: volume
      real(DP), dimension(NDIM) :: Ip_avg
      real(DP) ::  mremaining, mtot, mcumul, G, mass_noise
      real(DP), dimension(:), allocatable :: mass
      real(DP)  :: beta 
      integer(I4B), parameter :: MASS_NOISE_FACTOR = 5  !! The number of digits of random noise that get added to the minimum mass value to prevent identical masses from being generated in a single run 
      integer(I4B), parameter :: NFRAGMAX_UNSCALED = 100  !! Maximum number of fragments that can be generated
 !! Minimum number of fragments that can be generated 
      integer(I4B), dimension(:), allocatable :: ind
      logical :: flipper
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      call ieee_set_halting_mode(IEEE_ALL,.false.)
     
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
            ! The remainder from the third bin will be distributed among nfrag-2 bodies. 

            !Add a small amount of noise to the last digits of the minimum mass value so that multiple fragments don't get generated with identical mass values
            call random_number(mass_noise)
            mass_noise = 1.0_DP + mass_noise * epsilon(1.0_DP) * 10**(MASS_NOISE_FACTOR)
            min_mfrag = (param%min_GMfrag / param%GU) * mass_noise
            
            mremaining = impactors%mass_dist(iMrem)
            nfrag = iMrem - 1 
            Mslr = impactors%mass_dist(iMslr)
            nfragmax = ceiling(NFRAGMAX_UNSCALED / param%nfrag_reduction)
            Mrat = max((mremaining + Mslr) / min_mfrag, 1.0_DP)
            nfrag = max(ceiling(nfragmax * (log(Mrat) + 1.0_DP)), NFRAGMIN)

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
            call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)
            return
         case default
            write(*,*) "fraggle_util_set_mass_dist_fragments error: Unrecognized regime code",impactors%regime
         end select

         associate(fragments => self%fragments)
            fragments%mtot = mtot
            allocate(mass, mold=fragments%mass)

            ! Make the first two bins the same as the Mlr and Mslr values that came from collision_regime
            mass(1) = impactors%mass_dist(iMlr) 
            mass(2) = impactors%mass_dist(iMslr)

            ! Recompute the slope parameter beta so that we span the complete size range
            if (abs(Mslr - min_mfrag) < epsilon(min_mfrag) * min_mfrag) Mslr = Mslr + impactors%mass_dist(iMrem) / nfrag
            mremaining = impactors%mass_dist(iMrem)

            ! The mass will be distributed in a power law where N>M=(M/Mslr)**(-beta/3)
            ! Use Brent's method solver to get the logspace slope of the mass function
            Mrat = (mremaining + Mslr) / Mslr
            beta = 3.0_DP
            call solve_roots(fraggle_util_sfd_function,beta)

            do i = iMrem,nfrag
               mass(i) = Mslr * (i-1)**(-beta/3.0_DP)
            end do

            ! There may still be some small residual due to round-off error. If so, simply add it to the last bin of the mass distribution.
            mremaining = fragments%mtot - sum(mass(1:nfrag))
            if (mremaining < 0.0_DP) then
               mass(iMlr) = mass(iMlr) + mremaining
            else
               mass(iMslr) = mass(iMslr) + mremaining
            end if

            ! Sort the distribution in descending order by mass so that the largest fragment is always the first
            call util_sort(-mass, ind)
            call util_sort_rearrange(mass, ind, nfrag)
            call move_alloc(mass, fragments%mass)

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
#ifdef DOCONLOC
            do concurrent(i = istart:nfrag) shared(fragments,Ip_avg)
#else
            do concurrent(i = istart:nfrag)
#endif
               fragments%Ip(:, i) = Ip_avg(:)
            end do

            ! For catastrophic impacts, we will assign each of the n>2 fragments to one of the two original bodies so that the fragment cloud occupies 
            ! roughly the same space as both original bodies. For all other disruption cases, we use body 2 as the center of the cloud.
            fragments%origin_body(1) = 1_I1B
            fragments%origin_body(2) = 2_I1B
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

      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)
      return

   end subroutine fraggle_util_set_mass_dist

   module function fraggle_util_sfd_function(x) result(y)
      !! author: David A. Minton
      !!
      !! Function to be used in the Brent's method solver to find the slope of the mass distribution
      implicit none
      ! Arguments
      real(DP), intent(in) :: x
      ! Result
      real(DP)             :: y
      ! Internals
      integer(I4B) :: i

      y = Mrat - 1.0_DP
      do i = iMrem, nfrag
         y = y - (i-1)**(-x/3.0_DP)
      end do
   end function fraggle_util_sfd_function

end submodule s_fraggle_util