!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(symba) s_symba_gr
   use swiftest
contains

   pure module subroutine symba_gr_p4_pl(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_symba_p4.f90
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Step size

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         select type(nbody_system)
         class is (symba_nbody_system)
            pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE .and. pl%levelg(1:npl) == nbody_system%irec
            call helio_gr_p4_pl(pl, nbody_system, param, dt)
            pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE 
         end select
      end associate

   return
   end subroutine symba_gr_p4_pl


   pure module subroutine symba_gr_p4_tp(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_symba_p4.f90
      implicit none
      ! Arguments
      class(symba_tp),              intent(inout) :: self   !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Step size

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         select type(nbody_system)
         class is (symba_nbody_system)
            tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE .and. tp%levelg(1:ntp) == nbody_system%irec
            call helio_gr_p4_tp(tp, nbody_system, param, dt)
            tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE 
         end select
      end associate

   return
   end subroutine symba_gr_p4_tp

end submodule s_symba_gr
