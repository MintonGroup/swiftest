! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(helio) s_helio_gr
   use swiftest
contains

   pure module subroutine helio_gr_kick_getacch_pl(self, param) 
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of massive bodies
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_kick_getacch.f90
      implicit none
      ! Arguments
      class(helio_pl),            intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      
      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         call swiftest_gr_kick_getacch(pl%mu, pl%rh, pl%lmask, npl, param%inv_c2, pl%agr) 
         pl%ah(:,1:npl) = pl%ah(:,1:npl) + pl%agr(:,1:npl)
      end associate

      return
   end subroutine helio_gr_kick_getacch_pl


   pure module subroutine helio_gr_kick_getacch_tp(self, param)
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of test particles
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_helio_kick_getacch.f90
      implicit none
      ! Arguments
      class(helio_tp),            intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
      
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         call swiftest_gr_kick_getacch(tp%mu, tp%rh, tp%lmask, ntp, param%inv_c2, tp%agr) 
         tp%ah(:,1:ntp) = tp%ah(:,1:ntp) + tp%agr(:,1:ntp)
      end associate

      return
   end subroutine helio_gr_kick_getacch_tp
   

   pure module subroutine helio_gr_p4_pl(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_helio_p4.f90
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Swiftest particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Step size
      ! Internals
      integer(I4B) :: i, npl

      if (self%nbody == 0) return

      associate(lmask => self%lmask, rh => self%rh, vb => self%vb, inv_c2 => param%inv_c2)
         npl = self%nbody
#ifdef DOCONLOC
         do concurrent(i = 1:npl, lmask(i)) shared(inv_c2, lmask, rh, vb, dt)
#else
         do concurrent(i = 1:npl, lmask(i))
#endif
            call swiftest_gr_p4_pos_kick(inv_c2, rh(1,i), rh(2,i), rh(3,i), vb(1,i), vb(2,i), vb(3,i), dt)
         end do
      end associate
 
     return
   end subroutine helio_gr_p4_pl


   pure module subroutine helio_gr_p4_tp(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_helio_p4.f90
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self   !! Swiftest particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Step size
      ! Internals
      integer(I4B) :: i, ntp

      if (self%nbody == 0) return

      associate(rh => self%rh, vb => self%vb, lmask => self%lmask, inv_c2 => param%inv_c2)
         ntp = self%nbody
#ifdef DOCONLOC
         do concurrent(i = 1:ntp, lmask(i)) shared(inv_c2, lmask, rh, vb, dt)
#else
         do concurrent(i = 1:ntp, lmask(i))
#endif
            call swiftest_gr_p4_pos_kick(inv_c2, rh(1,i), rh(2,i), rh(3,i), vb(1,i), vb(2,i), vb(3,i), dt)
         end do
      end associate
 
     return
   end subroutine helio_gr_p4_tp

end submodule s_helio_gr