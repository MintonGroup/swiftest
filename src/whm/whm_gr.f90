!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(whm) s_whm_gr
   use swiftest
contains

   pure module subroutine whm_gr_kick_getacch_pl(self, param) 
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of massive bodies
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_kick_getacch.f90
      implicit none
      ! Arguments
      class(whm_pl),              intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: suma
      
      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody, inv_c2 => param%inv_c2)
         call swiftest_gr_kick_getacch(pl%muj, pl%xj, pl%lmask, npl, param%inv_c2, pl%agr) 
         suma(:) = 0.0_DP
         pl%ah(:, 1) = pl%ah(:, 1) + pl%agr(:, 1)
         do i = 2, npl
            suma(:) = suma(:) + pl%Gmass(i) * pl%agr(:, i) / pl%eta(i)
            pl%ah(:, i) = pl%ah(:, i) + pl%agr(:, i) + suma(:)
         end do
      end associate

      return
   end subroutine whm_gr_kick_getacch_pl


   pure module subroutine whm_gr_kick_getacch_tp(self, param)
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of test particles
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_kick_getacch.f90
      implicit none
      ! Arguments
      class(whm_tp),              intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
      ! Internals
      
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, inv_c2 => param%inv_c2)
         call swiftest_gr_kick_getacch(tp%mu, tp%rh, tp%lmask, ntp, param%inv_c2, tp%agr) 
         tp%ah(:,1:ntp) = tp%ah(:,1:ntp) + tp%agr(:,1:ntp)
      end associate

      return
   end subroutine whm_gr_kick_getacch_tp
   

   pure module subroutine whm_gr_p4_pl(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      implicit none
      ! Arguments
      class(whm_pl),                intent(inout) :: self   !! Swiftest particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Step size
      ! Internals
      integer(I4B) :: i, npl

      if (self%nbody == 0) return

      associate(xj => self%xj, vj => self%vj, lmask => self%lmask, inv_c2 => param%inv_c2)
         npl = self%nbody
#ifdef DOCONLOC
         do concurrent(i = 1:npl, lmask(i)) shared(lmask, inv_c2, xj, vj,dt)
#else
         do concurrent(i = 1:npl, lmask(i))
#endif
            call swiftest_gr_p4_pos_kick(inv_c2, xj(1,i), xj(2,i), xj(3,i), vj(1,i), vj(2,i), vj(3,i), dt)
         end do
      end associate
 
     return
   end subroutine whm_gr_p4_pl


   pure module subroutine whm_gr_p4_tp(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      implicit none
      ! Arguments
      class(whm_tp),                intent(inout) :: self  !! Swiftest particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt    !! Step size
      ! Internals
      integer(I4B) :: i, ntp

      associate(rh => self%rh, vh => self%vh, lmask => self%lmask, inv_c2 => param%inv_c2)
         ntp = self%nbody
         if (ntp == 0) return
#ifdef DOCONLOC
         do concurrent(i = 1:ntp, lmask(i)) shared(lmask, rh, vh, inv_c2, dt)
#else
         do concurrent(i = 1:ntp, lmask(i))
#endif
            call swiftest_gr_p4_pos_kick(inv_c2, rh(1,i), rh(2,i), rh(3,i), vh(1,i), vh(2,i), vh(3,i), dt)
         end do
      end associate
 
     return
   end subroutine whm_gr_p4_tp

end submodule s_whm_gr