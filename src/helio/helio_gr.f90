submodule(helio_classes) s_helio_gr
   use swiftest
contains

   module pure subroutine helio_gr_kick_getacch_pl(self, param) 
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
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: suma
      real(DP), dimension(:, :), allocatable       :: aj
      real(DP)                                     :: beta, rjmag4
      
      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         call gr_kick_getacch(pl%mu, pl%xh, pl%lmask, npl, param%inv_c2, pl%agr) 
         pl%ah(:,1:npl) = pl%ah(:,1:npl) + pl%agr(:,1:npl)
      end associate

      return
   end subroutine helio_gr_kick_getacch_pl


   module pure subroutine helio_gr_kick_getacch_tp(self, param)
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
      ! Internals
      integer(I4B)                                 :: i
      real(DP)                                     :: rjmag4, beta
      
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         call gr_kick_getacch(tp%mu, tp%xh, tp%lmask, ntp, param%inv_c2, tp%agr) 
         tp%ah(:,1:ntp) = tp%ah(:,1:ntp) + tp%agr(:,1:ntp)
      end associate

      return
   end subroutine helio_gr_kick_getacch_tp
   

   module pure subroutine helio_gr_p4_pl(self, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_helio_p4.f90
      implicit none
      ! Arguments
      class(helio_pl),            intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: dt     !! Step size
      ! Internals
      integer(I4B)                                 :: i

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         do concurrent(i = 1:npl, pl%lmask(i))
            call gr_p4_pos_kick(param, pl%xh(:, i), pl%vb(:, i), dt)
         end do
      end associate
 
     return
   end subroutine helio_gr_p4_pl

   module pure subroutine helio_gr_p4_tp(self, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_helio_p4.f90
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self  !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: dt    !! Step size
      ! Internals
      integer(I4B)                              :: i

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         do concurrent(i = 1:ntp, tp%lmask(i))
            call gr_p4_pos_kick(param, tp%xh(:, i), tp%vb(:, i), dt)
         end do
      end associate
 
     return
   end subroutine helio_gr_p4_tp

end submodule s_helio_gr