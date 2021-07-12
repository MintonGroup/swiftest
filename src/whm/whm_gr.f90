submodule(whm_classes) s_whm_gr
   use swiftest
contains
   module subroutine whm_gr_getacch_pl(self, param) !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of massive bodies
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_getacch.f90
      implicit none
      ! Arguments
      class(whm_pl),              intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: suma
      real(DP), dimension(:, :), allocatable       :: aj
      real(DP)                                     :: beta, rjmag4
      
      associate(n => self%nbody, mu => self%muj, c2 => param%inv_c2, &
         ah => self%ah, xj => self%xj, GMpl => self%Gmass, eta => self%eta)
         if (n == 0) return
         allocate(aj, mold = ah)
         do i = 1, n
            rjmag4 = (dot_product(xj(:, i), xj(:, i)))**2
            beta   = - mu(i)**2 * c2 
            aj(:, i)  = 2 * beta * xj(:, i) / rjmag4
         end do
         suma(:) = 0.0_DP
         ah(:, 1) = ah(:, 1) + aj(:, 1)
         do i = 2, n
            suma(:) = suma(:) + GMpl(i) * aj(:, i) / eta(i)
            ah(:, i) = ah(:, i) + aj(:, i) + suma(:)
         end do
      end associate
      return
   end subroutine whm_gr_getacch_pl

   module subroutine whm_gr_getacch_tp(self, param)
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of test particles
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_getacch.f90
      implicit none
      ! Arguments
      class(whm_tp),              intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP)                                     :: rjmag4, beta
      
      associate(n => self%nbody, mu => self%mu,& 
         c2 => param%inv_c2, ah => self%ah, xh => self%xh, status => self%status)
         if (n == 0) return
         do i = 1, n
            rjmag4 = (dot_product(xh(:, i), xh(:, i)))**2
            beta = - mu(i)**2 * c2 
            ah(:, i) = ah(:, i) + beta * xh(:, i) / rjmag4
         end do
      end associate
      return
   end subroutine whm_gr_getacch_tp

   module pure subroutine whm_gr_p4_pl(self, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      implicit none
      ! Arguments
      class(whm_pl),              intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of on parameters 
      real(DP),                   intent(in)    :: dt     !! Step size
      ! Internals
      integer(I4B)                                 :: i

      associate(pl => self, npl => self%nbody)
         if (npl == 0) return
         do i = 1, npl
            call gr_p4_pos_kick(param, pl%xj(:, i), pl%vj(:, i), dt)
         end do
      end associate
 
     return
   end subroutine whm_gr_p4_pl

   module pure subroutine whm_gr_p4_tp(self, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      implicit none
      ! Arguments
      class(whm_tp),              intent(inout) :: self  !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of on parameters 
      real(DP),                   intent(in)    :: dt    !! Step size
      ! Internals
      integer(I4B)                              :: i

      associate(tp => self, ntp => self%nbody)
         if (ntp == 0) return
         do i = 1, ntp
            call gr_p4_pos_kick(param, tp%xh(:, i), tp%vh(:, i), dt)
         end do
      end associate
 
     return
   end subroutine whm_gr_p4_tp

   module pure subroutine whm_gr_pv2vh_pl(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for massive bodies
      !! in a WHM object
      implicit none
      ! Arguments
      class(whm_pl),              intent(inout) :: self  !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of on parameters 
      ! Internals
      integer(I4B)                              :: i
      real(DP), dimension(:,:), allocatable     :: vh    !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(pl => self, npl => self%nbody)
         if (npl == 0) return
         allocate(vh, mold = pl%vh)
         do i = 1, npl
            call gr_pseudovel2vel(param, pl%mu(i), pl%xh(:, i), pl%vh(:, i), vh(:, i))
         end do
         call move_alloc(vh, pl%vh)
      end associate
      return
   end subroutine whm_gr_pv2vh_pl

   module pure subroutine whm_gr_pv2vh_tp(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from pseudovelocity to heliocentric velocity for test particles bodies
      !! in a WHM object
      implicit none
      ! Arguments
      class(whm_tp),              intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of on parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(:,:), allocatable        :: vh !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(tp => self, ntp => self%nbody)
         if (ntp == 0) return
         allocate(vh, mold = tp%vh)
         do i = 1, ntp
            call gr_pseudovel2vel(param, tp%mu(i), tp%xh(:, i), tp%vh(:, i), vh(:, i))
         end do
         call move_alloc(vh, tp%vh)
      end associate
      return
   end subroutine whm_gr_pv2vh_tp

   module pure subroutine whm_gr_vh2pv_pl(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for massive bodies
      !! in a WHM object
      implicit none
      ! Arguments
      class(whm_pl),              intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of on parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(:,:), allocatable        :: pv !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(pl => self, npl => self%nbody)
         if (npl == 0) return
         allocate(pv, mold = pl%vh)
         do i = 1, npl
            call gr_vel2pseudovel(param, pl%mu(i), pl%xh(:, i), pl%vh(:, i), pv(:, i))
         end do
         call move_alloc(pv, pl%vh)
      end associate
      return
   end subroutine whm_gr_vh2pv_pl

   module pure subroutine whm_gr_vh2pv_tp(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper function that converts from heliocentric velocity to pseudovelocity for teset particles
      !! in a WHM object
      implicit none
      ! Arguments
      class(whm_tp),              intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of on parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(:,:), allocatable        :: pv !! Temporary holder of pseudovelocity for in-place conversion
      
      associate(tp => self, ntp => self%nbody)
         if (ntp == 0) return
         allocate(pv, mold = tp%vh)
         do i = 1, ntp
            call gr_vel2pseudovel(param, tp%mu(i), tp%xh(:, i), tp%vh(:, i), pv(:, i))
         end do
         call move_alloc(pv, tp%vh)
      end associate
      return
   end subroutine whm_gr_vh2pv_tp

end submodule s_whm_gr