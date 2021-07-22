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
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: suma
      real(DP), dimension(:, :), allocatable       :: aj
      real(DP)                                     :: beta, rjmag4
      
      associate(pl => self, npl => self%nbody, inv_c2 => param%inv_c2)
         if (npl == 0) return
         allocate(aj, mold = pl%ah)
         do i = 1, npl
            rjmag4 = (dot_product(pl%xj(:, i), pl%xj(:, i)))**2
            beta   = -pl%muj(i)**2 * inv_c2 
            aj(:, i)  = 2 * beta * pl%xj(:, i) / rjmag4
         end do
         suma(:) = 0.0_DP
         pl%ah(:, 1) = pl%ah(:, 1) + aj(:, 1)
         do i = 2, npl
            suma(:) = suma(:) + pl%Gmass(i) * aj(:, i) / pl%eta(i)
            pl%ah(:, i) = pl%ah(:, i) + aj(:, i) + suma(:)
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
      
      associate(tp => self, ntp => self%nbody, inv_c2 => param%inv_c2)
         if (ntp == 0) return
         do i = 1, ntp
            rjmag4 = (dot_product(tp%xh(:, i), tp%xh(:, i)))**2
            beta = - tp%mu(i)**2 * inv_c2 
            tp%ah(:, i) = tp%ah(:, i) + beta * tp%xh(:, i) / rjmag4
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
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
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
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
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

end submodule s_whm_gr