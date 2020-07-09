submodule(whm_classes) whm_gr_implementations
contains
   module procedure whm_gr_getacch_pl 
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of massive bodies
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_getacch.f90
      use swiftest
      implicit none

      integer(I4B)                   :: i
      real(DP), dimension(NDIM)      :: suma
      real(DP), dimension(:, :), allocatable :: aj
      real(DP), dimension(:), allocatable :: rjmag4, beta
      
      associate(n => self%nbody, msun => cb%Gmass, mu => self%muj, c2 => config%inv_c2, &
         ah => self%ah, xj => self%xj, GMpl => self%Gmass, eta => self%eta)
         if (n == 0) return
         allocate(rjmag4(n))
         allocate(beta(n))
         allocate(aj, mold = ah)
         rjmag4(:) = (xj(1:n, :) .dot. xj(1:n, :))**2
         beta(:)   = - mu(1:n)**2 * c2 
         do concurrent (i = 1:n)
            aj(i, :) = beta(i) * xj(i,:) / rjmag4(i)
         end do
         suma(:) = 0.0_DP
         ah(1, :) = ah(1, :) + aj(1, :)
         do i = 2, n
            suma(:) = suma(:) + GMpl(i) * aj(i, :) / eta(i)
            ah(i, :) = ah(i, :) + aj(i, :) + suma(:)
         end do
      end associate
      return
   end procedure whm_gr_getacch_pl

   module procedure whm_gr_getacch_tp
      !! author: David A. Minton
      !!
      !! Compute relativisitic accelerations of test particles
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_getacch.f90
      use swiftest
      implicit none

      integer(I4B)                   :: i
      real(DP), dimension(NDIM)      :: suma
      real(DP), dimension(:, :), allocatable :: aj
      real(DP), dimension(:), allocatable :: rjmag4, beta
      
      associate(n => self%nbody, msun => cb%Gmass, mu => self%mu,& 
         c2 => config%inv_c2, ah => self%ah, xh => self%xh, status => self%status)
         if (n == 0) return
         allocate(rjmag4(n))
         allocate(beta(n))
         rjmag4(:) = (xh(1:n, :) .dot. xh(1:n, :))**2
         beta(:) = - mu(1:n)**2 * c2 
         ah(1, :) = ah(1, :) + aj(2, :)
         do concurrent(i = 2:n, status(i) == ACTIVE)
            ah(i, :) = ah(i, :) + aj(i, :) 
         end do
      end associate
      return
   end procedure whm_gr_getacch_tp

   module procedure whm_gr_p4_pl
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      use swiftest
      implicit none
      integer(I4B) :: i

      associate(n => self%nbody, xj => self%xh, vj => self%vj, status => self%status, c2 => config%inv_c2)
         if (n == 0) return
         do concurrent (i = 1:n, status(i) == ACTIVE)
            call p4_func(xj(i, :), vj(i, :), dt, c2)
         end do
      end associate
 
     return
   end procedure whm_gr_p4_pl

   module procedure whm_gr_p4_tp
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      use swiftest
      implicit none
      integer(I4B) :: i

      associate(n => self%nbody, xh => self%xh, vh => self%vh, status => self%status, c2 => config%inv_c2)
         if (n == 0) return
         do concurrent (i = 1:n, status(i) == ACTIVE)
            call p4_func(xh(i, :), vh(i, :), dt, c2)
         end do
      end associate
 
     return
   end procedure whm_gr_p4_tp

   pure subroutine p4_func(x, v, dt, c2)
      use swiftest
      implicit none
      real(DP), dimension(:), intent(inout) :: x
      real(DP), dimension(:), intent(in)    :: v
      real(DP),               intent(in)    :: dt, c2 
      real(DP), dimension(NDIM)             :: dr
      real(DP)                              :: vmag2
      integer(I4B) :: i

      vmag2 = v .dot. v 
      dr(:) = -c2 * vmag2 * v(:)
      x(:) = x(:) + dr(:) * dt

      return
   end subroutine p4_func  

end submodule whm_gr_implementations