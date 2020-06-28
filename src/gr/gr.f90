submodule(swiftest_classes) gr_implementations
contains
   module procedure gr_p4_body
      !! author: David A. Minton
      !!
      !! Position kick due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_whm_p4.f90
      use swiftest
      implicit none

      associate(n => self%nbody)
         select type(self)
         class is (whm_pl)
            associate(xj => self%xj, vj => self%vj)
               call p4_func(xj(1:n, :), vj(1:n, :))
            end associate
         class is (whm_tp)
            associate(xh => self%xh, vh => self%vh)
               call p4_func(xh(1:n, :), vh(1:n, :))
            end associate
         end select
      end associate
 
     return
     contains
         pure subroutine p4_func(x,v)
            implicit none
            real(DP), dimension(:,:), intent(inout) :: x
            real(DP), dimension(:,:), intent(in)    :: v
            real(DP), dimension(:,:), allocatable :: dr, vjmag2
            real(DP), dimension(:),   allocatable :: vmag2
            integer(I4B) :: i

            vmag2(:) = v(:, :) .dot. v(:, :)
            do concurrent (i = 1:NDIM)
               dr(:, i) = - config%inv_c2 * vmag2(:) * v(:,i)
               x(:, i)  =   x(:, i) + dr(:, i) * dt
            end do

            return
         end subroutine p4_func  

   end procedure gr_p4_body

end submodule gr_implementations