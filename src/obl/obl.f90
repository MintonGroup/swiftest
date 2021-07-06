submodule (swiftest_classes) s_obl
   use swiftest
contains
   module subroutine obl_acc_body(self, cb)
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      !!      Returned values do not include monopole term or terms higher than J4

      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self !! Swiftest generic body object
      class(swiftest_cb),   intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B) :: i
      real(DP)     :: r2, irh, rinv2, t0, t1, t2, t3, fac1, fac2

      associate(n => self%nbody, xh => self%xh, vh => self%vh, ah => self%ah)
         do i = 1, n 
            r2 = dot_product(self%xh(:, i), self%xh(:, i))
            irh = 1.0_DP / sqrt(r2)
            rinv2 = irh**2
            t0 = -cb%Gmass * rinv2 * rinv2 * irh
            t1 = 1.5_DP * cb%j2rp2
            t2 = self%xh(3, i) * self%xh(3, i) * rinv2
            t3 = 1.875_DP * cb%j4rp4 * rinv2
            fac1 = t0 * (t1 - t3 - (5 * t1 - (14.0_DP - 21.0_DP * t2) * t3) * t2)
            fac2 = 2 * t0 * (t1 - (2.0_DP - (14.0_DP * t2 / 3.0_DP)) * t3)
            self%aobl(:, i) = fac1 * self%xh(:, i)
            self%aobl(3, i) = fac2 * self%xh(3, i) + self%aobl(3, i)
         end do
         select type(self)
         class is (swiftest_pl)
            do i = 1, NDIM
               cb%aobl(i) = -sum(self%Gmass(1:n) * self%aobl(i, 1:n)) / cb%Gmass
            end do
         end select

         do i = 1, NDIM
            self%ah(i, 1:n) = self%ah(i, 1:n) + self%aobl(i, 1:n) - cb%aobl(i)
         end do
      end associate
      return

   end subroutine obl_acc_body

end submodule s_obl
