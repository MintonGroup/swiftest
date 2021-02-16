submodule (swiftest_classes) s_obl
contains
   module procedure obl_acc_body 
      !! author: David A. Minton
      !!
      !! Compute the barycentric accelerations of bodies due to the oblateness of the central body
      !!      Returned values do not include monopole term or terms higher than J4

      !! Adapted from David E. Kaufmann's Swifter routine: obl_acc.f90 and obl_acc_tp.f90
      !! Adapted from Hal Levison's Swift routine obl_acc.f and obl_acc_tp.f
      use swiftest
      implicit none
      integer(I4B) :: i
      real(DP)     :: irh, rinv2, t0, t1, t2, t3, fac1, fac2

      associate(n => self%nbody, aobl => self%aobl, xh => self%xh, j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, &
                msun => cb%Gmass, aoblcb => cb%aobl, ah => self%ah)
         do i = 1, n 
            irh = 1.0_DP / norm2(xh(:,i))
            rinv2 = irh**2
            t0 = -msun * rinv2**2 * irh
            t1 = 1.5_DP * j2rp2
            t2 = xh(3, i) * xh(3, i) * rinv2
            t3 = 1.875_DP * j4rp4 * rinv2
            fac1 = t0 * (t1 - t3 - (5.0_DP * t1 - (14.0_DP - 21.0_DP * t2) * t3) * t2)
            fac2 = 2.0_DP * t0 * (t1 - (2.0_DP - (14.0_DP * t2 / 3.0_DP)) * t3)
            aobl(:, i) = fac1 * xh(:, i)
            aobl(3, i) = fac2 * xh(3, i) + aobl(3, i)
         end do
         select type(self)
         class is (swiftest_pl)
            associate(Mpl => self%Gmass)
               do i = 1, NDIM
                  aoblcb(i) = -sum(Mpl(1:n) * aobl(i, 1:n)) / msun
               end do
            end associate
         end select

         do i = 1, NDIM
            ah(i, 1:n) = ah(i, 1:n) + aobl(i, 1:n) - aoblcb(i)
         end do
      end associate
      return

   end procedure obl_acc_body

end submodule s_obl
