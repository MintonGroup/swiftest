submodule (swiftest_classes) obl_implementations
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
      real(DP)     :: rinv2, t0, t1, t2, t3, fac1, fac2

      associate(n => self%nbody, aobl => self%aobl, xh => self%xh, j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, &
                msun => cb%Gmass, aoblcb => cb%aobl)
         do concurrent (i = 1:n) 
            rinv2 = irh(i)**2
            t0 = -msun * rinv2 * rinv2 * irh(i)
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
      end associate
      return

   end procedure obl_acc_body

   ! module procedure obl_pot_pl ! (self, cb, irh
   !    !! author: David A. Minton
   !    !!
   !    !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
   !    !!
   !    !!      Returned value does not include monopole term or terms higher than J4
   !    !!      Reference: MacMillan, W. D. 1958. The Theory of the Potential, (Dover Publications), 363.
   !    !!
   !    !! Adapted from David E. Kaufmann's Swifter routine: obl_pot.f90
   !    !! Adapted from Hal Levison's Swift routine obl_pot.f
   !    use swiftest
   !    implicit none
   !    integer(I4B) :: i
   !    real(DP)     :: rinv2, t0, t1, t2, t3, p2, p4

   !    associate(pl => self, npl => self%nbody, Mpl => self%Gmass, xh => self%xh, &
   !             j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, Mcb => cb%Gmass)
   !       oblpot = 0.0_DP
   !       !do concurrent (i = 1:npl) !shared(npl, oblpot, j2rp2, j4rp4, Mcb, Mpl, xh, irh) &
   !                                 !local(rinv2, t0, t1, t2, t3, p2, p4)
   !       !$omp simd
   !       do i = 1, npl
   !          rinv2 = irh(i)**2
   !          t0 = Mcb * Mpl(i) * rinv2 * irh(i)
   !          t1 = j2rp2
   !          t2 = xh(3, i)**2 * rinv2
   !          t3 = j4rp4 * rinv2
   !          p2 = 0.5_DP * (3 * t2 - 1.0_DP)
   !          p4 = 0.125_DP * ((35 * t2 - 30.0_DP) * t2 + 3.0_DP)
   !          oblpot = oblpot + t0 * (t1 * p2 + t3 * p4)
   !       end do
   !    end associate 
   !    return

   ! end procedure obl_pot_pl

end submodule obl_implementations
