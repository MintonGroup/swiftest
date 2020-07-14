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
      integer(I4B) :: i, npl
      real(DP)     :: rinv2, t0, t1, t2, t3, fac1, fac2, msun

      associate(bd => self, n => self%nbody, aobl => self%aobl, &
                j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, msun => cb%Gmass)
         do concurrent (i = 1:n)
            rinv2 = irh(i)**2
            t0 = -msun * rinv2**2 * irh(i)
            t1 = 1.5_DP * j2rp2
            t2 = xh(3, i)**2 * rinv2
            t3 = 1.875_DP * j4rp4 * rinv2
            fac1 = t0 * (t1 - t3 - (5 * t1 - (14.0_DP - 21 * t2) * t3) * t2)
            fac2 = 2 * t0 * (t1 - (2.0_DP - (14 * t2 / 3.0_DP)) * t3)
            aobl(:, i) = fac1 * xh(:, i)
            aobl(3, i) = fac2 * xh(3, i) + aobl(3, i)
         end do
      end associate
      select type(pl => self)
      class is (swiftest_pl)
         do concurrent (i = 1:NDIM)
            cb%aobl(i) = -sum(pl%aobl(i, 1:pl%nbody) * pl%Gmass(1:pl%nbody)) / cb%Gmass
         end do
      end select
      return

   end procedure obl_acc_body

   module procedure obl_pot_pl ! (self, cb, irh
      !! author: David A. Minton
      !!
      !! Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
      !!
      !!      Returned value does not include monopole term or terms higher than J4
      !!      Reference: MacMillan, W. D. 1958. The Theory of the Potential, (Dover Publications), 363.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: obl_pot.f90
      !! Adapted from Hal Levison's Swift routine obl_pot.f
      use swiftest
      implicit none
      integer(I4B) :: i
      real(DP)     :: rinv2, t0, t1, t2, t3, p2, p4

      associate(pl => self, npl => self%nbody, mupl => self%Gmass, xh => self%xh, &
               j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, mu => cb%Gmass)
         oblpot = 0.0_DP
         do concurrent (i = 1:npl)
            rinv2 = irh(i)**2
            t0 = mu * mupl(i) * rinv2 * irh(i)
            t1 = j2rp2
            t2 = xh(3, i)**2 * rinv2
            t3 = j4rp4 * rinv2
            p2 = 0.5_DP * (3 * t2 - 1.0_DP)
            p4 = 0.125_DP * ((35 * t2 - 30.0_DP) * t2 + 3.0_DP)
            oblpot = oblpot + t0 * (t1 * p2 + t3 * p4)
         end do
      end associate 
      return

   end procedure obl_pot_pl

end submodule obl_implementations
