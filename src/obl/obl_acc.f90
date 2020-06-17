submodule (swiftest_data_structures) s_obl_acc
contains
   module procedure obl_acc
   !! author: David A. Minton
   !!
   !! Compute the barycentric accelerations of planets due to the oblateness of the central body
   !!      Returned values do not include monopole term or terms higher than J4

   !! Adapted from David E. Kaufmann's Swifter modules: obl_acc.f90
   !! Adapted from Hal Levison's Swift routine obl_acc.f
use swiftest
implicit none
   integer(I4B)          :: i, npl
   real(DP)            :: rinv2, t0, t1, t2, t3, fac1, fac2, msun

! executable code
   msun = swiftest_pla%mass(1)
   npl = swiftest_pla%nbody
   do concurrent (i = 2:npl)
      rinv2 = irh(i)**2
      t0 = -msun * rinv2 * rinv2 *irh(i)
      t1 = 1.5_DP * j2rp2
      t2 = xh(3, i) * xh(3, i) * rinv2
      t3 = 1.875_DP * j4rp4 * rinv2
      fac1 = t0 * (t1 - t3 - (5 * t1 - (14.0_DP - 21 * t2) * t3) * t2)
      fac2 = 2 * t0 * (t1 - (2.0_DP - (14 * t2 / 3.0_DP)) * t3)
      aobl(:, i) = fac1*xh(:, i)
      aobl(3, i) = fac2*xh(3, i) + aobl(3, i)
   end do
   aobl(:, 1) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   do concurrent (i = 2:npl)
      aobl(:, 1) = aobl(:, 1) - swiftest_pla%mass(i) * aobl(:, i) / msun
   end do

   return

   end procedure obl_acc
end submodule s_obl_acc
