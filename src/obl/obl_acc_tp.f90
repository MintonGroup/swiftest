submodule (swiftest_classes) s_obl_acc_tp
contains
   module procedure obl_acc_tp
   !! author: David A. Minton
   !!
   !! Compute the barycentric accelerations of test particles due to the oblateness of the central body
   !!      Returned values do not include monopole term or terms higher than J4
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: obl_acc_tp.f90
   !! Adapted from Hal Levison's Swift routine obl_acc_tp.f
   use swiftest
   implicit none
   integer(I4B) :: i
   real(DP)   :: rinv2, t0, t1, t2, t3, fac1, fac2

   do concurrent (i = 1:ntp)
      rinv2 = irht(i)**2
      t0 = -msun * rinv2 * rinv2 * irht(i)
      t1 = 1.5_DP * config%j2rp2
      t2 = xht(3, i) * xht(3, i) * rinv2
      t3 = 1.875_DP * config%j4rp4 * rinv2
      fac1 = t0 * (t1 - t3 - (5 * t1 - (14.0_DP - 21 * t2) * t3) * t2)
      fac2 = 2 * t0 * (t1 - (2.0_DP - (14 * t2 / 3.0_DP)) * t3)
      aoblt(:, i) = fac1 * xht(:, i)
      aoblt(3, i) = aoblt(3, i) + fac2*xht(3, i)
   end do

   return

   end procedure obl_acc_tp
end submodule s_obl_acc_tp
