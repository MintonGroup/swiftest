submodule (swiftest_data_structures) s_drift_kepu_stumpff
contains
   module procedure drift_kepu_stumpff
   !! author: David A. Minton
   !!
   !! Compute Stumpff functions needed for Kepler drift in universal variables
   !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 171 - 172.
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_kepu_stumpff.f90
   !! Adapted from Hal Levison's Swift routine drift_kepu_stumpff.f
   use swiftest
   implicit none
   integer(i4b) :: i, n
   real(DP)   :: xm

! executable code
   n = 0
   xm = 0.1_DP
   do while (abs(x) >= xm)
      n = n + 1
      x = x / 4.0_DP
   end do
   c2 = (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * &
           (1.0_DP - x / 182.0_DP) / 132.0_DP) / 90.0_DP) / 56.0_DP) /        &
            30.0_DP) / 12.0_DP) / 2.0_DP
   c3 = (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * &
           (1.0_DP - x / 210.0_DP) / 156.0_DP) / 110.0_DP) / 72.0_DP) /       &
            42.0_DP) / 20.0_DP ) / 6.0_DP
   c1 = 1.0_DP - x*c3
   c0 = 1.0_DP - x*c2
   if (n /= 0) then
      do i = n, 1, -1
         c3 = (c2 + c0 * c3) / 4.0_DP
         c2 = c1 * c1 / 2.0_DP
         c1 = c0 * c1
         c0 = 2 * c0 * c0 - 1.0_DP
         x = x * 4
      end do
   end if

   return

   end procedure drift_kepu_stumpff
end submodule s_drift_kepu_stumpff
