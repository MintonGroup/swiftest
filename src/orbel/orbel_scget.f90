submodule (swiftest_data_structures) s_orbel_scget
contains
   module procedure orbel_scget
   !! author: David A. Minton
   !!
   !! Efficiently compute the sine and cosine of an input angle
   !!      Input angle must be in radians
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: orbel_scget.f90
   !! Adapted from Hal Levison's Swift routine orbel_scget.f
use swiftest
implicit none
   integer(i4b) :: nper
   real(DP)   :: x

! executable code
   nper = angle/twopi
   x = angle - nper*twopi
   if (x < 0.0_DP) x = x + twopi
   sx = sin(x)
   cx = sqrt(1.0_DP - sx*sx)
   if ((x > piby2) .and. (x < pi3by2)) cx = -cx

   return

   end procedure orbel_scget
end submodule s_orbel_scget
