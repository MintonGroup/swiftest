submodule (swiftest_classes) s_orbel_scget
contains
   module procedure orbel_scget
   !! author: David A. Minton
   !!
   !! Efficiently compute the sine and cosine of an input angle
   !!      Input angle must be in radians
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: orbel_scget.f90
   !! Adapted from Hal Levison's Swift routine orbel_scget.f
   use swiftest
   implicit none
   integer(I4B) :: nper
   real(DP)   :: x

! executable code
   nper = angle / TWOPI
   x = angle - nper * TWOPI
   if (x < 0.0_DP) x = x + TWOPI
   sx = sin(x)
   cx = sqrt(1.0_DP - sx * sx)
   if ((x > PIBY2) .and. (x < PI3BY2)) cx = -cx

   return

   end procedure orbel_scget
end submodule s_orbel_scget
