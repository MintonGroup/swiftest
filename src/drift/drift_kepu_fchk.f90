submodule (swiftest_classes) s_drift_kepu_fchk
contains
   module procedure drift_kepu_fchk
   !! author: David A. Minton
   !!
   !! Computes the value of f, the function whose root we are trying to find in universal variables
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_kepu_fchk.f90
   !! Adapted from Martin Duncan's Swift routine drift_kepu_fchk.f
   use swiftest
   implicit none
   real(DP) :: x, c0, c1, c2, c3

! executable code
   x = s * s * alpha
   call drift_kepu_stumpff(x, c0, c1, c2, c3)
   c1 = c1 * s
   c2 = c2 * s * s
   c3 = c3 * s * s * s
   f = r0 * c1 + u * c2 + mu * c3 - dt

   return

   end procedure drift_kepu_fchk
end submodule s_drift_kepu_fchk

