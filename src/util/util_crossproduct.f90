submodule (util) s_util_crossproduct
contains
   module procedure util_crossproduct
   !! author: David A. Minton, Carlisle Wishard, and Jennifer L. L. Pouplin
   !!
   !! Compute Hill sphere radii of planets
   !!***************************************************************************
   use swiftest
   implicit none

   ans(1) = ar1(2) * ar2(3) - ar1(3) * ar2(2)
   ans(2) = ar1(3) * ar2(1) - ar1(1) * ar2(3)
   ans(3) = ar1(1) * ar2(2) - ar1(2) * ar2(1)
   
   return

   end procedure util_crossproduct
end submodule s_util_crossproduct
