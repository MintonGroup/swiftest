submodule (util) s_util_hills
   use swiftest
contains
   module procedure util_hills
   !! author: David A. Minton
   !!
   !! Compute Hill sphere radii of planets
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: util_hills.f90
   !! Adapted from Hal Levison's Swift routine util_hills.f
   integer(I4B)  :: i
   real(DP)      :: msun, mp, mu, energy, ap, r, v2

   msun = swiftest_plA%mass(1)
   do i = 2, npl
      mp = swiftest_plA%mass(i)
      if (mp > 0.0_DP) then
         mu = msun + mp
         r = norm2(swiftest_plA%xh(:, i))
         v2 = dot_product(swiftest_plA%vh(:, i), swiftest_plA%vh(:, i))
         energy = 0.5_DP*v2 - mu/r
         ap = -0.5_DP*mu/energy
         swiftest_plA%rhill(i) = ap*(((mp/mu)/3.0_DP)**(1.0_DP/3.0_DP))
      else
         swiftest_plA%rhill(i) = 0.0_DP
      end if
   end do

   return

   end procedure util_hills
end submodule s_util_hills
