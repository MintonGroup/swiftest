submodule (swiftest_classes) s_obl_pot
contains
   module procedure obl_pot
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
   integer(I4B)          :: i, npl
   real(DP)            :: rinv2, t0, t1, t2, t3, p2, p4, mu

   oblpot = 0.0_DP
   mu = swiftest_plA%mass(1)
   npl = swiftest_plA%nbody
   do concurrent (i = 2:npl)
      rinv2 = irh(i)**2
      t0 = mu * swiftest_plA%mass(i) * rinv2 * irh(i)
      t1 = config%j2rp2
      t2 = xh(3, i) * xh(3, i) * rinv2
      t3 = config%j4rp4 * rinv2
      p2 = 0.5_DP * (3 * t2 - 1.0_DP)
      p4 = 0.125_DP * ((35 * t2 - 30.0_DP) * t2 + 3.0_DP)
      oblpot = oblpot + t0 * (t1 * p2 + t3 * p4)
   end do

   return

   end procedure obl_pot
end submodule s_obl_pot
