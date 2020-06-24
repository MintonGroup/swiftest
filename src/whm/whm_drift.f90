submodule(whm) s_whm_drift
contains
   module procedure whm_drift
   !! author: David A. Minton
   !!
   !! Loop through planets and call Danby drift routine
   !!
   !! Adapted from Hal Levison's Swift routine drift.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_drift.f90
   use swiftest
   implicit none
   integer(I4B)          :: i, iflag
   real(DP)            :: etajm1, etaj, mu, msun
   type(swifter_pl), pointer :: swifter_plp
   type(whm_pl), pointer   :: whm_plp
   real(DP)            :: dtp, energy, vmag2, rmag

! executable code
   whm_plp => whm_pl1p
   swifter_plp => whm_plp%swifter
   msun = swifter_plp%mass
   etajm1 = msun
   do i = 2, npl
      whm_plp => whm_plp%nextp
      swifter_plp => whm_plp%swifter
      etaj = etajm1 + swifter_plp%mass
      mu = msun*etaj/etajm1
      rmag = sqrt(dot_product(whm_plp%xj(:), whm_plp%xj(:)))
      vmag2 = dot_product(whm_plp%vj(:), whm_plp%vj(:))
      energy = 0.5_DP*vmag2 - mu/rmag
      dtp = dt * (1.0_DP + 3 * c2 * energy)
      call drift_one(mu, whm_plp%xj(:), whm_plp%vj(:), dtp, iflag)
      if (iflag /= 0) then
         write(*, *) " planet ", swifter_plp%id, " is lost!!!!!!!!!!"
         write(*, *) mu, dt
         write(*, *) whm_plp%xj(:)
         write(*, *) whm_plp%vj(:)
         write(*, *) " stopping "
         call util_exit(failure)
      end if
      etajm1 = etaj
   end do

   return

   end procedure whm_drift
end submodule s_whm_drift
