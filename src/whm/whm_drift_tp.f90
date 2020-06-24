submodule(whm) s_whm_drift_tp
contains
   module procedure whm_drift_tp(ntp, whm_tp1p, mu, dt, c2)
   !! author: David A. Minton
   !!
   !! Loop through test particles and call Danby drift routine
   !!
   !! Adapted from Hal Levison's Swift routine drift_tp.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f90
   use swiftest
   implicit none
   integer(I4B)          :: i, iflag
   type(swifter_tp), pointer :: swifter_tpp
   type(whm_tp), pointer   :: whm_tpp
   real(DP)            :: dtp, energy, vmag2, rmag

! executable code
   whm_tpp => whm_tp1p
   do i = 1, ntp
      swifter_tpp => whm_tpp%swifter
      if (swifter_tpp%status == active) then
         rmag = sqrt(dot_product(swifter_tpp%xh(:), swifter_tpp%xh(:)))
         vmag2 = dot_product(swifter_tpp%vh(:), swifter_tpp%vh(:))
         energy = 0.5_DP*vmag2 - mu/rmag
         dtp = dt * (1.0_DP + 3 * c2 * energy)
         call drift_one(mu, swifter_tpp%xh(:), swifter_tpp%vh(:), dt, iflag)
         if (iflag /= 0) then
            swifter_tpp%status = discarded_drifterr
            write(*, *) "particle ", swifter_tpp%id, " lost due to error in danby drift"
         end if
      end if
      whm_tpp => whm_tpp%nextp
   end do

   return

   end procedure whm_drift_tp
end submodule s_whm_drift_tp
