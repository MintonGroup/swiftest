submodule (symba) s_symba_helio_drift_tp
contains
   module procedure symba_helio_drift_tp
   !! author: David A. Minton
   !!
   !! Loop through test particles and call Danby drift routine
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_helio_drift_tp.f90
   !! Adapted from Hal Levison's Swift routine symba5_helio_drift.f
   use swiftest
   implicit none
   integer(I4B)          :: i, iflag

   associate(symba_tpA%xh => xh, symba_tpA%vb => vb)
   do i = 1, ntp
      if ((symba_tpA%levelg(i) == irec) .and. (symba_tpA%status(i) == ACTIVE)) then
         call drift_one(mu, xh(1,i), xh(2,i), xh(3,i), vb(1,i), vb(2,i), vb(3,i), dt, iflag)
         if (iflag /= 0) then
            symba_tpA%status(i) = DISCARDED_DRIFTERR
            write(*, *) "particle ", symba_tpA%name(i), " lost due to error in danby drift"
         end if
      end if
   end do

   return

   end procedure symba_helio_drift_tp
end submodule s_symba_helio_drift_tp
