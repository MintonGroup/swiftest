submodule (helio) s_helio_lindrift_tp
contains
module procedure helio_lindrift_tp
   !! author: David A. Minton
   !!
   !! Perform linear drift of test particles due to barycentric momentum of Sun
   !! New vectorized version included
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift_tp.f90
   !! Adapted from Hal Levison's Swift routine helio_lindrift_tp.f
   use swiftest
   integer(I4B)          :: i

! executable code
   do i = 1, ntp
      if (swiftest_tpA%status(i) == ACTIVE) then 
         swiftest_tpA%xh(:,i) = swiftest_tpA%xh(:,i) + pt(:) * dt
      end if
   end do

   return

   end procedure helio_lindrift_tp
end submodule s_helio_lindrift_tp