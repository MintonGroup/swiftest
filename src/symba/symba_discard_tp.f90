submodule (symba) s_symba_discard_tp
contains
   module procedure symba_discard_tp
   !! author: David A. Minton
   !!
   !! Call discard routine to determine spilled test particles, then remove them from active list
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_discard_tp.f90
use swiftest
implicit none
   logical(lgt)          :: lclosel = .false.
   integer(I4B)          :: i

! executable code
   call discard_tpt, dt, npl, ntp, symba_plA, symba_tpA, rmin, rmax, config%rmaxu, qmin, &
      qmin_alo, qmin_ahi, qmin_coord, lclosel,  &
      lrhill_present)
   return

   end procedure symba_discard_tp
end submodule s_symba_discard_tp
