submodule (helio) s_helio_discard
contains
   module procedure helio_discard
   !! author: David A. Minton
   !!
   !! Call discard routine to determine spilled test particles, then remove them from ACTIVE list
   !!
   !! Adapted from David E. Kaufmann's Swifter helio_discard.f90
   use swiftest
   implicit none

   integer(I4B)    :: i

   call discard(helio_plA, helio_tpA, config, t, dt)
   call helio_tpA%spill(helio_tp_discard)

   return
   end procedure helio_discard
end submodule s_helio_discard