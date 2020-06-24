submodule (symba) s_symba_step_helio
contains
   module procedure symba_step_helio
   !! author: David A. Minton
   !!
   !! Step planets and test particles ahead in democratic heliocentric coordinates
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_step_helio.f90
   use swiftest
   implicit none
   logical                      :: lfirsttp
   real(DP), dimension(NDIM)              :: ptb, pte
   real(DP), dimension(NDIM, nplm)          :: xbeg, xend

! executable code
   lfirsttp = lfirst
   call symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, config%nplmax, helio_plA, config%j2rp2, config%j4rp4, dt, xbeg, xend, ptb, pte)
   if (ntp > 0) call helio_step_tp(lfirsttp, lextra_force, t, nplm, config%nplmax, ntp, config%ntpmax, helio_plA, helio_tpA, config%j2rp2, config%j4rp4,  &
      dt, xbeg, xend, ptb, pte)

   return

   end procedure symba_step_helio
end submodule s_symba_step_helio
