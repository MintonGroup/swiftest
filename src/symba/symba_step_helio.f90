submodule (symba) s_symba_step_helio
contains
   module procedure symba_step_helio
   !! author: David A. Minton
   !!
   !! Step planets and test particles ahead in democratic heliocentric coordinates
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_step_helio.f90
   use swiftest
   implicit none
   logical(lgt)                     :: lfirsttp
   real(DP), dimension(ndim)              :: ptb, pte
   real(DP), dimension(ndim, nplm)          :: xbeg, xend

! executable code
   lfirsttp = lfirst
   call symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, nplmax, helio_pla, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
   if (ntp > 0) call helio_step_tp(lfirsttp, lextra_force, t, nplm, nplmax, ntp, ntpmax, helio_pla, helio_tpa, j2rp2, j4rp4,  &
      dt, xbeg, xend, ptb, pte)

   return

   end procedure symba_step_helio
end submodule s_symba_step_helio
