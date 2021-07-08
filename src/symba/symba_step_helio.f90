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
   real(DP), dimension(NDIM)              :: ptbeg, ptend
   real(DP), dimension(npl, NDIMm)          :: xbeg, xend

! executable code
   lfirsttp = lfirst
   call symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, param%nplmax, helio_plA, param%j2rp2, param%j4rp4, dt, xbeg, xend, ptbeg, ptend)
   if (ntp > 0) call helio_step_tp(lfirsttp, lextra_force, t, nplm, param%nplmax, ntp, param%ntpmax, helio_plA, helio_tpA, param%j2rp2, param%j4rp4,  &
      dt, xbeg, xend, ptbeg, ptend)

   return

   end procedure symba_step_helio
end submodule s_symba_step_helio
