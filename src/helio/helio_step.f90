submodule (helio) s_helio_step
contains
   module procedure helio_step
   !! author: David A. Minton
   !!
   !! Step planets and active test particles ahead in heliocentric coordinatee
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
   !! Adapted from Hal Levison's Swift routine helio_step.f
   use swiftest
   logical(lgt)                     :: lfirsttp
   logical(lgt), save                 :: lmalloc = .true.
   real(DP), dimension(ndim)            :: ptb, pte
   real(DP), dimension(:, :), allocatable, save :: xbeg, xend

   if (lmalloc) then
      allocate(xbeg(ndim, nplmax), xend(ndim, nplmax))
      lmalloc = .false.
   end if
   lfirsttp = lfirst
   call helio_step_pl(helio_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)
   if (ntp > 0) call helio_step_tp(helio_tpA, helio_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)

   return

   end procedure helio_step
end submodule s_helio_step