submodule (helio) s_helio_step
contains
   module procedure helio_step
   !! author: David A. Minton
   !!
   !! Step massive bodies and active test particles ahead in heliocentric coordinates
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
   !! Adapted from Hal Levison's Swift routine helio_step.f
   use swiftest
   logical                                :: lfirsttp
   logical, save                          :: lmalloc = .true.
   real(DP), dimension(NDIM)              :: ptb, pte
   real(DP), dimension(:, :), allocatable, save :: xbeg, xend

   if (lmalloc) then
      allocate(xbeg(NDIM, config%nplmax), xend(NDIM, config%nplmax))
      lmalloc = .false.
   end if
   lfirsttp = lfirst
   call helio_plA%step(helio_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)
   if (helio_tpA%nbody > 0) call helio_tpA%step(helio_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)

   return

   end procedure helio_step
end submodule s_helio_step