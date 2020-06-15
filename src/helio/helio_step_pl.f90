submodule (helio) s_helio_step_pl
contains
module procedure helio_step_pl
   !! author: David A. Minton
   !!
   !! Step planets ahead Democratic Heliocentric method
   !!
   !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
   !! Adapted from Hal Levison's Swift routine helio_step_pl.f
   use swiftest
   logical(lgt)     :: lflag
   integer(I4B)     :: i
   real(DP)         :: dth, msys

   dth = 0.5_DP*dt
   lflag = lfirst
   if (lfirst) then
      call coord_vh2vb(npl, helio_pla%swiftest, msys)
      lfirst = .false.
   end if
   call helio_lindrift(npl, helio_pla%swiftest, dth, ptb)
   call helio_getacch(lflag, lextra_force, t, npl, nplmax, helio_pla, j2rp2, j4rp4)
   lflag = .true.
   call helio_kickvb(npl, helio_pla, dth)
   do i = 2, npl
      xbeg(:, i) = helio_pla%swiftest%xh(:,i)
   end do
   call helio_drift(npl, helio_pla%swiftest, dt)
   do i = 2, npl
      xend(:, i) = helio_pla%swiftest%xh(:,i)
   end do
   call helio_getacch(lflag, lextra_force, t+dt, npl, nplmax, helio_pla, j2rp2, j4rp4)
   call helio_kickvb(npl, helio_pla, dth)
   call helio_lindrift(npl, helio_pla%swiftest, dth, pte)
   call coord_vb2vh(npl, helio_pla%swiftest)

   return

   end procedure helio_step_pl
end submodule s_helio_step_pl