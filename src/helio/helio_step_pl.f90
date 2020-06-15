submodule (helio) s_helio_step_pl
contains
module procedure helio_step_pl
   !! author: David A. Minton
   !!
   !! Step plAnets ahead Democratic Heliocentric method
   !!
   !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
   !! Adapted from Hal Levison's Swift routine helio_step_pl.f
   use swiftest
   logical(lgt)     :: lflag
   integer(I4B)     :: i,npl
   real(DP)         :: dth, msys

   npl = helio_plA%nbody
   dth = 0.5_DP * dt
   lflag = lfirst
   if (lfirst) then
      call coord_vh2vb(helio_plA, msys)
      lfirst = .false.
   end if
   call helio_lindrift_pl(helio_plA, dth, ptb)
   call helio_getacch_pl(helio_plA, config, t, lflag) 
   lflag = .true.
   call helio_kickvb_pl(helio_plA, dth)
   do i = 2, npl
      xbeg(:, i) = helio_plA%xh(:,i)
   end do
   call helio_drift_pl(helio_plA, dt)
   do i = 2, npl
      xend(:, i) = helio_plA%xh(:,i)
   end do
   call helio_getacch_pl(helio_plA, config, t + dt, lflag) 
   call helio_kickvb_pl(helio_plA, dth)
   call helio_lindrift_pl(helio_plA, dth, pte)
   call coord_vb2vh(helio_plA)

   return

   end procedure helio_step_pl
end submodule s_helio_step_pl