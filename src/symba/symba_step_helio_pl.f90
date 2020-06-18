submodule (symba) s_symba_step_helio_pl
contains
   module procedure symba_step_helio_pl
   !! author: David A. Minton
   !!
   !! Step planets ahead in democratic heliocentric coordinates
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_step_helio_pl.f90
   !! Adapted from Hal Levison's Swift routines symba5_step_helio.f and helio_step_pl.f
use swiftest
implicit none
   logical(lgt)          :: lflag
   integer(I4B)          :: i
   real(DP)            :: dth, msys

! executable code

  
   dth = 0.5_DP*dt
   lflag = lfirst
   if (lfirst) then
      call coord_vh2vb(npl, helio_plA%swiftest, msys) 
      lfirst = .false.
   end if
   
   call helio_lindrift(npl, helio_plA%swiftest, dth, ptb)

   call symba_helio_getacch(lflag, lextra_force, t, npl, nplm, config%nplmax, helio_plA, config%j2rp2, config%j4rp4) 
   lflag = .true.

   call helio_kickvb(npl, helio_plA, dth)

   do i = 2, nplm
      xbeg(:, i) = helio_plA%swiftest%xh(:,i)
   end do
   call helio_drift(npl, helio_plA%swiftest, dt) 

   do i = 2, nplm
      xend(:, i) = helio_plA%swiftest%xh(:,i)
   end do
   call symba_helio_getacch(lflag, lextra_force, t+dt, npl, nplm, config%nplmax, helio_plA, config%j2rp2, config%j4rp4) 

   call helio_kickvb(npl, helio_plA, dth)

   call helio_lindrift(npl, helio_plA%swiftest, dth, pte)

   call coord_vb2vh(npl, helio_plA%swiftest)

   return

   end procedure symba_step_helio_pl
end submodule s_symba_step_helio_pl
