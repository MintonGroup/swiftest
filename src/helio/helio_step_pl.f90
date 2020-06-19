submodule (helio) s_helio_step_pl
contains
module procedure helio_step_pl
   !! author: David A. Minton
   !!
   !! Step massive bodies ahead Democratic Heliocentric method
   !!
   !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
   !! Adapted from Hal Levison's Swift routine helio_step_pl.f
   use swiftest
   implicit none

   logical      :: lflag
   integer(I4B)     :: i,npl
   real(DP)         :: dth, msys

   npl = self%nbody
   dth = 0.5_DP * dt
   lflag = lfirst
   if (lfirst) then
      call self%vh2vb()
      lfirst = .false.
   end if
   call self%lindrift(dth, ptb)
   call self%getacch(config, t, lflag) 
   lflag = .true.
   call self%kick(dth)
   xbeg(:, 2:npl) = self%xh(:, 2:npl)
   call self%drift(dt)
   xend(:, 2:npl) = self%xh(:, 2:npl)
   call self%getacch(config, t + dt, lflag) 
   call self%kick(dth)
   call self%lindrift(dth, pte)
   call self%vb2vh()

   return

   end procedure helio_step_pl
end submodule s_helio_step_pl
