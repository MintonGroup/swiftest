submodule (helio) s_helio_step_tp
contains
module procedure helio_step_tp
   !! author: David A. Minton
   !!
   !! Step active test particles ahead using Democratic Heliocentric method
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_step_tp.f90
   !! Adapted from Hal Levison's Swift routine helio_step_tp.f
   use swiftest
   implicit none

   logical  :: lflag !! Flag to indicate that this is the first call
   real(DP) :: dth   !! Half step size
   real(DP) :: mu    !! Central mass term

! executable code
   dth = 0.5_DP * dt
   lflag = lfirst
   mu = helio_plA%mass(1)
   if (lfirst) then
      call self%vh2vb(vs = -ptb)
      lfirst = .false.
   end if
   call self%lindrift(dth, ptb)
   call self%getacch(config, t, lflag, helio_plA, xbeg)
   lflag = .true.
   call self%kick(dth)
   call self%drift(mu, dt)
   call self%getacch(config, t + dt, lflag, helio_plA, xend)
   call self%kick(dth)
   call self%lindrift(dth, pte)
   call self%vb2vh(vs = -pte)

   return

   end procedure helio_step_tp
end submodule s_helio_step_tp