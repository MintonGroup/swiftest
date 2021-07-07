submodule(helio_classes) s_helio_step
   use swiftest
contains
   module subroutine helio_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
      implicit none
      ! Arguments
      class(helio_nbody_system),  intent(inout) :: self   !! Helio nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      real(DP),                   intent(in)    :: t      !! Simulation time
      real(DP),                   intent(in)    :: dt     !! Current stepsize

      associate(system => self, cb => self%cb, pl => self%pl, tp => self%tp)
         call pl%set_rhill(cb)
         tp%lfirst = pl%lfirst
         call pl%step(system, param, t, dt)
         call tp%step(system, param, t, dt)
      end associate
      return
   end subroutine helio_step_system 

   module subroutine helio_step_pl(self, system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies ahead Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
      !! Adapted from Hal Levison's Swift routine helio_step_pl.f
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nboody system
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals 
      integer(I4B)     :: i
      real(DP)         :: dth, msys
 
      select type(system)
      class is (helio_nbody_system)
         associate(pl => self, cb => system%cb, ptb => system%ptb, pte => system%pte)
            dth = 0.5_DP * dt
            if (pl%lfirst) then
               call pl%vh2vb(cb)
               pl%lfirst = .false.
            end if
            call pl%lindrift(system, dth, ptb)
            call pl%get_accel(system, param, t)
            call pl%kick(dth)
            call system%set_beg_end(xbeg = pl%xh)
            call pl%drift(system, param, dt)
            call system%set_beg_end(xend = pl%xh)
            call pl%get_accel(system, param, t + dt)
            call pl%kick(dth)
            call pl%lindrift(system, dth, pte)
            call pl%vb2vh(cb)
         end associate
      end select
   
      return
   
   end subroutine helio_step_pl

   module subroutine helio_step_tp(self, system, param, t, dt)

      !! author: David A. Minton
      !!
      !! Step active test particles ahead using Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_step_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_step_tp.f
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self    !! Helio test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nboody system
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      real(DP),                     intent(in)    :: dt     !! Stepsiz
      ! Internals
      real(DP) :: dth   !! Half step size
   
      select type(system)
      class is (helio_nbody_system)
         associate(tp => self, cb => system%cb, pl => system%pl, xbeg => system%xbeg, xend => system%xend, ptb => system%ptb, pte => system%pte)
            dth = 0.5_DP * dt
            if (tp%lfirst) then
               call tp%vh2vb(vbcb = -ptb)
               tp%lfirst = .false.
            end if
            call tp%lindrift(system, dth, ptb)
            call tp%get_accel(system, param, t, xbeg)
            call tp%kick(dth)
            call tp%drift(system, param, dt)
            call tp%get_accel(system, param, t + dt, xend)
            call tp%kick(dth)
            call tp%lindrift(system, dth, pte)
            call tp%vb2vh(vbcb = -pte)
         end associate
      end select
   
      return
   
   end subroutine helio_step_tp

end submodule s_helio_step
