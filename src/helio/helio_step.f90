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
         tp%lfirst = pl%lfirst
         call pl%set_rhill(cb)
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

      if (self%nbody == 0) return
      associate(pl => self)
         select type(cb => system%cb)
         class is (helio_cb)
            dth = 0.5_DP * dt
            if (pl%lfirst) then
               call pl%vh2vb(cb)
               pl%lfirst = .false.
            end if
            call pl%lindrift(cb, dth, lbeg=.true.)
            call pl%accel(system, param, t)
            call pl%kick(dth)
            call pl%set_beg_end(xbeg = pl%xh)
            call pl%drift(system, param, dt)
            call pl%set_beg_end(xend = pl%xh)
            call pl%accel(system, param, t + dt)
            call pl%kick(dth)
            call pl%lindrift(cb, dth, lbeg=.false.)
            call pl%vb2vh(cb)
         end select
      end associate
   
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
   
      if (self%nbody == 0) return

      associate(tp => self)
         select type(cb => system%cb)
         class is (helio_cb)
            dth = 0.5_DP * dt
            if (tp%lfirst) then
               call tp%vh2vb(vbcb = -cb%ptbeg)
               tp%lfirst = .false.
            end if
            call tp%lindrift(cb, dth, lbeg=.true.)
            call tp%accel(system, param, t, lbeg=.true.)
            call tp%kick(dth)
            call tp%drift(system, param, dt)
            call tp%accel(system, param, t + dt, lbeg=.false.)
            call tp%kick(dth)
            call tp%lindrift(cb, dth, lbeg=.false.)
            call tp%vb2vh(vbcb = -cb%ptend)
         end select
      end associate
   
      return
   
   end subroutine helio_step_tp

end submodule s_helio_step
