submodule(helio_classes) s_helio_step
   use swiftest
contains

   module subroutine helio_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates.
      !!
      !! Currently there's no difference between this and the WHM system stepper, so this is just
      !! a wrapper function to keep the method calls consistent for inherited types.
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
      implicit none
      ! Arguments
      class(helio_nbody_system),  intent(inout) :: self   !! Helio nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      real(DP),                   intent(in)    :: t      !! Simulation time
      real(DP),                   intent(in)    :: dt     !! Current stepsize

      call whm_step_system(self, param, t, dt)

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
      real(DP) :: dth   !! Half step size 

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
            call pl%kick(system, param, t, dth, lbeg=.true.)
            call pl%drift(system, param, dt)
            call pl%kick(system, param, t + dt, dth, lbeg=.false.)
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
      real(DP),                     intent(in)    :: dt     !! Stepsize
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
            call tp%kick(system, param, t, dth, lbeg=.true.)
            call tp%drift(system, param, dt)
            call tp%kick(system, param, t + dt, dth, lbeg=.false.)
            call tp%lindrift(cb, dth, lbeg=.false.)
            call tp%vb2vh(vbcb = -cb%ptend)
         end select
      end associate
   
      return
   end subroutine helio_step_tp

end submodule s_helio_step
