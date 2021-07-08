submodule(whm_classes) s_whm_step
   use swiftest
contains

   module subroutine whm_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step.f90
      implicit none
      ! Arguments
      class(whm_nbody_system),    intent(inout) :: self  !! WHM nbody system object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters of on parameters 
      real(DP),                   intent(in)    :: t     !! Current simulation time
      real(DP),                   intent(in)    :: dt    !! Current stepsize

      associate(system => self, cb => self%cb,  pl => self%pl, tp => self%tp, ntp => self%tp%nbody)
         call pl%set_rhill(cb)
         call pl%set_beg_end(xbeg = pl%xh)
         call pl%step(system, param, t, dt)
         if (ntp > 0) then
            call pl%set_beg_end(xend = pl%xh)
            call tp%step(system, param, t, dt)
         end if
      end associate
      return
   end subroutine whm_step_system 

   module subroutine whm_step_pl(self, system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step planets ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_pl.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_pl.f90
      !logical, save :: lfirst = .true.
      implicit none
      ! Arguments
      class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t     !! Current simulation time
      real(DP),                     intent(in)    :: dt    !! Current stepsize
      ! Internals
      real(DP)                                    :: dth
      
      associate(pl => self, cb => system%cb)
         dth = 0.5_DP * dt
         if (pl%lfirst) then
            call pl%h2j(cb)
            call pl%get_accel(system, param, t)
            pl%lfirst = .false.
         end if
         call pl%kick(dth)
         call pl%vh2vj(cb) 
         !If GR enabled, calculate the p4 term before and after each drift
         if (param%lgr) call pl%gr_p4(param, dth)
         call pl%drift(system, param, dt)
         if (param%lgr) call pl%gr_p4(param, dth)
         call pl%j2h(cb)
         call pl%get_accel(system, param, t + dt)
         call pl%kick(dth)
      end associate
      return
   end subroutine whm_step_pl

   module subroutine whm_step_tp(self, system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t     !! Current simulation time
      real(DP),                     intent(in)    :: dt    !! Current stepsize
      ! Internals
      real(DP)                                    :: dth

      select type(system)
      class is (whm_nbody_system)
         associate(tp => self, cb => system%cb, pl => system%pl)
            dth = 0.5_DP * dt
            if (tp%lfirst) then
               call tp%get_accel(system, param, t, pl%xbeg)
               tp%lfirst = .false.
            end if
            call tp%kick(dth)
            if (param%lgr) call tp%gr_p4(param, dth)
            call tp%drift(system, param, dt)
            if (param%lgr) call tp%gr_p4(param, dth)
            call tp%get_accel(system, param, t + dt, pl%xend)
            call tp%kick(dth)
         end associate
      end select
      return
   end subroutine whm_step_tp   

end submodule s_whm_step
