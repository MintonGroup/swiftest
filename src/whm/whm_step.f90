submodule(whm_classes) s_whm_step
   use swiftest
contains

   module subroutine whm_step_system(self, param, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step.f90
      implicit none
      ! Arguments
      class(whm_nbody_system),    intent(inout) :: self    !! WHM nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters of on parameters 
      integer(I4B),               intent(in)    :: dt     !! Current stepsize

      select type (system => self)
      class is (whm_nbody_system)
      select type(cb => self%cb)
      class is (whm_cb)
      select type(pl => self%pl)
      class is (whm_pl)
      select type(tp => self%tp)
      class is (whm_tp)
      associate(ntp => tp%nbody, npl => pl%nbody, t => param%t, dt => param%dt)
         call pl%set_rhill(cb)
         call tp%set_beg_end(xbeg = pl%xh)
         call pl%step(system, param, dt)
         if (ntp > 0) then
            call tp%set_beg_end(xend = pl%xh)
            call tp%step(system, param, dt)
         end if
      end associate
      end select
      end select
      end select
      end select
      return
   end subroutine whm_step_system 

   module subroutine whm_step_pl(self, system, param, dt)
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
      integer(I4B),                 intent(in)    :: dt     !! Current stepsize
      ! Internals
      real(DP)                                     :: dth
      
      associate(cb => system%cb, t => param%t)
         dth = 0.5_DP * dt
         if (pl%lfirst) then
            call pl%h2j(cb)
            call pl%getacch(system, param, t)
            pl%lfirst = .false.
         end if

         call pl%kickvh(dth)
         call pl%vh2vj(cb) 
         !If GR enabled, calculate the p4 term before and after each drift
         if (param%lgr) call pl%gr_p4(param, dth)
         call pl%drift(cb, param, dt)
         if (param%lgr) call pl%gr_p4(param, dth)
         call pl%j2h(cb)
         call pl%getacch(system, param, t + dt)
         call pl%kickvh(dth)
      end associate
      return
   end subroutine whm_step_pl

   module subroutine whm_step_tp(self, system, param, dt)
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
      integer(I4B),                 intent(in)    :: dt     !! Current stepsize
      ! Internals
      real(DP)                                     :: dth

      associate(cb => system%cb, pl => system%pl, t => param%t, xbeg => self%xbeg, xend => self%xend)
         dth = 0.5_DP * dt
         if (tp%lfirst) then
            call tp%getacch(system, param, t, xbeg)
            tp%lfirst = .false.
         end if
         call tp%kickvh(dth)
         !If GR enabled, calculate the p4 term before and after each drift
         if (param%lgr) call tp%gr_p4(param, dth)
         call tp%drift(cb, param, dt)
         if (param%lgr) call tp%gr_p4(param, dth)
         call tp%getacch(system, param, t + dt, xend)
         call tp%kickvh(dth)
      end associate
      return
   end subroutine whm_step_tp   

end submodule s_whm_step
