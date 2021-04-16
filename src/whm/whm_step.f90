submodule(whm_classes) s_whm_step
   use swiftest
contains


   module subroutine whm_step_system(self, config)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step.f90
      implicit none
      ! Arguments
      class(whm_nbody_system),      intent(inout) :: self    !! WHM nbody system object
      class(swiftest_configuration), intent(in)    :: config  !! Input collection of on parameters 

      select type(cb => self%cb)
      class is (whm_cb)
      select type(pl => self%pl)
      class is (whm_pl)
      select type(tp => self%tp)
      class is (whm_tp)
      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt)
         call pl%set_rhill(cb)
         call tp%set_beg_end(xbeg = pl%xh)
         call pl%step(cb, config, t, dt)
         if (ntp > 0) then
            call tp%set_beg_end(xend = pl%xh)
            call tp%step(cb, pl, config, t, dt)
         end if
      end associate
      end select
      end select
      end select
      return
   end subroutine whm_step_system 

   module subroutine whm_step_pl(self, cb, config, t, dt)
      !! author: David A. Minton
      !!
      !! Step planets ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_pl.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_pl.f90
      !logical, save :: lfirst = .true.
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structure
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals
      real(DP)                                     :: dth
      
      associate(pl => self, xh => self%xh, vh => self%vh, ah => self%ah, &
               xj => self%xj, vj => self%vj)
         dth = 0.5_DP * dt
         if (pl%lfirst) then
            call pl%h2j(cb)
            call pl%getacch(cb, config, t)
            pl%lfirst = .false.
         end if

         call pl%kickvh(dth)
         call pl%vh2vj(cb) 
         !If GR enabled, calculate the p4 term before and after each drift
         if (config%lgr) call pl%gr_p4(config, dth)
         call pl%drift(cb, config, dt)
         if (config%lgr) call pl%gr_p4(config, dth)
         call pl%j2h(cb)
         call pl%getacch(cb, config, t + dt)
         call pl%kickvh(dth)
      end associate
      return
   end subroutine whm_step_pl

   module subroutine whm_step_tp(self, cb, pl, config, t, dt)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                 intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structure
      class(whm_pl),                 intent(inout) :: pl     !! WHM massive body data structure
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals
      real(DP)                                     :: dth

      associate(tp => self, xht => self%xh, vht => self%vh, aht => self%ah, &
         xbeg => self%xbeg, xend => self%xend)
         dth = 0.5_DP * dt
         if (tp%lfirst) then
            call tp%getacch(cb, pl, config, t, xbeg)
            tp%lfirst = .false.
         end if
         call tp%kickvh(dth)
         !If GR enabled, calculate the p4 term before and after each drift
         if (config%lgr) call tp%gr_p4(config, dth)
         call tp%drift(cb, config, dt)
         if (config%lgr) call tp%gr_p4(config, dth)
         call tp%getacch(cb, pl, config, t + dt, xend)
         call tp%kickvh(dth)
      end associate
      return
   end subroutine whm_step_tp   

end submodule s_whm_step
