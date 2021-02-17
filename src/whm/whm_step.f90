submodule(whm_classes) s_whm_step
contains
   module procedure whm_step_system
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step.f90
      use swiftest
      implicit none

      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt)
         if (allocated(tp%xbeg)) deallocate(tp%xbeg)
         allocate(tp%xbeg, source=pl%xh)
         call pl%step(cb, config, t, dt)
         if (ntp > 0) then
            if (allocated(tp%xend)) deallocate(tp%xend)
            allocate(tp%xend, source=pl%xh)
            call tp%step(cb, pl, config, t, dt)
         end if
      end associate
   end procedure whm_step_system 

   module procedure whm_step_pl
      !! author: David A. Minton
      !!
      !! Step planets ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_pl.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_pl.f90
      !logical, save :: lfirst = .true.
      real(DP) :: dth
      
      dth = 0.5_DP * dt
      if (self%lfirst) then
         call self%h2j(cb)
         call self%getacch(cb, config, t)
         self%lfirst = .false.
      end if

      call self%kickvh(dth)
      call self%vh2vj(cb) 
      !If GR enabled, calculate the p4 term before and after each drift
      if (config%lgr) call self%gr_p4(config, dth)
      call self%drift(cb, config, dt)
      if (config%lgr) call self%gr_p4(config, dth)
      call self%j2h(cb)
      call self%getacch(cb, config, t + dt)
      call self%kickvh(dth)

   end procedure whm_step_pl

   module procedure whm_step_tp
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_tp.f90
      real(DP) :: dth
      dth = 0.5_DP * dt
      if (self%lfirst) then
         call self%getacch(cb, pl, config, t, self%xbeg)
         self%lfirst = .false.
      end if
      call self%kickvh(dth)
      !If GR enabled, calculate the p4 term before and after each drift
      if (config%lgr) call self%gr_p4(config, dth)
      call self%drift(cb, config, dt)
      if (config%lgr) call self%gr_p4(config, dth)
      call self%getacch(cb, pl, config, t + dt, self%xend)
      call self%kickvh(dth)
      return
   end procedure whm_step_tp   

end submodule s_whm_step
