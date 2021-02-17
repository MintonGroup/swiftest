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
         call tp%set_beg_end(xbeg = pl%xh)
         call pl%step(cb, config, t, dt)
         if (ntp > 0) then
            call tp%set_beg_end(xend = pl%xh)
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
      
      associate(pl => self, xh => self%xh, vh => self%vh, ah => self%ah)
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

   end procedure whm_step_pl

   module procedure whm_step_tp
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_tp.f90
      real(DP) :: dth

      associate(tp => self, xht => self%xh, vht => self%vh, aht => self%ah)
         dth = 0.5_DP * dt
         if (tp%lfirst) then
            call tp%getacch(cb, pl, config, t, tp%xbeg)
            tp%lfirst = .false.
         end if
         call tp%kickvh(dth)
         !If GR enabled, calculate the p4 term before and after each drift
         if (config%lgr) call tp%gr_p4(config, dth)
         call tp%drift(cb, config, dt)
         if (config%lgr) call tp%gr_p4(config, dth)
         call tp%getacch(cb, pl, config, t + dt, tp%xend)
         call tp%kickvh(dth)
      end associate
      return
   end procedure whm_step_tp   

end submodule s_whm_step
