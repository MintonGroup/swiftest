submodule(whm_classes) whm_step
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
      !real(DP), dimension(:, :), allocatable, save :: xbeg, xend
      logical, save :: lfirst = .true.
      real(DP) :: dth
  
      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt, &
         xh => pl%xh, vh => pl%vh, xj => pl%xj, vj => pl%vj, ah => pl%ah,  eta => pl%eta, & ! These two lines of associations aid in debugging with gdb
         xht => tp%xh, vht => tp%vh, aht => tp%ah) 
         dth = 0.5_DP * dt 
         if (lfirst) then
            call pl%h2j(cb)
            call pl%getacch(cb, config, t)
            call tp%getacch(cb, pl, config, t)
            lfirst = .false.
         end if

         ! ****** Kick  ******
         call pl%kickvh(dth)
         call tp%kickvh(dth)
         call pl%vh2vj(cb) 
         ! *******************

         !If GR enabled, calculate the p4 term before and after each drift
         if (config%lgr) then 
            call pl%gr_p4(config, dth)
            call tp%gr_p4(config, dth)
         end if

         ! ****** Drift ******
         call pl%drift(cb, config, dt)
         call tp%drift(cb, config, dt)
         ! *******************

         if (config%lgr) then
            call pl%gr_p4(config, dth)
            call tp%gr_p4(config, dth)
         end if
         call pl%j2h(cb)

         call pl%getacch(cb, config, t + dt)
         call tp%getacch(cb, pl, config, t + dt)

         ! ****** Kick  ******
         call pl%kickvh(dth)
         call tp%kickvh(dth)
         ! *******************
      end associate

   end procedure whm_step_system 

end submodule whm_step
