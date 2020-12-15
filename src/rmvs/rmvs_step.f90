submodule(rmvs_classes) rmvs_step
contains
   module procedure rmvs_step_system
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step.f90
      use swiftest
      implicit none
      logical, save :: lfirst = .true.
      real(DP) :: dth
      real(DP), dimension(:,:), allocatable :: xtmp, vtmp, atmp
  
      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt, &
         xh => pl%xh, vh => pl%vh, xj => pl%xj, vj => pl%vj, ah => pl%ah,  eta => pl%eta, & ! These two lines of associations aid in debugging with gdb
         xht => tp%xh, vht => tp%vh, aht => tp%ah, irij3 => tp%irij3) 
         dth = 0.5_DP * dt 
         if (lfirst) then
            if (ntp > 0) then
               allocate(xtmp, source = xht)
               allocate(vtmp, source = vht)
               allocate(atmp, source = aht)
            end if
            call pl%h2j(cb)
            call pl%getacch(cb, config, t)
            call tp%getacch(cb, pl, config, t)
            if (ntp > 0) atmp = aht
            lfirst = .false.
         end if

         ! ****** Check for close encounters ***** !
         !call rmvs_chk

         ! ****** Kick  ******
         call pl%kickvh(dth)
         call tp%kickvh(dth)
         if (ntp > 0) vtmp = vht
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
         if (ntp > 0) then
            xtmp = xht
            vtmp = vht
         end if 
         ! *******************

         if (config%lgr) then
            call pl%gr_p4(config, dth)
            call tp%gr_p4(config, dth)
         end if
         call pl%j2h(cb)

         call pl%getacch(cb, config, t + dt)
         call tp%getacch(cb, pl, config, t + dt)
         if (ntp > 0) atmp = aht

         ! ****** Kick  ******
         call pl%kickvh(dth)
         call tp%kickvh(dth)
         if (ntp > 0) vtmp = vht
         ! *******************
      end associate

   end procedure rmvs_step_system 

end submodule rmvs_step
