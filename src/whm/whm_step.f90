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
  
      associate(ntp => self%tp%nbody, npl => self%pl%nbody, &
                t => config%t, dt => config%dt, is_tp => self%tp%nbody > 0, is_pl => self%pl%nbody > 0)
         dth = 0.5_DP * dt 
         !> Note: The nested select statements serve to make this method a "bridge" between the polymorphic swiftest_nbody_system class
         !> in which the cb, pl, and tp components are allocatable abstract classes and the actual integrator-specific methods that are
         !> called internally. Before this point, the actual types of cb, pl, and tp are ambiguous. The select type constructs remove the 
         !> ambiguity. - D. Minton
         select type(config => config)
            class is (whm_configuration)
            select type(cb => self%cb)
               class is (whm_cb)
               select type(pl => self%pl)
               class is (whm_pl)
               associate(xh => pl%xh, vh => pl%vh, xj => pl%xj, vj => pl%vj, ah => pl%ah,  eta => pl%eta) ! These associations aid in debugging with gdb
                  select type(tp => self%tp)
                  class is (whm_tp)
                  associate(xht => tp%xh, vht => tp%vh, aht => tp%ah) ! These associations aid in debugging with gdb
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
                  end select
               end associate
               end select
            end select
         end select
      end associate

   end procedure whm_step_system 

end submodule whm_step
