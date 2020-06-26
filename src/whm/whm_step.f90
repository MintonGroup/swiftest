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
      real(DP), dimension(:, :), allocatable, save :: xbeg, xend
  
      associate(ntp => self%tp%nbody, npl => self%pl%nbody)
         !> Note: The nesterd select statements serve to make this method a "bridge" between the polymorphic swiftest_nbody_system class
         !> in which the cb, pl, and tp components are allocatable abstract classes and the actual integrator-specific methods that are
         !> called internally. Before this point, the actual types of cb, pl, and tp are ambiguous. The select type constructs remove the 
         !> ambiguity. - D. Minton
         select type(cb => self%cb) 
         class is (whm_central_body)
            select type(pl => self%pl)
            class is (whm_pl)
               if (ntp > 0) xbeg(:, :) = pl%xh(:, 1:npl)
               call pl%step(cb, config, t, dt)
               if (ntp > 0) xend(:, :) = pl%xh(:, 1:npl)
               select type(tp => self%tp)
               class is (whm_tp)
                  call tp%step(cb, pl, config, t, dt, xbeg, xend)
               end select
            end select
         end select
      end associate
      self%lintegrate = .true.
      return
   end procedure whm_step_system



end submodule whm_step
