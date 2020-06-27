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

   module procedure whm_step_pl
      !! author: David A. Minton
      !!
      !! Step planets ahead using kick-drift-kick algorithm
      !!
      !! Adapted from Hal Levison's Swift routine step_kdk_pl.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_pl.f90
      use swiftest
      implicit none
      real(dp) :: dth
   
      dth = 0.5_DP * dt
      if (lfirst) then
         call self%h2j(npl, whm_pl1p)
         call whm_getacch(lextra_force, t, npl, nplmax, whm_pl1p, j2rp2, j4rp4, c2)
         lfirst = .false.
      end if
      call whm_kickvh(npl, whm_pl1p, dth)
      call coord_vh2vj(npl, whm_pl1p)
      call gr_whm_p4(npl, whm_pl1p, dth, c2)
      call whm_drift(npl, whm_pl1p, dt, c2)
      call gr_whm_p4(npl, whm_pl1p, dth, c2)
      call coord_j2h(npl, whm_pl1p)
      call whm_getacch(lextra_force, t+dt, npl, nplmax, whm_pl1p, j2rp2, j4rp4, c2)
      call whm_kickvh(npl, whm_pl1p, dth)
   
      return
   
   end procedure whm_step_pl

   module procedure whm_step_tp
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using kick-drift-kick algorithm
      !!
      !! Adapted from Hal Levison's Swift routine step_kdk_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_tp.f90
      use swiftest
      implicit none
      real(DP) :: dth
      logical :: lfirsttp = .true.
   
      dth = 0.5_dp*dt
      if (lfirsttp) then
         call whm_getacch_tp(lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, xbeg, j2rp2, j4rp4)
         firsttp = .false.
      end if
      call whm_kickvh_tp(ntp, whm_tp1p, dth)
      call whm_drift_tp(ntp, whm_tp1p, whm_pl1p%swifter%mass, dt, c2)
      call whm_getacch_tp(lextra_force, t+dt, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, xend, j2rp2, j4rp4)
      call whm_kickvh_tp(ntp, whm_tp1p, dth)
   
      return
   
   end procedure whm_step_tp



end submodule whm_step
