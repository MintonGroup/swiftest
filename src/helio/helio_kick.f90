! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(helio) s_helio_kick
   use swiftest
contains

   module subroutine helio_kick_getacch_pl(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_getacch.f90
      !! Adapted from Hal Levison's Swift routine helio_kick_getacch.f
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step

      if (self%nbody == 0) return

      associate(cb => nbody_system%cb, pl => self, npl => self%nbody)
         call pl%accel_int(param)
         if (param%lnon_spherical_cb) then 
            call pl%accel_non_spherical_cb(nbody_system)
            if (lbeg) then
               cb%aoblbeg = cb%aobl
            else
               cb%aoblend = cb%aobl
            end if
            ! TODO: Implement tides
            ! if (param%ltides) then
            !    call pl%accel_tides(nbody_system)
            !    if (lbeg) then
            !       cb%atidebeg = cb%atide
            !    else
            !       cb%atideend = cb%atide
            !    end if
            ! end if
         end if
         if (param%lextra_force) call pl%accel_user(nbody_system, param, t, lbeg)
         if (param%lgr) call pl%accel_gr(param)
      end associate

      return
   end subroutine helio_kick_getacch_pl


   module subroutine helio_kick_getacch_tp(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_getacch_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_kick_getacch_tp.f
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self   !! Helio test particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
   
      if (self%nbody == 0) return

      associate(tp => self, cb => nbody_system%cb, pl => nbody_system%pl, npl => nbody_system%pl%nbody)
         nbody_system%lbeg = lbeg
         if (npl > 0) then
            if (nbody_system%lbeg) then
               call tp%accel_int(param, pl%Gmass(1:npl), pl%rbeg(:,1:npl), npl)
            else
               call tp%accel_int(param, pl%Gmass(1:npl), pl%rend(:,1:npl), npl)
            end if
         end if
         if (param%lnon_spherical_cb) call tp%accel_non_spherical_cb(nbody_system)
         if (param%lextra_force) call tp%accel_user(nbody_system, param, t, lbeg)
         if (param%lgr) call tp%accel_gr(param)
      end associate

      return
   end subroutine helio_kick_getacch_tp


   module subroutine helio_kick_vb_pl(self, nbody_system, param, t, dt, lbeg)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f 
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_vb.f90
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Swiftest generic body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      ! Internals
      integer(I4B) :: i, npl

      if (self%nbody == 0) return

      associate(pl => self)
         npl = self%nbody
         pl%ah(:, 1:npl) = 0.0_DP
         call pl%accel(nbody_system, param, t, lbeg)
         if (lbeg) then
            call pl%set_beg_end(rbeg = pl%rh)
         else
            call pl%set_beg_end(rend = pl%rh)
         end if
#ifdef DOCONLOC
         do concurrent(i = 1:npl, pl%lmask(i)) shared(pl,dt)
#else
         do concurrent(i = 1:npl, pl%lmask(i)) 
#endif
            pl%vb(1, i) = pl%vb(1, i) + pl%ah(1, i) * dt
            pl%vb(2, i) = pl%vb(2, i) + pl%ah(2, i) * dt
            pl%vb(3, i) = pl%vb(3, i) + pl%ah(3, i) * dt
         end do
      end associate
   
      return
   end subroutine helio_kick_vb_pl


   module subroutine helio_kick_vb_tp(self, nbody_system, param, t, dt, lbeg)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_vb_tp.f90
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self !! Swiftest generic body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t     !! Current time
      real(DP),                     intent(in)    :: dt    !! Stepsize
      logical,                      intent(in)    :: lbeg  !! Logical flag indicating whether this is the beginning of the half step or not. 
      ! Internals
      integer(I4B) :: i, ntp

      if (self%nbody == 0) return

      associate(tp => self)
         ntp = self%nbody
         tp%ah(:, 1:ntp) = 0.0_DP
         call tp%accel(nbody_system, param, t, lbeg)
#ifdef DOCONLOC
         do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp,dt)
#else
         do concurrent(i = 1:ntp, tp%lmask(i)) 
#endif
            tp%vb(:, i) = tp%vb(:, i) + tp%ah(:, i) * dt
         end do
      end associate
   
      return
   end subroutine helio_kick_vb_tp

end submodule s_helio_kick