submodule(helio_classes) s_helio_kick
   use swiftest
contains
module subroutine helio_kick_getacch_pl(self, system, param, t, lbeg)
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of massive bodies
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_kick_getacch.f90
   !! Adapted from Hal Levison's Swift routine helio_kick_getacch.f
   implicit none
   ! Arguments
   class(helio_pl),              intent(inout) :: self   !! Helio massive body particle data structure
   class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
   class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
   real(DP),                     intent(in)    :: t      !! Current simulation time
   logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step

   associate(cb => system%cb, pl => self, npl => self%nbody)
      call pl%accel_int()
      if (param%loblatecb) then 
         cb%aoblbeg = cb%aobl
         call pl%accel_obl(system)
         cb%aoblend = cb%aobl
         if (param%ltides) then
            cb%atidebeg = cb%atide
            call pl%accel_tides(system)
            cb%atideend = cb%atide
         end if
      end if
      if (param%lextra_force) call pl%accel_user(system, param, t)
      !if (param%lgr) call pl%gr_accel(param)
   end associate

   return
   end subroutine helio_kick_getacch_pl

   module subroutine helio_kick_getacch_tp(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_getacch_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_kick_getacch_tp.f
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self   !! Helio test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
   
      associate(tp => self, cb => system%cb, pl => system%pl, npl => system%pl%nbody)
         if (present(lbeg)) system%lbeg = lbeg
         if (system%lbeg) then
            call tp%accel_int(pl%Gmass(:), pl%xbeg(:,:), npl)
         else
            call tp%accel_int(pl%Gmass(:), pl%xend(:,:), npl)
         end if
         if (param%loblatecb) call tp%accel_obl(system)
         if (param%lextra_force) call tp%accel_user(system, param, t)
         !if (param%lgr) call tp%gr_accel(param)
      end associate
      return
   end subroutine helio_kick_getacch_tp

   module subroutine helio_kick_vb_pl(self, system, param, t, dt, mask, lbeg)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f 
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_vb.f90
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Swiftest generic body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical, dimension(:),        intent(in)    :: mask   !! Mask that determines which bodies to kick
      logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      ! Internals
      integer(I4B) :: i

      associate(pl => self, npl => self%nbody)
         if (npl ==0) return
         pl%ah(:,:) = 0.0_DP
         call pl%accel(system, param, t)
         if (lbeg) then
            call pl%set_beg_end(xbeg = pl%xh)
         else
            call pl%set_beg_end(xend = pl%xh)
         end if
         do concurrent(i = 1:npl, mask(i)) 
            pl%vb(:, i) = pl%vb(:, i) + pl%ah(:, i) * dt
         end do
      end associate
   
      return
   
   end subroutine helio_kick_vb_pl

   module subroutine helio_kick_vb_tp(self, system, param, t, dt, mask, lbeg)
      !! author: David A. Minton
      !!
      !! Kick barycentric velocities of bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_kick_vb_tp.f90
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self !! Swiftest generic body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t     !! Current time
      real(DP),                     intent(in)    :: dt    !! Stepsize
      logical, dimension(:),        intent(in)    :: mask  !! Mask that determines which bodies to kick
      logical,                      intent(in)    :: lbeg  !! Logical flag indicating whether this is the beginning of the half step or not. 
      ! Internals
      integer(I4B) :: i

      associate(tp => self, ntp => self%nbody)
         if (ntp ==0) return
         tp%ah(:,:) = 0.0_DP
         call tp%accel(system, param, t, lbeg)
         do concurrent(i = 1:ntp, mask(i)) 
            tp%vb(:, i) = tp%vb(:, i) + tp%ah(:, i) * dt
         end do
      end associate
   
      return
   
   end subroutine helio_kick_vb_tp
end submodule s_helio_kick