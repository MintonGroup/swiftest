! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(whm) s_whm_kick
   use swiftest
contains

   module subroutine whm_kick_getacch_pl(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch.f90
      implicit none
      ! Arguments
      class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest central body particle data structure
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t       !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)                                :: i
      real(DP), dimension(NDIM)                   :: ah0

      if (self%nbody == 0) return

      associate(cb => nbody_system%cb, pl => self, npl => self%nbody)
         call pl%set_ir3()

         ah0(:) = whm_kick_getacch_ah0(pl%Gmass(2:npl), pl%rh(:,2:npl), npl-1)
         do i = 1, npl
            pl%ah(:, i) = pl%ah(:, i) + ah0(:)
         end do

         call whm_kick_getacch_ah1(cb, pl) 
         call whm_kick_getacch_ah2(cb, pl) 
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
            !    cb%atidebeg = cb%aobl
            !    call pl%accel_tides(nbody_system)
            !    cb%atideend = cb%atide
            ! end if
         end if

         if (param%lgr) call pl%accel_gr(param) 
         if (param%lextra_force) call pl%accel_user(nbody_system, param, t, lbeg)
      end associate

      return
   end subroutine whm_kick_getacch_pl


   module subroutine whm_kick_getacch_tp(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest central body particle data structure
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)                                :: i, npl, ntp
      real(DP), dimension(NDIM)                   :: ah0
   
      associate(tp => self, pl => nbody_system%pl, cb => nbody_system%cb)
         npl = nbody_system%pl%nbody
         ntp = self%nbody
         if (ntp == 0) return
         nbody_system%lbeg = lbeg

         if(npl > 0) then
            if (lbeg) then
               ah0(:) = whm_kick_getacch_ah0(pl%Gmass(1:npl), pl%rbeg(:, 1:npl), npl)
#ifdef DOCONLOC
               do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp,ah0)
#else
               do concurrent(i = 1:ntp, tp%lmask(i))
#endif
                  tp%ah(:, i) = tp%ah(:, i) + ah0(:)
               end do
               call tp%accel_int(param, pl%Gmass(1:npl), pl%rbeg(:, 1:npl), npl)
            else
               ah0(:) = whm_kick_getacch_ah0(pl%Gmass(1:npl), pl%rend(:, 1:npl), npl)
#ifdef DOCONLOC
               do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp,ah0)
#else
               do concurrent(i = 1:ntp, tp%lmask(i))
#endif
                  tp%ah(:, i) = tp%ah(:, i) + ah0(:)
               end do
               call tp%accel_int(param, pl%Gmass(1:npl), pl%rend(:, 1:npl), npl)
            end if
         end if

         if (param%lnon_spherical_cb) call tp%accel_non_spherical_cb(nbody_system)
         if (param%lextra_force) call tp%accel_user(nbody_system, param, t, lbeg)
         if (param%lgr) call tp%accel_gr(param) 
      end associate

      return
   end subroutine whm_kick_getacch_tp


   function whm_kick_getacch_ah0(mu, rhp, n) result(ah0)
      !! author: David A. Minton
      !!
      !! Compute zeroth term heliocentric accelerations of planets 
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)         :: mu
      real(DP), dimension(:,:), intent(in)         :: rhp
      integer(I4B),             intent(in)         :: n
      ! Result
      real(DP), dimension(NDIM)                    :: ah0
      ! Internals
      real(DP)                                     :: fac, r2, ir3h
      integer(I4B)                                 :: i

      ah0(:) = 0.0_DP
      do i = 1, n
         r2 = dot_product(rhp(:, i), rhp(:, i))
         ir3h = 1.0_DP / (r2 * sqrt(r2)) 
         fac = mu(i) * ir3h 
         ah0(:) = ah0(:) - fac * rhp(:, i)
      end do

      return
   end function whm_kick_getacch_ah0


   pure subroutine whm_kick_getacch_ah1(cb, pl)
      !! author: David A. Minton
      !!
      !! Compute first term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah1.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah1.f90
      implicit none
      ! Arguments
      class(swiftest_cb), intent(in)    :: cb !! WHM central body object
      class(whm_pl),      intent(inout) :: pl !! WHM massive body object
      ! Internals
      integer(I4B)                 :: i, npl
      real(DP), dimension(NDIM)    :: ah1h, ah1j

      npl = pl%nbody
#ifdef DOCONLOC
      do concurrent (i = 2:npl, pl%lmask(i)) shared(pl,cb) local(ah1j,ah1h)
#else
      do concurrent (i = 2:npl, pl%lmask(i))
#endif
         ah1j(:) = pl%xj(:, i) * pl%ir3j(i)
         ah1h(:) = pl%rh(:, i) * pl%ir3h(i)
         pl%ah(:, i) = pl%ah(:, i) + cb%Gmass * (ah1j(:) - ah1h(:))
      end do
   
      return
   end subroutine whm_kick_getacch_ah1


   pure subroutine whm_kick_getacch_ah2(cb, pl)
      !! author: David A. Minton
      !!
      !! Compute second term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah2.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kick_getacch_ah2.f90
      implicit none
      ! Arguments
      class(swiftest_cb), intent(in)    :: cb !! Swiftest central body object
      class(whm_pl),      intent(inout) :: pl !! WHM massive body object
      ! Internals
      integer(I4B)                 :: i, npl
      real(DP)                     :: etaj, fac
      real(DP), dimension(NDIM)    :: ah2, ah2o
   
      npl = pl%nbody
      ah2(:) = 0.0_DP
      ah2o(:) = 0.0_DP
      etaj = cb%Gmass
#ifdef DOCONLOC
      do concurrent(i = 2:npl, pl%lmask(i)) shared(pl,cb,ah2,ah2o) local(etaj,fac)
#else
      do concurrent(i = 2:npl, pl%lmask(i))
#endif
         etaj = etaj + pl%Gmass(i - 1)
         fac = pl%Gmass(i) * cb%Gmass * pl%ir3j(i) / etaj
         ah2(:) = ah2o(:) + fac * pl%xj(:, i)
         pl%ah(:,i) = pl%ah(:, i) + ah2(:)
         ah2o(:) = ah2(:)
      end do
   
      return
   end subroutine whm_kick_getacch_ah2


   module subroutine whm_kick_vh_pl(self, nbody_system, param, t, dt, lbeg)
      !! author: David A. Minton
      !!
      !! Kick heliocentric velocities of massive bodies
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f 
      !! Adapted from David E. Kaufmann's Swifter routine whm_kickvh.f90 
      implicit none
      ! Arguments
      class(whm_pl),                intent(inout) :: self  !! WHM massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      ! Internals
      integer(I4B) :: i, npl

      associate(pl => self, cb => nbody_system%cb)
         npl = self%nbody
         if (npl == 0) return
         if (lbeg) then
            if (pl%lfirst) then
               call pl%h2j(cb)
               pl%ah(:, 1:npl) = 0.0_DP
               call pl%accel(nbody_system, param, t, lbeg)
               pl%lfirst = .false.
            end if
            call pl%set_beg_end(rbeg = pl%rh)
         else
            pl%ah(:, 1:npl) = 0.0_DP
            call pl%accel(nbody_system, param, t, lbeg)
            call pl%set_beg_end(rend = pl%rh)
         end if
#ifdef DOCONLOC
         do concurrent(i = 1:npl, pl%lmask(i)) shared(pl,dt)
#else
         do concurrent(i = 1:npl, pl%lmask(i))
#endif
            pl%vh(:, i) = pl%vh(:, i) + pl%ah(:, i) * dt
         end do
      end associate

      return
   end subroutine whm_kick_vh_pl


   module subroutine whm_kick_vh_tp(self, nbody_system, param, t, dt, lbeg)
      !! author: David A. Minton
      !!
      !! Kick heliocentric velocities of test particles
      !!
      !! Adapted from Martin Duncan and Hal Levison's Swift routine kickvh_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_kickvh_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                intent(inout) :: self   !! WHM massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      ! Internals
      integer(I4B) :: i, ntp

      if (self%nbody == 0) return

      associate(tp => self)
         ntp = self%nbody
         if (tp%lfirst) then
#ifdef DOCONLOC
            do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp)
#else
            do concurrent(i = 1:ntp, tp%lmask(i))
#endif
               tp%ah(:, i) = 0.0_DP
            end do
            call tp%accel(nbody_system, param, t, lbeg=.true.)
            tp%lfirst = .false.
         end if
         if (.not.lbeg) then
#ifdef DOCONLOC
            do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp)
#else
            do concurrent(i = 1:ntp, tp%lmask(i))
#endif
               tp%ah(:, i) = 0.0_DP
            end do
            call tp%accel(nbody_system, param, t, lbeg)
         end if
#ifdef DOCONLOC
         do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp)
#else
         do concurrent(i = 1:ntp, tp%lmask(i))
#endif
            tp%vh(:, i) = tp%vh(:, i) + tp%ah(:, i) * dt
         end do
      end associate

      return
   end subroutine whm_kick_vh_tp

end submodule s_whm_kick
