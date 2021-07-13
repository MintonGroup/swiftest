submodule (helio_classes) s_helio_getacch
   use swiftest
contains
   module subroutine helio_getacch_pl(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_getacch.f90
      !! Adapted from Hal Levison's Swift routine helio_getacch.f
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step

      associate(cb => system%cb, pl => self, npl => self%nbody)
         call helio_getacch_int_pl(pl, t)
         if (param%loblatecb) then 
            cb%aoblbeg = cb%aobl
            call pl%accel_obl(system)
            cb%aoblend = cb%aobl
         end if
         if (param%lextra_force) call pl%accel_user(system, param, t)
         if (param%ltides) call pl%accel_tides(system)
         !if (param%lgr) call pl%gr_accel(param)
      end associate

      return
      end subroutine helio_getacch_pl

      module subroutine helio_getacch_tp(self, system, param, t, lbeg)
         !! author: David A. Minton
         !!
         !! Compute heliocentric accelerations of test particles
         !!
         !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_tp.f90
         !! Adapted from Hal Levison's Swift routine helio_getacch_tp.f
         implicit none
         ! Arguments
         class(helio_tp),              intent(inout) :: self   !! WHM test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
         real(DP),                     intent(in)    :: t      !! Current time
         logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
         ! Internals
         logical, save                               :: lmalloc = .true.
         integer(I4B)                                :: i
         real(DP)                                    :: r2, mu
         real(DP), dimension(:), allocatable, save   :: irh, irht
      
         associate(tp => self, ntp => self%nbody, cb => system%cb,  npl => system%pl%nbody)
            if (present(lbeg)) system%lbeg = lbeg
            call helio_getacch_int_tp(tp, system, param, t)
            if (param%loblatecb) call tp%accel_obl(system)
            if (param%lextra_force) call tp%accel_user(system, param, t)
            !if (param%lgr) call tp%gr_accel(param)
         end associate
         return
      end subroutine helio_getacch_tp

      subroutine helio_getacch_int_pl(pl, t)
         !! author: David A. Minton
         !!
         !! Compute direct cross term heliocentric accelerations of massive bodiese
         !!
         !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int.f90
         !! Adapted from Hal Levison's Swift routine getacch_ah3.f
         implicit none
         ! Arguments
         class(helio_pl), intent(inout) :: pl     !! Helio massive body particle data structure
         real(DP),        intent(in)    :: t        !! Current time
         ! Internals
         integer(I4B)              :: i, j
         real(DP)                  :: rji2, irij3, faci, facj
         real(DP), dimension(NDIM) :: dx
      
         associate(npl => pl%nbody)
            pl%ah(:,:) = 0.0_DP
            do i = 1, npl - 1
               do j = i + 1, npl
                  dx(:) = pl%xh(:,j) - pl%xh(:,i)
                  rji2 = dot_product(dx(:), dx(:))
                  irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                  faci = pl%Gmass(i) * irij3
                  facj = pl%Gmass(j) * irij3
                  pl%ah(:,i) = pl%ah(:,i) + facj * dx(:)
                  pl%ah(:,j) = pl%ah(:,j) - faci * dx(:)
               end do
            end do
         end associate
      
         return
      end subroutine helio_getacch_int_pl

      subroutine helio_getacch_int_tp(tp, system, param, t)
         !! author: David A. Minton
         !!
         !! Compute direct cross term heliocentric accelerations of test particles
         !!
         !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int_tp.f90
         !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
         implicit none
         ! Arguments
         class(helio_tp),              intent(inout) :: tp     !! Helio test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
         real(DP),                     intent(in)    :: t      !! Current times
         ! Internals
         integer(I4B)                                :: i, j
         real(DP)                                    :: r2, fac
         real(DP), dimension(NDIM)                   :: dx
         real(DP), dimension(:, :), allocatable      :: xhp

         associate(ntp => tp%nbody, pl => system%pl, npl => system%pl%nbody, lbeg => system%lbeg)
            if (lbeg) then
               allocate(xhp, source=pl%xbeg)
            else
               allocate(xhp, source=pl%xend)
            end if
      
            tp%ah(:,:) = 0.0_DP
            do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               do j = 1, npl
                  dx(:) = tp%xh(:,i) - xhp(:,j)
                  r2 = dot_product(dx(:), dx(:))
                  fac = pl%Gmass(j) / (r2 * sqrt(r2))
                  tp%ah(:,i) = tp%ah(:,i) - fac * dx(:)
               end do
            end if
            end do
         end associate 
         return
      end subroutine helio_getacch_int_tp
   
end submodule s_helio_getacch
