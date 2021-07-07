submodule (helio_classes) s_helio_getacch
   use swiftest
contains
   module subroutine helio_getacch_pl(self, system, param, t)
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
      ! Internals

      associate(cb => system%cb, pl => self, npl => self%nbody)
         call helio_getacch_int_pl(pl, t)
         pl%ah(:,:) = pl%ahi(:,:)
         if (param%loblatecb) call pl%obl_acc(cb)
         if (param%lextra_force) call pl%user_getacch(system, param, t)
         !if (param%lgr) call pl%gr_get_accel(param)
      end associate

      return
      end subroutine helio_getacch_pl

      module subroutine helio_getacch_tp(self, system, param, t, xhp)
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
         real(DP), dimension(:,:),     intent(in)    :: xhp    !! Heliocentric positions of planets at the current substep
         ! Internals
         logical, save                                :: lmalloc = .true.
         integer(I4B)                                 :: i
         real(DP)                                     :: r2, mu
         real(DP), dimension(:), allocatable, save    :: irh, irht
      
         associate(tp => self, ntp => self%nbody, cb => system%cb, npl => system%pl%nbody)
            select type(pl => system%pl)
            class is (helio_pl)
               call helio_getacch_int_tp(tp, pl, t, xhp)
               tp%ah(:,:) = tp%ahi(:,:)
               if (param%loblatecb) call tp%obl_acc(cb) 
               if (param%lextra_force) call tp%user_getacch(system, param, t)
               !if (param%lgr) call tp%gr_get_accel(param)
            end select
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
            pl%ahi(:,:) = 0.0_DP
            do i = 1, npl - 1
               do j = i + 1, npl
                  dx(:) = pl%xh(:,j) - pl%xh(:,i)
                  rji2 = dot_product(dx(:), dx(:))
                  irij3 = 1.0_DP / (rji2 * sqrt(rji2))
                  faci = pl%Gmass(i) * irij3
                  facj = pl%Gmass(j) * irij3
                  pl%ahi(:,i) = pl%ahi(:,i) + facj * dx(:)
                  pl%ahi(:,j) = pl%ahi(:,j) - faci * dx(:)
               end do
            end do
         end associate
      
         return
      end subroutine helio_getacch_int_pl

      subroutine helio_getacch_int_tp(tp, pl, t, xhp)
         !! author: David A. Minton
         !!
         !! Compute direct cross term heliocentric accelerations of test particles
         !!
         !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int_tp.f90
         !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
         implicit none
         ! Arguments
         class(helio_tp),          intent(inout) :: tp     !! Helio test particle data structure
         class(swiftest_pl),       intent(inout) :: pl       !! Helio massive body particle data structure
         real(DP),                 intent(in)    :: t        !! Current time
         real(DP), dimension(:,:), intent(in)    :: xhp     !! Heliocentric positions of planets
         ! Internals
         integer(I4B)              :: i, j
         real(DP)                  :: r2, fac
         real(DP), dimension(NDIM) :: dx
      
         associate(ntp => tp%nbody, npl => pl%nbody)
            tp%ahi(:,:) = 0.0_DP
            do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               do j = 1, npl
                  dx(:) = tp%xh(:,i) - xhp(:,j)
                  r2 = dot_product(dx(:), dx(:))
                  fac = pl%Gmass(j) / (r2 * sqrt(r2))
                  tp%ahi(:,i) = tp%ahi(:,i) - fac * dx(:)
               end do
            end if
            end do
         end associate 
         return
      end subroutine helio_getacch_int_tp
   
end submodule s_helio_getacch
