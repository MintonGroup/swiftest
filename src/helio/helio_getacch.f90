submodule (helio_classes) s_helio_getacch
contains
module subroutine helio_getacch_pl(self, cb, config, t)
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of massive bodies
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_getacch.f90
   !! Adapted from Hal Levison's Swift routine helio_getacch.f
   use swiftest
   implicit none
   ! Arguments
   class(helio_pl),               intent(inout) :: self     !! Helio massive body particle data structure
   class(swiftest_cb),            intent(inout) :: cb       !! Swiftest central body particle data structure
   class(swiftest_configuration), intent(in)    :: config   !! Input collection of 
   real(DP),                      intent(in)    :: t        !! Current time
   ! Internals
   logical, save                    :: lmalloc = .true.
   integer(I4B)                     :: i
   real(DP)                       :: r2
   real(DP), dimension(:), allocatable, save    :: irh
   real(DP), dimension(:, :), allocatable, save :: xh_loc, aobl

   associate(npl => self%nbody)
      !if (lflag) then
         self%ahi(:,2:npl) = 0.0_DP
         call helio_getacch_int_pl(self, t)
      !end if
      !if (config%loblatecb) call self%obl_acc(cb) TODO: Fix this
      !else
      self%ah(:,:) = self%ahi(:,:)
      !end if
      if (config%lextra_force) call self%user_getacch(cb, config, t)
   end associate

   return
   end subroutine helio_getacch_pl

   module subroutine helio_getacch_tp(self, cb, pl, config, t, xh)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_getacch_tp.f
      use swiftest
      implicit none
      ! Arguments
      class(helio_tp),               intent(inout) :: self   !! Helio test particle data structure
      class(swiftest_cb),                 intent(inout) :: cb     !! Swiftest central body particle data structuree 
      class(whm_pl),                 intent(inout) :: pl     !! WHM massive body particle data structure. 
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      ! Internals
      logical, save                                :: lmalloc = .true.
      integer(I4B)                                 :: i
      real(DP)                                     :: r2, mu
      real(DP), dimension(:), allocatable, save    :: irh, irht
      real(DP), dimension(:, :), allocatable, save :: aobl, xht, aoblt
   
   ! executable code
      associate(ntp => self%nbody, npl => pl%nbody)
         !if (lflag) then
            self%ahi(:,:) = 0.0_DP
            call helio_getacch_int_tp(self, pl, t, xh)
         !end if
         !if (config%loblatecb) call self%obl_acc(cb) TODO: Fix this
         self%ah(:,:) = self%ahi(:,:)
         if (config%lextra_force) call self%user_getacch(cb, config, t)
         if (config%lgr) call self%gr_getacch(cb, config)
      end associate
      return
   end subroutine helio_getacch_tp

   module subroutine helio_getacch_int_pl(self, t)
      !! author: David A. Minton
      !!
      !! Compute direct cross term heliocentric accelerations of massive bodiese
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int.f90
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      use swiftest
      implicit none
      ! Arguments
      class(helio_pl),               intent(inout) :: self     !! Helio massive body particle data structure
      real(DP),                      intent(in)    :: t        !! Current time
      ! Internals
      integer(I4B)              :: i, j
      real(DP)                  :: rji2, irij3, faci, facj
      real(DP), dimension(NDIM) :: dx
   
      associate(npl => self%nbody)
         do i = 2, npl - 1
            do j = i + 1, npl
               dx(:) = self%xh(:,j) - self%xh(:,i)
               rji2 = dot_product(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               faci = self%Gmass(i) * irij3
               facj = self%Gmass(j) * irij3
               self%ahi(:,i) = self%ahi(:,i) + facj * dx(:)
               self%ahi(:,i) = self%ahi(:,j) - faci * dx(:)
            end do
         end do
      end associate
   
      return
   end subroutine helio_getacch_int_pl

   module subroutine helio_getacch_int_tp(self, pl, t, xh)
      !! author: David A. Minton
      !!
      !! Compute direct cross term heliocentric accelerations of test particles
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_int_tp.f90
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      use swiftest
      implicit none
      ! Arguments
      class(helio_tp),               intent(inout) :: self     !! Helio test particle data structure
      class(whm_pl),               intent(inout) :: pl       !! Helio massive body particle data structure
      real(DP),                      intent(in)    :: t        !! Current time
      real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      ! Internals
      integer(I4B)              :: i, j
      real(DP)                  :: r2, fac
      real(DP), dimension(NDIM) :: dx
   
      associate(ntp => self%nbody, npl => pl%nbody)
         do i = 1, ntp
         if (self%status(i) == ACTIVE) then
            do j = 2, npl
               dx(:) = self%xh(:,i) - xh(:,j)
               r2 = dot_product(dx(:), dx(:))
               fac = pl%Gmass(j) / (r2 * sqrt(r2))
               self%ahi(:,i) = self%ahi(:,i) - fac * dx(:)
            end do
         end if
         end do
      end associate 
      return
   end subroutine helio_getacch_int_tp
   
end submodule s_helio_getacch
