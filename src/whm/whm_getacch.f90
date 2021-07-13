submodule(whm_classes) s_whm_getacch
   use swiftest
contains
   module subroutine whm_getacch_pl(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch.f90
      implicit none
      ! Arguments
      class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest central body particle data structure
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
      real(DP),                     intent(in)    :: t       !! Current time
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)                                :: i
      real(DP), dimension(NDIM)                   :: ah0

      associate(cb => system%cb, pl => self, npl => self%nbody)
         if (npl == 0) return
         call pl%set_ir3()

         ah0 = whm_getacch_ah0(pl%Gmass(2:npl), pl%xh(:,2:npl), npl-1)
         do i = 1, npl
            pl%ah(:, i) = ah0(:)
         end do

         call whm_getacch_ah1(cb, pl) 
         call whm_getacch_ah2(cb, pl) 
         call whm_getacch_ah3(pl)

         if (param%loblatecb) then
            cb%aoblbeg = cb%aobl
            call pl%accel_obl(system)
            cb%aoblend = cb%aobl
            if (param%ltides) then
               cb%atidebeg = cb%aobl
               call pl%accel_tides(system)
               cb%atideend = cb%atide
            end if
         end if

         if (param%lgr) call pl%accel_gr(param) 

         if (param%lextra_force) call pl%accel_user(system, param, t)
      end associate
      return
   end subroutine whm_getacch_pl

   module subroutine whm_getacch_tp(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest central body particle data structure
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
      real(DP),                     intent(in)    :: t      !! Current time
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      ! Internals
      integer(I4B)                                :: i
      real(DP), dimension(NDIM)                   :: ah0
      real(DP), dimension(:,:), allocatable       :: xhp
   
      associate(tp => self, ntp => self%nbody, pl => system%pl, cb => system%cb, npl => system%pl%nbody)
         if (ntp == 0 .or. npl == 0) return
         if (present(lbeg)) system%lbeg = lbeg

         if (system%lbeg) then
            allocate(xhp, source=pl%xbeg)
         else
            allocate(xhp, source=pl%xend)
         end if

         ah0(:) = whm_getacch_ah0(pl%Gmass(:), xhp(:,:), npl)
         do i = 1, ntp
            tp%ah(:, i) = ah0(:)
         end do
         call whm_getacch_ah3_tp(system, xhp)
         if (param%loblatecb) call tp%accel_obl(system)
         if (param%lextra_force) call tp%accel_user(system, param, t)
         if (param%lgr) call tp%accel_gr(param) 
      end associate
      return
   end subroutine whm_getacch_tp

   function whm_getacch_ah0(mu, xhp, n) result(ah0)
      !! author: David A. Minton
      !!
      !! Compute zeroth term heliocentric accelerations of planets 
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)         :: mu
      real(DP), dimension(:,:), intent(in)         :: xhp
      integer(I4B),             intent(in)         :: n
      ! Result
      real(DP), dimension(NDIM)                    :: ah0
      ! Internals
      real(DP)                                     :: fac, r2, ir3h, irh
      integer(I4B)                                 :: i

      ah0(:) = 0.0_DP
      do i = 1, n
         r2 = dot_product(xhp(:, i), xhp(:, i))
         irh = 1.0_DP / sqrt(r2)
         ir3h = irh / r2
         fac = mu(i) * ir3h 
         ah0(:) = ah0(:) - fac * xhp(:, i)
      end do

      return
   end function whm_getacch_ah0

   pure subroutine whm_getacch_ah1(cb, pl)
      !! author: David A. Minton
      !!
      !! Compute first term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah1.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah1.f90
      implicit none
      ! Arguments
      class(swiftest_cb), intent(in)    :: cb !! WHM central body object
      class(whm_pl),      intent(inout) :: pl !! WHM massive body object
      ! Internals
      integer(I4B)                 :: i
      real(DP), dimension(NDIM)    :: ah1h, ah1j

      associate(npl => pl%nbody)
         do i = 2, npl
            ah1j(:) = pl%xj(:, i) * pl%ir3j(i)
            ah1h(:) = pl%xh(:, i) * pl%ir3h(i)
            pl%ah(:, i) = pl%ah(:, i) + cb%Gmass * (ah1j(:) - ah1h(:))
         end do
      end associate
   
      return
   
   end subroutine whm_getacch_ah1

   pure subroutine whm_getacch_ah2(cb, pl)
      !! author: David A. Minton
      !!
      !! Compute second term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah2.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah2.f90
      implicit none
      ! Arguments
      class(swiftest_cb), intent(in)    :: cb !! Swiftest central body object
      class(whm_pl),      intent(inout) :: pl !! WHM massive body object
      ! Internals
      integer(I4B)                 :: i
      real(DP)                     :: etaj, fac
      real(DP), dimension(NDIM)    :: ah2, ah2o
   
      associate(npl => pl%nbody)
         ah2(:) = 0.0_DP
         ah2o(:) = 0.0_DP
         etaj = cb%Gmass
         do i = 2, npl
            etaj = etaj + pl%Gmass(i - 1)
            fac = pl%Gmass(i) * cb%Gmass * pl%ir3j(i) / etaj
            ah2(:) = ah2o + fac * pl%xj(:, i)
            pl%ah(:,i) = pl%ah(:, i) + ah2(:)
            ah2o(:) = ah2(:)
         end do
      end associate
   
      return
   end subroutine whm_getacch_ah2

   pure subroutine whm_getacch_ah3(pl)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
      implicit none

      class(whm_pl),           intent(inout) :: pl
      integer(I4B)                           :: i, j
      real(DP)                               :: rji2, irij3, faci, facj
      real(DP), dimension(NDIM)              :: dx
      real(DP), dimension(:,:), allocatable  :: ah3
   
      associate(npl => pl%nbody)
         allocate(ah3, mold=pl%ah)
         ah3(:, :) = 0.0_DP

         do i = 1, npl - 1
            do j = i + 1, npl
               dx(:) = pl%xh(:, j) - pl%xh(:, i)
               rji2  = dot_product(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               faci = pl%Gmass(i) * irij3
               facj = pl%Gmass(j) * irij3
               ah3(:, i) = ah3(:, i) + facj * dx(:)
               ah3(:, j) = ah3(:, j) - faci * dx(:)
            end do
         end do
         do i = 1, NDIM
            pl%ah(i, 1:npl) = pl%ah(i, 1:npl) + ah3(i, 1:npl)
         end do
         deallocate(ah3)
      end associate
   
      return
   end subroutine whm_getacch_ah3

   pure subroutine whm_getacch_ah3_tp(system, xhp) 
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
      implicit none
      ! Arguments
      class(swiftest_nbody_system),  intent(inout) :: system !! WHM nbody system object
      real(DP), dimension(:,:),      intent(in)    :: xhp    !! Heliocentric positions of planets at the current substep
      ! Internals
      integer(I4B)                         :: i, j
      real(DP)                             :: rji2, irij3, fac
      real(DP), dimension(NDIM)            :: dx, acc

      associate(ntp => system%tp%nbody, npl => system%pl%nbody, tp => system%tp, pl => system%pl) 
         if (ntp == 0) return
         do i = 1, ntp
            acc(:) = 0.0_DP
            do j = 1, npl
               dx(:) = tp%xh(:, i) - xhp(:, j)
               rji2 = dot_product(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               fac = pl%Gmass(j) * irij3
               acc(:) = acc(:) - fac * dx(:)
            end do
            tp%ah(:, i) = tp%ah(:, i) + acc(:)
         end do
      end associate
      return
   end subroutine whm_getacch_ah3_tp
end submodule s_whm_getacch
