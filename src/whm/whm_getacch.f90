submodule(whm_classes) s_whm_getacch
   use swiftest
contains
   module subroutine whm_getacch_pl(self, cb, config, t)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch.f90
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self     !! WHM massive body particle data structure
      class(swiftest_cb),            intent(inout) :: cb  !! Swiftest central body particle data structure
      class(swiftest_configuration), intent(in)    :: config   !! Input collection of 
      real(DP),                      intent(in)    :: t        !! Current time
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: ah0

      associate(pl => self, npl => self%nbody, j2rp2 => cb%j2rp2, &
         ah => self%ah, xh => self%xh, xj => self%xj, vh => self%vh, vj => self%vj)
         if (npl == 0) return
         call pl%set_ir3()

         ah0 = whm_getacch_ah0(pl%Gmass(2:npl), pl%xh(:,2:npl))
         do i = 1, NDIM
            pl%ah(i, 1:npl) = ah0(i)
         end do
         call whm_getacch_ah1(cb, pl) 
         call whm_getacch_ah2(cb, pl) 
         call whm_getacch_ah3(pl)

         if (config%loblatecb) call pl%obl_acc(cb)
         if (config%lextra_force) call pl%user_getacch(cb, config, t)
         if (config%lgr) call pl%gr_getacch(cb, config) 

      end associate
      return
   end subroutine whm_getacch_pl

   module subroutine whm_getacch_tp(self, cb, pl, config, t, xh)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                 intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Generic Swiftest central body particle data structuree 
      class(whm_pl),                 intent(inout) :: pl     !! Generic Swiftest massive body particle data structure. 
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      ! Internals
      integer(I4B)                                 :: i
      real(DP), dimension(NDIM)                    :: ah0
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody, j2rp2 => cb%j2rp2, aht => self%ah, &
               ir3h => pl%ir3h, GMpl => pl%Gmass)
         if (ntp == 0 .or. npl == 0) return

         ah0 = whm_getacch_ah0(pl%Gmass(:), xh(:,:))
         do i = 1, ntp
            tp%ah(:, i) = ah0(:)
         end do
         call whm_getacch_ah3_tp(cb, pl, tp, xh)
         if (config%loblatecb) call tp%obl_acc(cb)
         if (config%lextra_force) call tp%user_getacch(cb, config, t)
         if (config%lgr) call tp%gr_getacch(cb, config) 
      end associate
      return
   end subroutine whm_getacch_tp

   pure function whm_getacch_ah0(mu, xh) result(ah0)
      !! author: David A. Minton
      !!
      !! Compute zeroth term heliocentric accelerations of planets 
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(in)           :: mu
      real(DP), dimension(:,:), intent(in)         :: xh
      ! Result
      real(DP), dimension(NDIM)                    :: ah0
      ! Internals
      real(DP)                                     :: fac, r2, ir3h
      integer(I4B)                                 :: i, n

      n = size(mu)

      ah0(:) = 0.0_DP
      do i = 1, n
         r2 = dot_product(xh(:, i), xh(:, i))
         ir3h = 1.0_DP / (r2 * sqrt(r2))
         fac = mu(i) * ir3h 
         ah0(:) = ah0(:) - fac * xh(:, i)
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
      class(swiftest_cb), intent(in)  :: cb !! Swiftest central body object
      class(whm_pl), intent(inout)    :: pl !! WHM massive body object
      ! Internals
      integer(I4B)                    :: i
      real(DP), dimension(NDIM)       :: ah1h, ah1j

      associate(npl => pl%nbody, msun => cb%Gmass, xh => pl%xh, xj => pl%xj, ir3j => pl%ir3j, ir3h => pl%ir3h )
         do i = 2, npl
            ah1j(:) = xj(:, i) * ir3j(i)
            ah1h(:) = xh(:, i) * ir3h(i)
            pl%ah(:, i) = pl%ah(:, i) + msun * (ah1j(:) - ah1h(:))
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
      class(swiftest_cb), intent(in)     :: cb !! Swiftest central body object
      class(whm_pl),      intent(inout)  :: pl !! WHM massive body object
      ! Internals
      integer(I4B)                                 :: i
      real(DP)                                     :: etaj, fac
      real(DP), dimension(NDIM)                    :: ah2, ah2o
   
      associate(npl => pl%nbody, Gmsun => cb%Gmass, xh => pl%xh, xj => pl%xj, Gmpl => pl%Gmass, ir3j => pl%ir3j)
         ah2(:) = 0.0_DP
         ah2o(:) = 0.0_DP
         etaj = Gmsun
         do i = 2, npl
            etaj = etaj + Gmpl(i - 1)
            fac = Gmpl(i) * Gmsun * ir3j(i) / etaj
            ah2(:) = ah2o + fac * xj(:, i)
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

      class(whm_pl),           intent(inout)       :: pl
      integer(I4B)                                 :: i, j
      real(DP)                                     :: rji2, irij3, faci, facj
      real(DP), dimension(NDIM)                    :: dx
      real(DP), dimension(:,:), allocatable        :: ah3
   
      associate(npl => pl%nbody, xh => pl%xh, Gmpl => pl%Gmass) 
         allocate(ah3, mold=pl%ah)
         ah3(:, 1:npl) = 0.0_DP

         do i = 1, npl - 1
            do j = i + 1, npl
               dx(:) = xh(:, j) - xh(:, i)
               rji2  = dot_product(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               faci = Gmpl(i) * irij3
               facj = Gmpl(j) * irij3
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

   pure subroutine whm_getacch_ah3_tp(cb, pl, tp, xh) 
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
      implicit none
      ! Arguments
      class(swiftest_cb), intent(in)               :: cb  !! Swiftest central body object
      class(whm_pl), intent(in)                    :: pl  !! WHM massive body object
      class(whm_tp), intent(inout)                 :: tp  !! WHM test particle object
      real(DP), dimension(:,:), intent(in)         :: xh  !! Position vector of massive bodies at required point in step
      ! Internals
      integer(I4B)                                 :: i, j
      real(DP)                                     :: rji2, irij3, fac
      real(DP), dimension(NDIM)                    :: dx
      real(DP), dimension(:,:), allocatable        :: aht

      associate(ntp => tp%nbody, npl => pl%nbody, msun => cb%Gmass,  Gmpl => pl%Gmass, &
                  xht => tp%xh)
  
         allocate(aht, source=tp%ah)
         if (ntp == 0) return
         do j = 1, npl
            !$omp simd private(dx,rji2,irij3,fac) reduction(-:aht)
            do i = 1, tp%nbody
               dx(:) = tp%xh(:, i) - xh(:, j)
               rji2 = dot_product(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               fac = pl%Gmass(j) * irij3
               aht(:, i) = aht(:, i) - fac * dx(:)
            end do
         end do
         call move_alloc(aht, tp%ah)
      end associate
      return
   end subroutine whm_getacch_ah3_tp
end submodule s_whm_getacch
