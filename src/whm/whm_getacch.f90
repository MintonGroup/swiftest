submodule(whm_classes) whm_getacch_implementations
contains
   module procedure whm_getacch_pl  
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of planets
   !!
   !! Adapted from Hal Levison's Swift routine getacch.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch.f90
   use swiftest
   implicit none
   integer(I4B)                     :: i
   real(DP), dimension(:), allocatable, save    :: fac
   real(DP), dimension(NDIM) :: ah0
   real(DP) :: r2

   associate(pl => self, npl => self%nbody, j2rp2 => cb%j2rp2, &
      ah1 => self%ah1, ah2 => self%ah2, ah3 => self%ah3, ah => self%ah, &
      xh => self%xh, xj => self%xj, vh => self%vh, vj => self%vj)
      if (npl == 0) return
      call pl%set_ir3()

      ah0 = whm_getacch_ah0(pl%Gmass(2:npl), pl%xh(:,2:npl))
      call whm_getacch_ah1(cb, pl, pl%ir3h, pl%ir3j)
      call whm_getacch_ah2(cb, pl, pl%ir3j)
      call whm_getacch_ah3(pl)
      do i = 1, NDIM
         pl%ah(i, 1:npl) = ah0(i) + pl%ah1(i, 1:npl) + pl%ah2(i, 1:npl) + pl%ah3(i, 1:npl)
      end do

      if (j2rp2 /= 0.0_DP) call self%obl_acc(cb)
      if (config%lextra_force) call pl%user_getacch(cb, config, t)
      if (config%lgr) call pl%gr_getacch(cb, config) 

   end associate
   return
   end procedure whm_getacch_pl

   module procedure whm_getacch_tp 
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_tp.f90
      use swiftest
      implicit none
      integer(I4B)                                 :: i
      real(DP), dimension(:), allocatable, save    :: fac
      real(DP), dimension(NDIM)                    :: ah0
      real(DP)                                     :: r2
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody, j2rp2 => cb%j2rp2, aht => self%ah, &
               ir3h => pl%ir3h, GMpl => pl%Gmass)
         if (ntp == 0 .or. npl == 0) return

         ah0 = whm_getacch_ah0(pl%Gmass(:), xh(:,:))
         call whm_getacch_ah3_tp(cb, pl, tp, xh)
         do i = 1, ntp
            tp%ah(:, i) = tp%ah(:, i) + ah0(:)
         end do
         if (j2rp2 /= 0.0_DP) call tp%obl_acc(cb)
         if (config%lextra_force) call tp%user_getacch(cb, config, t)
         if (config%lgr) call tp%gr_getacch(cb, config) 
      end associate
      return
   end procedure whm_getacch_tp

   pure function whm_getacch_ah0(mu, xh) result(ah0)
      !! author: David A. Minton
      !!
      !! Compute zeroth term heliocentric accelerations of planets 
      use swiftest
      implicit none
      real(DP), dimension(:), intent(in)    :: mu
      real(DP), dimension(:,:), intent(in)  :: xh
      real(DP) :: fac, r2, irh, ir3h
      real(DP), dimension(NDIM) :: ah0
      integer(I4B) :: i, n

      n = size(mu)

      ah0(:) = 0.0_DP
      do i = 1, n
         r2 = dot_product(xh(:, i), xh(:, i))
         irh = 1.0_DP / sqrt(r2) 
         ir3h = irh / r2
         fac = mu(i) * ir3h !/ (norm2(xh(:, i)))**3
         ah0(:) = ah0(:) - fac * xh(:, i)
      end do

      return
   end function whm_getacch_ah0

   pure subroutine whm_getacch_ah1(cb, pl, ir3h, ir3j)
      !! author: David A. Minton
      !!
      !! Compute first term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah1.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah1.f90
      use swiftest
      implicit none
      class(whm_cb), intent(in) :: cb
      class(whm_pl), intent(inout)        :: pl
      real(DP), dimension(:), intent(in) :: ir3h, ir3j

      integer(I4B)              :: i
      real(DP), dimension(NDIM) :: ah1h, ah1j

      associate(npl => pl%nbody, msun => cb%Gmass, xh => pl%xh, xj => pl%xj, ah1 => pl%ah1)
         ah1(:,:) = 0.0_DP
         do i = 2, npl
            ah1j(:) = xj(:, i) * ir3j(i)
            ah1h(:) = xh(:, i) * ir3h(i)
            ah1(:, i) = msun * (ah1j(:) - ah1h(:))
         end do
      end associate
   
      return
   
   end subroutine whm_getacch_ah1


   pure subroutine whm_getacch_ah2(cb, pl, ir3j) 
      !! author: David A. Minton
      !!
      !! Compute second term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah2.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah2.f90
      use swiftest
      implicit none

      class(whm_cb), intent(in)    :: cb
      class(whm_pl),           intent(inout) :: pl
      real(DP), dimension(:),  intent(in)    :: ir3j
      integer(I4B)                           :: i
      real(DP)                               :: etaj, fac
   
      associate(npl => pl%nbody, Gmsun => cb%Gmass, xh => pl%xh, xj => pl%xj, ah2 => pl%ah2, Gmpl => pl%Gmass)
         ah2(:, :) = 0.0_DP
         etaj = Gmsun
         do i = 2, npl
            etaj = etaj + Gmpl(i - 1)
            fac = Gmpl(i) * Gmsun * ir3j(i) / etaj
            ah2(:, i) = ah2(:, i - 1) + fac * xj(:, i)
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
      use swiftest
      implicit none

      class(whm_pl),           intent(inout) :: pl
      integer(I4B)                           :: i, j
      real(DP)                               :: rji2, irij3, faci, facj
      real(DP), dimension(NDIM)              :: dx
   
      associate(npl => pl%nbody, xh => pl%xh, ah3 => pl%ah3, Gmpl => pl%Gmass) 
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
      use swiftest
      implicit none
      class(whm_cb), intent(in) :: cb 
      class(whm_pl), intent(in) :: pl 
      class(whm_tp), intent(inout) :: tp
      real(DP), dimension(:,:), intent(in) :: xh
      integer(I4B)          :: i, j
      real(DP)            :: rji2, irij3, fac
      real(DP), dimension(NDIM) :: dx, acc

      associate(ntp => tp%nbody, npl => pl%nbody, msun => cb%Gmass,  Gmpl => pl%Gmass, &
                  xht => tp%xh, aht => tp%ah)
   
         if (ntp == 0) return
         aht(:,:) = 0.0_DP
         do i = 1, ntp
            do j = 1, npl
               dx(:) = xht(:, i) - xh(:, j)
               rji2 = dot_product(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2 * sqrt(rji2))
               fac = Gmpl(j) * irij3
               aht(:, i) = aht(:, i) - fac * dx(:)
            end do
         end do
      end associate
      return
   end subroutine whm_getacch_ah3_tp
end submodule whm_getacch_implementations
