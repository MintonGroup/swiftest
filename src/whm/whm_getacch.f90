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
   logical, save :: lmalloc = .false.  
   real(DP), dimension(:), allocatable, save    :: irh, irj, ir3h, ir3j, fac
   real(DP), dimension(NDIM) :: ah0
   real(DP) :: r2

   associate(pl => self, npl => self%nbody, Gmpl => self%Gmass, xj => self%xj, xh => self%xh, &
             ah => self%ah, ah1 => self%ah1, ah2 => self%ah2, ah3 => self%ah3, &
             j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, aobl => self%aobl, aobl0 => cb%aobl)
      if (npl == 0) return
      if (.not. allocated(irj)) allocate(irj(npl))
      if (.not. allocated(ir3j)) allocate(ir3j(npl))
      if (.not. allocated(irh)) allocate(irh(npl))
      if (.not. allocated(ir3h)) allocate(ir3h(npl))
      if (.not. allocated(fac)) allocate(fac(npl))
      !do concurrent(i = 1:npl)
      do i = 1, npl
         r2 = dot_product(xj(:, i), xj(:, i))
         irj(i)= 1.0_DP / sqrt(r2)
         ir3j(i) = irj(i) / r2

         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i)= 1.0_DP / sqrt(r2)
         ir3h(i) = irh(i) / r2
      end do

      ah0(1) = 0.0_DP
      fac(2:npl) = Gmpl(2:npl) * ir3h(2:npl)
      !do concurrent (i = 1:NDIM)
      do i = 1, NDIM
         ah0(i) = -sum(fac(2:npl) * xh(i, 2:npl))
      end do
      call whm_getacch_ah1(cb, pl, ir3h, ir3j)
      call whm_getacch_ah2(cb, pl, ir3j)
      call whm_getacch_ah3(pl)
      
      !do concurrent (i = 1:NDIM)
      do i = 1, NDIM
         ah(i, 1:npl) = ah0(i) + ah1(i, 1:npl) + ah2(i, 1:npl) + ah3(i, 1:npl)
      end do

      if (j2rp2 /= 0.0_DP) then
         call  self%obl_acc(cb, irh)

         !do concurrent (i = 1:NDIM)
         do i = 1, NDIM
            ah(i, 1:npl) = ah(i, 1:npl) + aobl(i, 1:npl) - aobl0(i)
         end do
      end if
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
      real(DP), dimension(:), allocatable, save    :: irh, ir3h
      real(DP), dimension(:), allocatable, save    :: irht, r2p, r2t
      real(DP), dimension(NDIM)                    :: ah0
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody, aht => self%ah, &
                status => self%status, xht => self%xh, Gmpl => pl%Gmass, xh => pl%xh,&
                j2rp2 => cb%j2rp2, j4rp4 => cb%j4rp4, aobl => self%aobl, aobl0 => cb%aobl)
         if (ntp == 0 .or. npl == 0) return
         if(.not.allocated(irh)) allocate(irh(npl))
         if (.not.allocated(ir3h)) allocate(ir3h(npl))
         if (.not. allocated(irht)) allocate(irht(ntp))
         if (.not. allocated(r2p)) allocate(r2p(npl))
         if (.not. allocated(r2t)) allocate(r2t(ntp))

         !$omp simd
         do i = 1, npl
            r2p(i) = dot_product(xh(:, i), xh(:, i))
         end do
         irh(:)= 1.0_DP / sqrt(r2p(:))
         ir3h(:) = irh(:) / r2p(:)

         !$omp simd
         do i = 1, ntp
            r2t(i) = dot_product(xht(:, i), xht(:, i))
         end do
         irht(:) = 1.0_DP / sqrt(r2t(:))

         ah0(:) = 0.0_DP
         if(.not.allocated(fac)) allocate(fac(npl))
         fac(:) = Gmpl(:) * ir3h(:)

         do concurrent (i = 1:NDIM)
            ah0(i) = -sum(fac(:) * xh(i, :))
         end do
         call whm_getacch_ah3_tp(cb, pl, tp)
         do concurrent (i = 1:ntp)
            aht(:, i) = aht(:, i) + ah0(:)
         end do
         if (j2rp2 /= 0.0_DP) then
            call tp%obl_acc(cb, irht)
            do concurrent (i = 1:ntp)
               aht(:, i) = aht(:, i) + tp%aobl(:, i) - aobl0(:)
            end do
         end if
         if (config%lextra_force) call tp%user_getacch(cb, config, t)
         if (config%lgr) call tp%gr_getacch(cb, config) 
      end associate
      return
   end procedure whm_getacch_tp

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
         do concurrent (i = 2:pl%nbody) !local(ah1j, ah1h) shared(pl, cb, xj, xh, ir3h, ir3j, msun, ah1)
         !do i = 2, npl
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
         ah2(:, 1) = 0.0_DP
         etaj = Gmsun
         do i = 2, npl
            etaj = etaj + Gmpl(i - 1)
            fac = Gmpl(i) * Gmsun * ir3j(i) / etaj
            ah2(:, i) = ah2(:, i - 1) + fac * xj(:, i)
         end do
      end associate
   
      return
   end subroutine whm_getacch_ah2

   subroutine whm_getacch_ah3(pl)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
      use swiftest
      implicit none

      class(whm_pl),           intent(inout) :: pl
      real(DP) :: rji2
   
      associate(npl => pl%nbody, xh => pl%xh, ah3 => pl%ah3, Gmpl => pl%Gmass, &
         nk => pl%num_comparisons, k_plpl => pl%k_eucl, irij3 => pl%irij3) 
         ah3(:, 1:npl) = 0.0_DP

         call pl%eucl_irij3()
         call ah3p(nk, k_plpl, xh, irij3, Gmpl, ah3)

      end associate
   
      return
   end subroutine whm_getacch_ah3

   subroutine ah3p(nk, k_plpl, xh, irij3, Gmpl, ah3)
      use swiftest
      implicit none
      integer(I4B), intent(in) :: nk
      integer(I4B), dimension(:,:), intent(in) :: k_plpl
      real(DP), dimension(:,:), intent(in) :: xh
      real(DP), dimension(:), intent(in) :: irij3, Gmpl
      real(DP), dimension(:,:), intent(out) :: ah3

      integer(I4B)                           :: i, j, k
      real(DP)                               :: rji2, faci, facj
      real(DP), dimension(NDIM)              :: dx

      !!$omp parallel do schedule(static) default(private) &
      !!$omp shared (nk, k_plpl, xh, irij3, Gmpl) &
      !!$omp reduction(+:ah3)
      !$omp simd
      do k = 1, nk
         i = k_plpl(1, k)
         j = k_plpl(2, k)
         dx(:) = xh(:, j) - xh(:, i)
         faci = Gmpl(i) * irij3(k)
         facj = Gmpl(j) * irij3(k)
         ah3(:, i) = ah3(:, i) + facj * dx(:)
         ah3(:, j) = ah3(:, j) - faci * dx(:)
      end do

      return
   end subroutine ah3p

   subroutine whm_getacch_ah3_tp(cb, pl, tp) 
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
      integer(I4B)          :: i, j
      real(DP)            :: irij3, rji2, fac
      real(DP), dimension(NDIM) :: dx, acc

      associate(ntp => tp%nbody, npl => pl%nbody, msun => cb%Gmass,  Gmpl => pl%Gmass, &
                  xht => tp%xh, xhp => pl%xh, aht => tp%ah)
   
         if (ntp == 0) return

         !dir$ parallel private(acc, dx, rji2, irij3, fac)
         do j = 1, ntp
            acc(:) = 0.0_DP
            !$omp simd 
            do i = 1, npl
               dx(:) = xht(:, j) - xhp(:, i)
               rji2  = dot_product(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2 * sqrt(rji2)) 
               fac = Gmpl(i) * irij3
               acc(:) = acc(:) - fac * dx(:) 
            end do
            aht(:, j) = acc(:)
         end do

      end associate
      return
   end subroutine whm_getacch_ah3_tp

end submodule whm_getacch_implementations
