submodule(whm_classes) s_whm_getacch
contains
   module procedure whm_getacch_pl  ! (self, cb, config, t)
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of planets
   !!
   !! Adapted from Hal Levison's Swift routine getacch.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_getacch.f90
   use swiftest
   implicit none
   logical(lgt), save                 :: lmalloc = .true.
   integer(I4B)                     :: i
   real(DP)                       ::  fac
   real(DP), dimension(:), allocatable :: r2, fac
   real(DP), dimension(:), allocatable, save    :: irh, irj, ir3h, ir3j
   real(DP), dimension(:, :), allocatable, save :: xh, aobl
   real(DP), dimension(NDIM) :: ah0

   associate(n => self%nbody, xj => self%xj, xh => self%xh)
      !if (lmalloc) then
      !   allocate(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax), irj(nplmax), ir3h(nplmax), ir3j(nplmax))
      !   lmalloc = .false.
      !end if
      !swifter_pl1p => whm_pl1p%swifter
      !whm_plp => whm_pl1p
      allocate(r2(n))
      r2(1:n) = xj(1:n, 1:NDIM) .dot. xj(1:n, 1:NDIM)
      irj(1:n)= 1.0_DP / sqrt(r2(1:n))
      ir3j(1:n) = irj(1:n) / r2(1:n)

      r2(1:n) = xh(1:n, 1:NDIM) .dot. xh(1:n, 1:NDIM)
      irh(1:n)= 1.0_DP / sqrt(r2(1:n))
      ir3h(1:n) = irh(1:n) / r2(1:n)

      ah0(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)

      fac(1:n) = self%Gmass(1:n) * ir3h(1:n)
      do concurrent (i = 1:NDIM)
         ah0(i) - sum(fac(1:n) * xh(1:n, i))
      end do
      call whm_getacch_ah1(npl, whm_pl1p, ir3h, ir3j)
      call whm_getacch_ah2(npl, whm_pl1p, ir3j)
      call whm_getacch_ah3(npl, whm_pl1p)
      whm_plp => whm_pl1p
      whm_plp%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      do i = 2, npl
         whm_plp => whm_plp%nextp
         whm_plp%ah(:) = ah0(:) + whm_plp%ah1(:) + whm_plp%ah2(:) + whm_plp%ah3(:)
      end do
      if (j2rp2 /= 0.0_DP) then
         swifter_plp => swifter_pl1p
         do i = 2, npl
            swifter_plp => swifter_plp%nextp
            xh(:, i) = swifter_plp%xh(:)
         end do
         call obl_acc(npl, swifter_pl1p, j2rp2, j4rp4, xh, irh, aobl)
         whm_plp => whm_pl1p
         do i = 2, npl
            whm_plp => whm_plp%nextp
            whm_plp%ah(:) = whm_plp%ah(:) + aobl(:, i) - aobl(:, 1)
         end do
      end if
      if (lextra_force) call whm_user_getacch(t, npl, whm_pl1p)
      call gr_whm_getacch(npl, whm_pl1p, c2)

      deallocate(r2)
   end associate
   return

   end procedure whm_getacch_pl

   module procedure whm_getacch_tp !(cb, pl, config, t, dt, xh)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_tp.f90
      use swiftest
      implicit none
      logical(lgt), save                 :: lmalloc = .true.
      integer(I4B)                     :: i
      real(DP)                       :: r2, fac, mu
      real(DP), dimension(:), allocatable, save    :: irh, ir3h
      real(DP), dimension(:), allocatable, save    :: irht
      real(DP), dimension(:, :), allocatable, save :: aobl
      real(DP), dimension(:, :), allocatable, save :: xht, aoblt
      type(swifter_pl), pointer            :: swifter_pl1p, swifter_plp
      type(swifter_tp), pointer            :: swifter_tp1p, swifter_tpp
      type(whm_pl), pointer                :: whm_plp
      type(whm_tp), pointer                :: whm_tpp
   
   ! executable code
      if (lmalloc) then
         allocate(aobl(ndim, nplmax), irh(nplmax), ir3h(nplmax), xht(ndim, ntpmax), aoblt(ndim, ntpmax), irht(ntpmax))
         lmalloc = .false.
      end if
      swifter_pl1p => whm_pl1p%swifter
      swifter_tp1p => whm_tp1p%swifter
      do i = 2, npl
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP/sqrt(r2)
         ir3h(i) = irh(i)/r2
      end do
      swifter_tpp => swifter_tp1p
      do i = 1, ntp
         xht(:, i) = swifter_tpp%xh(:)
         r2 = dot_product(xht(:, i), xht(:, i))
         irht(i) = 1.0_DP/sqrt(r2)
         swifter_tpp => swifter_tpp%nextp
      end do
      ah0(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      swifter_plp => swifter_pl1p
      do i = 2, npl
         swifter_plp => swifter_plp%nextp
         fac = swifter_plp%mass*ir3h(i)
         ah0(:) = ah0(:) - fac*xh(:, i)
      end do
      call whm_getacch_ah3_tp(npl, ntp, whm_pl1p, whm_tp1p, xh)
      whm_tpp => whm_tp1p
      do i = 1, ntp
         whm_tpp%ah(:) = whm_tpp%ah(:) + ah0(:)
         whm_tpp => whm_tpp%nextp
      end do
      if (j2rp2 /= 0.0_DP) then
         call obl_acc(npl, swifter_pl1p, j2rp2, j4rp4, xh, irh, aobl)
         mu = whm_pl1p%swifter%mass
         call obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
         whm_tpp => whm_tp1p
         do i = 1, ntp
            whm_tpp%ah(:) = whm_tpp%ah(:) + aoblt(:, i) - aobl(:, 1)
            whm_tpp => whm_tpp%nextp
         end do
      end if
      if (lextra_force) call whm_user_getacch_tp(t, ntp, whm_tp1p)
   
      return
   
   end procedure whm_getacch_tp

   module procedure whm_getacch_ah1 !(npl, whm_pl1p, ir3h, ir3j)
      !! author: David A. Minton
      !!
      !! Compute first term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah1.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah1.f90
      use swiftest
      implicit none
      integer(I4B)          :: i
      real(DP)            :: msun
      real(DP), dimension(ndim) :: ah1h, ah1j
      type(whm_pl), pointer   :: whm_plp
   
   ! executable code
      msun = whm_pl1p%swifter%mass
      whm_pl1p%ah1(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      if (npl > 1) then
         whm_plp => whm_pl1p%nextp
         whm_plp%ah1(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      end if
      do i = 3, npl
         whm_plp => whm_plp%nextp
         ah1j(:) = whm_plp%xj(:)*ir3j(i)
         ah1h(:) = whm_plp%swifter%xh(:)*ir3h(i)
         whm_plp%ah1(:) = msun*(ah1j(:) - ah1h(:))
      end do
   
      return
   
   end procedure whm_getacch_ah1

   module procedure whm_getacch_ah2 !(npl, whm_pl1p, ir3j)
      !! author: David A. Minton
      !!
      !! Compute second term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah2.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah2.f90
      use swiftest
      implicit none
      integer(I4B)      :: i
      real(DP)          :: etaj, fac, msun
      type(whm_pl), pointer :: whm_plp, whm_plop
   
   ! executable code
      msun = whm_pl1p%swifter%mass
      whm_pl1p%ah2(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      if (npl > 1) then
         whm_plp => whm_pl1p%nextp
         whm_plp%ah2(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
         etaj = msun
      endif
      do i = 3, npl
         whm_plop => whm_plp
         whm_plp => whm_plp%nextp
         etaj = etaj + whm_plop%swifter%mass
         fac = whm_plp%swifter%mass*msun*ir3j(i)/etaj
         whm_plp%ah2(:) = whm_plop%ah2(:) + fac*whm_plp%xj(:)
      end do
   
      return
   
   end procedure whm_getacch_ah2

   module procedure whm_getacch_ah3 !(npl, whm_pl1p)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of planets
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
      use swiftest
      implicit none
      integer(I4B)          :: i, j
      real(DP)            :: rji2, irij3, faci, facj
      real(DP), dimension(ndim) :: dx
      type(whm_pl), pointer   :: whm_plip, whm_pljp
   
   ! executable code
      whm_plip => whm_pl1p
      do i = 1, npl
         whm_plip%ah3(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
         whm_plip => whm_plip%nextp
      end do
      whm_plip => whm_pl1p
      do i = 2, npl - 1
         whm_plip => whm_plip%nextp
         whm_pljp => whm_plip
         do j = i + 1, npl
            whm_pljp => whm_pljp%nextp
            dx(:) = whm_pljp%swifter%xh(:) - whm_plip%swifter%xh(:)
            rji2 = dot_product(dx(:), dx(:))
            irij3 = 1.0_DP/(rji2*sqrt(rji2))
            faci = whm_plip%swifter%mass*irij3
            facj = whm_pljp%swifter%mass*irij3
            whm_plip%ah3(:) = whm_plip%ah3(:) + facj*dx(:)
            whm_pljp%ah3(:) = whm_pljp%ah3(:) - faci*dx(:)
         end do
      end do
   
      return
   
   end procedure whm_getacch_ah3


   module procedure whm_getacch_ah3_tp !(npl, ntp, whm_pl1p, whm_tp1p, xh)
      !! author: David A. Minton
      !!
      !! Compute direct cross (third) term heliocentric accelerations of test particles
      !!
      !! Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_getacch_ah3.f90
      use swiftest
      implicit none
      integer(I4B)          :: i, j
      real(DP)            :: rji2, irij3, fac
      real(DP), dimension(ndim) :: dx, acc, xht
      type(whm_pl), pointer   :: whm_plp
      type(whm_tp), pointer   :: whm_tpp
   
   ! executable code
      whm_tpp => whm_tp1p
      do i = 1, ntp
         acc(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
         xht(:) = whm_tpp%swifter%xh(:)
         whm_plp => whm_pl1p
         do j = 2, npl
            whm_plp => whm_plp%nextp
            dx(:) = xht(:) - xh(:, j)
            rji2 = dot_product(dx(:), dx(:))
            irij3 = 1.0_DP/(rji2*sqrt(rji2))
            fac = whm_plp%swifter%mass*irij3
            acc(:) = acc(:) - fac*dx(:)
         end do
         whm_tpp%ah(:) = acc(:)
         whm_tpp => whm_tpp%nextp
      end do
   
      return
   
   end procedure whm_getacch_ah3_tp
end submodule s_whm_getacch
