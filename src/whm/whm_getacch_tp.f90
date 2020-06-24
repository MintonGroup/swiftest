submodule(whm) s_whm_getacch_tp
contains
   module procedure whm_getacch_tp(lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, xh, j2rp2, j4rp4)
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
end submodule s_whm_getacch_tp
