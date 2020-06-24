submodule(whm) s_whm_getacch
contains
   module procedure whm_getacch(lextra_force, t, npl, nplmax, whm_pl1p, j2rp2, j4rp4, c2)
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
   real(DP)                       :: r2, fac
   real(DP), dimension(:), allocatable, save    :: irh, irj, ir3h, ir3j
   real(DP), dimension(:, :), allocatable, save :: xh, aobl
   type(swifter_pl), pointer            :: swifter_pl1p, swifter_plp
   type(whm_pl), pointer                :: whm_plp

! executable code
   if (lmalloc) then
      allocate(xh(ndim, nplmax), aobl(ndim, nplmax), irh(nplmax), irj(nplmax), ir3h(nplmax), ir3j(nplmax))
      lmalloc = .false.
   end if
   swifter_pl1p => whm_pl1p%swifter
   whm_plp => whm_pl1p
   do i = 2, npl
      whm_plp => whm_plp%nextp
      r2 = dot_product(whm_plp%xj(:), whm_plp%xj(:))
      irj(i) = 1.0_DP/sqrt(r2)
      ir3j(i) = irj(i)/r2
   end do
   swifter_plp => swifter_pl1p
   do i = 2, npl
      swifter_plp => swifter_plp%nextp
      r2 = dot_product(swifter_plp%xh(:), swifter_plp%xh(:))
      irh(i) = 1.0_DP/sqrt(r2)
      ir3h(i) = irh(i)/r2
   end do
   ah0(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   swifter_plp => swifter_pl1p%nextp
   do i = 3, npl
      swifter_plp => swifter_plp%nextp
      fac = swifter_plp%mass*ir3h(i)
      ah0(:) = ah0(:) - fac*swifter_plp%xh(:)
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

   return

   end procedure whm_getacch
end submodule s_whm_getacch
