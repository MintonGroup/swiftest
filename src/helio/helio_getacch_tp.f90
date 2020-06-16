submodule (helio) s_helio_getacch_tp
contains
   module procedure helio_getacch_tp
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of test particles
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_tp.f90
   !! Adapted from Hal Levison's Swift routine helio_getacch_tp.f
   use swiftest
   logical, save                 :: lmalloc = .true.
   integer(I4B)                     :: i
   real(DP)                       :: r2, mu
   real(DP), dimension(:), allocatable, save    :: irh, irht
   real(DP), dimension(:, :), allocatable, save :: aobl, xht, aoblt

! executable code
   if (lflag) then
      helio_tpA%ahi(:,:) = 0.0_DP
      call helio_getacch_int_tp(helio_plA, helio_tpA)
   end if
   if (config%j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(aobl(NDIM, config%nplmax), irh(config%nplmax), xht(NDIM, config%ntpmax), &
                     aoblt(NDIM, config%ntpmax), irht(config%ntpmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP / sqrt(r2)
      end do
      call obl_acc(npl, helio_plA, config%j2rp2, config%j4rp4, xh, irh, aobl)
      mu = helio_plA%mass(1)
      do i = 1, ntp
         xht(:, i) = helio_tpA%xh(:,i)
         r2 = dot_product(xht(:, i), xht(:, i))
         irht(i) = 1.0_DP / sqrt(r2)
      end do
      call obl_acc_tp(ntp, xht, config%j2rp2, config%j4rp4, irht, aoblt, mu)
      helio_tpA%ah(:,:) = helio_tpA%ahi(:,:) + aoblt(:, :) - aobl(:, 1)
   else
      helio_tpA%ah(:,:) = helio_tpA%ahi(:,:)
   end if
   if (config%lextra_force) call helio_user_getacch_tp(helio_tpA, t)

   return

   end procedure helio_getacch_tp
end submodule s_helio_getacch_tp