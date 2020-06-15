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
      if (config%vectorize) then
         helio_tpA%ahi(:,:) = 0.0_DP
      else
         do i = 1, ntp
            helio_tpA%ahi(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
         end do
      end if
      call helio_getacch_int_tp(npl, ntp, helio_plA%swiftest, helio_tpA)
   end if
   if (j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(aobl(ndim, nplmax), irh(nplmax), xht(ndim, ntpmax), aoblt(ndim, ntpmax), irht(ntpmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP / sqrt(r2)
      end do
      call obl_acc(npl, helio_plA%swiftest, j2rp2, j4rp4, xh, irh, aobl)
      mu = helio_plA%swiftest%mass(1)
      do i = 1, ntp
         xht(:, i) = helio_tpA%swiftest%xh(:,i)
         r2 = dot_product(xht(:, i), xht(:, i))
         irht(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
      if (config%lvectorize) then
         helio_tpA%ah(:,:) = helio_tpA%ahi(:,:) + aoblt(:, :) - aobl(:, 1)
      else
         do i = 1, ntp
            helio_tpA%ah(:,i) = helio_tpA%ahi(:,i) + aoblt(:, i) - aobl(:, 1)
         end do
      end if
   else
      if (config%vectorize) then
         helio_tpA%ah(:,:) = helio_tpA%ahi(:,:)
      else
         do i = 1, ntp
            helio_tpA%ah(:,i) = helio_tpA%ahi(:,i)
         end do
      end if
   end if
   if (lextra_force) call helio_user_getacch_tp(t, ntp, helio_tpA)

   return

   end procedure helio_getacch_tp
end submodule s_helio_getacch_tp