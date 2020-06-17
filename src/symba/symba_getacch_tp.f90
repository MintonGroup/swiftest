submodule (symba) s_symba_getacch_tp
contains
   module procedure symba_getacch_tp
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of test particles
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_getacch_tp.f90
   !! Adapted from Hal Levison's Swift routine symba5_getacch.f
   use swiftest
   implicit none
   logical(lgt), save                 :: lmalloc = .true.
   integer(I4B)                     :: i, j, index_pl, index_tp
   real(DP)                       :: rji2, irij3, faci, facj, r2, fac, mu
   real(DP), dimension(ndim)            :: dx
   real(DP), dimension(:), allocatable, save    :: irh, irht
   real(DP), dimension(:, :), allocatable, save :: aobl, xht, aoblt

! executable code
   !removed by d. minton
   !helio_tpp => symba_tp1p
   !^^^^^^^^^^^^^^^^^^^^
   ! openmp parallelization added by d. minton
   do i = 1, ntp
      !added by d. minton
      !helio_tpp => symba_tp1p%symba_tppa(i)%thisp
      !^^^^^^^^^^^^^^^^^^
      symba_tpa%ah(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      if (symba_tpa%status(i) == active) then
         !swifter_plp => swifter_pl1p
         !do j = 2, nplm
         do j = 2, nplm
            !swifter_plp => swifter_plp%nextp
            dx(:) = symba_tpa%xh(:,i) - xh(:, j)
            r2 = dot_product(dx(:), dx(:))
            fac = symba_pla%mass(j)/(r2*sqrt(r2))
            symba_tpa%ah(:,i) = symba_tpa%ah(:,i) - fac*dx(:)
         end do
      end if
      !removed by d. minton
      !helio_tpp => helio_tpp%nextp
      !^^^^^^^^^^^^^^^^^^^^
   end do
   do i = 1, npltpenc
      !swifter_plp => pltpenc_list(i)%plp%swifter
      !helio_tpp => pltpenc_list(i)%tpp
      index_pl = pltpenc_list%indexpl(i)
      index_tp = pltpenc_list%indextp(i)
      if (symba_tpa%status(index_tp) == active) then
         dx(:) = symba_tpa%xh(:,index_tp) - symba_pla%xh(:,index_pl)
         r2 = dot_product(dx(:), dx(:))
         fac = symba_pla%mass(index_pl)/(r2*sqrt(r2))
         symba_tpa%ah(:,index_tp) = symba_tpa%ah(:,index_tp) + fac*dx(:)
      end if
   end do
   if (j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(aobl(ndim, nplmax), irh(nplmax), xht(ndim, ntpmax), aoblt(ndim, ntpmax), irht(ntpmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc(symba_pla, j2rp2, j4rp4, symba_pla%xh(:,:), irh, aobl)
      mu = symba_pla%mass(1)
      do i = 1, ntp
         xht(:, i) = symba_tpa%xh(:,i) !optimize
         r2 = dot_product(xht(:, i), xht(:, i))
         irht(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
      do i = 1, ntp
         if (symba_tpa%status(i) == active) &
         symba_tpa%ah(:,i) = symba_tpa%ah(:,i) + aoblt(:, i) - aobl(:, 1)
      end do
   end if
   if (lextra_force) call symba_user_getacch_tp(t, ntp, symba_tpa)

   return

   end procedure symba_getacch_tp
end submodule s_symba_getacch_tp
