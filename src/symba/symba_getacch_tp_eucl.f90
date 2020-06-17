submodule (symba) s_symba_getacch_tp_eucl
contains
   module procedure symba_getacch_tp_eucl
   !! author: Jacob R. Elliott
   !!
   !! Compute heliocentric accelerations of test particles.
   !!      Accelerations in an encounter are not included here
   !! Uses the single-loop blocking to evaluate the Euclidean distance matrix
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_getacch_tp_eucl.f90
   !! Adapted from Hal Levison's Swift routine symba5_getacch.f
use swiftest
implicit none
   logical(lgt), save                 :: lmalloc = .true.
   integer(I4B)                     :: i, j, k, index_pl, index_tp
   real(DP)                       :: rji2, irij3, faci, facj, r2, fac, mu
   real(DP), dimension(ndim)            :: dx
   real(DP), dimension(:), allocatable, save    :: irh, irht
   real(DP), dimension(:, :), allocatable, save :: aobl, xht, aoblt
   real(DP), dimension(ndim, ntp)         :: ah

! executable code

   ah(:,1:ntp) = 0.0_DP

   ! call util_dist_eucl_pltp(npl, ntp, symba_pla%helio%swiftest%xh, symba_tpa%helio%swiftest%xh, &
   !    num_pltp_comparisons, k_pltp, dist_pltp_array)

!$omp parallel do default(none) schedule(static) &
!$omp shared(num_pltp_comparisons, symba_pla, symba_tpa, k_pltp) &
!$omp private(k, i, j, dx, r2, fac) &
!$omp reduction(+:ah)
   do k = 1,num_pltp_comparisons
      j = k_pltp(2,k)
      if (symba_tpa%helio%swiftest%status(j) == active) then
         i = k_pltp(1,k)
         dx(:) = symba_tpa%helio%swiftest%xh(:,k_pltp(2,k)) - symba_pla%helio%swiftest%xh(:,k_pltp(1,k))
         r2 = dot_product(dx(:), dx(:))
         fac = symba_pla%helio%swiftest%mass(i)/(r2*sqrt(r2))
         ah(:,j) = ah(:,j) - fac*dx(:)
      endif
   enddo
!$omp end parallel do

   symba_tpa%helio%ah(:,1:ntp) = ah(:,1:ntp)

   !removed by d. minton
   !helio_tpp => symba_tp1p%helio
   !^^^^^^^^^^^^^^^^^^^^
   ! openmp parallelization added by d. minton
   ! $omp parallel do schedule(static) default(none) &
   ! $omp private(i,helio_tpp,swifter_tpp,swifter_plp,dx,r2,fac) &
   ! $omp shared(ntp,npl,symba_tp1p,swifter_pl1p,xh) 
   ! do i = 1, ntp
   !    !added by d. minton
   !    !helio_tpp => symba_tp1p%symba_tppa(i)%thisp%helio
   !    !^^^^^^^^^^^^^^^^^^
   !    symba_tpa%helio%ah(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   !    if (symba_tpa%helio%swiftest%status(i) == active) then
   !       !swifter_plp => swifter_pl1p
   !       !do j = 2, nplm
   !       do j = 2, nplm
   !          !swifter_plp => swifter_plp%nextp
   !          dx(:) = symba_tpa%helio%swiftest%xh(:,i) - xh(:, j)
   !          r2 = dot_product(dx(:), dx(:))
   !          fac = symba_pla%helio%swiftest%mass(j)/(r2*sqrt(r2))
   !          symba_tpa%helio%ah(:,i) = symba_tpa%helio%ah(:,i) - fac*dx(:)
   !       end do
   !    end if
   !    !removed by d. minton
   !    !helio_tpp => helio_tpp%nextp
   !    !^^^^^^^^^^^^^^^^^^^^
   ! end do
   ! $omp end parallel do
   ! openmp parallelization added by d. minton
   ! $omp parallel do schedule (static) default(none) &
   ! $omp private(i,swifter_plp,helio_tpp,dx,r2,fac) &
   ! $omp shared(pltpenc_list,npltpenc)
   do i = 1, npltpenc
      !swifter_plp => pltpenc_list(i)%plp%helio%swifter
      !helio_tpp => pltpenc_list(i)%tpp%helio
      index_pl = pltpenc_list%indexpl(i)
      index_tp = pltpenc_list%indextp(i)
      if (symba_tpa%helio%swiftest%status(index_tp) == active) then
         dx(:) = symba_tpa%helio%swiftest%xh(:,index_tp) - symba_pla%helio%swiftest%xh(:,index_pl)
         r2 = dot_product(dx(:), dx(:))
         fac = symba_pla%helio%swiftest%mass(index_pl)/(r2*sqrt(r2))
         symba_tpa%helio%ah(:,index_tp) = symba_tpa%helio%ah(:,index_tp) + fac*dx(:)
      end if
   end do
   ! $omp end parallel do
   if (j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(aobl(ndim, nplmax), irh(nplmax), xht(ndim, ntpmax), aoblt(ndim, ntpmax), irht(ntpmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc(symba_pla%helio%swiftest, j2rp2, j4rp4, symba_pla%helio%swiftest%xh(:,:), irh, aobl)
      mu = symba_pla%helio%swiftest%mass(1)
      do i = 1, ntp
         xht(:, i) = symba_tpa%helio%swiftest%xh(:,i) !optimize
         r2 = dot_product(xht(:, i), xht(:, i))
         irht(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
      do i = 1, ntp
         if (symba_tpa%helio%swiftest%status(i) == active) &
         symba_tpa%helio%ah(:,i) = symba_tpa%helio%ah(:,i) + aoblt(:, i) - aobl(:, 1)
      end do
   end if
   if (lextra_force) call symba_user_getacch_tp(t, ntp, symba_tpa)

   return

   end procedure symba_getacch_tp_eucl
end submodule s_symba_getacch_tp_eucl
