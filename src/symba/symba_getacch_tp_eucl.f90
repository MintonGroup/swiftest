submodule (symba) s_symba_getacch_tp_eucl
contains
   module procedure symba_getacch_tp_eucl
   !! author: Jacob R. Elliott
   !!
   !! Compute heliocentric accelerations of test particles.
   !!      Accelerations in an encounter are not included here
   !! Uses the single-loop blocking to evaluate the Euclidean distance matrix
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_getacch_tp_eucl.f90
   !! Adapted from Hal Levison's Swift routine symba5_getacch.f
use swiftest
implicit none
   logical , save                 :: lmalloc = .true.
   integer(I4B)                     :: i, j, k, index_pl, index_tp
   real(DP)                       :: rji2, irij3, faci, facj, r2, fac, mu
   real(DP), dimension(NDIM)            :: dx
   real(DP), dimension(:), allocatable, save    :: irh, irht
   real(DP), dimension(:, :), allocatable, save :: aobl, xht, aoblt
   real(DP), dimension(NDIM, ntp)         :: ah

! executable code

   ah(:,1:ntp) = 0.0_DP

   ! call util_dist_eucl_pltp(npl, ntp, symba_plA%xh, symba_tpA%xh, &
   !    num_pltp_comparisons, k_pltp, dist_pltp_array)

!$omp parallel do default(none) schedule(static) &
!$omp shared(num_pltp_comparisons, symba_plA, symba_tpA, k_pltp) &
!$omp private(k, i, j, dx, r2, fac) &
!$omp reduction(+:ah)
   do k = 1,num_pltp_comparisons
      j = k_pltp(2,k)
      if (symba_tpA%status(j) == ACTIVE) then
         i = k_pltp(1,k)
         dx(:) = symba_tpA%xh(:,k_pltp(2,k)) - symba_plA%xh(:,k_pltp(1,k))
         r2 = dot_product(dx(:), dx(:))
         fac = symba_plA%mass(i)/(r2*sqrt(r2))
         ah(:,j) = ah(:,j) - fac*dx(:)
      endif
   enddo
!$omp end parallel do

   symba_tpA%ah(:,1:ntp) = ah(:,1:ntp)

   do i = 1, npltpenc
      index_pl = pltpenc_list%indexpl(i)
      index_tp = pltpenc_list%indextp(i)
      if (symba_tpA%status(index_tp) == ACTIVE) then
         dx(:) = symba_tpA%xh(:,index_tp) - symba_plA%xh(:,index_pl)
         r2 = dot_product(dx(:), dx(:))
         fac = symba_plA%mass(index_pl)/(r2*sqrt(r2))
         symba_tpA%ah(:,index_tp) = symba_tpA%ah(:,index_tp) + fac*dx(:)
      end if
   end do
   ! $omp end parallel do
   if (config%j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(aobl(NDIM, config%nplmax), irh(config%nplmax), xht(NDIM, config%ntpmax), aoblt(NDIM, config%ntpmax), irht(config%ntpmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc(symba_plA, config%j2rp2, config%j4rp4, symba_plA%xh(:,:), irh, aobl)
      mu = symba_plA%mass(1)
      do i = 1, ntp
         xht(:, i) = symba_tpA%xh(:,i) !optimize
         r2 = dot_product(xht(:, i), xht(:, i))
         irht(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc_tp(ntp, xht, config%j2rp2, config%j4rp4, irht, aoblt, mu)
      do i = 1, ntp
         if (symba_tpA%status(i) == ACTIVE) &
         symba_tpA%ah(:,i) = symba_tpA%ah(:,i) + aoblt(:, i) - aobl(:, 1)
      end do
   end if
   if (lextra_force) call symba_user_getacch_tp(t, ntp, symba_tpA)

   return

   end procedure symba_getacch_tp_eucl
end submodule s_symba_getacch_tp_eucl
