submodule (symba) s_symba_getacch_eucl
contains
   module procedure symba_getacch_eucl
   !! author: Jacob R. Elliott
   !!
   !! Same as symba_getacch but now uses the single-loop blocking to evaluate the Euclidean distance matrix
   !!      Accelerations in an encounter are not included here
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_getacch.f90
   !! Adapted from Hal Levison's Swift routine symba5_getacch.f
use swiftest
implicit none
   logical(lgt), save                 :: lmalloc = .true.
   integer(I4B)                     :: i, j, index_i, index_j, k, counter
   real(DP)                       :: rji2, irij3, faci, facj, r2, fac
   real(DP), dimension(ndim)            :: dx
   real(DP), dimension(ndim, npl)         :: ah
   real(DP), dimension(:), allocatable, save    :: irh
   real(DP), dimension(:, :), allocatable, save :: xh, aobl
   ! real(DP), allocatable, dimension(:,:) :: dist_plpl_array


!executable code
 
   symba_pla%helio%ah(:,2:npl) = 0.0_DP
   ah(:,2:npl) = 0.0_DP
   
   ! call util_dist_eucl_plpl(npl,symba_pla%helio%swiftest%xh, num_plpl_comparisons, k_plpl, dist_plpl_array) ! does not care about mtiny

! there is floating point arithmetic round off error in this loop
! for now, we will keep it in the serial operation, so we can easily compare
! it to the older swifter versions

!$omp parallel do default(none) schedule(static) &
!$omp num_threads(min(omp_get_max_threads(),ceiling(num_plpl_comparisons/10000.))) &
!$omp shared (num_plpl_comparisons, k_plpl, symba_pla) &
!$omp private (i, j, k, dx, rji2, irij3, faci, facj) &
!$omp reduction(+:ah)
   do k = 1, num_plpl_comparisons
      i = k_plpl(1,k)
      j = k_plpl(2,k)
      
      if ((.not. symba_pla%lmerged(i) .or. (.not. symba_pla%lmerged(j)) .or. &
         (symba_pla%index_parent(i) /= symba_pla%index_parent(j)))) then
         
         dx(:) = symba_pla%helio%swiftest%xh(:,k_plpl(2,k)) - symba_pla%helio%swiftest%xh(:,k_plpl(1,k))
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP/(rji2*sqrt(rji2))
         faci = symba_pla%helio%swiftest%mass(i)*irij3
         facj = symba_pla%helio%swiftest%mass(j)*irij3
         ah(:,i) = ah(:,i) + facj*dx(:)
         ah(:,j) = ah(:,j) - faci*dx(:)

      endif
   end do
!$omp end parallel do

   symba_pla%helio%ah(:,2:npl) = ah(:,2:npl)

   do i = 1, nplplenc
      index_i = plplenc_list%index1(i)
      index_j = plplenc_list%index2(i)
      if ((.not. symba_pla%lmerged(index_i)) .or. (.not. symba_pla%lmerged(index_j))  &
          .or. (symba_pla%index_parent(index_i) /= symba_pla%index_parent(index_j))) then !need to update parent/children
         dx(:) = symba_pla%helio%swiftest%xh(:,index_j) - symba_pla%helio%swiftest%xh(:,index_i)
         rji2 = dot_product(dx(:), dx(:))
         irij3 = 1.0_DP/(rji2*sqrt(rji2))
         faci = symba_pla%helio%swiftest%mass(index_i)*irij3
         facj = symba_pla%helio%swiftest%mass(index_j)*irij3
         symba_pla%helio%ah(:,index_i) = symba_pla%helio%ah(:,index_i) - facj*dx(:)
         symba_pla%helio%ah(:,index_j) = symba_pla%helio%ah(:,index_j) + faci*dx(:)
      end if
   end do

   if (j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(xh(ndim, npl), aobl(ndim, npl), irh(npl))
         lmalloc = .false.
      end if
      do i = 2, npl
         
         r2 = dot_product(symba_pla%helio%swiftest%xh(:,i), symba_pla%helio%swiftest%xh(:,i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_acc(symba_pla%helio%swiftest, j2rp2, j4rp4, symba_pla%helio%swiftest%xh(:,:), irh, aobl)
      do i = 2, npl
         symba_pla%helio%ah(:,i) = symba_pla%helio%ah(:,i) + aobl(:, i) - aobl(:, 1)
      end do
   end if

   if (lextra_force) call symba_user_getacch(t, npl, symba_pla)

   return

   end procedure symba_getacch_eucl
end submodule s_symba_getacch_eucl
