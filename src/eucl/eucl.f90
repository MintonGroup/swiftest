submodule (swiftest_classes) eucl_implementations
contains

   module procedure eucl_dist_index_plpl
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      use swiftest
      implicit none

      integer(I4B) :: i, j, k, kp, p

      associate(npl => self%nbody, num_comparisons => self%num_comparisons)
         num_comparisons = (npl * (npl - 1) / 2) ! number of entries in a strict lower triangle, nplm x npl, minus first column
         if (allocated(self%k_eucl)) deallocate(self%k_eucl) ! Reset the index array if it's been set previously
         if (allocated(self%irij3)) deallocate(self%irij3)  
         allocate(self%k_eucl(2, num_comparisons))
         allocate(self%irij3(num_comparisons))
        !do concurrent(k = 1:num_comparisons)
         do k = 1, num_comparisons
            kp = npl * (npl - 1) / 2 - k
            p = floor((sqrt(1._DP + 8 * kp) - 1._DP) / 2._DP)
            i = k - (npl - 1) * (npl - 2) / 2 + p * (p + 1) / 2 + 1
            j = npl - 1 - p 
            self%k_eucl(1, k) = i 
            self%k_eucl(2, k) = j
         end do
      end associate
      return

   end procedure eucl_dist_index_plpl

   module procedure eucl_dist_index_pltp
      !! author: Jacob Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix
      use swiftest
      integer(I4B)          :: i, j
  
      if (allocated(self%irij3)) deallocate(self%irij3)  
      allocate(self%irij3(self%nbody, pl%nbody))

      return
    end procedure eucl_dist_index_pltp

   module procedure eucl_irij3_plpl
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Efficient parallel loop-blocking algrorithm for evaluating the Euclidean distance matrix for planet-planet
      use swiftest
      implicit none
      integer(I4B) :: k, i, j
      real(DP), dimension(NDIM) :: dx
      real(DP) :: rji2

      !!$omp num_threads(min(omp_get_max_threads(),ceiling(num_comparisons/10000.))) &

      !!$omp parallel do schedule(static) default(none) &
      !!$omp shared (self) &
      !!$omp private(k, i, j, dx, rji2)
      do k = 1, self%num_comparisons
      !do concurrent (k = 1:self%num_comparisons)
         i = self%k_eucl(1, k)
         j = self%k_eucl(2, k)
         dx(:) = self%xh(:, j) - self%xh(:, i)
         rji2  = dot_product(dx(:), dx(:))
         self%irij3(k) = 1.0_DP / (rji2 * sqrt(rji2))
      end do
      !!$omp end parallel do

      return
   end procedure eucl_irij3_plpl

   module procedure eucl_irij3_pltp
      !! author: David A. Minton
      !!
      !! Calculates the term 1.0_DP / (rji2 * sqrt(rji2)) where rji2 is the square of the Euclidean distance between pl-tp pairs
      use swiftest
      implicit none

      call irij3_pltp(self%nbody, pl%nbody, self%xh, pl%xh, self%irij3)

   end procedure eucl_irij3_pltp

   subroutine irij3_pltp(ntp, npl, xht, xhp, irij3)
      use swiftest
      implicit none
      integer(I4B), intent(in) :: ntp, npl 
      real(DP), dimension(:,:), intent(in) :: xht, xhp
      real(DP), dimension(:,:), intent(out) :: irij3

      integer(I4B) :: i, j
      real(DP), dimension(NDIM) :: dx
      real(DP) :: rji2

      !$omp parallel do schedule(static) default(private) &
      !$omp shared (ntp, npl, irij3, xhp, xht) 
      do j = 1, npl
         !$omp simd
         do i = 1, ntp
            dx(:) = xht(:, i) - xhp(:, j)
            rji2  = dot_product(dx(:), dx(:))
            irij3(i, j) = 1.0_DP / (rji2 * sqrt(rji2))
         end do
      end do
      !$omp end parallel do

      return
   end subroutine irij3_pltp

   !    integer(I4B)          :: k
     ! associate(k_eucl => self%k_eucl, ntp => self%nbody, num_pltp_comparisons => self%num_comparisons, &
     !    test_particles => intp, planets => inpl)
     !           
     !    !$omp parallel do default(none) schedule(static) &
     !    !$omp shared (num_pltp_comparisons, test_particles, planets, outvar, k_eucl) &
     !    !$omp private (k)
     !    do k = 1, num_pltp_comparisons
     !       outvar(:, k) = test_particles(:, k_eucl(2, k)) - planets(:, k_eucl(1, k))
     !    end do
     !    !$omp end parallel do
     ! end associate



         !   nb = 10 ! number of blocks
   !   np = (nplm-1)/nb ! number of planets per block
   !   nt = ntp/nb ! number of test particles per block

   !   np = 1000

   !   counter = 1

   !   do i = 2,nplm,np
   !      ii_end = min(i+np-1, nplm)
   !        do j = 1,ntp,np
   !          jj_end = min(j+np-1, ntp)
   !             do ii = i, ii_end
   !                  do jj = j, jj_end
   !                       k_eucl(1,counter) = ii
   !                       k_eucl(2,counter) = jj
   !                       counter = counter + 1
   !                  enddo
   !             enddo
   !        enddo
!    !   enddo
!       return
   

   ! module procedure eucl_dist_index_plplm
   !    !! author: Jacob R. Elliott
   !    !!
   !    !! Turns i,j indices into k index for use in the Euclidean distance matrix (for SyMBA objects with MTINY)
   !    use swiftest
   !    integer(I4B)          :: i,j,counter

   !    associate(k_eucl => self%k_eucl, npl => self%nbody, num_comparisons => self%num_comparisons)
   !       num_comparisons = ((npl - 1) * (npl - 2) / 2) - ( (npl-nplm-1) * ((npl-nplm-1)+1)/2 )! number of entries in a strict lower triangle, nplm x npl, minus first column
   !       allocate(k_eucl(2,num_comparisons))
   !       ! this is a 'fancier' code, but so far i think it runs slower
   !       ! so leaving it in, but commenting it out
   !       ! i think it's because of the 'mod' call, but i haven't profiled it yet
   !       ! don't forget to uncomment the 'k' declaration up top!
   !       ! allocate(k(num_comparisons))
   
   !       ! m = ceiling(sqrt(2. * num_comparisons))
   
   !       ! k = (/(i, i=1,num_comparisons, 1)/)
   
   !       ! ik_eucl = m - nint( sqrt( dble(2) * (dble(1) + num_comparisons - k))) + 1
   !       ! jk_eucl = mod(k + (ik_eucl - 1) * ik_eucl / 2 - 1, m) + 2
   
   !       ! brute force the index creation
   
   !       !$omp parallel do default(none) schedule(dynamic) &
   !       !$omp shared (k_eucl, npl, nplm) &
   !       !$omp private (i, j, counter)
   !       do i = 2, nplm
   !          counter = (i - 2) * npl - i * (i - 1) / 2 + 2
   !       do j = i + 1, npl
   !             k_eucl(1,counter) = i
   !             k_eucl(2,counter) = j
   !             counter = counter + 1
   !          end do
   !       end do
   !       !$omp end parallel do
   
   !       return
   !    end associate
   ! end procedure eucl_dist_index_plplm



 end submodule eucl_implementations
