submodule (swiftest_classes) eucl_implementations
contains
   module procedure eucl_dist_index_plpl
      !! author: Jacob R. Elliott
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix
      use swiftest
      integer(I4B)          :: i,j,counter

      associate(npl => self%nbody, num_comparisons => self%num_comparisons)
         num_comparisons = ((npl - 1) * (npl - 2) / 2) ! number of entries in a strict lower triangle, nplm x npl, minus first column
         allocate(self%k_plpl(2, num_comparisons))
         ! this is a 'fancier' code, but so far i think it runs slower
         ! so leaving it in, but commenting it out
         ! i think it's because of the 'mod' call, but i haven't profiled it yet
         ! don't forget to uncomment the 'k' declaration up top!
         ! allocate(k(num_comparisons))
   
         ! m = ceiling(sqrt(2. * num_comparisons))
   
         ! k = (/(i, i=1,num_comparisons, 1)/)
   
         ! ik_plpl = m - nint( sqrt( dble(2) * (dble(1) + num_comparisons - k))) + 1
         ! jk_plpl = mod(k + (ik_plpl - 1) * ik_plpl / 2 - 1, m) + 2
   
         ! brute force the index creation
   
         !$omp parallel do default(none) schedule(dynamic) &
         !$omp shared (k_plpl, npl) &
         !$omp private (i, j, counter)
         do i = 1, npl - 1
            do concurrent(j = i + 1:npl)
               counter = (i - 1) * npl - i * (i - 1) / 2 + (j - i)
               self%k_plpl(1, counter) = i
               self%k_plpl(2, counter) = j
            end do
         end do
         !$omp end parallel do
   
         return
      end associate
   end procedure eucl_dist_index_plpl

   module procedure eucl_dist_index_pltp
      !! author: Jacob Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix
      use swiftest
      integer(I4B)          :: i,j,ii,jj,nb,np,nt,counter,ii_end,jj_end
  
      associate(k_pltp => self%k_pltp, ntp => self%nbody, num_comparisons => self%num_comparisons, &
                npl => pl%nbody, nplm => pl%nbody)
         num_comparisons = (nplm - 1) * ntp ! number of entries in our distance array
   
         allocate(self%k_pltp(2, num_comparisons))
   
         ! !$omp parallel do schedule(static) default(none) &
         ! !$omp shared(k_pltp, nplm, ntp) &
         ! !$omp private(i, j, counter)
         !    do i = 2,nplm
         !       counter = (i-2) * ntp + 1
         !       do j = 1,ntp
         !          k_pltp(1,counter) = i
         !          k_pltp(2,counter) = j
         !          counter = counter + 1
         !       enddo
         !    enddo
         ! !$omp end parallel do
      
         np = 500
      
         counter = 1
      
         do i = 2, nplm, np
            ii_end = min(i + np - 1, nplm)
            do j = 1, ntp, np
               jj_end = min(j + np - 1, ntp)
               do ii = i, ii_end
                  do jj = j, jj_end
                     k_pltp(1,counter) = ii
                     k_pltp(2,counter) = jj
                     counter = counter + 1
                  end do
               end do
            end do
         end do
 
      end associate
      return
   
    end procedure eucl_dist_index_pltp

   ! module procedure eucl_dist_index_plplm
   !    !! author: Jacob R. Elliott
   !    !!
   !    !! Turns i,j indices into k index for use in the Euclidean distance matrix (for SyMBA objects with MTINY)
   !    use swiftest
   !    integer(I4B)          :: i,j,counter

   !    associate(k_plpl => self%k_plpl, npl => self%nbody, num_comparisons => self%num_comparisons)
   !       num_comparisons = ((npl - 1) * (npl - 2) / 2) - ( (npl-nplm-1) * ((npl-nplm-1)+1)/2 )! number of entries in a strict lower triangle, nplm x npl, minus first column
   !       allocate(k_plpl(2,num_comparisons))
   !       ! this is a 'fancier' code, but so far i think it runs slower
   !       ! so leaving it in, but commenting it out
   !       ! i think it's because of the 'mod' call, but i haven't profiled it yet
   !       ! don't forget to uncomment the 'k' declaration up top!
   !       ! allocate(k(num_comparisons))
   
   !       ! m = ceiling(sqrt(2. * num_comparisons))
   
   !       ! k = (/(i, i=1,num_comparisons, 1)/)
   
   !       ! ik_plpl = m - nint( sqrt( dble(2) * (dble(1) + num_comparisons - k))) + 1
   !       ! jk_plpl = mod(k + (ik_plpl - 1) * ik_plpl / 2 - 1, m) + 2
   
   !       ! brute force the index creation
   
   !       !$omp parallel do default(none) schedule(dynamic) &
   !       !$omp shared (k_plpl, npl, nplm) &
   !       !$omp private (i, j, counter)
   !       do i = 2, nplm
   !          counter = (i - 2) * npl - i * (i - 1) / 2 + 2
   !       do j = i + 1, npl
   !             k_plpl(1,counter) = i
   !             k_plpl(2,counter) = j
   !             counter = counter + 1
   !          end do
   !       end do
   !       !$omp end parallel do
   
   !       return
   !    end associate
   ! end procedure eucl_dist_index_plplm


   module procedure eucl_dist_plpl
      !! author: Jacob R. Elliott
      !!
      !! Efficient parallel loop-blocking algrorithm for evaluating the Euclidean distance matrix for planet-planet
      use swiftest
      implicit none
      integer(I4B) :: k

      associate(k_plpl => self%k_plpl, num_comparisons => self%num_comparisons)
      
         !$omp parallel do schedule(static) default(none) &
         !$omp num_threads(min(omp_get_max_threads(),ceiling(num_comparisons/10000.))) &
         !$omp shared (outvar, invar, num_comparisons, k_plpl) &
         !$omp private(k)
         do k = 1,num_comparisons
            outvar(:, k) = invar(:, k_plpl(2, k)) - invar(:, k_plpl(1, k))
         end do
         !$omp end parallel do
      end associate

      return
   
   end procedure eucl_dist_plpl

   module procedure eucl_dist_pltp
      !! author: Jacob R. Elliott
      !!
      !! Efficient parallel loop-blocking algrorithm for evaluating the Euclidean distance matrix for planet-test particles
      use swiftest
      implicit none
   
      integer(I4B)          :: k
      associate(k_pltp => self%k_pltp, ntp => self%nbody, num_pltp_comparisons => self%num_comparisons, &
         test_particles => intp, planets => inpl)
                
         !$omp parallel do default(none) schedule(static) &
         !$omp shared (num_pltp_comparisons, test_particles, planets, outvar, k_pltp) &
         !$omp private (k)
         do k = 1, num_pltp_comparisons
            outvar(:, k) = test_particles(:, k_pltp(2, k)) - planets(:, k_pltp(1, k))
         end do
         !$omp end parallel do
      end associate
      return
   
   end procedure eucl_dist_pltp


end submodule eucl_implementations
