submodule (util) s_util_dist_eucl_pltp
contains
   module procedure util_dist_eucl_pltp
   !! author: Jacob R. Elliott
   !!
   !! Efficient parallel loop-blocking algrorithm for evaluating the Euclidean distance matrix for planet-test particles
   use swiftest
   implicit none

   integer(I4B)          :: k
   !$omp parallel do default(none) schedule(static) &
   !$omp shared (num_pltp_comparisons, test_particles, massive bodies, outvar, k_pltp) &
   !$omp private (k)
   do k = 1,num_pltp_comparisons
      outvar(:,k) = test_particles(:,k_pltp(2,k)) - massive bodies(:,k_pltp(1,k))
   end do
   !$omp end parallel do
   return

   end procedure util_dist_eucl_pltp
end submodule s_util_dist_eucl_pltp
