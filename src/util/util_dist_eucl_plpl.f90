submodule (util) s_util_dist_eucl_plpl
contains
   module procedure util_dist_eucl_plpl
   !! author: Jacob R. Elliott
   !!
   !! Efficient parallel loop-blocking algrorithm for evaluating the Euclidean distance matrix for planet-planet
   use swiftest
   implicit none
   integer(I4B) :: k
   
   !$omp parallel do schedule(static) default(none) &
   !$omp num_threads(min(omp_get_max_threads(),ceiling(num_comparisons/10000.))) &
   !$omp shared (outvar, invar, num_comparisons, k_plpl) &
   !$omp private(k)
   do k = 1,num_comparisons
      outvar(:,k) = invar(:,k_plpl(2,k)) - invar(:,k_plpl(1,k))
   end do
   !$omp end parallel do

   return

   end procedure util_dist_eucl_plpl
end submodule s_util_dist_eucl_plpl
