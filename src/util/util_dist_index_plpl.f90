submodule (util) s_util_dist_index_plpl
contains
   module procedure util_dist_index_plpl
   !! author: Jacob R. Elliott
   !!
   !! Turns i,j indices into k index for use in the Euclidean distance matrix
   use swiftest
   integer(I4B)          :: i,j,counter

   num_comparisons = ((npl - 1) * (npl - 2) / 2) - ( (npl-nplm-1) * ((npl-nplm-1)+1)/2 )! number of entries in a strict lower triangle, nplm x npl, minus first column
   allocate(k_plpl(2,num_comparisons))
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
   !$omp shared (k_plpl, npl, nplm) &
   !$omp private (i, j, counter)
   do i = 2,nplm
      counter = (i - 2) * npl - i*(i-1)/2 + 2
      do j = i+1,npl
         k_plpl(1,counter) = i
         k_plpl(2,counter) = j
         counter = counter + 1
      end do
   end do
  !$omp end parallel do

   return

   end procedure util_dist_index_plpl
end submodule s_util_dist_index_plpl
