submodule (util) s_util_dist_index_pltp
contains
   module procedure util_dist_index_pltp
   !! author: Jacob Elliott and David A. Minton
   !!
   !! Turns i,j indices into k index for use in the Euclidean distance matrix
   use swiftest
   integer(I4B)          :: i,j,ii,jj,nb,np,nt,counter,ii_end,jj_end

   num_comparisons = (nplm -1) * ntp ! number of entries in our distance array

   allocate(k_pltp(2,num_comparisons))

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

   do i = 2,nplm,np
      ii_end = min(i+np-1, nplm)
      do j = 1,ntp,np
         jj_end = min(j+np-1, ntp)
         do ii = i, ii_end
            do jj = j, jj_end
               k_pltp(1,counter) = ii
               k_pltp(2,counter) = jj
               counter = counter + 1
            end do
         end do
      end do
   end do

   return

   end procedure util_dist_index_pltp
end submodule s_util_dist_index_pltp
