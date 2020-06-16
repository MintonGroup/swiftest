submodule (helio) s_helio_drift_pl
contains
   module procedure helio_drift_pl     
   !! author: David A. Minton
   !!
   !! Loop through planets and call Danby drift routine
   !! New vectorized version included
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_drift.f90
   !! Adapted from Hal Levison's Swift routine drift.f
   use swiftest
   implicit none

   integer(I4B) :: npl 
   integer(I4B), dimension(:),allocatable :: iflag
   integer(I4B) :: i

   npl = helio_plA%nbody
   allocate(iflag(npl))
   call drift_one(helio_plA%mu_vec(2:npl), helio_plA%xh(1,2:npl),&
                                           helio_plA%xh(2,2:npl),& 
                                           helio_plA%xh(3,2:npl),& 
                                           helio_plA%vb(1,2:npl),& 
                                           helio_plA%vb(2,2:npl),& 
                                           helio_plA%vb(3,2:npl),&
                                           helio_plA%dt_vec(2:npl), iflag(2:npl))
   if (any(iflag(2:npl) /= 0)) then
      do i = 1,npl
         if (iflag(i) /= 0) then
            write(*, *) " Planet ", helio_plA%name(i), " is lost!!!!!!!!!!"
            write(*, *) mu, dt
            write(*, *) helio_plA%xh(:,i)
            write(*, *) helio_plA%vb(:,i)
            write(*, *) " Stopping "
         end if
      end do
      deallocate(iflag)
      call util_exit(failure)
   end if
   deallocate(iflag)
   
   return


   end procedure helio_drift_pl
end submodule s_helio_drift_pl
