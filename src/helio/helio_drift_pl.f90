submodule (helio) s_helio_drift_pl
contains
   module procedure helio_drift_pl     
   !! author: David A. Minton
   !!
   !! Loop through massive bodies and call Danby drift routine
   !! New vectorized version included
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_drift.f90
   !! Adapted from Hal Levison's Swift routine drift.f
   use swiftest
   implicit none

   integer(I4B), dimension(:),allocatable :: iflag !! Vectorized error code flag
   integer(I4B) :: i !! Loop counter

   associate(npl => self%nbody)
      allocate(iflag(npl))
      call drift_one(self%mu(2:npl), self%xh(1,2:npl), self%xh(2,2:npl), self%xh(3,2:npl),& 
                                         self%vb(1,2:npl), self%vb(2,2:npl), self%vb(3,2:npl),&
                                         iflag(2:npl))
      if (any(iflag(2:npl) /= 0)) then
         do i = 1,npl
            if (iflag(i) /= 0) then
               write(*, *) " Massive body ", self%name(i), " is lost!!!!!!!!!!"
               write(*, *) mu, dt
               write(*, *) self%xh(:,i)
               write(*, *) self%vb(:,i)
               write(*, *) " Stopping "
            end if
         end do
         deallocate(iflag)
         call util_exit(FAILURE)
      end if
      deallocate(iflag)
   end associate
   
   return
   end procedure helio_drift_pl
end submodule s_helio_drift_pl
