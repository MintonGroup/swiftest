submodule (helio) s_helio_drift_tp
contains
   module procedure helio_drift_tp     
   !! author: David A. Minton
   !!
   !! Loop through test particles and call Danby drift routine
   !! New vectorized version included
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_drift_tp.f90
   !! Adapted from Hal Levison's Swift routine drift_tp.f
   use swiftest
   implicit none

   integer(I4B), dimension(:),allocatable :: iflag !! Vectorized error code flag
   integer(I4B) :: i !! Loop counter

   associate(ntp => self%nbody)
      allocate(iflag(ntp))
      call drift_one(self%mu(:), self%xh(:, 1), self%xh(:, 2), self%xh(:, 3),& 
                                     self%vb(:, 1), self%vb(:, 2), self%vb(:, 3),&
                                     iflag(:))
      if (any(iflag(:) /= 0 )) then
         do i = 1,ntp
            if (iflag(i) /= 0) then
               write(*, *) "Particle ", self%name(i), " lost due to error in Danby drift"
            end if
         end do
         deallocate(iflag)
      end if
      deallocate(iflag)
   end associate

   return
   end procedure helio_drift_tp
end submodule s_helio_drift_tp
