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

   integer(I4B) :: npl 
   integer(I4B), dimension(:),allocatable :: iflag
   integer(I4B) :: i

   npl = self%nbody
   allocate(iflag(npl))
   call drift_one(self%mu_vec(2:npl), self%xh(1,2:npl), self%xh(2,2:npl), self%xh(3,2:npl),& 
                                      self%vb(1,2:npl), self%vb(2,2:npl), self%vb(3,2:npl),&
                                      self%dt_vec(2:npl), iflag(2:npl))
   if (any(iflag(2:npl) /= 0 )) then
      do i = 1,npl
         if (iflag(i) /= 0) then
            write(*, *) "Particle ", self%name(i), " lost due to error in Danby drift"
         end if
      end do
      deallocate(iflag)
   end if
   deallocate(iflag)
   

   return

   end procedure helio_drift_tp
end submodule s_helio_drift_tp