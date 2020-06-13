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
   integer(I4B), dimension(:),allocatable :: iflag_vec
   integer(I4B) :: i, iflag

   npl = helio_tpA%nbody
   if (lvectorize) then
      allocate(iflag_vec(npl))
      call drift_one_vec(helio_tpA%mu_vec(2:npl), helio_tpA%xh(1,2:npl),&
                                                  helio_tpA%xh(2,2:npl),& 
                                                  helio_tpA%xh(3,2:npl),& 
                                                  helio_tpA%vb(1,2:npl),& 
                                                  helio_tpA%vb(2,2:npl),& 
                                                  helio_tpA%vb(3,2:npl),&
                                                  helio_tpA%dt_vec(2:npl), iflag_vec(2:npl))
      if (any(iflag_vec(2:npl) /= 0)) then
         do i = 1,npl
            if (iflag_vec(i) /= 0) then
               write(*, *) "Particle ", helio_tpA%name(i), " lost due to error in Danby drift"
            end if
         end do
         deallocate(iflag_vec)
      end if
      deallocate(iflag_vec)
   else
      do i = 2, npl
         call drift_one(mu, helio_tpA%xh(:,i), helio_tpA%vb(:,i), dt, iflag)
         if (iflag /= 0) then
               write(*, *) "Particle ", helio_tpA%name(i), " lost due to error in Danby drift"
         end if
      end do
   end if
   

   return

   end procedure helio_drift_tp
end submodule s_helio_drift_tp