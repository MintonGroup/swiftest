submodule(rmvs_classes) s_rmvs_obl
   use swiftest
contains  
   module subroutine rmvs_obl_acc_in(self, cb)
      !! author: David A. Minton
      !!
      !! Compute the oblateness acceleration in the inner encounter region with planets 
      !! 
      implicit none
      ! Arguments
      class(rmvs_pl),            intent(inout) :: self !! RMVS massive body object
      class(swiftest_cb),        intent(inout) :: cb   !! Swiftest central body object
      ! Internals
      integer(I4B)                             :: i
      real(DP), dimension(:, :), allocatable   :: xh_original 

      associate(pl => self)
         allocate(xh_original, source=pl%xh)
         do i = 0, NTPHENC
            pl%xh(:,:) = pl%xin(:,:,i) ! Temporarily replace heliocentric position with inner substep values to calculate the oblateness terms
            call pl%obl_acc(cb)
            pl%aoblin(:,:,i) = pl%aobl(:,:) ! Save the oblateness acceleration on the planet for this substep
         end do
         ! Put back the original heliocentric position for the planets
         pl%xh(:,:) = xh_original(:,:)
         deallocate(xh_original)
      end associate

      return
   end subroutine rmvs_obl_acc_in
end submodule s_rmvs_obl