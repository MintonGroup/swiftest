submodule(rmvs_classes) s_rmvs_obl
contains  

   module procedure rmvs_obl_acc_in
      !! author: David A. Minton
      !!
      !! Compute the oblateness acceleration in the inner encounter region with planets 
      !! 

      use swiftest
      implicit none
      integer(I4B) :: i
      real(DP), dimension(:, :), allocatable       :: xh_original 

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
   end procedure rmvs_obl_acc_in
end submodule s_rmvs_obl