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

   module procedure rmvs_obl_acc_tp
      !! author: David A. Minton
      !!
      !! During a close encounter, the encountering planet and Sun are swapped. This subroutine ensures that the coorrect
      !! solar oblateness acceleration term is computed
      !! 
      implicit none
      real(DP), dimension(:,:), allocatable :: xh

      associate(tp => self)
         if (tp%lenc) then ! This is a close encounter body, so it will switch back to the heliocentric frame for the oblateness calculation
            allocate(xh, source=tp%xh)
            tp%xh(:,:) = tp%xheliocen(:,:)
            call obl_acc_body(tp, tp%cb)
            tp%xh(:,:) = xh(:,:)
            deallocate(xh)
         else
            call obl_acc_body(tp, cb)
         end if
      end associate

      return
   end procedure rmvs_obl_acc_tp
end submodule s_rmvs_obl