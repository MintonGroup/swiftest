submodule(whm_classes) s_whm_discard_spill 
contains
   module procedure whm_discard_spill
   !! author: David A. Minton
   !!
   !! Move spilled (discarded) WHM test particle structure from active list to discard list
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
   use swiftest
   implicit none
   integer(I4B)          :: i

   associate(ntp => keeps%nbody)
      do i = 1, NDIM
         discards%ah(:, i) = pack(keeps%ah(1:ntp, i),       lspill_list(1:ntp))
         keeps%ah(:, i)    = pack(keeps%ah(1:ntp, i), .not. lspill_list(1:ntp))
      end do
   end associate

   return

   end procedure whm_discard_spill
end submodule s_whm_discard_spill
