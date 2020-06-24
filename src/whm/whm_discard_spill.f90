submodule(whm) s_whm_discard_spill 
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


   return

   end procedure whm_discard_spill
end submodule s_whm_discard_spill
