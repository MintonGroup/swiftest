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
   integer(I4B) :: i

   select type(keeps)
   class is (whm_pl)
      select type(discards)
      class is (whm_pl)
         associate(npl => keeps%nbody)
            discards%eta(:) = pack(keeps%eta(1:npl),       lspill_list(1:npl))
            keeps%eta(:)    = pack(keeps%eta(1:npl), .not. lspill_list(1:npl))

            discards%muj(:) = pack(keeps%muj(1:npl),       lspill_list(1:npl))
            keeps%muj(:)    = pack(keeps%muj(1:npl), .not. lspill_list(1:npl))

            do i = 1, NDIM
               discards%xj(i, :) = pack(keeps%xj(i, 1:npl),       lspill_list(1:npl))
               keeps%xj(i, :)    = pack(keeps%xj(i, 1:npl), .not. lspill_list(1:npl))
   
               discards%vj(i, :) = pack(keeps%vj(i, 1:npl),       lspill_list(1:npl))
               keeps%vj(i, :)    = pack(keeps%vj(i, 1:npl), .not. lspill_list(1:npl))
            end do
         end associate
      end select
   !class is (whm_tp)
   end select
   call discard_spill_body(keeps, discards, lspill_list)
   return

   end procedure whm_discard_spill
end submodule s_whm_discard_spill
