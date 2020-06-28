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
      class is (whm_tp)
         associate(npl => self%nbody)
            do concurrent (i = 1:NDIM)
               discards%xj(:, i) = pack(keeps%xj(1:npl, i),       lspill_list(1:npl))
               keeps%xj(:, i)    = pack(keeps%xj(1:npl, i), .not. lspill_list(1:npl))
   
               discards%vj(:, i) = pack(keeps%vj(1:npl, i),       lspill_list(1:npl))
               keeps%vj(:, i)    = pack(keeps%vj(1:npl, i), .not. lspill_list(1:npl))
               
               discards%ah1(:, i) = pack(keeps%ah1(1:npl, i),       lspill_list(1:npl))
               keeps%ah1(:, i)    = pack(keeps%ah1(1:npl, i), .not. lspill_list(1:npl))
   
               discards%ah2(:, i) = pack(keeps%ah2(1:npl, i),       lspill_list(1:npl))
               keeps%ah2(:, i)    = pack(keeps%ah2(1:npl, i), .not. lspill_list(1:npl))
   
               discards%ah3(:, i) = pack(keeps%ah3(1:npl, i),       lspill_list(1:npl))
               keeps%ah3(:, i)    = pack(keeps%ah3(1:npl, i), .not. lspill_list(1:npl))
            end do
         end associate
      end select
   !class is (whm_tp)
   end select
   call discard_spill_body(keeps, discards, lspill_list)
   return

   end procedure whm_discard_spill
end submodule s_whm_discard_spill
