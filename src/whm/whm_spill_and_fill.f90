submodule(whm_classes) s_whm_spill_and_fill
contains
   module procedure whm_spill_pl
   !! author: David A. Minton
   !!
   !! Move spilled (discarded) WHM test particle structure from active list to discard list
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
   use swiftest
   implicit none

   integer(I4B) :: i

   associate(keeps => self, npl => self%nbody)
      select type(discards)
      class is (whm_pl)
         discards%eta(:) = pack(keeps%eta(1:npl),       lspill_list(1:npl))
         keeps%eta(:)    = pack(keeps%eta(1:npl), .not. lspill_list(1:npl))

         discards%muj(:) = pack(keeps%muj(1:npl),       lspill_list(1:npl))
         keeps%muj(:)    = pack(keeps%muj(1:npl), .not. lspill_list(1:npl))

         discards%ir3j(:) = pack(keeps%ir3j(1:npl),       lspill_list(1:npl))
         keeps%ir3j(:)    = pack(keeps%ir3j(1:npl), .not. lspill_list(1:npl))

         do i = 1, NDIM
            discards%xj(i, :) = pack(keeps%xj(i, 1:npl),       lspill_list(1:npl))
            keeps%xj(i, :)    = pack(keeps%xj(i, 1:npl), .not. lspill_list(1:npl))

            discards%vj(i, :) = pack(keeps%vj(i, 1:npl),       lspill_list(1:npl))
            keeps%vj(i, :)    = pack(keeps%vj(i, 1:npl), .not. lspill_list(1:npl))
         end do
         call util_spill_pl(keeps, discards, lspill_list)
      class default
         write(*,*) 'Error! spill method called for incompatible return type on whm_pl'
      end select
   end associate

   return

   end procedure whm_spill_pl

   module procedure whm_fill_pl
      !! author: David A. Minton
      !!
      !! Insert new WHM test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      use swiftest
      implicit none
   
      integer(I4B) :: i
   
      associate(keeps => self)
         select type(inserts)
         class is (whm_pl)
            keeps%eta(:)  = merge(inserts%eta(:),  keeps%eta(:),  lfill_list(:))
            keeps%muj(:)  = merge(inserts%muj(:),  keeps%muj(:),  lfill_list(:))
            keeps%ir3j(:) = merge(inserts%ir3j(:), keeps%ir3j(:), lfill_list(:))
   
            do i = 1, NDIM
               keeps%xj(i, :) = merge(inserts%xj(i, :), keeps%xj(i, :), lfill_list(:))
               keeps%vj(i, :) = merge(inserts%vj(i, :), keeps%vj(i, :), lfill_list(:))
            end do
            call util_fill_pl(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on whm_pl'
         end select
      end associate
   
      return
   
      end procedure whm_fill_pl

end submodule s_whm_spill_and_fill
