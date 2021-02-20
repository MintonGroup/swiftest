submodule(rmvs_classes) s_rmvs_spill_and_fill
contains
module procedure rmvs_spill_pl
   !! author: David A. Minton
   !!
   !! Move spilled (discarded) RMVS test particle structure from active list to discard list
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine discard_discard_spill.f90
   use swiftest
   implicit none

   integer(I4B) :: i

   associate(keeps => self, npl => self%nbody)
      select type(discards)
      class is (rmvs_pl)

         discards%nenc(:)    = pack(keeps%nenc(1:npl),         lspill_list(1:npl))
         keeps%nenc(:)       = pack(keeps%nenc(1:npl),   .not. lspill_list(1:npl))

         discards%tpenc1P(:) = pack(keeps%tpenc1P(1:npl),       lspill_list(1:npl))
         keeps%tpenc1P(:)    = pack(keeps%tpenc1P(1:npl), .not. lspill_list(1:npl))

         call whm_spill_pl(keeps, discards, lspill_list)
      class default
         write(*,*) 'Error! spill method called for incompatible return type on rmvs_pl'
      end select
   end associate

   return

   end procedure rmvs_spill_pl

   module procedure rmvs_fill_pl
      !! author: David A. Minton
      !!
      !! Insert new RMVS massive body structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      use swiftest
      implicit none
   
      integer(I4B) :: i
   
      associate(keeps => self)
         select type(inserts)
         class is (rmvs_pl)
   
            keeps%nenc(:)    = merge(inserts%nenc(:),    keeps%nenc(:),    lfill_list(:))
            keeps%tpenc1P(:) = merge(inserts%tpenc1P(:), keeps%tpenc1P(:), lfill_list(:))
   
            call whm_fill_pl(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on rmvs_pl'
         end select
      end associate
   
      return
   
   end procedure rmvs_fill_pl

   module procedure rmvs_spill_tp
   !! author: David A. Minton
   !!
   !! Move spilled (discarded) RMVS test particle structure from active list to discard list
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
   use swiftest
   implicit none

   integer(I4B) :: i

   associate(keeps => self, ntp => self%nbody)
      select type(discards)
      class is (rmvs_tp)
         discards%lperi(:)  = pack(keeps%lperi(1:ntp),       lspill_list(1:ntp))
         keeps%lperi(:)     = pack(keeps%lperi(1:ntp), .not. lspill_list(1:ntp))

         discards%plperP(:) = pack(keeps%plperP(1:ntp),       lspill_list(1:ntp))
         keeps%plperP(:)    = pack(keeps%plperP(1:ntp), .not. lspill_list(1:ntp))

         discards%plencP(:) = pack(keeps%plencP(1:ntp),       lspill_list(1:ntp))
         keeps%plencP(:)    = pack(keeps%plencP(1:ntp), .not. lspill_list(1:ntp))

         discards%tpencP(:) = pack(keeps%tpencP(1:ntp),       lspill_list(1:ntp))
         keeps%tpencP(:)    = pack(keeps%tpencP(1:ntp), .not. lspill_list(1:ntp))

         call util_spill_tp(keeps, discards, lspill_list)
      class default
         write(*,*) 'Error! spill method called for incompatible return type on rmvs_tp'
      end select
   end associate

   return

   end procedure rmvs_spill_tp

   module procedure rmvs_fill_tp
      !! author: David A. Minton
      !!
      !! Insert new RMVS test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      use swiftest
      implicit none
   
      associate(keeps => self)
         select type(inserts)
         class is (rmvs_tp)

            keeps%lperi(:)  = merge(inserts%lperi(:),  keeps%lperi(:), lfill_list(:)) 
            keeps%plperP(:) = merge(inserts%plperP(:), keeps%plperP(:), lfill_list(:))
            keeps%plencP(:) = merge(inserts%plencP(:), keeps%plencP(:), lfill_list(:))
            keeps%tpencP(:) = merge(inserts%tpencP(:), keeps%tpencP(:), lfill_list(:))
   
            call util_fill_tp(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on rmvs_tp'
         end select
      end associate
   
      return
   
      end procedure rmvs_fill_tp

end submodule s_rmvs_spill_and_fill
