submodule(rmvs_classes) s_rmvs_spill_and_fill
   use swiftest
contains
   module subroutine rmvs_spill_pl(self, discards, lspill_list)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) RMVS test particle structure from active list to discard list
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine discard_discard_spill.f90
      implicit none
      ! Arguments
      class(rmvs_pl),                        intent(inout) :: self      !! Swiftest massive body body object
      class(swiftest_body),                  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:),                 intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      ! Internals
      integer(I4B)                                         :: i

      associate(keeps => self)
         select type(discards)
         class is (rmvs_pl)

            discards%nenc(:)    = pack(keeps%nenc(:),         lspill_list(:))
            discards%tpenc1P(:) = pack(keeps%tpenc1P(:),       lspill_list(:))
            if (count(.not.lspill_list(:))  > 0) then
               keeps%nenc(:)       = pack(keeps%nenc(:),   .not. lspill_list(:))
               keeps%tpenc1P(:)    = pack(keeps%tpenc1P(:), .not. lspill_list(:))
            end if
            call whm_spill_pl(keeps, discards, lspill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on rmvs_pl'
         end select
      end associate

      return

      end subroutine rmvs_spill_pl

      module subroutine rmvs_fill_pl(self, inserts, lfill_list)
         !! author: David A. Minton
         !!
         !! Insert new RMVS massive body structure into an old one. 
         !! This is the inverse of a fill operation.
         !! 
         implicit none
         ! Arguments
         class(rmvs_pl),        intent(inout) :: self       !! RMVS massive body object
         class(swiftest_body),  intent(inout) :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
         ! Internals
         integer(I4B)                                      :: i
      
         associate(keeps => self)
            select type(inserts)
            class is (rmvs_pl)
      
               keeps%nenc(:)    = unpack(keeps%nenc(:),    .not.lfill_list(:), keeps%nenc(:))
               keeps%nenc(:)    = unpack(inserts%nenc(:),    lfill_list(:), keeps%nenc(:))
               
               keeps%tpenc1P(:) = unpack(keeps%tpenc1P(:), .not.lfill_list(:), keeps%tpenc1P(:))
               keeps%tpenc1P(:) = unpack(inserts%tpenc1P(:), lfill_list(:), keeps%tpenc1P(:))
               
               call whm_fill_pl(keeps, inserts, lfill_list)
            class default
               write(*,*) 'Error! spill method called for incompatible return type on rmvs_pl'
            end select
         end associate
      
         return
      
   end subroutine rmvs_fill_pl

   module subroutine rmvs_spill_tp(self, discards, lspill_list)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) RMVS test particle structure from active list to discard list
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(rmvs_tp),                        intent(inout) :: self       !! RMVS test particle object
      class(swiftest_body),                  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:),                 intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      ! Internals
      integer(I4B)                                         :: i

      associate(keeps => self)
         select type(discards)
         class is (rmvs_tp)
            discards%lperi(:)  = pack(keeps%lperi(:),       lspill_list(:))
            discards%plperP(:) = pack(keeps%plperP(:),       lspill_list(:))
            discards%plencP(:) = pack(keeps%plencP(:),       lspill_list(:))
            discards%tpencP(:) = pack(keeps%tpencP(:),       lspill_list(:))
            if (count(.not.lspill_list(:))  > 0) then
               keeps%lperi(:)     = pack(keeps%lperi(:), .not. lspill_list(:))
               keeps%plperP(:)    = pack(keeps%plperP(:), .not. lspill_list(:))
               keeps%plencP(:)    = pack(keeps%plencP(:), .not. lspill_list(:))
               keeps%tpencP(:)    = pack(keeps%tpencP(:), .not. lspill_list(:))
            end if

            call util_spill_tp(keeps, discards, lspill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on rmvs_tp'
         end select
      end associate

      return

   end subroutine rmvs_spill_tp

   module subroutine rmvs_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new RMVS test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(rmvs_tp),                     intent(inout) :: self        !! RMVS massive body object
      class(swiftest_body),               intent(inout) :: inserts     !!  Inserted object 
      logical, dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
   
      associate(keeps => self)
         select type(inserts)
         class is (rmvs_tp)

            keeps%lperi(:)  = unpack(keeps%lperi(:),  .not.lfill_list(:), keeps%lperi(:))
            keeps%lperi(:)  = unpack(inserts%lperi(:),  lfill_list(:), keeps%lperi(:))
            
            keeps%plperP(:) = unpack(keeps%plperP(:), .not.lfill_list(:), keeps%plperP(:))
            keeps%plperP(:) = unpack(inserts%plperP(:), lfill_list(:), keeps%plperP(:))
            
            keeps%plencP(:) = unpack(keeps%plencP(:), .not.lfill_list(:), keeps%plencP(:))
            keeps%plencP(:) = unpack(inserts%plencP(:), lfill_list(:), keeps%plencP(:))
            
            keeps%tpencP(:) = unpack(keeps%tpencP(:), .not.lfill_list(:), keeps%tpencP(:))
            keeps%tpencP(:) = unpack(inserts%tpencP(:), lfill_list(:), keeps%tpencP(:))
            
            call util_fill_tp(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on rmvs_tp'
         end select
      end associate
   
      return

   end subroutine rmvs_fill_tp

end submodule s_rmvs_spill_and_fill
