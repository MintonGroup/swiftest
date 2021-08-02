submodule(rmvs_classes) s_rmvs_util
   use swiftest
contains

   module subroutine rmvs_util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new RMVS massive body structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(rmvs_pl),        intent(inout) :: self       !! RMVS massive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B) :: i

      associate(keeps => self)
         select type(inserts)
         class is (rmvs_pl)

            call util_fill(keeps%nenc, inserts%nenc, lfill_list)
            call util_fill(keeps%tpenc1P, inserts%tpenc1P, lfill_list)
            call util_fill(keeps%plind, inserts%plind, lfill_list)

            call whm_util_fill_pl(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on rmvs_pl'
         end select
      end associate

      return
   end subroutine rmvs_util_fill_pl


   module subroutine rmvs_util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new RMVS test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(rmvs_tp),        intent(inout) :: self       !! RMVS test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (rmvs_tp)

            call util_fill(keeps%lperi, inserts%lperi, lfill_list)
            call util_fill(keeps%plperP, inserts%plperP, lfill_list)
            call util_fill(keeps%plencP, inserts%plencP, lfill_list)
            
            call util_fill_tp(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on rmvs_tp'
         end select
      end associate

      return
   end subroutine rmvs_util_fill_tp

   module subroutine rmvs_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a RMVS massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(rmvs_pl), intent(inout) :: self       !! RMVS massive body object
      character(*),   intent(in)     :: sortby    !! Sorting attribute
      logical,        intent(in)     :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(self%nbody) :: ind
      integer(I4B) :: direction

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("nenc")
            call util_sort(direction * pl%nenc(1:npl), ind(1:npl))
         case("tpenc1P")
            call util_sort(direction * pl%tpenc1P(1:npl), ind(1:npl))
         case("plind")
            call util_sort(direction * pl%plind(1:npl), ind(1:npl))
         case("outer", "inner", "planetocentric", "lplanetocentric")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call whm_util_sort_pl(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)

      end associate
      return
   end subroutine rmvs_util_sort_pl


   module subroutine rmvs_util_sort_tp(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a RMVS test particle object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(rmvs_tp), intent(inout) :: self      !! RMVS test particle object
      character(*),   intent(in)    :: sortby    !! Sorting attribute
      logical,        intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(self%nbody) :: ind
      integer(I4B)                        :: direction

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("plperP")
            call util_sort(direction * tp%plperP(1:ntp), ind(1:ntp))
         case("plencP")
            call util_sort(direction * tp%plencP(1:ntp), ind(1:ntp))
         case("lperi", "cb_heliocentric", "xheliocentric", "index", "ipleP", "lplanetocentric")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class (*NOTE whm_tp does not need its own sort method, so we go straight to the swiftest_tp method)
            call util_sort_tp(tp, sortby, ascending)
            return
         end select

         call tp%rearrange(ind)

      end associate
      return
   end subroutine rmvs_util_sort_tp

   module subroutine rmvs_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange RMVS massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(rmvs_pl),               intent(inout) :: self !! RMVS massive body object
      integer(I4B),   dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      ! Internals
      class(rmvs_pl), allocatable :: pl_sorted  !! Temporary holder for sorted body
      integer(I4B) :: i

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange_pl(pl,ind)
         allocate(pl_sorted, source=self)
         pl%eta(1:npl) = pl_sorted%eta(ind(1:npl))
         pl%xj(:,1:npl) = pl_sorted%xj(:,ind(1:npl))
         pl%vj(:,1:npl) = pl_sorted%vj(:,ind(1:npl))
         pl%muj(1:npl) = pl_sorted%muj(ind(1:npl))
         pl%ir3j(1:npl) = pl_sorted%ir3j(ind(1:npl))
         deallocate(pl_sorted)
      end associate

      return
   end subroutine rmvs_util_sort_rearrange_pl


   module subroutine rmvs_util_sort_rearrange_tp(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange RMVS test particle object in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(rmvs_tp),                intent(inout) :: self !! RMVS test particle object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      ! Internals
      class(rmvs_tp), allocatable :: tp_sorted  !! Temporary holder for sorted body

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange_tp(tp,ind)
         allocate(tp_sorted, source=self)
         tp%lperi(1:ntp) = tp_sorted%lperi(ind(1:ntp))
         tp%plperP(1:ntp) = tp_sorted%plperP(ind(1:ntp))
         tp%plencP(1:ntp) = tp_sorted%plencP(ind(1:ntp))
         tp%xheliocentric(:,1:ntp) = tp_sorted%xheliocentric(:,ind(1:ntp))
         deallocate(tp_sorted)
      end associate

      return
   end subroutine rmvs_util_sort_rearrange_tp
   

   module subroutine rmvs_util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) RMVS test particle structure from active list to discard list
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine discard_discard_spill.f90
      implicit none
      ! Arguments
      class(rmvs_pl),        intent(inout) :: self         !! RMVS massive body body object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: i

      associate(keeps => self)
         select type(discards)
         class is (rmvs_pl)
            call util_spill(keeps%nenc, discards%nenc, lspill_list, ldestructive)
            call util_spill(keeps%tpenc1P, discards%tpenc1P, lspill_list, ldestructive)
            call util_spill(keeps%plind, discards%plind, lspill_list, ldestructive)
            call whm_util_spill_pl(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on rmvs_pl'
         end select
      end associate

      return
   end subroutine rmvs_util_spill_pl

   
   module subroutine rmvs_util_spill_tp(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) RMVS test particle structure from active list to discard list
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(rmvs_tp),        intent(inout) :: self        !! RMVS test particle object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B) :: i

      associate(keeps => self)
         select type(discards)
         class is (rmvs_tp)
            call util_spill(keeps%lperi, discards%lperi, lspill_list, ldestructive)
            call util_spill(keeps%plperP, discards%plperP, lspill_list, ldestructive)
            call util_spill(keeps%plencP, discards%plencP, lspill_list, ldestructive)

            call util_spill_tp(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) 'Error! spill method called for incompatible return type on rmvs_tp'
         end select
      end associate

      return
   end subroutine rmvs_util_spill_tp

end submodule s_rmvs_util
