submodule(whm_classes) s_whm_util
   use swiftest
contains

   module subroutine whm_util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one massive body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(whm_pl),                   intent(inout) :: self         !! WHM massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (whm_pl)
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%eta, source%eta, nold, nsrc, lsource_mask)
            call util_append(self%muj, source%muj, nold, nsrc, lsource_mask)
            call util_append(self%ir3j, source%ir3j, nold, nsrc, lsource_mask)
            call util_append(self%xj, source%xj, nold, nsrc, lsource_mask)
            call util_append(self%vj, source%vj, nold, nsrc, lsource_mask)

            call util_append_pl(self, source, lsource_mask)
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class whm_pl or its descendents"
         call util_exit(FAILURE)
      end select

      return
   end subroutine whm_util_append_pl

   module subroutine whm_util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new WHM test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(whm_pl),                      intent(inout) :: self       !! WHM massive body object
      class(swiftest_body),               intent(in)    :: inserts    !! inserted object 
      logical, dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B)                                      :: i
   
      associate(keeps => self)
         select type(inserts)
         class is (whm_pl)
            call util_fill(keeps%eta, inserts%eta, lfill_list)
            call util_fill(keeps%muj, inserts%muj, lfill_list)
            call util_fill(keeps%ir3j, inserts%ir3j, lfill_list)
            call util_fill(keeps%xj, inserts%xj, lfill_list)
            call util_fill(keeps%vj, inserts%vj, lfill_list)

            call util_fill_pl(keeps, inserts, lfill_list)
         class default
            write(*,*) "Invalid object passed to the fill method. Inserts must be of class whm_pl or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
   
      return
   end subroutine whm_util_fill_pl


   module subroutine whm_util_index_eucl_plpl(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper for the indexing method for WHM massive bodies. Sorts the massive bodies by heliocentric distance and then flattens the pl-pl upper triangular matrix
      implicit none
      ! Arguments
      class(whm_pl),              intent(inout) :: self  !! WHM massive body object
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters

      call self%sort("ir3h", ascending=.false.)
      call util_index_eucl_plpl(self, param)

      return
   end subroutine whm_util_index_eucl_plpl


   module subroutine whm_util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a massive body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(whm_pl), intent(inout) :: self  !! WHM massive body object
      integer(I4B),  intent(in)    :: nnew  !! New size neded

      call util_resize(self%eta, nnew)
      call util_resize(self%xj, nnew)
      call util_resize(self%vj, nnew)
      call util_resize(self%muj, nnew)
      call util_resize(self%ir3j, nnew)

      call util_resize_pl(self, nnew)

      return
   end subroutine whm_util_resize_pl


   module subroutine whm_util_set_ir3j(self)
      !! author: David A. Minton
      !!
      !! Sets the inverse Jacobi and heliocentric radii cubed (1/rj**3 and 1/rh**3)
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self    !! WHM massive body object
      ! Internals
      integer(I4B)                                 :: i
      real(DP)                                     :: r2, ir

      if (self%nbody > 0) then
         do i = 1, self%nbody
            r2 = dot_product(self%xh(:, i), self%xh(:, i))
            ir = 1.0_DP / sqrt(r2)
            self%ir3h(i) = ir / r2
            r2 = dot_product(self%xj(:, i), self%xj(:, i))
            ir = 1.0_DP / sqrt(r2)
            self%ir3j(i) = ir / r2
         end do
      end if

      return
   end subroutine whm_util_set_ir3j


   module subroutine whm_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a WHM massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(whm_pl), intent(inout) :: self      !! WHM massive body object
      character(*),  intent(in)    :: sortby    !! Sorting attribute
      logical,       intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
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
         case("eta")
            call util_sort(direction * pl%eta(1:npl), ind(1:npl))
         case("muj")
            call util_sort(direction * pl%muj(1:npl), ind(1:npl))
         case("ir3j")
            call util_sort(direction * pl%ir3j(1:npl), ind(1:npl))
         case("xj", "vj")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default
            call util_sort_pl(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)
      end associate

      return
   end subroutine whm_util_sort_pl


   module subroutine whm_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange WHM massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(whm_pl),               intent(inout) :: self !! WHM massive body object
      integer(I4B),  dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange(pl%eta,  ind, npl)
         call util_sort_rearrange(pl%xj,   ind, npl)
         call util_sort_rearrange(pl%vj,   ind, npl)
         call util_sort_rearrange(pl%muj,  ind, npl)
         call util_sort_rearrange(pl%ir3j, ind, npl)

         call util_sort_rearrange_pl(pl,ind)
      end associate

      return
   end subroutine whm_util_sort_rearrange_pl


   module subroutine whm_util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) WHM test particle structure from active list to discard list
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(whm_pl),                         intent(inout) :: self        !! WHM massive body object
      class(swiftest_body),                  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:),                 intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      ! Internals
      integer(I4B)                                         :: i
      associate(keeps => self)
         select type(discards)
         class is (whm_pl)
            call util_spill(keeps%eta, discards%eta, lspill_list, ldestructive)
            call util_spill(keeps%muj, discards%muj, lspill_list, ldestructive)
            call util_spill(keeps%ir3j, discards%ir3j, lspill_list, ldestructive)
            call util_spill(keeps%xj, discards%xj, lspill_list, ldestructive)
            call util_spill(keeps%vj, discards%vj, lspill_list, ldestructive)

            call util_spill_pl(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class whm_pl or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine whm_util_spill_pl
   
end submodule s_whm_util
