submodule(whm_classes) s_whm_util
   use swiftest
contains
   module subroutine whm_util_spill_pl(self, discards, lspill_list)
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
   ! Internals
   integer(I4B)                                         :: i
   associate(keeps => self)
      select type(discards)
      class is (whm_pl)
         discards%eta(:) = pack(keeps%eta(:),       lspill_list(:))
         discards%muj(:) = pack(keeps%muj(:),       lspill_list(:))
         discards%ir3j(:) = pack(keeps%ir3j(:),       lspill_list(:))
         do i = 1, NDIM
            discards%xj(i, :) = pack(keeps%xj(i, :),       lspill_list(:))
            discards%vj(i, :) = pack(keeps%vj(i, :),       lspill_list(:))
         end do

         if (count(.not.lspill_list(:))  > 0) then 
            keeps%eta(:)    = pack(keeps%eta(:), .not. lspill_list(:))
            keeps%muj(:)    = pack(keeps%muj(:), .not. lspill_list(:))
            keeps%ir3j(:)    = pack(keeps%ir3j(:), .not. lspill_list(:))
            do i = 1, NDIM
               keeps%xj(i, :)    = pack(keeps%xj(i, :), .not. lspill_list(:))
               keeps%vj(i, :)    = pack(keeps%vj(i, :), .not. lspill_list(:))
            end do
         end if
         call util_spill_pl(keeps, discards, lspill_list)
      class default
         write(*,*) 'Error! spill method called for incompatible return type on whm_pl'
      end select
   end associate

   return

   end subroutine whm_util_spill_pl

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
      class(swiftest_body),               intent(inout) :: inserts    !! inserted object 
      logical, dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      ! Internals
      integer(I4B)                                      :: i
   
      associate(keeps => self)
         select type(inserts)
         class is (whm_pl)
            keeps%eta(:)  = unpack(keeps%eta(:),  .not.lfill_list(:), keeps%eta(:))
            keeps%eta(:)  = unpack(inserts%eta(:),  lfill_list(:), keeps%eta(:))
            
            keeps%muj(:)  = unpack(keeps%muj(:),  .not.lfill_list(:), keeps%muj(:))
            keeps%muj(:)  = unpack(inserts%muj(:),  lfill_list(:), keeps%muj(:))
            
            keeps%ir3j(:) = unpack(keeps%ir3j(:), .not.lfill_list(:), keeps%ir3j(:))
            keeps%ir3j(:) = unpack(inserts%ir3j(:), lfill_list(:), keeps%ir3j(:))
            
   
            do i = 1, NDIM
               keeps%xj(i, :) = unpack(keeps%xj(i, :), .not.lfill_list(:), keeps%xj(i, :))
               keeps%xj(i, :) = unpack(inserts%xj(i, :), lfill_list(:), keeps%xj(i, :))
            
               keeps%vj(i, :) = unpack(keeps%vj(i, :), .not.lfill_list(:), keeps%vj(i, :))
               keeps%vj(i, :) = unpack(inserts%vj(i, :), lfill_list(:), keeps%vj(i, :))
            end do
            call util_fill_pl(keeps, inserts, lfill_list)
         class default
            write(*,*) 'Error! fill method called for incompatible return type on whm_pl'
         end select
      end associate
   
      return
   
      end subroutine whm_util_fill_pl

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
         ! Internals
         class(whm_pl), allocatable :: pl_sorted  !! Temporary holder for sorted body
         integer(I4B) :: i
   
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
      end subroutine whm_util_sort_rearrange_pl
   
end submodule s_whm_util
