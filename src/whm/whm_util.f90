! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(whm) s_whm_util
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
         call util_append(self%eta, source%eta, lsource_mask=lsource_mask)
         call util_append(self%muj, source%muj, lsource_mask=lsource_mask)
         call util_append(self%ir3j, source%ir3j, lsource_mask=lsource_mask)
         call util_append(self%xj, source%xj, lsource_mask=lsource_mask)
         call util_append(self%vj, source%vj, lsource_mask=lsource_mask)

         call swiftest_util_append_pl(self, source, lsource_mask)
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class whm_pl or its descendents"
         call base_util_exit(FAILURE)
      end select

      return
   end subroutine whm_util_append_pl


   module subroutine whm_util_dealloc_pl(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(whm_pl),  intent(inout) :: self !! WHM massive body object

      if (allocated(self%eta)) deallocate(self%eta)
      if (allocated(self%muj)) deallocate(self%muj)
      if (allocated(self%xj)) deallocate(self%xj)
      if (allocated(self%vj)) deallocate(self%vj)
      if (allocated(self%ir3j)) deallocate(self%ir3j)

      call swiftest_util_dealloc_pl(self)

      return
   end subroutine whm_util_dealloc_pl


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
   
      associate(keeps => self)
         select type(inserts)
         class is (whm_pl)
            call util_fill(keeps%eta, inserts%eta, lfill_list)
            call util_fill(keeps%muj, inserts%muj, lfill_list)
            call util_fill(keeps%ir3j, inserts%ir3j, lfill_list)
            call util_fill(keeps%xj, inserts%xj, lfill_list)
            call util_fill(keeps%vj, inserts%vj, lfill_list)

            call swiftest_util_fill_pl(keeps, inserts, lfill_list)
         class default
            write(*,*) "Invalid object passed to the fill method. Inserts must be of class whm_pl or its descendents!"
            call base_util_exit(FAILURE)
         end select
      end associate
   
      return
   end subroutine whm_util_fill_pl


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

      call swiftest_util_resize_pl(self, nnew)

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
            r2 = dot_product(self%rh(:, i), self%rh(:, i))
            ir = 1.0_DP / sqrt(r2)
            self%ir3h(i) = ir / r2
            r2 = dot_product(self%xj(:, i), self%xj(:, i))
            ir = 1.0_DP / sqrt(r2)
            self%ir3j(i) = ir / r2
         end do
      end if

      return
   end subroutine whm_util_set_ir3j


   module subroutine whm_util_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate WHM planet structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_util_setup.f90
      implicit none
      ! Arguments
      class(whm_pl),             intent(inout) :: self  !! Swiftest test particle object
      integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class
      call swiftest_util_setup_pl(self, n, param) 
      if (n == 0) return

      allocate(self%eta(n))
      allocate(self%muj(n))
      allocate(self%xj(NDIM, n))
      allocate(self%vj(NDIM, n))
      allocate(self%ir3j(n))

      self%eta(:)   = 0.0_DP
      self%muj(:)   = 0.0_DP
      self%xj(:,:)  = 0.0_DP
      self%vj(:,:)  = 0.0_DP
      self%ir3j(:) = 0.0_DP

      return
   end subroutine whm_util_setup_pl 


   module subroutine whm_util_set_mu_eta_pl(self, cb)
      !! author: David A. Minton
      !!
      !! Sets the Jacobi mass value eta for all massive bodies
      implicit none
      ! Arguments
      class(whm_pl),      intent(inout) :: self   !! WHM nbody_system object
      class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body object
      ! Internals
      integer(I4B)                                 :: i

      associate(pl => self, npl => self%nbody)
         if (npl == 0) return
         call swiftest_util_set_mu_pl(pl, cb)
         pl%eta(1) = cb%Gmass + pl%Gmass(1)
         pl%muj(1) = pl%eta(1)
         do i = 2, npl
            pl%eta(i) = pl%eta(i - 1) + pl%Gmass(i)
            pl%muj(i) = cb%Gmass * pl%eta(i) / pl%eta(i - 1)
         end do
      end associate

      return
   end subroutine whm_util_set_mu_eta_pl


   module subroutine whm_util_setup_initialize_system(self, system_history, param)
      !! author: David A. Minton
      !!
      !! Initialize a WHM nbody system from files
      !!
      implicit none
      ! Arguments
      class(whm_nbody_system),                 intent(inout) :: self            !! WHM nbody system object
      class(swiftest_storage),    allocatable, intent(inout) :: system_history  !! Stores the system history between output dumps
      class(swiftest_parameters),              intent(inout) :: param           !! Current run configuration parameters 

      call swiftest_util_setup_initialize_system(self, system_history, param)
      ! First we need to make sure that the massive bodies are sorted by heliocentric distance before computing jacobies
      call swiftest_util_set_ir3h(self%pl)
      call self%pl%sort("ir3h", ascending=.false.)
      call self%pl%flatten(param)

      ! Make sure that the discard list gets allocated initially
      call self%pl_discards%setup(0, param)
      call self%tp_discards%setup(0, param)
      call self%pl%set_mu(self%cb)
      call self%tp%set_mu(self%cb)
      if (param%lgr .and. param%in_type == "ASCII") then !! pseudovelocity conversion for NetCDF input files is handled by NetCDF routines
         call self%pl%v2pv(param)
         call self%tp%v2pv(param)
      end if

      return
   end subroutine whm_util_setup_initialize_system


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
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B) :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("eta")
            call util_sort(direction * pl%eta(1:npl), ind)
         case("muj")
            call util_sort(direction * pl%muj(1:npl), ind)
         case("ir3j")
            call util_sort(direction * pl%ir3j(1:npl), ind)
         case("xj", "vj")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default
            call swiftest_util_sort_pl(pl, sortby, ascending)
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

         call swiftest_util_sort_rearrange_pl(pl,ind)
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

      associate(keeps => self)
         select type(discards)
         class is (whm_pl)
            call util_spill(keeps%eta, discards%eta, lspill_list, ldestructive)
            call util_spill(keeps%muj, discards%muj, lspill_list, ldestructive)
            call util_spill(keeps%ir3j, discards%ir3j, lspill_list, ldestructive)
            call util_spill(keeps%xj, discards%xj, lspill_list, ldestructive)
            call util_spill(keeps%vj, discards%vj, lspill_list, ldestructive)

            call swiftest_util_spill_pl(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class whm_pl or its descendents!"
            call base_util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine whm_util_spill_pl
   
end submodule s_whm_util
