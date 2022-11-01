!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_sort
   use swiftest
contains

   module subroutine util_sort_body(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest body structure in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self      !! Swiftest body object
      character(*),         intent(in)    :: sortby    !! Sorting attribute
      logical,              intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(body => self, n => self%nbody)
         select case(sortby)
         case("id")
            call util_sort(direction * body%id(1:n), ind)
         case("status")
            call util_sort(direction * body%status(1:n), ind)
         case("ir3h")
            call util_sort(direction * body%ir3h(1:n), ind)
         case("a")
            call util_sort(direction * body%a(1:n), ind)
         case("e")
            call util_sort(direction * body%e(1:n), ind)
         case("inc")
            call util_sort(direction * body%inc(1:n), ind)
         case("capom")
            call util_sort(direction * body%capom(1:n), ind)
         case("mu")
            call util_sort(direction * body%mu(1:n), ind)
         case("lfirst", "nbody", "ldiscard", "xh", "vh", "xb", "vb", "ah", "aobl", "atide", "agr")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not found!'
            return
         end select

         call body%rearrange(ind)

      end associate

      return
   end subroutine util_sort_body


   pure module subroutine util_sort_dp(arr)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array in place into ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(inout) :: arr

      call qsort_DP(arr)

      return
   end subroutine util_sort_dp


   pure module subroutine util_sort_index_dp(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quick sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here).
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine allocates it.
      !!
      implicit none
      ! Arguments
      real(DP),     dimension(:),              intent(in)    :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I4B) :: n, i
      real(DP), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_DP(tmparr, ind)
   
      return
   end subroutine util_sort_index_dp


   recursive pure subroutine qsort_DP(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quicksort sort.
      !!
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(inout)           :: arr
      integer(I4B),dimension(:),intent(out), optional :: ind
      !! Internals
      integer :: iq

      if (size(arr) > 1) then
         if (present(ind)) then
            call partition_DP(arr, iq, ind)
            call qsort_DP(arr(:iq-1),ind(:iq-1))
            call qsort_DP(arr(iq:),  ind(iq:))
         else
            call partition_DP(arr, iq)
            call qsort_DP(arr(:iq-1))
            call qsort_DP(arr(iq:))
         end if
      end if

      return
   end subroutine qsort_DP

 
   pure subroutine partition_DP(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on DP type
      !!
      implicit none
      ! Arguments
      real(DP),     intent(inout), dimension(:)           :: arr
      integer(I4B), intent(inout), dimension(:), optional :: ind
      integer(I4B), intent(out)                           :: marker
      ! Internals
      integer(I4B) :: i, j, itmp, narr, ipiv
      real(DP) :: temp
      real(DP) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2
      x = arr(ipiv)
      i = 0
      j = narr + 1
   
      do
         j = j - 1
         do
            if (arr(j) <= x) exit
            j = j - 1
         end do
         i = i + 1
         do
            if (arr(i) >= x) exit
            i = i + 1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine partition_DP
 

   pure module subroutine util_sort_i4b(arr)
      !! author: David A. Minton
      !!
      !! Sort input integer array in place into ascending numerical order using quick sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(inout) :: arr

      call qsort_I4B(arr)

      return
   end subroutine util_sort_i4b


   pure module subroutine util_sort_index_I4B(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input integer array by index in ascending numerical order using quicksort.
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine allocates it.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:),              intent(in)  :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I4B) :: n, i
      integer(I4B), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_I4B(tmparr, ind)

      return
   end subroutine util_sort_index_I4B


   pure module subroutine util_sort_index_I4B_I8Bind(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input integer array by index in ascending numerical order using quicksort.
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine allocates it.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:),              intent(in)  :: arr
      integer(I8B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I8B) :: n, i
      integer(I4B), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1_I8B, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_I4B_I8Bind(tmparr, ind)

      return
   end subroutine util_sort_index_I4B_I8Bind


   recursive pure subroutine qsort_I4B(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input I4B array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(inout)          :: arr
      integer(I4B), dimension(:), intent(out),  optional :: ind
      ! Internals
      integer(I4B) :: iq

      if (size(arr) > 1) then
         if (present(ind)) then
            call partition_I4B(arr, iq, ind)
            call qsort_I4B(arr(:iq-1),ind(:iq-1))
            call qsort_I4B(arr(iq:),  ind(iq:))
         else
            call partition_I4B(arr, iq)
            call qsort_I4B(arr(:iq-1))
            call qsort_I4B(arr(iq:))
         end if
      end if

      return
   end subroutine qsort_I4B

   recursive pure subroutine qsort_I4B_I8Bind(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input I4B array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(inout)          :: arr
      integer(I8B), dimension(:), intent(out),  optional :: ind
      ! Internals
      integer(I8B) :: iq

      if (size(arr) > 1_I8B) then
         if (present(ind)) then
            call partition_I4B_I8Bind(arr, iq, ind)
            call qsort_I4B_I8Bind(arr(:iq-1_I8B),ind(:iq-1_I8B))
            call qsort_I4B_I8Bind(arr(iq:),  ind(iq:))
         else
            call partition_I4B_I8Bind(arr, iq)
            call qsort_I4B_I8Bind(arr(:iq-1_I8B))
            call qsort_I4B_I8Bind(arr(iq:))
         end if
      end if

      return
   end subroutine qsort_I4B_I8Bind


   recursive pure subroutine qsort_I8B_I8Bind(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input I8B array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      integer(I8B), dimension(:), intent(inout)          :: arr
      integer(I8B), dimension(:), intent(out),  optional :: ind
      ! Internals
      integer(I8B) :: iq

      if (size(arr) > 1_I8B) then
         if (present(ind)) then
            call partition_I8B_I8Bind(arr, iq, ind)
            call qsort_I8B_I8Bind(arr(:iq-1_I8B),ind(:iq-1_I8B))
            call qsort_I8B_I8Bind(arr(iq:),  ind(iq:))
         else
            call partition_I8B_I8Bind(arr, iq)
            call qsort_I8B_I8Bind(arr(:iq-1_I8B))
            call qsort_I8B_I8Bind(arr(iq:))
         end if
      end if

      return
   end subroutine qsort_I8B_I8Bind

 
   pure subroutine partition_I4B(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on I4B type
      !!
      implicit none
      ! Arguments
      integer(I4B), intent(inout), dimension(:)           :: arr
      integer(I4B), intent(inout), dimension(:), optional :: ind
      integer(I4B), intent(out)                           :: marker
      ! Internals
      integer(I4B) :: i, j, itmp, narr, ipiv
      integer(I4B) :: temp
      integer(I4B) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2
      x = arr(ipiv)
      i = 0
      j = narr + 1
   
      do
         j = j - 1
         do
            if (arr(j) <= x) exit
            j = j - 1
         end do
         i = i + 1
         do
            if (arr(i) >= x) exit
            i = i + 1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine partition_I4B

   pure subroutine partition_I4B_I8Bind(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on I4B type
      !!
      implicit none
      ! Arguments
      integer(I4B), intent(inout), dimension(:)           :: arr
      integer(I8B), intent(inout), dimension(:), optional :: ind
      integer(I8B), intent(out)                           :: marker
      ! Internals
      integer(I8B) :: i, j, itmp, narr, ipiv
      integer(I4B) :: temp
      integer(I8B) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2_I8B
      x = arr(ipiv)
      i = 0_I8B
      j = narr + 1_I8B
   
      do
         j = j - 1_I8B
         do
            if (arr(j) <= x) exit
            j = j - 1_I8B
         end do
         i = i + 1_I8B
         do
            if (arr(i) >= x) exit
            i = i + 1_I8B
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1_I8B
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine partition_I4B_I8Bind

   pure subroutine partition_I8B_I8Bind(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on I8B type with I8B index
      !!
      implicit none
      ! Arguments
      integer(I8B), intent(inout), dimension(:)           :: arr
      integer(I8B), intent(inout), dimension(:), optional :: ind
      integer(I8B), intent(out)                           :: marker
      ! Internals
      integer(I8B) :: i, j, itmp, narr, ipiv
      integer(I8B) :: temp
      integer(I8B) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2_I8B
      x = arr(ipiv)
      i = 0_I8B
      j = narr + 1_I8B
   
      do
         j = j - 1_I8B
         do
            if (arr(j) <= x) exit
            j = j - 1_I8B
         end do
         i = i + 1_I8B
         do
            if (arr(i) >= x) exit
            i = i + 1_I8B
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1_I8B
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine partition_I8B_I8Bind


   pure module subroutine util_sort_sp(arr)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array in place into ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(inout) :: arr

      call qsort_SP(arr)

      return
   end subroutine util_sort_sp


   pure module subroutine util_sort_index_sp(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quicksort.
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine allocates it.
      !!
      implicit none
      ! Arguments
      real(SP),     dimension(:),              intent(in)    :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I4B) :: n, i
      real(SP), dimension(:), allocatable :: tmparr

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      allocate(tmparr, mold=arr)
      tmparr(:) = arr(ind(:))
      call qsort_SP(tmparr, ind)
   
      return
   end subroutine util_sort_index_sp


   recursive pure subroutine qsort_SP(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input DP precision array by index in ascending numerical order using quicksort.
      !!
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(inout)           :: arr
      integer(I4B),dimension(:),intent(out), optional :: ind
      !! Internals
      integer :: iq

      if (size(arr) > 1) then
         if (present(ind)) then
            call partition_SP(arr, iq, ind)
            call qsort_SP(arr(:iq-1),ind(:iq-1))
            call qsort_SP(arr(iq:),  ind(iq:))
         else
            call partition_SP(arr, iq)
            call qsort_SP(arr(:iq-1))
            call qsort_SP(arr(iq:))
         end if
      end if

      return
   end subroutine qsort_SP


   pure subroutine partition_SP(arr, marker, ind)
      !! author: David A. Minton
      !!
      !! Partition function for quicksort on SP type
      !!
      implicit none
      ! Arguments
      real(SP),     intent(inout), dimension(:)           :: arr
      integer(I4B), intent(inout), dimension(:), optional :: ind
      integer(I4B), intent(out)                           :: marker
      ! Internals
      integer(I4B) :: i, j, itmp, narr, ipiv
      real(SP) :: temp
      real(SP) :: x   ! pivot point

      narr = size(arr)

      ! Get center as pivot, as this is likely partially sorted
      ipiv = narr / 2
      x = arr(ipiv)
      i = 0
      j = narr + 1
   
      do
         j = j - 1
         do
            if (arr(j) <= x) exit
            j = j - 1
         end do
         i = i + 1
         do
            if (arr(i) >= x) exit
            i = i + 1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            if (present(ind)) then
               itmp = ind(i)
               ind(i) = ind(j)
               ind(j) = itmp
            end if
         else if (i == j) then
            marker = i + 1
            return
         else
            marker = i
            return
         endif
      end do
  
      return
   end subroutine partition_SP


   module subroutine util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self      !! Swiftest massive body object
      character(*),       intent(in)    :: sortby    !! Sorting attribute
      logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("Gmass","mass")
            call util_sort(direction * pl%Gmass(1:npl), ind)
         case("rhill")
            call util_sort(direction * pl%rhill(1:npl), ind)
         case("renc")
            call util_sort(direction * pl%renc(1:npl), ind)
         case("radius")
            call util_sort(direction * pl%radius(1:npl), ind)
         case("density")
            call util_sort(direction * pl%density(1:npl), ind)
         case("k2")
            call util_sort(direction * pl%k2(1:npl), ind)
         case("Q")
            call util_sort(direction * pl%Q(1:npl), ind)
         case("tlag")
            call util_sort(direction * pl%tlag(1:npl), ind)
         case("xbeg", "xend", "vbeg", "Ip", "rot", "k_plpl", "nplpl")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call util_sort_body(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)

      end associate

      return
   end subroutine util_sort_pl


   module subroutine util_sort_tp(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest test particle object  in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self      !! Swiftest test particle object
      character(*),       intent(in)    :: sortby    !! Sorting attribute
      logical,            intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("peri")
            call util_sort(direction * tp%peri(1:ntp), ind)
         case("atp")
            call util_sort(direction * tp%atp(1:ntp), ind)
         case("isperi")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call util_sort_body(tp, sortby, ascending)
            return
         end select

         call tp%rearrange(ind)

      end associate

      return
   end subroutine util_sort_tp


   module subroutine util_sort_rearrange_body(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(swiftest_body),               intent(inout) :: self !! Swiftest body object
      integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(n => self%nbody)
         call util_sort_rearrange(self%id,       ind, n)
         call util_sort_rearrange(self%info,     ind, n)
         call util_sort_rearrange(self%status,   ind, n)
         call util_sort_rearrange(self%ldiscard, ind, n)
         call util_sort_rearrange(self%xh,       ind, n)
         call util_sort_rearrange(self%vh,       ind, n)
         call util_sort_rearrange(self%xb,       ind, n)
         call util_sort_rearrange(self%vb,       ind, n)
         call util_sort_rearrange(self%ah,       ind, n)
         call util_sort_rearrange(self%ir3h,     ind, n)
         call util_sort_rearrange(self%mu,       ind, n)
         call util_sort_rearrange(self%lmask,    ind, n)
         call util_sort_rearrange(self%a,        ind, n)
         call util_sort_rearrange(self%e,        ind, n)
         call util_sort_rearrange(self%inc,      ind, n)
         call util_sort_rearrange(self%capom,    ind, n)
         call util_sort_rearrange(self%omega,    ind, n)
         call util_sort_rearrange(self%capm,     ind, n)
         call util_sort_rearrange(self%aobl,     ind, n)
         call util_sort_rearrange(self%atide,    ind, n)
         call util_sort_rearrange(self%agr,      ind, n)
      end associate

      return
   end subroutine util_sort_rearrange_body


   pure module subroutine util_sort_rearrange_arr_char_string(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of character string in-place from an index list.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B),          dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                                     intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      character(len=STRMAX), dimension(:), allocatable                :: tmp !! Temporary copy of arry used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_char_string


   pure module subroutine util_sort_rearrange_arr_DP(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of DP type in-place from an index list.
      implicit none
      ! Arguments
      real(DP),     dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),              intent(in)  :: ind !! Index to rearrange against
      integer(I4B),                            intent(in)  :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      real(DP), dimension(:), allocatable :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_DP


   pure module subroutine util_sort_rearrange_arr_DPvec(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of (NDIM,n) DP-type vectors in-place from an index list.
      implicit none
      ! Arguments
      real(DP),     dimension(:,:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),                intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                              intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      real(DP), dimension(:,:), allocatable :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(:,1:n) = arr(:, ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_DPvec


   pure module subroutine util_sort_rearrange_arr_I4B(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of integers in-place from an index list.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                             intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      integer(I4B), dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_I4B

   pure module subroutine util_sort_rearrange_arr_I4B_I8Bind(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of integers in-place from an index list.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I8B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I8B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      integer(I4B), dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0_I8B) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_I4B_I8Bind


   pure module subroutine util_sort_rearrange_arr_logical(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of logicals in-place from an index list.
      implicit none
      ! Arguments
      logical,      dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      logical, dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_logical


   pure module subroutine util_sort_rearrange_arr_logical_I8Bind(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of logicals in-place from an index list.
      implicit none
      ! Arguments
      logical,      dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I8B), dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I8B),                            intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      logical, dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)
      tmp(1:n) = arr(ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_logical_I8Bind


   module subroutine util_sort_rearrange_arr_info(arr, ind, n)
      !! author: David A. Minton
      !!
      !! Rearrange a single array of particle information type in-place from an index list.
      implicit none
      ! Arguments
      type(swiftest_particle_info),  dimension(:), allocatable, intent(inout) :: arr !! Destination array 
      integer(I4B),                  dimension(:),              intent(in)    :: ind !! Index to rearrange against
      integer(I4B),                                             intent(in)    :: n   !! Number of elements in arr and ind to rearrange
      ! Internals
      type(swiftest_particle_info),  dimension(:), allocatable                :: tmp !! Temporary copy of array used during rearrange operation

      if (.not. allocated(arr) .or. n <= 0) return
      allocate(tmp, mold=arr)

      call util_copy_particle_info_arr(arr, tmp, ind)
      call move_alloc(tmp, arr)

      return
   end subroutine util_sort_rearrange_arr_info


   module subroutine util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      class(swiftest_pl),               intent(inout) :: self !! Swiftest massive body object
      integer(I4B),       dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange(pl%mass,    ind, npl)
         call util_sort_rearrange(pl%Gmass,   ind, npl)
         call util_sort_rearrange(pl%rhill,   ind, npl)
         call util_sort_rearrange(pl%xbeg,    ind, npl)
         call util_sort_rearrange(pl%vbeg,    ind, npl)
         call util_sort_rearrange(pl%radius,  ind, npl)
         call util_sort_rearrange(pl%density, ind, npl)
         call util_sort_rearrange(pl%Ip,      ind, npl)
         call util_sort_rearrange(pl%rot,     ind, npl)
         call util_sort_rearrange(pl%k2,      ind, npl)
         call util_sort_rearrange(pl%Q,       ind, npl)
         call util_sort_rearrange(pl%tlag,    ind, npl)

         if (allocated(pl%k_plpl)) deallocate(pl%k_plpl)

         call util_sort_rearrange_body(pl, ind)
      end associate

      return
   end subroutine util_sort_rearrange_pl


   module subroutine util_sort_rearrange_tp(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(swiftest_tp),                 intent(inout) :: self !! Swiftest test particle object
      integer(I4B),         dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange(tp%isperi, ind, ntp)
         call util_sort_rearrange(tp%peri,   ind, ntp)
         call util_sort_rearrange(tp%atp,    ind, ntp)

         call util_sort_rearrange_body(tp, ind)
      end associate

      return
   end subroutine util_sort_rearrange_tp

end submodule s_util_sort
