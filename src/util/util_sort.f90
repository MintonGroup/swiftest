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


   module subroutine util_sort_dp(arr)
      !! author: David A. Minton
      !!
      !! Sort input double precision array in place into ascending numerical order using insertion sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !!
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(inout) :: arr
      ! Internals
      real(DP)                :: tmp
      integer(I4B)            :: n, i, j

      n = size(arr)
      do i = 2, n
         tmp = arr(i)
         do j = i - 1, 1, -1
            if (arr(j) <= tmp) exit
            arr(j + 1) = arr(j)
         end do
         arr(j + 1) = tmp
      end do

      return
   end subroutine util_sort_dp


   module subroutine util_sort_index_dp(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input double precision array by index in ascending numerical order using insertion sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here).
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine allocates it.
      !!
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(in)  :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      real(DP) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      do i = 2, n
         tmp = arr(ind(i))
         do j = i - 1, 1, -1
            if (arr(ind(j)) <= tmp) exit
            ind(j + 1) = ind(j)
         end do
         ind(j + 1) = i
      end do

      return
   end subroutine util_sort_index_dp


   module subroutine util_sort_i4b(arr)
      !! author: David A. Minton
      !!
      !! Sort input integer array in place into ascending numerical order using insertion sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(inout) :: arr
      ! Internals
      integer(I4B) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      do i = 2, n
         tmp = arr(i)
         do j = i - 1, 1, -1
            if (arr(j) <= tmp) exit
            arr(j + 1) = arr(j)
         end do
         arr(j + 1) = tmp
      end do

      return
   end subroutine util_sort_i4b


   module subroutine util_sort_index_i4b(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input integer array by index in ascending numerical order using insertion sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here).
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine allocates it.
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(in)  :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      integer(I4B) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      do i = 2, n
         tmp = arr(ind(i))
         do j = i - 1, 1, -1
            if (arr(ind(j)) <= tmp) exit
            ind(j + 1) = ind(j)
         end do
         ind(j + 1) = i
      end do

      return
   end subroutine util_sort_index_i4b


   module subroutine util_sort_sp(arr)
      !! author: David A. Minton
      !!
      !! Sort input single precision array in place into ascending numerical order using insertion sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(inout) :: arr
      ! Internals
      real(SP) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      do i = 2, n
         tmp = arr(i)
         do j = i - 1, 1, -1
            if (arr(j) <= tmp) exit
            arr(j + 1) = arr(j)
         end do
         arr(j + 1) = tmp
      end do

      return
   end subroutine util_sort_sp


   module subroutine util_sort_index_sp(arr, ind)
      !! author: David A. Minton
      !!
      !! Sort input single precision array by index in ascending numerical order using insertion sort.
      !! This algorithm works well for partially sorted arrays (which is usually the case here).
      !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously
      !! sorted array). If it is not allocated, this subroutine allocates it.
      !!
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(in)  :: arr
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind
      ! Internals
      real(SP) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      if (.not.allocated(ind)) then
         allocate(ind(n))
         ind = [(i, i=1, n)]
      end if
      do i = 2, n
         tmp = arr(ind(i))
         do j = i - 1, 1, -1
            if (arr(ind(j)) <= tmp) exit
            ind(j + 1) = ind(j)
         end do
         ind(j + 1) = i
      end do

      return
   end subroutine util_sort_index_sp


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


   module subroutine util_sort_rearrange_arr_char_string(arr, ind, n)
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


   module subroutine util_sort_rearrange_arr_DP(arr, ind, n)
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


   module subroutine util_sort_rearrange_arr_DPvec(arr, ind, n)
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


   module subroutine util_sort_rearrange_arr_I4B(arr, ind, n)
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


   module subroutine util_sort_rearrange_arr_logical(arr, ind, n)
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
      integer(I4B) :: i


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
