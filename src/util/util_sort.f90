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
      integer(I4B), dimension(self%nbody) :: ind
      integer(I4B)                        :: direction

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(body => self, n => self%nbody)
         select case(sortby)
         case("id")
            call util_sort(direction * body%id(1:n), ind(1:n))
         case("status")
            call util_sort(direction * body%status(1:n), ind(1:n))
         case("ir3h")
            call util_sort(direction * body%ir3h(1:n), ind(1:n))
         case("a")
            call util_sort(direction * body%a(1:n), ind(1:n))
         case("e")
            call util_sort(direction * body%e(1:n), ind(1:n))
         case("inc")
            call util_sort(direction * body%inc(1:n), ind(1:n))
         case("capom")
            call util_sort(direction * body%capom(1:n), ind(1:n))
         case("mu")
            call util_sort(direction * body%mu(1:n), ind(1:n))
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
      integer(I4B), dimension(self%nbody) :: ind
      integer(I4B)                        :: direction

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("Gmass","mass")
            call util_sort(direction * pl%Gmass(1:npl), ind(1:npl))
         case("rhill")
            call util_sort(direction * pl%rhill(1:npl), ind(1:npl))
         case("radius")
            call util_sort(direction * pl%radius(1:npl), ind(1:npl))
         case("density")
            call util_sort(direction * pl%density(1:npl), ind(1:npl))
         case("k2")
            call util_sort(direction * pl%k2(1:npl), ind(1:npl))
         case("Q")
            call util_sort(direction * pl%Q(1:npl), ind(1:npl))
         case("tlag")
            call util_sort(direction * pl%tlag(1:npl), ind(1:npl))
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
      integer(I4B), dimension(self%nbody) :: ind
      integer(I4B)                        :: direction

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("peri")
            call util_sort(direction * tp%peri(1:ntp), ind(1:ntp))
         case("atp")
            call util_sort(direction * tp%atp(1:ntp), ind(1:ntp))
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
      ! Internals
      class(swiftest_body), allocatable :: body_sorted  !! Temporary holder for sorted body

      associate(n => self%nbody)
         allocate(body_sorted, source=self)
         self%id(1:n) = body_sorted%id(ind(1:n))
         self%name(1:n) = body_sorted%name(ind(1:n))
         self%status(1:n) = body_sorted%status(ind(1:n))
         self%ldiscard(1:n) = body_sorted%ldiscard(ind(1:n))
         self%xh(:,1:n) = body_sorted%xh(:,ind(1:n))
         self%vh(:,1:n) = body_sorted%vh(:,ind(1:n))
         self%xb(:,1:n) = body_sorted%xb(:,ind(1:n))
         self%vb(:,1:n) = body_sorted%vb(:,ind(1:n))
         self%ah(:,1:n) = body_sorted%ah(:,ind(1:n))
         self%ir3h(1:n) = body_sorted%ir3h(ind(1:n))
         self%a(1:n) = body_sorted%a(ind(1:n))
         self%e(1:n) = body_sorted%e(ind(1:n))
         self%inc(1:n) = body_sorted%inc(ind(1:n))
         self%capom(1:n) = body_sorted%capom(ind(1:n))
         self%omega(1:n) = body_sorted%omega(ind(1:n))
         self%capm(1:n) = body_sorted%capm(ind(1:n))
         self%mu(1:n) = body_sorted%mu(ind(1:n))
         if (allocated(self%aobl))  self%aobl(:,1:n) = body_sorted%aobl(:,ind(1:n))
         if (allocated(self%atide)) self%atide(:,1:n) = body_sorted%atide(:,ind(1:n))
         if (allocated(self%agr))   self%agr(:,1:n) = body_sorted%agr(:,ind(1:n))
         deallocate(body_sorted)
      end associate

      return
   end subroutine util_sort_rearrange_body


   module subroutine util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange Swiftest massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      class(swiftest_pl),               intent(inout) :: self !! Swiftest massive body object
      integer(I4B),       dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      ! Internals
      class(swiftest_pl), allocatable :: pl_sorted  !! Temporary holder for sorted body

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange_body(pl,ind)
         allocate(pl_sorted, source=self)
         pl%mass(1:npl) = pl_sorted%mass(ind(1:npl))
         pl%Gmass(1:npl) = pl_sorted%Gmass(ind(1:npl))
         pl%rhill(1:npl) = pl_sorted%rhill(ind(1:npl))
         pl%xbeg(:,1:npl) = pl_sorted%xbeg(:,ind(1:npl))
         pl%xend(:,1:npl) = pl_sorted%xend(:,ind(1:npl))
         pl%vbeg(:,1:npl) = pl_sorted%vbeg(:,ind(1:npl))
         if (allocated(pl%radius))  pl%radius(1:npl) = pl_sorted%radius(ind(1:npl))
         if (allocated(pl%density)) pl%density(1:npl) = pl_sorted%density(ind(1:npl))
         if (allocated(pl%Ip))      pl%Ip(:,1:npl) = pl_sorted%Ip(:,ind(1:npl))
         if (allocated(pl%rot))     pl%rot(:,1:npl) = pl_sorted%rot(:,ind(1:npl))
         if (allocated(pl%k2))      pl%k2(1:npl) = pl_sorted%k2(ind(1:npl))
         if (allocated(pl%Q))       pl%Q(1:npl) = pl_sorted%Q(ind(1:npl))
         if (allocated(pl%tlag))    pl%tlag(1:npl) = pl_sorted%tlag(ind(1:npl))

         deallocate(pl_sorted)
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
      ! Internals
      class(swiftest_tp), allocatable :: tp_sorted  !! Temporary holder for sorted body

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange_body(tp,ind)
         allocate(tp_sorted, source=self)
         tp%isperi(1:ntp) = tp_sorted%isperi(ind(1:ntp))
         tp%peri(1:ntp) = tp_sorted%peri(ind(1:ntp))
         tp%atp(1:ntp) = tp_sorted%atp(ind(1:ntp))
         deallocate(tp_sorted)
      end associate

      return
   end subroutine util_sort_rearrange_tp


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
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !!
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(in)  :: arr
      integer(I4B), dimension(:), intent(out) :: ind
      ! Internals
      real(DP) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      ind = [(i, i=1, n)]
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
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !!
      implicit none
      ! Arguments
      integer(I4B), dimension(:), intent(in)  :: arr
      integer(I4B), dimension(:), intent(out) :: ind
      ! Internals
      integer(I4B) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      ind = [(i, i=1, n)]
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
      !! This algorithm works well for partially sorted arrays (which is usually the case here)
      !!
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(in)  :: arr
      integer(I4B), dimension(:), intent(out) :: ind
      ! Internals
      real(SP) :: tmp
      integer(I4B) :: n, i, j

      n = size(arr)
      ind = [(i, i=1, n)]
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

end submodule s_util_sort
