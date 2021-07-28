submodule (swiftest_classes) s_util_sort
   use swiftest
contains
   module subroutine util_sort_body(self, sortby, reverse)
      !! author: David A. Minton
      !!
      !! Sort a Swiftest body structure in-place. 
      !! sortby is a string. The only valid input the body class takes is "id," which is also the default value. 
      !! Sort order is ascending order by default. Set reverse=.true. to sort in descending order.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self    !! Swiftest body object
      character(*),         intent(in)    :: sortby  !! Sorting attribute
      logical,              intent(in)    :: reverse !! Logical flag indicating whether or not the sorting should be in reverse (descending order)
      ! Internals
      class(swiftest_body),       allocatable :: body_sorted  !! Temporary holder for sorted body
      integer(I4B), dimension(:), allocatable :: ind

      associate(n => self%nbody)
         allocate(body_sorted, source=self)
         allocate(ind(n))
         select case(sortby)
         case("id")
            if (reverse) then
               call util_sort(-self%id(1:n), ind(1:n))
            else
               call util_sort(self%id(1:n), ind(1:n))
            end if
         end select

         self%id(1:n) = body_sorted%id(ind(1:n))
         self%name(1:n) = body_sorted%name(ind(1:n))
         self%status(1:n) = body_sorted%status(ind(1:n))
         self%ldiscard(1:n) = body_sorted%ldiscard(ind(1:n))
         self%xh(:,1:n) = body_sorted%xh(:,ind(1:n))
         self%vh(:,1:n) = body_sorted%vh(:,ind(1:n))
         self%xb(:,1:n) = body_sorted%xb(:,ind(1:n))
         self%vb(:,1:n) = body_sorted%vb(:,ind(1:n))
         self%ah(:,1:n) = body_sorted%ah(:,ind(1:n))
         self%aobl(:,1:n) = body_sorted%aobl(:,ind(1:n))
         self%atide(:,1:n) = body_sorted%atide(:,ind(1:n))
         self%agr(:,1:n) = body_sorted%agr(:,ind(1:n))
         self%ir3h(1:n) = body_sorted%ir3h(ind(1:n))
         self%a(1:n) = body_sorted%a(ind(1:n))
         self%e(1:n) = body_sorted%e(ind(1:n))
         self%inc(1:n) = body_sorted%inc(ind(1:n))
         self%capom(1:n) = body_sorted%capom(ind(1:n))
         self%omega(1:n) = body_sorted%omega(ind(1:n))
         self%capm(1:n) = body_sorted%capm(ind(1:n))
         self%mu(1:n) = body_sorted%mu(ind(1:n))
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
