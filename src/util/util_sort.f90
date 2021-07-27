submodule (swiftest_classes) s_util_sort
   use swiftest
contains
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
