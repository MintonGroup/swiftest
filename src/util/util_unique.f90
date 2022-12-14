!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_unique
   use swiftest
contains

   module subroutine util_unique_DP(input_array, output_array, index_map)
      !! author: David A. Minton
      !!
      !! Takes an input unsorted integer array and returns a new array of sorted, unique values (DP version)
      implicit none
      ! Arguments
      real(DP),     dimension(:),              intent(in)  :: input_array  !! Unsorted input array 
      real(DP),     dimension(:), allocatable, intent(out) :: output_array !! Sorted array of unique values 
      integer(I4B), dimension(:), allocatable, intent(out) :: index_map    !! An array of the same size as input_array that such that any for any index i, output_array(index_map(i)) = input_array(i)       
      ! Internals
      real(DP), dimension(:), allocatable :: unique_array
      integer(I4B) :: n
      real(DP) :: lo, hi

      allocate(unique_array, mold=input_array)
      allocate(index_map(size(input_array)))
      lo = minval(input_array) - 1
      hi = maxval(input_array)

      n = 0
      do 
         n = n + 1
         lo = minval(input_array(:), mask=input_array(:) > lo)
         unique_array(n) = lo
         where(input_array(:) == lo) index_map(:) = n
         if (lo >= hi) exit
      enddo
      allocate(output_array(n), source=unique_array(1:n)) 

      return
   end subroutine util_unique_DP


   module subroutine util_unique_I4B(input_array, output_array, index_map)
      !! author: David A. Minton
      !!
      !! Takes an input unsorted integer array and returns a new array of sorted, unique values (I4B version)
      implicit none
      ! Arguments
      integer(I4B), dimension(:),              intent(in)  :: input_array  !! Unsorted input array 
      integer(I4B), dimension(:), allocatable, intent(out) :: output_array !! Sorted array of unique values
      integer(I4B), dimension(:), allocatable, intent(out) :: index_map    !! An array of the same size as input_array that such that any for any index i, output_array(index_map(i)) = input_array(i)     
      ! Internals
      integer(I4B), dimension(:), allocatable :: unique_array
      integer(I4B) :: n, lo, hi

      allocate(unique_array, mold=input_array)
      allocate(index_map, mold=input_array)
      lo = minval(input_array) - 1
      hi = maxval(input_array)

      n = 0
      do 
         n = n + 1
         lo = minval(input_array(:), mask=input_array(:) > lo)
         unique_array(n) = lo
         where(input_array(:) == lo) index_map(:) = n
         if (lo >= hi) exit
      enddo
      allocate(output_array(n), source=unique_array(1:n)) 

      return
   end subroutine util_unique_I4B



end submodule s_util_unique