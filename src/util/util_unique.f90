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

   module subroutine util_unique(input_array, output_array)
      !! author: David A. Minton
      !!
      !! Takes an input unsorted integer array and returns a new array of sorted, unique values
      implicit none
      ! Arguments
      integer(I4B), dimension(:),              intent(in)  :: input_array
      integer(I4B), dimension(:), allocatable, intent(out) :: output_array
      ! Internals
      integer(I4B), dimension(:), allocatable :: unique_array
      integer(I4B) :: n, lo, hi

      allocate(unique_array, mold=input_array)
      lo = minval(input_array) - 1
      hi = maxval(input_array)

      n = 0
      do while (lo < hi)
         n = n + 1
         lo = minval(input_array, mask=input_array > lo)
         unique_array(n) = lo
      enddo
      allocate(output_array(n), source=unique_array(1:n)) 

      return
   end subroutine util_unique

end submodule s_util_unique