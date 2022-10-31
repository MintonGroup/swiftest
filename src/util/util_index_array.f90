!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_index_array
   use swiftest
contains

   module subroutine util_index_array(ind_arr, n)
      !! author: David A. Minton
      !!
      !! Creates or resizes an index array of size n where ind_arr = [1, 2, ... n].
      !! This subroutine assumes that if ind_arr is already allocated, it is a pre-existing index array of a different size.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: ind_arr !! Index array. Input is a pre-existing index array where n /= size(ind_arr). Output is a new index array ind_arr = [1, 2, ... n]
      integer(I4B),                            intent(in)    :: n       !! The new size of the index array
      ! Internals
      integer(I4B) :: nold, i
      integer(I4B), dimension(:), allocatable :: itmp

      if (allocated(ind_arr)) then
         nold = size(ind_arr)
         if (nold == n) return ! Nothing to do, so go home
      else
         nold = 0
      end if

      allocate(itmp(n))
      if (n >= nold) then
         if (nold > 0) itmp(1:nold) = ind_arr(1:nold)
         itmp(nold+1:n) = [(i, i = nold + 1, n)]
         call move_alloc(itmp, ind_arr)
      else
         itmp(1:n) = ind_arr(1:n)
         call move_alloc(itmp, ind_arr)
      end if

      return
   end subroutine util_index_array

end submodule s_util_index_array