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


   module subroutine util_get_idvalues_system(self, idvals)
      !! author: David A. Minton
      !!
      !! Returns an array of all id values saved in this snapshot
      implicit none
      ! Arguments
      class(swiftest_nbody_system),            intent(in)  :: self   !! Encounter snapshot object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      ! Internals
      integer(I4B) :: npl, ntp

      if (allocated(self%pl)) then
         npl = self%pl%nbody
      else
         npl = 0
      end if 
      if (allocated(self%tp)) then
         ntp = self%tp%nbody
      else
         ntp = 0
      end if

      allocate(idvals(1 + npl+ntp))

      idvals(1) = self%cb%id
      if (npl > 0) idvals(2:npl+1) = self%pl%id(:)
      if (ntp > 0) idvals(npl+2:npl+ntp+1) = self%tp%id(:)

      return

   end subroutine util_get_idvalues_system


   subroutine util_get_vals_storage(storage, idvals, tvals)
      !! author: David A. Minton
      !!
      !! Gets the id values in a storage object, regardless of whether it is encounter of collision
      ! Argument
      class(swiftest_storage(*)), intent(in)               :: storage !! Swiftest storage object
      integer(I4B), dimension(:), allocatable, intent(out) :: idvals  !! Array of all id values in all snapshots
      real(DP),     dimension(:), allocatable, intent(out) :: tvals   !! Array of all time values in all snapshots
      ! Internals
      integer(I4B) :: i, n, nlo, nhi, ntotal
      integer(I4B), dimension(:), allocatable :: itmp

      associate(nsnaps => storage%iframe)

         allocate(tvals(nsnaps))
         tvals(:) = 0.0_DP

         ! First pass to get total number of ids
         ntotal = 0
         do i = 1, nsnaps
            if (allocated(storage%frame(i)%item)) then
               select type(snapshot => storage%frame(i)%item)
               class is (swiftest_nbody_system)
                  tvals(i) = snapshot%t
                  call snapshot%get_idvals(itmp)
                  if (allocated(itmp)) then
                     n = size(itmp)
                     ntotal = ntotal + n
                  end if
               end select
            end if
         end do

         allocate(idvals(ntotal))
         nlo = 1
         ! Second pass to store all ids get all of the ids stored
         do i = 1, nsnaps
            if (allocated(storage%frame(i)%item)) then
               select type(snapshot => storage%frame(i)%item)
               class is (swiftest_nbody_system)
                  tvals(i) = snapshot%t
                  call snapshot%get_idvals(itmp)
                  if (allocated(itmp)) then
                     n = size(itmp)
                     nhi = nlo + n - 1
                     idvals(nlo:nhi) = itmp(1:n)
                     nlo = nhi + 1 
                  end if
               end select
            end if
         end do

      end associate 
      return
   end subroutine util_get_vals_storage


   module subroutine util_index_map_storage(self)
      !! author: David A. Minton
      !!
      !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      implicit none
      ! Arguments
      class(swiftest_storage(*)), intent(inout) :: self  !! Swiftest storage object
      ! Internals
      integer(I4B), dimension(:), allocatable :: idvals
      real(DP), dimension(:), allocatable :: tvals
 
      call util_get_vals_storage(self, idvals, tvals)

      call util_unique(idvals,self%idvals,self%idmap)
      self%nid = size(self%idvals)

      call util_unique(tvals,self%tvals,self%tmap)
      self%nt = size(self%tvals)

      return
   end subroutine util_index_map_storage

end submodule s_util_index_array