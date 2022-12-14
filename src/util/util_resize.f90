!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_resize
   use swiftest
contains

   module subroutine util_resize_arr_char_string(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. nnew = 0 will deallocate.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                     intent(in)    :: nnew !! New size
      ! Internals
      character(len=STRMAX), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(nnew))
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = ""
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = ""
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_char_string


   module subroutine util_resize_arr_DP(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of double precision type. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                        intent(in)    :: nnew !! New size
      ! Internals
      real(DP), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      real(DP), parameter :: init_val = 0.0_DP

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(nnew))
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = init_val
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = init_val
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_DP


   module subroutine util_resize_arr_DPvec(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of double precision vectors of size (NDIM, n). Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                          intent(in)    :: nnew !! New size
      ! Internals
      real(DP), dimension(:,:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      real(DP), dimension(NDIM), parameter :: init_val = 0.0_DP
      integer(I4B) :: i

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr, dim=2)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(NDIM, nnew))
      if (nold > 0) then
         if (nnew > nold) then
            tmp(:,1:nold) = arr(:,1:nold)
            do i = nold+1, nnew
               tmp(:,i) = init_val(:)
            end do
         else
            tmp(:,1:nnew) = arr(:,1:nnew)
         end if
      else
         do i = 1, nnew
            tmp(:, i) = init_val(:)
         end do
      end if
      call move_alloc(tmp, arr)

      return

      return
   end subroutine util_resize_arr_DPvec


   module subroutine util_resize_arr_I4B(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of integer type. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                            intent(in)    :: nnew !! New size
      ! Internals
      integer(I4B), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      integer(I4B), parameter :: init_val = -1

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(nnew))
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = init_val
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = init_val
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_I4B


   module subroutine util_resize_arr_info(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                         intent(in)    :: nnew !! New size
      ! Internals
      type(swiftest_particle_info), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(nnew))
      if (nnew > nold) then
         call util_copy_particle_info_arr(arr(1:nold), tmp(1:nold))
      else
         call util_copy_particle_info_arr(arr(1:nnew), tmp(1:nnew))
      end if

      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_info


   module subroutine util_resize_arr_logical(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of logical type. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                       intent(in)    :: nnew !! New size
      ! Internals
      logical, dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      logical, parameter :: init_val = .false.

      if (nnew < 0) return

      if (nnew == 0) then
         if (allocated(arr)) deallocate(arr)
         return
      end if
      
      if (allocated(arr)) then
         nold = size(arr)
      else
         nold = 0
      end if

      if (nnew == nold) return
      
      allocate(tmp(nnew))
      if (nold > 0) then
         if (nnew > nold) then
            tmp(1:nold) = arr(1:nold)
            tmp(nold+1:nnew) = init_val
         else
            tmp(1:nnew) = arr(1:nnew)
         end if
      else
         tmp(1:nnew) = init_val
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_logical


   module subroutine util_resize_body(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self  !! Swiftest body object
      integer(I4B),         intent(in)    :: nnew  !! New size neded

      call util_resize(self%info, nnew)
      call util_resize(self%id, nnew)
      call util_resize(self%status, nnew)
      call util_resize(self%ldiscard, nnew)
      call util_resize(self%lmask, nnew)
      call util_resize(self%mu, nnew)
      call util_resize(self%rh, nnew)
      call util_resize(self%vh, nnew)
      call util_resize(self%rb, nnew)
      call util_resize(self%vb, nnew)
      call util_resize(self%ah, nnew)
      call util_resize(self%aobl, nnew)
      call util_resize(self%atide, nnew)
      call util_resize(self%agr, nnew)
      call util_resize(self%ir3h, nnew)
      call util_resize(self%a, nnew)
      call util_resize(self%e, nnew)
      call util_resize(self%inc, nnew)
      call util_resize(self%capom, nnew)
      call util_resize(self%omega, nnew)
      call util_resize(self%capm, nnew)
      self%nbody = count(self%status(1:nnew) /= INACTIVE)

      return
   end subroutine util_resize_body


   module subroutine util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest massive body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self  !! Swiftest massive body object
      integer(I4B),       intent(in)    :: nnew  !! New size neded

      call util_resize_body(self, nnew)

      call util_resize(self%mass, nnew)
      call util_resize(self%Gmass, nnew)
      call util_resize(self%rhill, nnew)
      call util_resize(self%renc, nnew)
      call util_resize(self%radius, nnew)
      call util_resize(self%rbeg, nnew)
      call util_resize(self%xend, nnew)
      call util_resize(self%vbeg, nnew)
      call util_resize(self%density, nnew)
      call util_resize(self%Ip, nnew)
      call util_resize(self%rot, nnew)
      call util_resize(self%k2, nnew)
      call util_resize(self%Q, nnew)
      call util_resize(self%tlag, nnew)

      if (allocated(self%k_plpl)) deallocate(self%k_plpl)

      return
   end subroutine util_resize_pl


   module subroutine util_resize_tp(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest test particle against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self  !! Swiftest test particle object
      integer(I4B),       intent(in)    :: nnew  !! New size neded

      call util_resize_body(self, nnew)

      call util_resize(self%isperi, nnew)
      call util_resize(self%peri, nnew)
      call util_resize(self%atp, nnew)

      return
   end subroutine util_resize_tp


end submodule s_util_resize