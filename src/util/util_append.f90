!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_append
   use swiftest
contains

   module subroutine util_append_arr_char_string(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of character string type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      character(len=STRMAX), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                                     intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,               dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_char_string


   module subroutine util_append_arr_DP(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                        intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,  dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_DP


   module subroutine util_append_arr_DPvec(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision vector type of size (NDIM, n) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:,:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                          intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,  dimension(:),                intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(NDIM,nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(1, nold + 1:nold + nnew) = pack(source(1,1:nsrc), lsource_mask(1:nsrc))
      arr(2, nold + 1:nold + nnew) = pack(source(2,1:nsrc), lsource_mask(1:nsrc))
      arr(3, nold + 1:nold + nnew) = pack(source(3,1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_DPvec


   module subroutine util_append_arr_I4B(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of integer(I4B) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      integer(I4B), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                            intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,      dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_I4B


   module subroutine util_append_arr_info(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of particle information type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      type(swiftest_particle_info), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                                            intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,                       dimension(:),             intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew, i
      integer(I4B), dimension(:), allocatable :: idx

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      allocate(idx(nnew))

      idx = pack([(i, i = 1, nsrc)], lsource_mask(1:nsrc))

      call util_copy_particle_info_arr(source(1:nsrc), arr(nold+1:nold+nnew), idx)

      return
   end subroutine util_append_arr_info


   module subroutine util_append_arr_logical(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of logical type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      logical, dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                       intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical, dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_logical


   module subroutine util_append_body(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),  intent(inout) :: self         !! Swiftest body object
      class(swiftest_body),  intent(in)    :: source       !! Source object to append
      logical, dimension(:), intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nold, nsrc, nnew

      nold = self%nbody
      nsrc = source%nbody
      nnew = count(lsource_mask(1:nsrc))

      call util_append(self%info, source%info, nold, nsrc, lsource_mask)
      call util_append(self%id, source%id, nold, nsrc, lsource_mask)
      call util_append(self%status, source%status, nold, nsrc, lsource_mask)
      call util_append(self%ldiscard, source%ldiscard, nold, nsrc, lsource_mask)
      call util_append(self%lmask, source%lmask, nold, nsrc, lsource_mask)
      call util_append(self%mu, source%mu, nold, nsrc, lsource_mask)
      call util_append(self%xh, source%xh, nold, nsrc, lsource_mask)
      call util_append(self%vh, source%vh, nold, nsrc, lsource_mask)
      call util_append(self%xb, source%xb, nold, nsrc, lsource_mask)
      call util_append(self%vb, source%vb, nold, nsrc, lsource_mask)
      call util_append(self%ah, source%ah, nold, nsrc, lsource_mask)
      call util_append(self%aobl, source%aobl, nold, nsrc, lsource_mask)
      call util_append(self%atide, source%atide, nold, nsrc, lsource_mask)
      call util_append(self%agr, source%agr, nold, nsrc, lsource_mask)
      call util_append(self%ir3h, source%ir3h, nold, nsrc, lsource_mask)
      call util_append(self%a, source%a, nold, nsrc, lsource_mask)
      call util_append(self%e, source%e, nold, nsrc, lsource_mask)
      call util_append(self%inc, source%inc, nold, nsrc, lsource_mask)
      call util_append(self%capom, source%capom, nold, nsrc, lsource_mask)
      call util_append(self%omega, source%omega, nold, nsrc, lsource_mask)
      call util_append(self%capm, source%capm, nold, nsrc, lsource_mask)

      self%nbody = nold + nnew

      return
   end subroutine util_append_body


   module subroutine util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_pl),              intent(inout) :: self         !! Swiftest massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (swiftest_pl)
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%mass, source%mass, nold, nsrc, lsource_mask)
            call util_append(self%Gmass, source%Gmass, nold, nsrc, lsource_mask)
            call util_append(self%rhill, source%rhill, nold, nsrc, lsource_mask)
            call util_append(self%renc, source%renc, nold, nsrc, lsource_mask)
            call util_append(self%radius, source%radius, nold, nsrc, lsource_mask)
            call util_append(self%xbeg, source%xbeg, nold, nsrc, lsource_mask)
            call util_append(self%xend, source%xend, nold, nsrc, lsource_mask)
            call util_append(self%vbeg, source%vbeg, nold, nsrc, lsource_mask)
            call util_append(self%density, source%density, nold, nsrc, lsource_mask)
            call util_append(self%Ip, source%Ip, nold, nsrc, lsource_mask)
            call util_append(self%rot, source%rot, nold, nsrc, lsource_mask)
            call util_append(self%k2, source%k2, nold, nsrc, lsource_mask)
            call util_append(self%Q, source%Q, nold, nsrc, lsource_mask)
            call util_append(self%tlag, source%tlag, nold, nsrc, lsource_mask)

            if (allocated(self%k_plpl)) deallocate(self%k_plpl)

            call util_append_body(self, source, lsource_mask)
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_pl or its descendents"
         call util_exit(FAILURE)
      end select
   
      return
   end subroutine util_append_pl


   module subroutine util_append_tp(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_tp),              intent(inout) :: self         !! Swiftest test particle object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (swiftest_tp)
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%isperi, source%isperi, nold, nsrc, lsource_mask)
            call util_append(self%peri, source%peri, nold, nsrc, lsource_mask)
            call util_append(self%atp, source%atp, nold, nsrc, lsource_mask)

            call util_append_body(self, source, lsource_mask)
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_tp or its descendents"
         call util_exit(FAILURE)
      end select

      return
   end subroutine util_append_tp

end submodule s_util_append