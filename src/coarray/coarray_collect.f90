! Copyright 2024 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (coarray) s_coarray_collect
   use swiftest
contains

   module subroutine coarray_component_collect_DP_arr1D(var,dest_img)
      !! author: David A. Minton
      !!
      !! Collects components of a coarray derived type from all images and combines them into destination image component. 
      !! The default destination image is 1
      !! real(DP) 1D allocatable array version
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      real(DP), dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: i,img, ti, di, istart, iend
      integer(I4B), allocatable :: n[:],nmax[:]
      logical, allocatable :: isalloc[:]
      real(DP), dimension(:), allocatable :: vari1

      allocate(isalloc[*])
      allocate(n[*])
      allocate(nmax[*])

      if (present(dest_img)) then
         di = dest_img
      else
         di = 1
      end if

      isalloc = allocated(var)
      if (isalloc) then
         n = size(var)
      else
         n = 0
      end if

      nmax = n
      call co_max(nmax)

      allocate(tmp(nmax)[*])
      if (isalloc) tmp(1:n) = var(1:n)

      if (this_image() == di) then
         do img = 1, num_images()
            if (img /= di) then
               if (allocated(vari1)) deallocate(vari1)
               allocate(vari1, source=tmp(1:n[img])[img])
               call util_append(var, vari1)
               n = n + n[img]
            end if
         end do
         sync images(*)
      else
         sync images(di)
         if (allocated(var)) deallocate(var)
      end if

      deallocate(isalloc,n,tmp)

      return
   end subroutine coarray_component_collect_DP_arr1D


   module subroutine coarray_component_collect_DP_vec2D(var,dest_img)
      !! author: David A. Minton
      !!
      !! Collects components of a coarray derived type from all images and combines them into destination image component . The default destination image is 1
      !! real(DP) 2D allocatable array version
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      real(DP), dimension(:,:), codimension[:], allocatable :: tmp
      integer(I4B) :: i, img, ti, di, istart, iend
      integer(I4B), allocatable :: n[:], nmax[:]
      logical, allocatable :: isalloc[:]
      real(DP), dimension(:,:), allocatable :: vari1

      allocate(n[*])
      allocate(nmax[*])
      allocate(isalloc[*])

      if (present(dest_img)) then
         di = dest_img
      else
         di = 1
      end if

      isalloc = allocated(var)
      if (isalloc) then
         n = size(var,dim=2)
      else
         n = 0
      end if
      nmax = n
      call co_max(nmax)

      allocate(tmp(NDIM,nmax)[*])
      if (isalloc) tmp(:,1:n) = var(:,1:n)

      if (this_image() == di) then
         do img = 1, num_images()
            if ((img /= di).and.(isalloc[img])) then
               if (allocated(vari1)) deallocate(vari1)
               allocate(vari1,source=tmp(:,1:n[img])[img])
               call util_append(var, vari1)
               n = n + n[img]
            end if
         end do
         sync images(*)
      else
         sync images(di)
         if (allocated(var)) deallocate(var)
      end if

      deallocate(isalloc,n,nmax,tmp)

      return
   end subroutine coarray_component_collect_DP_vec2D


   module subroutine coarray_component_collect_I4B_arr1D(var,dest_img)
      !! author: David A. Minton
      !!
      !! Collects components of a coarray derived type from all images and combines them into destination image component. 
      !! The default destination image is 1
      !! integer(I4B) 1D allocatable array version
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      integer(I4B), dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: i,img, ti, di, istart, iend
      integer(I4B), allocatable, codimension[:] :: n, nmax
      logical, allocatable, codimension[:] :: isalloc[:]
      integer(I4B), dimension(:), allocatable :: vari1

      allocate(isalloc[*])
      allocate(n[*])
      allocate(nmax[*])

      if (present(dest_img)) then
         di = dest_img
      else
         di = 1
      end if

      isalloc = allocated(var)
      if (isalloc) then
         n = size(var)
      else
         n = 0
      end if
      nmax = n
      sync all
      call co_max(nmax)

      allocate(tmp(nmax)[*])
      if (isalloc) tmp(1:n) = var(1:n)

      if (this_image() == di) then
         do img = 1, num_images()
            if (img /= di) then
               if (allocated(vari1)) deallocate(vari1)
               allocate(vari1, source=tmp(1:n[img])[img])
               call util_append(var, vari1)
               n = n + n[img]
            end if
         end do
         sync images(*)
      else
         sync images(di)
         if (allocated(var)) deallocate(var)
      end if

      deallocate(isalloc,n,tmp)

      return
   end subroutine coarray_component_collect_I4B_arr1D


   module subroutine coarray_component_collect_lgt_arr1D(var,dest_img)
      !! author: David A. Minton
      !!
      !! Collects components of a coarray derived type from all images and combines them into destination image component. 
      !! The default destination image is 1
      !! logical 1D allocatable array version
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      logical, dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: i,img, ti, di, ntot, istart, iend
      integer(I4B), allocatable :: n[:], nmax[:]
      logical, allocatable :: isalloc[:]
      logical, dimension(:), allocatable :: vari1

      allocate(isalloc[*])
      allocate(n[*])
      allocate(nmax[*])

      if (present(dest_img)) then
         di = dest_img
      else
         di = 1
      end if

      isalloc = allocated(var)
      if (isalloc) then
         n = size(var)
      else
         n = 0
      end if
      nmax = n
      call co_max(nmax)
      allocate(tmp(nmax)[*])
      if (isalloc) tmp(1:n) = var(1:n)

      if (this_image() == di) then
         do img = 1, num_images()
            if (img /= di) then
               if (allocated(vari1)) deallocate(vari1)
               allocate(vari1, source=tmp(1:n[img])[img])
               call util_append(var, vari1)
               n = n + n[img]
            end if
         end do
         sync images(*)
      else
         sync images(di)
         if (allocated(var)) deallocate(var)
      end if

      deallocate(isalloc,n,tmp)

      return
   end subroutine coarray_component_collect_lgt_arr1D

end submodule s_coarray_collect