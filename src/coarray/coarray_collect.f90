!! Copyright 2023 - David Minton
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (coarray) s_coarray_collect
   use swiftest
contains

   module subroutine coarray_component_collect_DP_arr1D(var,dest_img)
      !! author: David A. Minton
      !!
      !! Collects components of a coarray derived type from all images and combines them into destination image component. The default destination image is 1
      !! real(DP) 1D allocatable array version
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      real(DP), dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: i,img, ti, di, istart, iend, nmax
      integer(I4B), allocatable :: n[:]
      logical, allocatable :: isalloc[:]
      real(DP), dimension(:), allocatable :: vari1

      allocate(isalloc[*])
      allocate(n[*])

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

      sync all

      nmax = 0
      do img = 1, num_images()
         if (n[img] > nmax) nmax = n[img]
      end do

      allocate(tmp(nmax)[*])
      if (isalloc) tmp(1:n) = var(1:n)

      if (this_image() == di) then
         do img = 1, num_images()
            if (img /= di) then
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


   module subroutine coarray_component_collect_DP_arr2D(var,dest_img)
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
      integer(I4B) :: i, img, ti, di, ntot, istart, iend, nmax
      integer(I4B), allocatable :: n1[:], n2[:]
      logical, allocatable :: isalloc[:]
      real(DP), dimension(:,:), allocatable :: vari1

      allocate(n1[*])
      allocate(n2[*])
      allocate(isalloc[*])

      if (present(dest_img)) then
         di = dest_img
      else
         di = 1
      end if

      isalloc = allocated(var)
      if (isalloc) then
         n1 = size(var,dim=1)
         n2 = size(var,dim=2)
      else
         n1 = 0
         n2 = 0
      end if
      sync all

      nmax = 0
      do img = 1, num_images()
         if (n2[img] > nmax) nmax = n2[img]
      end do

      allocate(tmp(NDIM,nmax)[*])
      if (isalloc) tmp(:,1:n2) = var(:,1:n2)

      if (this_image() == di) then
         do img = 1, num_images()
            if (img /= di) then
               allocate(vari1,source=tmp(:,1:n2[img])[img])
               call util_append(var, vari1)
               n2 = n2 + n2[img]
            end if
         end do
         sync images(*)
      else
         sync images(di)
         if (allocated(var)) deallocate(var)
      end if

      deallocate(isalloc,n1,n2,tmp)

      return
   end subroutine coarray_component_collect_DP_arr2D


   module subroutine coarray_component_collect_I4B(var,dest_img)
      !! author: David A. Minton
      !!
      !! Sums this component of a coarray derived type from all images and places the value in the destination image component. The default destination image is 1
      !! integer(I4B) version
      implicit none
      ! Arguments
      integer(I4B), intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      integer(I4B), allocatable :: tmp[:]
      integer(I4B) :: img, di

      if (present(dest_img)) then
         di = dest_img
      else
         di = 1
      end if

      allocate(tmp[*], source=var)
      
      if (this_image() == di) then
         var = 0
         do img = 1, num_images()
            var = var + tmp[img]
         end do
      else
         var = 0
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_collect_I4B


   module subroutine coarray_component_collect_I8B(var,dest_img)
      !! author: David A. Minton
      !!
      !! Sums this component of a coarray derived type from all images and places the value in the destination image component. The default destination image is 1
      !! integer(I8B) version
      implicit none
      ! Arguments
      integer(I8B), intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      integer(I8B), allocatable :: tmp[:]
      integer(I4B) :: img, di

      if (present(dest_img)) then
         di = dest_img
      else
         di = 1
      end if

      allocate(tmp[*], source=var)
      
      if (this_image() == di) then
         var = 0
         do img = 1, num_images()
            var = var + tmp[img]
         end do
      else
         var = 0
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_collect_I8B


   module subroutine coarray_component_collect_I4B_arr1D(var,dest_img)
      !! author: David A. Minton
      !!
      !! Collects components of a coarray derived type from all images and combines them into destination image component . The default destination image is 1
      !! integer(I4B) 1D allocatable array version
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      integer(I4B), dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: i,img, ti, di, istart, iend, nmax
      integer(I4B), allocatable :: n[:]
      logical, allocatable :: isalloc[:]
      integer(I4B), dimension(:), allocatable :: vari1

      allocate(isalloc[*])
      allocate(n[*])

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
      sync all
      nmax = 0
      do img = 1, num_images()
         if (n[img] > nmax) nmax = n[img]
      end do

      allocate(tmp(nmax)[*])
      if (isalloc) tmp(1:n) = var(1:n)

      if (this_image() == di) then
         do img = 1, num_images()
            if (img /= di) then
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
      !! Collects components of a coarray derived type from all images and combines them into destination image component . The default destination image is 1
      !! logical 1D allocatable array version
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: dest_img
      ! Internals
      logical, dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: i,img, ti, di, ntot, istart, iend, nmax
      integer(I4B), allocatable :: n[:]
      logical, allocatable :: isalloc[:]
      logical, dimension(:), allocatable :: vari1

      allocate(isalloc[*])
      allocate(n[*])

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
      sync all
      nmax = 0
      do img = 1, num_images()
         if (n[img] > nmax) nmax = n[img]
      end do

      allocate(tmp(nmax)[*])
      if (isalloc) tmp(1:n) = var(1:n)

      if (this_image() == di) then
         do img = 1, num_images()
            if (img /= di) then
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