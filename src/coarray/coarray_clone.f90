!! Copyright 2023 - David Minton
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (coarray) s_coarray_clone
   use swiftest
contains

   module subroutine coarray_component_clone_char(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! character scalar version
      implicit none
      ! Arguments
      character(len=*), intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      character(len=STRMAX),allocatable :: tmp[:]
      integer(I4B) :: img, si

      allocate(tmp[*])
      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      sync all
      if (this_image() == si) then
         do img = 1, num_images()
         tmp[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         var = tmp[si]
      end if
      deallocate(tmp)

      return
   end subroutine coarray_component_clone_char


   module subroutine coarray_component_clone_DP(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! real(DP) scalar version
      implicit none
      ! Arguments
      real(DP), intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      real(DP),allocatable :: tmp[:]
      integer(I4B) :: img, si

      allocate(tmp[*])

      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      sync all
      if (this_image() == si) then
         do img = 1, num_images()
            tmp[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         var = tmp[si]
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_clone_DP 


   module subroutine coarray_component_clone_DP_arr1D(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! real(DP) 1D allocatable array version
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      real(DP), dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: img, si
      integer(I4B), allocatable :: n[:]
      logical, allocatable :: isalloc[:]

      allocate(isalloc[*])
      allocate(n[*])

      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      isalloc = allocated(var)
      if (isalloc) n = size(var)
      sync all 
      if (.not. isalloc[si]) return

      allocate(tmp(n[si])[*])
      if (this_image() == si) then
         do img = 1, num_images()
            tmp(:)[img] = var 
         end do
         sync images(*)
      else
         sync images(si) 
         if (allocated(var)) deallocate(var)
         allocate(var, source=tmp)
      end if
      deallocate(isalloc,n,tmp)

      return
   end subroutine coarray_component_clone_DP_arr1D


   module subroutine coarray_component_clone_DP_arr2D(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! real(DP) 2D allocatable array version
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      real(DP), dimension(:,:), codimension[:], allocatable :: tmp
      integer(I4B) :: img, si
      integer(I4B), allocatable :: n1[:], n2[:]
      logical, allocatable :: isalloc[:]

      allocate(n1[*])
      allocate(n2[*])
      allocate(isalloc[*])

      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      isalloc = allocated(var)
      if (isalloc) then
         n1 = size(var,dim=1)
         n2 = size(var,dim=2)
      end if
      sync all 
      if (.not. isalloc[si]) return

      allocate(tmp(n1[si],n2[si])[*])
      if (this_image() == si) then
         do img = 1, num_images()
            tmp(:,:)[img] = var(:,:)
         end do
         sync images(*)
      else
         sync images(si)
         if (allocated(var)) deallocate(var)
         allocate(var, source=tmp)
      end if

      deallocate(isalloc,n1,n2,tmp)

      return
   end subroutine coarray_component_clone_DP_arr2D


   module subroutine coarray_component_clone_DP_vec1D(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! real(DP) 1D (NDIM) array version
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      real(DP), dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: img, si
      integer(I4B), allocatable :: n[:]

      allocate(n[*])
      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      allocate(tmp(NDIM)[*])
      if (this_image() == si) then
         do img = 1, num_images()
            tmp(:)[img] = var(:)
         end do
         sync images(*)
      else
         sync images(si) 
         var(:) = tmp(:)[si]
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_clone_DP_vec1D


   module subroutine coarray_component_clone_DP_vec2D(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! real(DP) 1D allocatable array version
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      real(DP), dimension(:,:), codimension[:], allocatable :: tmp
      integer(I4B) :: img, si
      integer(I4B), allocatable :: n[:]
      logical, allocatable :: isalloc[:]

      allocate(isalloc[*])
      allocate(n[*])

      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      isalloc = allocated(var)
      if (isalloc) n = size(var,dim=2)
      sync all 
      if (.not. isalloc[si]) return

      allocate(tmp(NDIM,n[si])[*])
      if (this_image() == si) then
         do img = 1, num_images()
            tmp(:,:)[img] = var(:,:)
         end do
         sync images(*)
      else
         sync images(si) 
         if (allocated(var)) deallocate(var)
         allocate(var, source=tmp)
      end if

      deallocate(isalloc,n,tmp)

      return
   end subroutine coarray_component_clone_DP_vec2D


   module subroutine coarray_component_clone_I4B(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! integer(I4B) scalar version
      implicit none
      ! Arguments
      integer(I4B), intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      integer(I4B),allocatable :: tmp[:]
      integer(I4B) :: img, si

      allocate(tmp[*])

      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      if (this_image() == si) then
         do img = 1, num_images()
            tmp[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         var = tmp[si]
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_clone_I4B


   module subroutine coarray_component_clone_I4B_arr1D(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! integer(I4B) 1D allocatable array version
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      integer(I4B), dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: img, si
      integer(I4B), allocatable :: n[:]
      logical, allocatable :: isalloc[:]

      allocate(isalloc[*])
      allocate(n[*])
      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      isalloc = allocated(var)
      if (isalloc) n = size(var)
      sync all 
      if (.not. isalloc[si]) return

      allocate(tmp(n[si])[*])
      if (this_image() == si) then
         do img = 1, num_images()
            tmp(:)[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         if (allocated(var)) deallocate(var)
         allocate(var, source=tmp)
      end if

      deallocate(isalloc,n,tmp)

      return
   end subroutine coarray_component_clone_I4B_arr1D


   module subroutine coarray_component_clone_I8B(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! integer(I4B) scalar version
      implicit none
      ! Arguments
      integer(I8B), intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      integer(I8B),allocatable :: tmp[:]
      integer(I4B) :: img, si

      allocate(tmp[*])
      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      if (this_image() == si) then
         do img = 1, num_images()
            tmp[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         var = tmp[si]
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_clone_I8B


   module subroutine coarray_component_clone_lgt(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! logical scalar version
      implicit none
      ! Arguments
      logical, intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      logical,allocatable :: tmp[:]
      integer(I4B) :: img, si

      allocate(tmp[*])

      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      if (this_image() == si) then
         do img = 1, num_images()
            tmp[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         var = tmp[si]
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_clone_lgt


   module subroutine coarray_component_clone_lgt_arr1D(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! logical 1D allocatable array version
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      logical, dimension(:), codimension[:], allocatable :: tmp
      integer(I4B) :: img, si
      integer(I4B), allocatable :: n[:]
      logical, allocatable :: isalloc[:]

      allocate(isalloc[*])
      allocate(n[*])
      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      isalloc = allocated(var)
      if (isalloc) n = size(var)
      sync all 
      if (.not. isalloc[si]) return

      allocate(tmp(n[si])[*])
      if (this_image() == si) then
         do img = 1, num_images()
            tmp(:)[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         if (allocated(var)) deallocate(var)
         allocate(var, source=tmp)
      end if

      deallocate(isalloc,n,tmp)

      return
   end subroutine coarray_component_clone_lgt_arr1D


#ifdef QUADPREC
   module subroutine coarray_component_clone_QP(var,src_img)
      !! author: David A. Minton
      !!
      !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
      !! real(DP) scalar version
      implicit none
      ! Arguments
      real(QP), intent(inout) :: var
      integer(I4B), intent(in),optional :: src_img
      ! Internals
      real(QP),allocatable :: tmp[:]
      integer(I4B) :: img, si

      allocate(tmp[*])
      if (present(src_img)) then
         si = src_img
      else
         si = 1
      end if

      if (this_image() == si) then
         do img = 1, num_images()
            tmp[img] = var 
         end do
         sync images(*)
      else
         sync images(si)
         var = tmp[si]
      end if

      deallocate(tmp)

      return
   end subroutine coarray_component_clone_QP 
#endif

end submodule s_coarray_clone