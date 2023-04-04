!! Copyright 2023 - David Minton
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module coarray
    !! author: David A. Minton
    !!
    !! Utilities that are used for coarray test particles
    !!
    use globals
    implicit none
    public

    interface cocopy
        module procedure coarray_component_copy_char
        module procedure coarray_component_copy_DP
        module procedure coarray_component_copy_DP_arr1D
        module procedure coarray_component_copy_I4B
        module procedure coarray_component_copy_I4B_arr1D
        module procedure coarray_component_copy_I8B
        module procedure coarray_component_copy_lgt
        module procedure coarray_component_copy_QP
    end interface

    contains

    subroutine coarray_component_copy_char(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! character scalar version
        implicit none
        ! Arguments
        character(len=*), intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        character(len=STRMAX),save :: tmp[*]
        integer(I4B) :: img, si
 
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
        end if
        sync all
        var = tmp[si]
  
        return
     end subroutine coarray_component_copy_char


    subroutine coarray_component_copy_DP(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! real(DP) scalar version
        implicit none
        ! Arguments
        real(DP), intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        real(DP),save :: tmp[*]
        integer(I4B) :: img, si
 
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
        end if
        sync all
        var = tmp[si]
  
        return
     end subroutine coarray_component_copy_DP 


     subroutine coarray_component_copy_DP_arr1D(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! real(DP) allocatable array version
        implicit none
        ! Arguments
        real(DP), dimension(:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        real(DP), dimension(:), codimension[:], allocatable :: tmp
        integer(I4B) :: img, si
        integer(I4B), save :: n[*]
 
        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        n = size(var)
        sync all
        allocate(tmp(n[si])[*])
        if (this_image() == si) then
            do img = 1, num_images()
                tmp(:)[img] = var(:)
            end do
        end if
        sync all
        if (this_image() /= si) then
            if (allocated(var)) deallocate(var)
            allocate(var, source=tmp)
        end if
  
        return
     end subroutine coarray_component_copy_DP_arr1D


     subroutine coarray_component_copy_DP_arr2D(var,src_img)
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
        integer(I4B), save :: n1[*], n2[*]
 
        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        n1 = size(var,dim=1)
        n2 = size(var,dim=2)
        sync all
        allocate(tmp(n1[si],n2[si])[*])
        if (this_image() == si) then
            do img = 1, num_images()
                tmp(:,:)[img] = var(:,:)
            end do
        end if
        sync all
        if (this_image() /= si) then
            if (allocated(var)) deallocate(var)
            allocate(var, source=tmp)
        end if
  
        return
     end subroutine coarray_component_copy_DP_arr2D


     subroutine coarray_component_copy_I4B(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! integer(I4B) scalar version
        implicit none
        ! Arguments
        integer(I4B), intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        integer(I4B),save :: tmp[*]
        integer(I4B) :: img, si
 
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
        end if
        sync all
        var = tmp[si]
  
        return
     end subroutine coarray_component_copy_I4B
 

    subroutine coarray_component_copy_I4B_arr1D(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! integer(I4B) allocatable array version
        implicit none
        ! Arguments
        integer(I4B), dimension(:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        integer(I4B), dimension(:), codimension[:], allocatable :: tmp
        integer(I4B) :: img, si
        integer(I4B), save :: n[*]
 
        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        n = size(var)
        sync all
        allocate(tmp(n[si])[*])
        if (this_image() == si) then
            do img = 1, num_images()
                tmp(:)[img] = var 
            end do
        end if
        sync all
        if (allocated(var)) deallocate(var)
        allocate(var, source=tmp)
  
        return
     end subroutine coarray_component_copy_I4B_arr1D


     subroutine coarray_component_copy_I8B(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! integer(I4B) scalar version
        implicit none
        ! Arguments
        integer(I8B), intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        integer(I8B),save :: tmp[*]
        integer(I4B) :: img, si
 
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
        end if
        sync all
        var = tmp[si]
  
        return
     end subroutine coarray_component_copy_I8B


     subroutine coarray_component_copy_lgt(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! logical scalar version
        implicit none
        ! Arguments
        logical, intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        logical,save :: tmp[*]
        integer(I4B) :: img, si
 
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
        end if
        sync all
        var = tmp[si]
  
        return
     end subroutine coarray_component_copy_lgt


     subroutine coarray_component_copy_QP(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! real(DP) scalar version
        implicit none
        ! Arguments
        real(QP), intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        real(QP),save :: tmp[*]
        integer(I4B) :: img, si
 
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
        end if
        sync all
        var = tmp[si]
  
        return
     end subroutine coarray_component_copy_QP 

end module coarray