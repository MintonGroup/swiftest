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

    interface coclone
        module procedure coarray_component_clone_char
        module procedure coarray_component_clone_DP
        module procedure coarray_component_clone_DP_arr1D
        module procedure coarray_component_clone_DP_arr2D
        module procedure coarray_component_clone_I4B
        module procedure coarray_component_clone_I4B_arr1D
        module procedure coarray_component_clone_I4B_arr2D
        module procedure coarray_component_clone_I8B
        module procedure coarray_component_clone_lgt
        module procedure coarray_component_clone_lgt_arr1D
        module procedure coarray_component_clone_QP
    end interface

    interface cocollect
        module procedure coarray_component_collect_DP_arr1D
        module procedure coarray_component_collect_DP_arr2D
        module procedure coarray_component_collect_I4B
        module procedure coarray_component_collect_I4B_arr1D
        module procedure coarray_component_collect_I4B_arr2D
        module procedure coarray_component_collect_I8B
        module procedure coarray_component_collect_lgt_arr1D
    end interface

    contains

    subroutine coarray_component_clone_char(var,src_img)
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
            sync images(*)
        else
            sync images(si)
            var = tmp[si]
        end if
  
        return
    end subroutine coarray_component_clone_char


    subroutine coarray_component_clone_DP(var,src_img)
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
            sync images(*)
        else
            sync images(si)
            var = tmp[si]
        end if

  
        return
    end subroutine coarray_component_clone_DP 


    subroutine coarray_component_clone_DP_arr1D(var,src_img)
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
        integer(I4B), save :: n[*]
        logical, save :: isalloc[*]

        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        sync all
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
  
        return
    end subroutine coarray_component_clone_DP_arr1D


    subroutine coarray_component_clone_DP_arr2D(var,src_img)
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
        logical, save :: isalloc[*]

        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        sync all
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

        return
    end subroutine coarray_component_clone_DP_arr2D


    subroutine coarray_component_clone_I4B(var,src_img)
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
            sync images(*)
        else
            sync images(si)
            var = tmp[si]
        end if
  
        return
    end subroutine coarray_component_clone_I4B


    subroutine coarray_component_clone_I4B_arr1D(var,src_img)
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
        integer(I4B), save :: n[*]
        logical, save :: isalloc[*]

        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        sync all
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

        return
    end subroutine coarray_component_clone_I4B_arr1D


    subroutine coarray_component_clone_I4B_arr2D(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! integer(I4B) 2D allocatable array version
        implicit none
        ! Arguments
        integer(I4B), dimension(:,:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        integer(I4B), dimension(:,:), codimension[:], allocatable :: tmp
        integer(I4B) :: img, si
        integer(I4B), save :: n1[*], n2[*]
        logical, save :: isalloc[*]

        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        sync all
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

        return
    end subroutine coarray_component_clone_I4B_arr2D


    subroutine coarray_component_clone_I8B(var,src_img)
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
            sync images(*)
        else
            sync images(si)
            var = tmp[si]
        end if
  
        return
    end subroutine coarray_component_clone_I8B


    subroutine coarray_component_clone_lgt(var,src_img)
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
            sync images(*)
        else
            sync images(si)
            var = tmp[si]
        end if
  
  
        return
    end subroutine coarray_component_clone_lgt


    subroutine coarray_component_clone_lgt_arr1D(var,src_img)
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
        integer(I4B), save :: n[*]
        logical, save :: isalloc[*]

        if (present(src_img)) then
            si = src_img
        else
            si = 1
        end if

        sync all
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
  
        return
    end subroutine coarray_component_clone_lgt_arr1D


    subroutine coarray_component_clone_QP(var,src_img)
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
            sync images(*)
        else
            sync images(si)
            var = tmp[si]
        end if
  
        return
    end subroutine coarray_component_clone_QP 


    subroutine coarray_component_collect_DP_arr1D(var,dest_img)
        !! author: David A. Minton
        !!
        !! Collects components of a coarray derived type from all images and combines them into destination image component . The default destination image is 1
        !! real(DP) 1D allocatable array version
        implicit none
        ! Arguments
        real(DP), dimension(:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: dest_img
        ! Internals
        real(DP), dimension(:), codimension[:], allocatable :: tmp
        integer(I4B) :: img, ti, di, ntot, istart, iend
        integer(I4B), save :: n[*]
        logical, save :: isalloc[*]

        if (present(dest_img)) then
            di = dest_img
        else
            di = 1
        end if

        sync all
        isalloc = allocated(var)
        if (isalloc) then
            n = size(var)
        else
            n = 0
        end if
        sync all
        ntot = 0
        do img = 1, num_images()
            ntot = ntot + n[img]
        end do

        allocate(tmp(ntot)[*])

        ti = this_image()

        istart = 1
        iend = n
        do img = 1, this_image() - 1 
            istart = istart + n[img]
            iend = iend + n[img]
        end do

        if (isalloc) then
            tmp(istart:iend) = var(:)
            deallocate(var)
        end if

        sync all
        if (this_image() == di) then
            allocate(var(ntot))
            istart = 1
            iend = n
            do img = 1, num_images()
                var(istart:iend) = tmp[img](istart:iend)
                istart = istart + n[img]
                iend = iend + n[img]
            end do
        end if

        return
    end subroutine coarray_component_collect_DP_arr1D


    subroutine coarray_component_collect_DP_arr2D(var,dest_img)
        !! author: David A. Minton
        !!
        !! Collects components of a coarray derived type from all images and combines them into destination image component . The default destination image is 1
        !! real(DP) 2D allocatable array version
        implicit none
        ! Arguments
        real(DP), dimension(:,:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: dest_img
        ! Internals
        integer(I4B), dimension(:,:), codimension[:], allocatable :: tmp
        integer(I4B) :: img, ti, di, ntot, istart, iend
        integer(I4B), save :: n1[*], n2[*]
        logical, save :: isalloc[*]

        if (present(dest_img)) then
            di = dest_img
        else
            di = 1
        end if

        sync all
        isalloc = allocated(var)
        if (isalloc) then
            n1 = size(var,dim=1)
            n2 = size(var,dim=2)
        else
            n1 = 0
            n2 = 0
        end if
        sync all
        ntot = 0
        do img = 1, num_images()
            ntot = ntot + n2[img]
        end do

        allocate(tmp(n1,ntot)[*])

        ti = this_image()

        istart = 1
        iend = n2
        do img = 1, this_image() - 1 
            istart = istart + n2[img]
            iend = iend + n2[img]
        end do

        if (isalloc) then
            tmp(:,istart:iend) = var(:,:)
            deallocate(var)
        end if

        sync all
        if (this_image() == di) then
            allocate(var(n1,ntot))

            istart = 1
            iend = n2
            do img = 1, num_images()
                var(:,istart:iend) = tmp[img](:,istart:iend)
                istart = istart + n2[img]
                iend = iend + n2[img]
            end do
        end if

        return
    end subroutine coarray_component_collect_DP_arr2D


    subroutine coarray_component_collect_I4B(var,dest_img)
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
        sync all

        if (this_image() == di) then
            var = 0
            do img = 1, num_images()
                var = var + tmp[img]
            end do
        else
            var = 0
        end if

        return
    end subroutine coarray_component_collect_I4B


    subroutine coarray_component_collect_I8B(var,dest_img)
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
        sync all

        if (this_image() == di) then
            var = 0
            do img = 1, num_images()
                var = var + tmp[img]
            end do
        else
            var = 0
        end if

        return
    end subroutine coarray_component_collect_I8B


    subroutine coarray_component_collect_I4B_arr1D(var,dest_img)
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
        integer(I4B) :: img, ti, di, ntot, istart, iend
        integer(I4B), save :: n[*]
        logical, save :: isalloc[*]

        if (present(dest_img)) then
            di = dest_img
        else
            di = 1
        end if

        sync all
        isalloc = allocated(var)
        if (isalloc) then
            n = size(var)
        else
            n = 0
        end if
        sync all
        ntot = 0
        do img = 1, num_images()
            ntot = ntot + n[img]
        end do

        allocate(tmp(ntot)[*])

        ti = this_image()

        istart = 1
        iend = n
        do img = 1, this_image() - 1 
            istart = istart + n[img]
            iend = iend + n[img]
        end do

        if (isalloc) then
            tmp(istart:iend) = var(:)
            deallocate(var)
        end if

        sync all
        if (this_image() == di) then
            allocate(var(ntot))
            istart = 1
            iend = n
            do img = 1, num_images()
                var(istart:iend) = tmp[img](istart:iend)
                istart = istart + n[img]
                iend = iend + n[img]
            end do
        end if

        return
    end subroutine coarray_component_collect_I4B_arr1D


    subroutine coarray_component_collect_I4B_arr2D(var,dest_img)
        !! author: David A. Minton
        !!
        !! Collects components of a coarray derived type from all images and combines them into destination image component . The default destination image is 1
        !! integer(I4B) 2D allocatable array version
        implicit none
        ! Arguments
        integer(I4B), dimension(:,:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: dest_img
        ! Internals
        integer(I4B), dimension(:,:), codimension[:], allocatable :: tmp
        integer(I4B) :: img, ti, di, ntot, istart, iend
        integer(I4B), save :: n1[*], n2[*]
        logical, save :: isalloc[*]

        if (present(dest_img)) then
            di = dest_img
        else
            di = 1
        end if

        sync all
        isalloc = allocated(var)
        if (isalloc) then
            n1 = size(var,dim=1)
            n2 = size(var,dim=2)
        else
            n1 = 0
            n2 = 0
        end if
        sync all
        ntot = 0
        do img = 1, num_images()
            ntot = ntot + n2[img]
        end do

        allocate(tmp(n1,ntot)[*])

        ti = this_image()

        istart = 1
        iend = n2
        do img = 1, this_image() - 1 
            istart = istart + n2[img]
            iend = iend + n2[img]
        end do

        if (isalloc) then
            tmp(:,istart:iend) = var(:,:)
            deallocate(var)
        end if

        sync all
        if (this_image() == di) then
            allocate(var(n1,ntot))

            istart = 1
            iend = n2
            do img = 1, num_images()
                var(:,istart:iend) = tmp[img](:,istart:iend)
                istart = istart + n2[img]
                iend = iend + n2[img]
            end do
        end if

        return
    end subroutine coarray_component_collect_I4B_arr2D


    subroutine coarray_component_collect_lgt_arr1D(var,dest_img)
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
        integer(I4B) :: img, ti, di, ntot, istart, iend
        integer(I4B), save :: n[*]
        logical, save :: isalloc[*]

        if (present(dest_img)) then
            di = dest_img
        else
            di = 1
        end if

        sync all
        isalloc = allocated(var)
        if (isalloc) then
            n = size(var)
        else
            n = 0
        end if
        sync all
        ntot = 0
        do img = 1, num_images()
            ntot = ntot + n[img]
        end do

        allocate(tmp(ntot)[*])

        ti = this_image()

        istart = 1
        iend = n
        do img = 1, this_image() - 1 
            istart = istart + n[img]
            iend = iend + n[img]
        end do

        if (isalloc) then
            tmp(istart:iend) = var(:)
            deallocate(var)
        end if

        sync all
        if (this_image() == di) then
            allocate(var(ntot))
            istart = 1
            iend = n
            do img = 1, num_images()
                var(istart:iend) = tmp[img](istart:iend)
                istart = istart + n[img]
                iend = iend + n[img]
            end do
        end if

        return
    end subroutine coarray_component_collect_lgt_arr1D


end module coarray