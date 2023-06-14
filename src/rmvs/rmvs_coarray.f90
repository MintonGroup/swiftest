!! Copyright 2023 - David Minton
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (rmvs) s_rmvs_coarray
use coarray
use swiftest
use whm
contains

    module subroutine rmvs_coarray_coclone_cb(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(rmvs_cb),intent(inout),codimension[*]  :: self  !! RMVS pl object

        call coclone(self%outer)
        call coclone(self%inner)
        call coclone(self%lplanetocentric)

        call swiftest_coarray_coclone_cb(self)

        return
    end subroutine rmvs_coarray_coclone_cb


    module subroutine rmvs_coarray_coclone_interp(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(rmvs_interp),intent(inout),codimension[*]  :: self  !! RMVS pl object

        call coclone(self%x)
        call coclone(self%v) 
        call coclone(self%aobl)
        call coclone(self%atide)

        return
    end subroutine rmvs_coarray_coclone_interp


    module subroutine rmvs_coarray_coclone_pl(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(rmvs_pl),intent(inout),codimension[*]  :: self  !! RMVS pl object

        call coclone(self%nenc)
        call coclone(self%tpenc1P)
        call coclone(self%plind)
        call coclone(self%outer)
        call coclone(self%lplanetocentric)

        call whm_coarray_coclone_pl(self)

        return
    end subroutine rmvs_coarray_coclone_pl


    module subroutine rmvs_coarray_coclone_system(self)
        !! author: David A. Minton
         !!
         !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(rmvs_nbody_system),intent(inout),codimension[*]  :: self  !! Swiftest body object
        ! Internals
        integer(I4B) :: i, img

        call coclone(self%lplanetocentric)
        call coclone(self%rts)
        call coclone(self%vbeg)

        call swiftest_coarray_coclone_system(self)

        return
    end subroutine rmvs_coarray_coclone_system


    module subroutine rmvs_coarray_coclone_tp(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(rmvs_tp),intent(inout),codimension[*]  :: self  !! RMVS pl object

        call coclone(self%lperi)
        call coclone(self%plperP) 
        call coclone(self%plencP)
        call coclone(self%index)
        call coclone(self%ipleP)
        call coclone(self%lplanetocentric)

        call swiftest_coarray_coclone_tp(self)

        return
    end subroutine rmvs_coarray_coclone_tp


    module subroutine rmvs_coarray_component_clone_interp_arr1D(var,src_img)
        implicit none
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! swiftest_particle_info scalar version
        ! Arguments
        type(rmvs_interp), dimension(:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        type(rmvs_interp), dimension(:), codimension[:], allocatable :: tmp
        integer(I4B) :: i,img, si
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
        do i = 1, n[si]
            call tmp(i)%coclone()
        end do
        if (this_image() /= si) then
            if (allocated(var)) deallocate(var)
            allocate(var, source=tmp)
        end if

        return
    end subroutine rmvs_coarray_component_clone_interp_arr1D


    module subroutine rmvs_coarray_cocollect_tp(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(rmvs_tp),intent(inout),codimension[*]  :: self  !! RMVS pl object

        call cocollect(self%lperi)
        call cocollect(self%plperP) 
        call cocollect(self%plencP)

        call swiftest_coarray_cocollect_tp(self)

        return
    end subroutine rmvs_coarray_cocollect_tp

end submodule s_rmvs_coarray