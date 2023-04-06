!! Copyright 2023 - David Minton
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_coarray
    use coarray
contains

    module subroutine swiftest_coarray_coclone_body(self)
        !! author: David A. Minton
         !!
         !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(swiftest_body),intent(inout),codimension[*]  :: self  !! Swiftest body object

        call coclone(self%lfirst)
        call coclone(self%nbody)
        call coclone(self%id)
        call coclone(self%info)
        call coclone(self%lmask)
        call coclone(self%status)
        call coclone(self%ldiscard)
        call coclone(self%lcollision)
        call coclone(self%lencounter)
        call coclone(self%mu)
        call coclone(self%rh)
        call coclone(self%vh)
        call coclone(self%rb)
        call coclone(self%vb)
        call coclone(self%ah)
        call coclone(self%aobl)
        call coclone(self%agr)
        call coclone(self%atide)
        call coclone(self%ir3h)
        call coclone(self%isperi)
        call coclone(self%peri)
        call coclone(self%atp)
        call coclone(self%a)
        call coclone(self%e)
        call coclone(self%inc)
        call coclone(self%capom)
        call coclone(self%omega)
        call coclone(self%capm)

        return
    end subroutine swiftest_coarray_coclone_body


    module subroutine swiftest_coarray_component_copy_info(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! swiftest_particle_info scalar version
        implicit none
        ! Arguments
        type(swiftest_particle_info), intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        type(swiftest_particle_info),save :: tmp[*]
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
    
    end subroutine swiftest_coarray_component_copy_info


    module subroutine swiftest_coarray_component_copy_info_arr1D(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! swiftest_particle_info 1D allocatable array version
        implicit none
        ! Arguments
        type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        type(swiftest_particle_info), dimension(:), codimension[:], allocatable :: tmp
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
     end subroutine swiftest_coarray_component_copy_info_arr1D


    module subroutine swiftest_coarray_collect_system(nbody_system)
       !! author: David A. Minton
       !!
       !! Collects all the test particles from other images into the image #1 test particle system
       implicit none
       ! Arguments
       class(swiftest_nbody_system), intent(inout) :: nbody_system[*]
       ! Internals
       integer(I4B) :: i,j
       integer(I4B), dimension(num_images()) :: ntp
       class(swiftest_tp), allocatable :: tp_img
 
       ! ntp(this_image()) = nbody_system%tp%nbody
       ! sync all
       ! if (this_image() == 1) then
       !    write(*,*) "Collecting test particles"
       !    write(*,*) "Image ",1," ntp: ",ntp(1)
       !    do i = 2, num_images()
       !       write(*,*) "Image ",i," ntp: ",ntp(i)
       !       allocate(tp_img, source=nbody_system[i]%tp)
       !       call nbody_system%tp%append(tp_img,lsource_mask=[(.true., j = 1, ntp(i))])
       !       deallocate(tp_img)
       !    end do
       !    write(*,*) "Total test particles: ",nbody_system%tp%nbody
       ! end if
 
       return
    end subroutine swiftest_coarray_collect_system
 
 
    module subroutine swiftest_coarray_distribute_system(nbody_system)
       !! author: David A. Minton
       !!
       !! Distributes test particles from image #1 out to all images.
       implicit none
       ! Arguments
       class(swiftest_nbody_system), intent(inout) :: nbody_system[*]
       ! Internals
       integer(I4B) :: i, istart, iend, ntot, num_per_image, ncopy
       class(swiftest_tp), allocatable :: tp_orig
       logical, dimension(:), allocatable :: lspill_list
       integer(I4B), codimension[*],save :: ntp
       class(swiftest_nbody_system), allocatable :: tmp_system
       class(swiftest_tp), allocatable :: tp
 
       ! ntp = nbody_system%tp%nbody
       ! sync all
 
       ! ntot = ntp[1]
       ! if (ntot == 0) return
 
       ! allocate(tp, mold=nbody_system%tp)
 
       ! write(*,*) "Image ",this_image(), "Distributing ",ntot
       ! allocate(lspill_list(ntot))
       ! num_per_image = ntot / num_images()
       ! istart = (this_image() - 1) * num_per_image + 1
       ! if (this_image() == num_images()) then
       !    iend = ntot
       ! else
       !    iend = this_image() * num_per_image
       ! end if
 
       ! if (this_image() == 1) then
       !    lspill_list(:) = .true.
       !    lspill_list(istart:iend) = .false.
       !    call nbody_system%tp%spill(tp,lspill_list(:), ldestructive=.true.)
       ! else
       !    lspill_list(:) = .false.
       !    lspill_list(istart:iend) = .true.
       !    tp%nbody = ntot
       !    call nbody_system%tp%spill(tp,lspill_list(:), ldestructive=.true.)
       ! end if
 
       ! write(*,*) "Image ",this_image(), "ntp: ",nbody_system%tp%nbody
 
       return
    end subroutine swiftest_coarray_distribute_system

end submodule s_swiftest_coarray