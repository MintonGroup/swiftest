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


    module subroutine swiftest_coarray_coclone_cb(self)
        !! author: David A. Minton
         !!
         !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(swiftest_cb),intent(inout),codimension[*]  :: self  !! Swiftest body object
        ! Internals
        integer(I4B) :: i

        call coclone(self%info)
        call coclone(self%id)
        call coclone(self%mass)
        call coclone(self%Gmass)
        call coclone(self%radius)
        call coclone(self%density)
        call coclone(self%j2rp2)
        call coclone(self%j4rp4)
        call coclone(self%k2)
        call coclone(self%Q)
        call coclone(self%tlag)
        call coclone(self%GM0)
        call coclone(self%dGM)
        call coclone(self%R0)
        call coclone(self%dR)

        call coclonevec(self%aobl)
        call coclonevec(self%atide)
        call coclonevec(self%aoblbeg)
        call coclonevec(self%aoblend)
        call coclonevec(self%atidebeg)
        call coclonevec(self%atideend)
        call coclonevec(self%rb)
        call coclonevec(self%vb)
        call coclonevec(self%agr)
        call coclonevec(self%Ip)
        call coclonevec(self%rot)
        call coclonevec(self%L0)
        call coclonevec(self%dL)

        return
    end subroutine swiftest_coarray_coclone_cb


    module subroutine swiftest_coarray_coclone_pl(self)
        !! author: David A. Minton
         !!
         !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(swiftest_pl),intent(inout),codimension[*]  :: self  !! Swiftest body object

        call coclone(self%mass)
        call coclone(self%Gmass)
        call coclone(self%rhill)
        call coclone(self%renc)
        call coclone(self%radius)
        call coclone(self%density)
        call coclone(self%rbeg)
        call coclone(self%rend)
        call coclone(self%vbeg)
        call coclone(self%Ip)
        call coclone(self%rot)
        call coclone(self%k2)
        call coclone(self%Q )
        call coclone(self%tlag)
        call coclone(self%kin)
        call coclone(self%lmtiny)
        call coclone(self%nplm)
        call coclone(self%nplplm)
        call coclone(self%nplenc)
        call coclone(self%ntpenc)

        call swiftest_coarray_coclone_body(self)

        return
    end subroutine swiftest_coarray_coclone_pl


    module subroutine swiftest_coarray_coclone_tp(self)
        !! author: David A. Minton
         !!
         !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(swiftest_tp),intent(inout),codimension[*]  :: self  !! Swiftest body object

        call coclone(self%nplenc)
        call swiftest_coarray_coclone_body(self)

        return
    end subroutine swiftest_coarray_coclone_tp


    module subroutine swiftest_coarray_coclone_system(self)
        !! author: David A. Minton
         !!
         !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(swiftest_nbody_system),intent(inout),codimension[*]  :: self  !! Swiftest body object
        ! Internals
        integer(I4B) :: i, img

        call self%cb%coclone()
        call self%pl%coclone()
        call self%tp%coclone()

        call coclone(self%maxid)
        call coclone(self%t)
        call coclone(self%GMtot)
        call coclone(self%ke_orbit)
        call coclone(self%ke_spin)
        call coclone(self%pe)
        call coclone(self%be)
        call coclone(self%te)
        call coclone(self%oblpot)
        do i = 1, NDIM
            call coclone(self%L_orbit(i))
            call coclone(self%L_spin(i))
            call coclone(self%L_total(i))
            call coclone(self%L_total_orig(i))
            call coclone(self%L_orbit_orig(i))
            call coclone(self%L_spin_orig(i))
            call coclone(self%L_escape(i))
        end do
        call coclone(self%ke_orbit_orig)
        call coclone(self%ke_spin_orig)
        call coclone(self%pe_orig)
        call coclone(self%be_orig)
        call coclone(self%te_orig)
        call coclone(self%be_cb)
        call coclone(self%E_orbit_orig)
        call coclone(self%GMtot_orig)
        call coclone(self%GMescape)
        call coclone(self%E_collisions)
        call coclone(self%E_untracked)
        call coclone(self%ke_orbit_error)
        call coclone(self%ke_spin_error)
        call coclone(self%pe_error)
        call coclone(self%be_error)
        call coclone(self%E_orbit_error)
        call coclone(self%Ecoll_error)
        call coclone(self%E_untracked_error)
        call coclone(self%te_error)
        call coclone(self%L_orbit_error)
        call coclone(self%L_spin_error)
        call coclone(self%L_escape_error)
        call coclone(self%L_total_error)
        call coclone(self%Mtot_error)
        call coclone(self%Mescape_error)
        call coclone(self%lbeg)

        return
    end subroutine swiftest_coarray_coclone_system

  
    module subroutine swiftest_coarray_component_clone_info(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! swiftest_particle_info scalar version
        implicit none
        ! Arguments
        type(swiftest_particle_info), intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        type(swiftest_particle_info),allocatable :: tmp[:]
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
    
        return
    end subroutine swiftest_coarray_component_clone_info
    
    
    module subroutine swiftest_coarray_component_clone_info_arr1D(var,src_img)
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
    
        return
    end subroutine swiftest_coarray_component_clone_info_arr1D


    module subroutine swiftest_coarray_component_clone_kin_arr1D(var,src_img)
        implicit none
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source image is 1
        !! swiftest_kinship allocatable array version
        ! Arguments
        type(swiftest_kinship), dimension(:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: src_img
        ! Internals
        type(swiftest_kinship), dimension(:), codimension[:], allocatable :: tmp
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
    end subroutine swiftest_coarray_component_clone_kin_arr1D


    module subroutine swiftest_coarray_component_collect_info_arr1D(var,dest_img)
        !! author: David A. Minton
        !!
        !! Collects components of a coarray derived type from all images and combines them into destination image component . The default destination image is 1
        !! swiftest_particle_info 1D allocatable array version
        implicit none
        ! Arguments
        type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: var
        integer(I4B), intent(in),optional :: dest_img
        ! Internals
        type(swiftest_particle_info), dimension(:), codimension[:], allocatable :: tmp
        integer(I4B) :: i,img, ti, di, ntot, istart, iend, nmax
        integer(I4B), allocatable :: n[:]
        logical, allocatable :: isalloc[:]
  
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
                 call util_append(var, tmp(1:n[img])[img])
                 n = n + n[img]
              end if
           end do
        end if
  
        return
    end subroutine swiftest_coarray_component_collect_info_arr1D

    
    module subroutine swiftest_coarray_cocollect_body(self)
        !! author: David A. Minton
        !!
        !! Collects all body object array components from all images and combines them into the image 1 body object
        implicit none
        ! Arguments
        class(swiftest_body),intent(inout), codimension[*] :: self !! Swiftest body object
        integer(I4B) :: i

        call cocollect(self%nbody)
        call cocollect(self%id)
        call cocollect(self%info)
        call cocollect(self%lmask)
        call cocollect(self%status)
        call cocollect(self%ldiscard)
        call cocollect(self%lcollision)
        call cocollect(self%lencounter)
        call cocollect(self%mu)
        call cocollect(self%rh)
        call cocollect(self%vh)
        call cocollect(self%rb)
        call cocollect(self%vb)
        call cocollect(self%ah)
        call cocollect(self%aobl)
        call cocollect(self%agr)
        call cocollect(self%atide)
        call cocollect(self%ir3h)
        call cocollect(self%isperi)
        call cocollect(self%peri)
        call cocollect(self%atp)
        call cocollect(self%a)
        call cocollect(self%e)
        call cocollect(self%inc)
        call cocollect(self%capom)
        call cocollect(self%omega)
        call cocollect(self%capm)

        return
    end subroutine swiftest_coarray_cocollect_body


    module subroutine swiftest_coarray_cocollect_tp(self)
        !! author: David A. Minton
        !!
        !! Collects all object array components from all images and combines them into the image 1 object
        implicit none
        ! Arguments
        class(swiftest_tp),intent(inout),codimension[*]  :: self  !! Swiftest body object

        call cocollect(self%npltp)
        call cocollect(self%nplenc)
        call swiftest_coarray_cocollect_body(self)

        return
    end subroutine swiftest_coarray_cocollect_tp


    module subroutine swiftest_coarray_collect_system(nbody_system, param)
        !! author: David A. Minton
        !!
        !! Collects all the test particles from other images into the image #1 test particle system
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
        ! Internals
        integer(I4B) :: i,j
        integer(I4B), codimension[*], save :: ntp
        class(swiftest_tp), allocatable, codimension[:] :: cotp
        character(len=NAMELEN) :: image_num_char

        if (.not.param%lcoarray) return

        sync all
        if (this_image() == 1) then
            write(image_num_char,*) num_images()
            write(param%display_unit,*) " Collecting test particles from " // trim(adjustl(image_num_char)) // " images."
        end if

        allocate(cotp[*], source=nbody_system%tp) 
        call cotp%cocollect()
        deallocate(nbody_system%tp)
        allocate(nbody_system%tp, source=cotp)

        deallocate(cotp)

        return
    end subroutine swiftest_coarray_collect_system
 
 
    module subroutine swiftest_coarray_distribute_system(nbody_system, param)
        !! author: David A. Minton
        !!
        !! Distributes test particles from image #1 out to all images.
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
        ! Internals
        integer(I4B) :: i, istart, iend, ntot, num_per_image, ncopy
        class(swiftest_tp), allocatable :: tp
        logical, dimension(:), allocatable :: lspill_list
        integer(I4B), codimension[*], save  :: ntp
        character(len=NAMELEN) :: image_num_char
        class(swiftest_tp), allocatable, codimension[:] :: cotp
        class(swiftest_tp), allocatable :: tmp

        if (.not.param%lcoarray) return
        sync all
        if (this_image() == 1) then
            write(image_num_char,*) num_images()
            write(param%display_unit,*) " Distributing test particles across " // trim(adjustl(image_num_char)) // " images."
        end if

        ntp = nbody_system%tp%nbody
        sync all
        ntot = ntp[1]
        if (ntot == 0) return
    
        allocate(lspill_list(ntot))
        num_per_image = ceiling(1.0_DP * ntot / num_images())
        istart = (this_image() - 1) * num_per_image + 1
        if (this_image() == num_images()) then
            iend = ntot
        else
            iend = this_image() * num_per_image
        end if
    
        lspill_list(:) = .true.
        lspill_list(istart:iend) = .false.

        allocate(cotp[*], source=nbody_system%tp)
        call cotp%coclone()
        if (this_image() /= 1) then
            deallocate(nbody_system%tp)
            allocate(nbody_system%tp, source=cotp)
        end if
        allocate(tmp, mold=nbody_system%tp)
        call nbody_system%tp%spill(tmp, lspill_list(:), ldestructive=.true.)

        deallocate(tmp, cotp)

        return
    end subroutine swiftest_coarray_distribute_system


    module subroutine swiftest_coarray_initialize_system(nbody_system, param)
        !! author: David A. Minton
        !!
        !! Distributes test particles from image #1 out to all images.
        implicit none
        ! Arguments
        class(swiftest_nbody_system), allocatable, intent(inout) :: nbody_system !! Swiftest nbody system 
        class(swiftest_parameters),                intent(inout) :: param        !! Current run configuration parameters 
        ! Internals
        class(swiftest_nbody_system), allocatable, codimension[:] :: tmp_system
        character(len=NAMELEN) :: image_num_char

        if (.not.param%lcoarray) return

        sync all
        if (this_image() == 1) then
            write(image_num_char,*) num_images()
            write(param%display_unit,*) " Cloning nbody system to " // trim(adjustl(image_num_char)) // " images."
        end if
        allocate(tmp_system[*], source=nbody_system)
        call tmp_system%coclone()
        if (this_image() /= 1) then
           if (allocated(nbody_system)) deallocate(nbody_system)
           allocate(nbody_system, source=tmp_system)
        end if
 
        return
    end subroutine swiftest_coarray_initialize_system


end submodule s_swiftest_coarray