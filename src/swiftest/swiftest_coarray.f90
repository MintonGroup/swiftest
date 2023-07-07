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


    module subroutine swiftest_coarray_balance_system(nbody_system, param)
        !! author: David A. Minton
        !!
        !! Checks whether or not the system needs to be rebalance. Rebalancing occurs when the difference between the number of test particles between the
        !! image with the smallest and largest number of test particles is larger than the number of images
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
        ! Internals
        integer(I4B), codimension[*], save :: ntp
        integer(I4B) :: img,ntp_min, ntp_max
        character(len=NAMELEN) :: min_str, max_str, diff_str, ni_str

        ntp = nbody_system%tp%nbody
        sync all
        write(param%display_unit,*) "Checking whether test particles need to be reblanced."
        ntp_min = huge(1)
        ntp_max = 0
        do img = 1, num_images()
            if (ntp[img] < ntp_min) ntp_min = ntp[img]
            if (ntp[img] > ntp_max) ntp_max = ntp[img]
        end do
        write(min_str,*) ntp_min
        write(max_str,*) ntp_max
        write(diff_str,*) ntp_max - ntp_min
        write(ni_str,*) num_images()
        write(param%display_unit,*) "ntp_min   : " // trim(adjustl(min_str))
        write(param%display_unit,*) "ntp_max   : " // trim(adjustl(max_str))
        write(param%display_unit,*) "difference: " // trim(adjustl(diff_str))
        flush(param%display_unit)
        sync all
        if (ntp_max - ntp_min >= num_images()) then
            write(param%display_unit,*) trim(adjustl(diff_str)) // ">=" // trim(adjustl(ni_str)) // ": Rebalancing"
            flush(param%display_unit)
            call nbody_system%coarray_collect(param)
            call nbody_system%coarray_distribute(param)
            write(param%display_unit,*) "Rebalancing complete"
        else
            write(param%display_unit,*) trim(adjustl(diff_str)) // "<" // trim(adjustl(ni_str)) // ": No rebalancing needed"
        end if
        flush(param%display_unit)
        return
    end subroutine swiftest_coarray_balance_system


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

    module subroutine swiftest_coarray_coclone_kin(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(swiftest_kinship),intent(inout),codimension[*]  :: self  !! Swiftest kinship object

        call coclone(self%parent)
        call coclone(self%nchild)
        call coclone(self%child)

        return
     end subroutine swiftest_coarray_coclone_kin

    module subroutine swiftest_coarray_coclone_nc(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(swiftest_netcdf_parameters),intent(inout),codimension[*]  :: self  !! Swiftest body object

        call coclone(self%file_name)
        call coclone(self%lfile_is_open)
        call coclone(self%out_type)
        call coclone(self%id)
        call coclone(self%tslot)
        call coclone(self%max_tslot)
        call coclone(self%idvals)
        call coclone(self%idslot)
        call coclone(self%max_idslot)
        call coclone(self%str_dimname)
        call coclone(self%str_dimid)
        call coclone(self%time_dimname)
        call coclone(self%time_dimid)
        call coclone(self%time_varid)
        call coclone(self%name_dimname)
        call coclone(self%name_dimid)
        call coclone(self%name_varid)
        call coclone(self%space_dimname)
        call coclone(self%space_dimid)
        call coclone(self%space_varid)
        call coclone(self%id_varname)
        call coclone(self%id_varid)
        call coclone(self%status_varname)
        call coclone(self%status_varid)
        call coclone(self%ptype_varname)
        call coclone(self%ptype_varid)
        call coclone(self%npl_varname)
        call coclone(self%npl_varid)
        call coclone(self%ntp_varname)
        call coclone(self%ntp_varid)
        call coclone(self%nplm_varname)
        call coclone(self%nplm_varid)
        call coclone(self%a_varname)
        call coclone(self%a_varid)
        call coclone(self%e_varname)
        call coclone(self%e_varid)
        call coclone(self%inc_varname)
        call coclone(self%inc_varid)
        call coclone(self%capom_varname)
        call coclone(self%capom_varid)
        call coclone(self%omega_varname)
        call coclone(self%omega_varid)
        call coclone(self%capm_varname)
        call coclone(self%capm_varid)
        call coclone(self%varpi_varname)
        call coclone(self%varpi_varid)
        call coclone(self%lam_varname)
        call coclone(self%lam_varid)
        call coclone(self%f_varname)
        call coclone(self%f_varid)
        call coclone(self%cape_varname)
        call coclone(self%cape_varid)
        call coclone(self%rh_varname)
        call coclone(self%rh_varid)
        call coclone(self%vh_varname)
        call coclone(self%vh_varid)
        call coclone(self%gr_pseudo_vh_varname)
        call coclone(self%gr_pseudo_vh_varid)
        call coclone(self%Gmass_varname)
        call coclone(self%Gmass_varid)
        call coclone(self%mass_varname)
        call coclone(self%mass_varid)
        call coclone(self%rhill_varname)
        call coclone(self%rhill_varid)
        call coclone(self%radius_varname)
        call coclone(self%radius_varid)
        call coclone(self%Ip_varname)
        call coclone(self%Ip_varid)
        call coclone(self%rot_varname)
        call coclone(self%rot_varid)
        call coclone(self%j2rp2_varname)
        call coclone(self%j2rp2_varid)
        call coclone(self%j4rp4_varname)
        call coclone(self%j4rp4_varid)
        call coclone(self%k2_varname)
        call coclone(self%k2_varid)
        call coclone(self%q_varname)
        call coclone(self%Q_varid)
        call coclone(self%ke_orb_varname)
        call coclone(self%KE_orb_varid)
        call coclone(self%ke_spin_varname)
        call coclone(self%KE_spin_varid)
        call coclone(self%pe_varname)
        call coclone(self%PE_varid)
        call coclone(self%be_varname)
        call coclone(self%BE_varid)
        call coclone(self%te_varname)
        call coclone(self%TE_varid)
        call coclone(self%L_orbit_varname)
        call coclone(self%L_orbit_varid)
        call coclone(self%L_spin_varname)
        call coclone(self%L_spin_varid)
        call coclone(self%L_escape_varname)
        call coclone(self%L_escape_varid)
        call coclone(self%E_collisions_varname)
        call coclone(self%E_collisions_varid)
        call coclone(self%E_untracked_varname)
        call coclone(self%E_untracked_varid)
        call coclone(self%GMescape_varname)
        call coclone(self%GMescape_varid)
        call coclone(self%origin_type_varname)
        call coclone(self%origin_type_varid)
        call coclone(self%origin_time_varname)
        call coclone(self%origin_time_varid)
        call coclone(self%collision_id_varname)
        call coclone(self%collision_id_varid)
        call coclone(self%origin_rh_varname)
        call coclone(self%origin_rh_varid)
        call coclone(self%origin_vh_varname)
        call coclone(self%origin_vh_varid)
        call coclone(self%discard_time_varname)
        call coclone(self%discard_time_varid)
        call coclone(self%discard_rh_varname)
        call coclone(self%discard_rh_varid)
        call coclone(self%discard_vh_varname)
        call coclone(self%discard_vh_varid)
        call coclone(self%discard_body_id_varname)
        call coclone(self%lpseudo_vel_exists)
        return
    end subroutine swiftest_coarray_coclone_nc

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
        integer(I4B) :: i

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

        deallocate(isalloc,n,tmp)
    
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
        integer(I4B) :: i, img, si
        integer(I4B), allocatable :: n[:]
        logical, allocatable :: isalloc[:]

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

        deallocate(isalloc,n,tmp)
  
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
        type(swiftest_particle_info), dimension(:), allocatable :: vari1
  
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
        class(swiftest_tp), allocatable, codimension[:] :: cotp
        character(len=NAMELEN) :: image_num_char

        if (.not.param%lcoarray) return

        if (this_image() == 1 .or. param%log_output) then
            write(image_num_char,*) num_images()
            write(param%display_unit,*) " Collecting test particles from " // trim(adjustl(image_num_char)) // " images."
            if (param%log_output) flush(param%display_unit)
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
        integer(I4B) :: istart, iend, ntot, num_per_image, ncopy
        logical, dimension(:), allocatable :: lspill_list
        integer(I4B), codimension[:], allocatable  :: ntp
        character(len=NAMELEN) :: image_num_char, ntp_num_char
        class(swiftest_tp), allocatable, codimension[:] :: cotp
        class(swiftest_tp), allocatable :: tmp

        if (.not.param%lcoarray) return

        allocate(ntp[*])
        ntp = nbody_system%tp%nbody
        sync all
        ntot = ntp[1]
        if (ntot == 0) return

        write(image_num_char,*) num_images()

        if (this_image() == 1 .or. param%log_output) then
            write(ntp_num_char,*) ntot
            write(param%display_unit,*) " Distributing " // trim(adjustl(ntp_num_char)) // " test particles across " // trim(adjustl(image_num_char)) // " images."
            if (param%log_output) flush(param%display_unit)
        end if

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

        write(image_num_char,*) this_image()
        write(ntp_num_char,*) nbody_system%tp%nbody
        write(param%display_unit,*) "Image " // trim(adjustl(image_num_char)) // " ntp: " // trim(adjustl(ntp_num_char))
        if (param%log_output) flush(param%display_unit)

        deallocate(ntp, lspill_list, tmp, cotp)

        return
    end subroutine swiftest_coarray_distribute_system


end submodule s_swiftest_coarray
