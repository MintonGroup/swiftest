! Copyight 2024 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

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


    module subroutine swiftest_coarray_coclone_kin(kin)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        type(swiftest_kinship),intent(inout),codimension[*]  :: kin  !! Swiftest kinship object

        call coclone(kin%parent)
        call coclone(kin%nchild)
        call coclone(kin%child)

        return
    end subroutine swiftest_coarray_coclone_kin


    module subroutine swiftest_coarray_coclone_nc(nc)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        type(swiftest_netcdf_parameters),intent(inout),codimension[*]  :: nc  !! Swiftest body object

        call coclone(nc%file_name)
        call coclone(nc%lfile_is_open)
        call coclone(nc%out_type)
        call coclone(nc%id)
        call coclone(nc%tslot)
        call coclone(nc%max_tslot)
        call coclone(nc%idvals)
        call coclone(nc%idslot)
        call coclone(nc%max_idslot)
        call coclone(nc%str_dimname)
        call coclone(nc%str_dimid)
        call coclone(nc%time_dimname)
        call coclone(nc%time_dimid)
        call coclone(nc%time_varid)
        call coclone(nc%name_dimname)
        call coclone(nc%name_dimid)
        call coclone(nc%name_varid)
        call coclone(nc%space_dimname)
        call coclone(nc%space_dimid)
        call coclone(nc%space_varid)
        call coclone(nc%id_varname)
        call coclone(nc%id_varid)
        call coclone(nc%status_varname)
        call coclone(nc%status_varid)
        call coclone(nc%ptype_varname)
        call coclone(nc%ptype_varid)
        call coclone(nc%npl_varname)
        call coclone(nc%npl_varid)
        call coclone(nc%ntp_varname)
        call coclone(nc%ntp_varid)
        call coclone(nc%nplm_varname)
        call coclone(nc%nplm_varid)
        call coclone(nc%a_varname)
        call coclone(nc%a_varid)
        call coclone(nc%e_varname)
        call coclone(nc%e_varid)
        call coclone(nc%inc_varname)
        call coclone(nc%inc_varid)
        call coclone(nc%capom_varname)
        call coclone(nc%capom_varid)
        call coclone(nc%omega_varname)
        call coclone(nc%omega_varid)
        call coclone(nc%capm_varname)
        call coclone(nc%capm_varid)
        call coclone(nc%varpi_varname)
        call coclone(nc%varpi_varid)
        call coclone(nc%lam_varname)
        call coclone(nc%lam_varid)
        call coclone(nc%f_varname)
        call coclone(nc%f_varid)
        call coclone(nc%cape_varname)
        call coclone(nc%cape_varid)
        call coclone(nc%rh_varname)
        call coclone(nc%rh_varid)
        call coclone(nc%vh_varname)
        call coclone(nc%vh_varid)
        call coclone(nc%gr_pseudo_vh_varname)
        call coclone(nc%gr_pseudo_vh_varid)
        call coclone(nc%Gmass_varname)
        call coclone(nc%Gmass_varid)
        call coclone(nc%mass_varname)
        call coclone(nc%mass_varid)
        call coclone(nc%rhill_varname)
        call coclone(nc%rhill_varid)
        call coclone(nc%radius_varname)
        call coclone(nc%radius_varid)
        call coclone(nc%Ip_varname)
        call coclone(nc%Ip_varid)
        call coclone(nc%rot_varname)
        call coclone(nc%rot_varid)
        call coclone(nc%j2rp2_varname)
        call coclone(nc%j2rp2_varid)
        call coclone(nc%j4rp4_varname)
        call coclone(nc%j4rp4_varid)
        call coclone(nc%k2_varname)
        call coclone(nc%k2_varid)
        call coclone(nc%q_varname)
        call coclone(nc%Q_varid)
        call coclone(nc%ke_orb_varname)
        call coclone(nc%KE_orb_varid)
        call coclone(nc%ke_spin_varname)
        call coclone(nc%KE_spin_varid)
        call coclone(nc%pe_varname)
        call coclone(nc%PE_varid)
        call coclone(nc%be_varname)
        call coclone(nc%BE_varid)
        call coclone(nc%te_varname)
        call coclone(nc%TE_varid)
        call coclone(nc%L_orbit_varname)
        call coclone(nc%L_orbit_varid)
        call coclone(nc%L_spin_varname)
        call coclone(nc%L_spin_varid)
        call coclone(nc%L_escape_varname)
        call coclone(nc%L_escape_varid)
        call coclone(nc%E_collisions_varname)
        call coclone(nc%E_collisions_varid)
        call coclone(nc%E_untracked_varname)
        call coclone(nc%E_untracked_varid)
        call coclone(nc%GMescape_varname)
        call coclone(nc%GMescape_varid)
        call coclone(nc%origin_type_varname)
        call coclone(nc%origin_type_varid)
        call coclone(nc%origin_time_varname)
        call coclone(nc%origin_time_varid)
        call coclone(nc%collision_id_varname)
        call coclone(nc%collision_id_varid)
        call coclone(nc%origin_rh_varname)
        call coclone(nc%origin_rh_varid)
        call coclone(nc%origin_vh_varname)
        call coclone(nc%origin_vh_varid)
        call coclone(nc%discard_time_varname)
        call coclone(nc%discard_time_varid)
        call coclone(nc%discard_rh_varname)
        call coclone(nc%discard_rh_varid)
        call coclone(nc%discard_vh_varname)
        call coclone(nc%discard_vh_varid)
        call coclone(nc%discard_body_id_varname)
        call coclone(nc%lpseudo_vel_exists)
        call coclone(nc%lc_lm_exists)
        return
    end subroutine swiftest_coarray_coclone_nc


    module subroutine swiftest_coarray_clone_param(param)
        !! author: David A. Minton 
        !! 
        !! Broadcasts the image 1 parameter to all other images in a parameter coarray 
        implicit none
        ! Arguments
        type(swiftest_parameters),intent(inout),codimension[*]  :: param  
        !! Collection of parameters 
        ! Internals
        call coclone(param%integrator)
        call coclone(param%param_file_name)
        call coclone(param%t0)
        call coclone(param%tstart)
        call coclone(param%tstop)
        call coclone(param%dt)
        call coclone(param%iloop)
        call coclone(param%nloops)
        call coclone(param%incbfile)
        call coclone(param%inplfile)
        call coclone(param%intpfile)
        call coclone(param%nc_in)
        call coclone(param%in_type)
        call coclone(param%in_form)
        call coclone(param%istep_out)
        call coclone(param%nstep_out)
        call coclone(param%fstep_out)
        call coclone(param%ltstretch)
        call coclone(param%outfile)
        call coclone(param%out_type)
        call coclone(param%out_form)
        call coclone(param%out_stat)
        call coclone(param%dump_cadence)
        call coclone(param%rmin)
        call coclone(param%rmax)
        call coclone(param%rmaxu)
        call coclone(param%qmin)
        call coclone(param%qmin_coord)
        call coclone(param%qmin_alo)
        call coclone(param%qmin_ahi)
        call coclone(param%MU2KG)
        call coclone(param%TU2S)
        call coclone(param%DU2M)
        call coclone(param%GU)
        call coclone(param%inv_c2)
        call coclone(param%GMTINY)
        call coclone(param%min_GMfrag)
        call coclone(param%nfrag_reduction)
        call coclone(param%lmtiny_pl)
        call coclone(param%collision_model)
        call coclone(param%encounter_save)
        call coclone(param%lenc_save_trajectory)
        call coclone(param%lenc_save_closest)
        call coclone(param%interaction_loops)
        call coclone(param%encounter_check_plpl)
        call coclone(param%encounter_check_pltp)
        call coclone(param%lflatten_interactions)
        call coclone(param%lencounter_sas_plpl)
        call coclone(param%lencounter_sas_pltp)
        call coclone(param%lrhill_present)
        call coclone(param%lextra_force)
        call coclone(param%lbig_discard)
        call coclone(param%lclose)
        call coclone(param%lenergy)
        call coclone(param%lnon_spherical_cb)
        call coclone(param%lrotation)
        call coclone(param%ltides)
        call coclone(param%E_orbit_orig)
        call coclone(param%GMtot_orig)
        call coclonevec(param%L_total_orig)
        call coclonevec(param%L_orbit_orig)
        call coclonevec(param%L_spin_orig)
        call coclonevec(param%L_escape)
        call coclone(param%GMescape)
        call coclone(param%E_collisions)
        call coclone(param%E_untracked)
        call coclone(param%lfirstenergy)
        call coclone(param%lfirstkick)
        call coclone(param%lrestart)
        call coclone(param%display_style)
        call coclone(param%display_unit)
        call coclone(param%log_output )
        call coclone(param%lgr)
        call coclone(param%lyarkovsky)
        call coclone(param%lyorp)
        call coclone(param%seed)
        call coclone(param%lcoarray)

        return
    end subroutine swiftest_coarray_clone_param 

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
        !! Copies a component of a coarray derived type from the specified source image to the current local one. The default source
        !! image is 1 swiftest_kinship allocatable array version
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
            call coclone(tmp(i))
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


    module subroutine swiftest_coarray_collect_system(nbody_system, param)
        !! author: David A. Minton
        !!
        !! Collects all the test particles from other images into the image #1 test particle system
        use whm
        use rmvs
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
        ! Internals
        integer(I4B) :: i,j
        type(whm_tp), allocatable, codimension[:] :: whm_cotp
        type(rmvs_tp), allocatable, codimension[:] :: rmvs_cotp
        character(len=NAMELEN) :: image_num_char

        if (.not.param%lcoarray) return

        if (this_image() == 1 .or. param%log_output) then
            write(image_num_char,*) num_images()
            write(param%display_unit,*) " Collecting test particles from " // trim(adjustl(image_num_char)) // " images."
            if (param%log_output) flush(param%display_unit)
        end if

        select type(tp => nbody_system%tp)
        class is (rmvs_tp)
            allocate(rmvs_cotp[*], source=tp)
            call rmvs_coarray_cocollect_tp(rmvs_cotp)
            if (this_image() /= 1) then
                deallocate(nbody_system%tp)
                allocate(nbody_system%tp, source=rmvs_cotp)
            end if
            deallocate(rmvs_cotp)
        class is (whm_tp)
            allocate(whm_cotp[*], source=tp)
            call whm_coarray_cocollect_tp(whm_cotp)
            if (this_image() /= 1) then
                deallocate(nbody_system%tp)
                allocate(nbody_system%tp, source=whm_cotp)
            end if
            deallocate(whm_cotp)
        end select

        return
    end subroutine swiftest_coarray_collect_system
 
 
    module subroutine swiftest_coarray_distribute_system(nbody_system, param)
        !! author: David A. Minton
        !!
        !! Distributes test particles from image #1 out to all images.
        use whm
        use rmvs
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param        !! Current run configuration parameters 
        ! Internals
        integer(I4B) :: istart, iend, ntot, num_per_image, ncopy
        logical, dimension(:), allocatable :: lspill_list
        integer(I4B), codimension[:], allocatable  :: ntp
        character(len=NAMELEN) :: image_num_char, ntp_num_char
        type(whm_tp), allocatable, codimension[:] :: whm_cotp
        type(rmvs_tp), allocatable, codimension[:] :: rmvs_cotp
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
            write(param%display_unit,*) " Distributing " // trim(adjustl(ntp_num_char)) // " test particles across " // &
                                                            trim(adjustl(image_num_char)) // " images."
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

        select type(tp => nbody_system%tp)
        class is (rmvs_tp)
            allocate(rmvs_cotp[*], source=tp)
            call coclone(rmvs_cotp)
            if (this_image() /= 1) then
                deallocate(nbody_system%tp)
                allocate(nbody_system%tp, source=rmvs_cotp)
            end if
            deallocate(rmvs_cotp)
        class is (whm_tp) 
            allocate(whm_cotp[*], source=tp)
            call coclone(whm_cotp)
            if (this_image() /= 1) then
                deallocate(nbody_system%tp)
                allocate(nbody_system%tp, source=whm_cotp)
            end if
            deallocate(whm_cotp)
        end select

        allocate(tmp, mold=nbody_system%tp)
        call nbody_system%tp%spill(tmp, lspill_list(:), ldestructive=.true.)

        write(image_num_char,*) this_image()
        write(ntp_num_char,*) nbody_system%tp%nbody
        write(param%display_unit,*) "Image " // trim(adjustl(image_num_char)) // " ntp: " // trim(adjustl(ntp_num_char))
        if (param%log_output) flush(param%display_unit)

        deallocate(ntp, lspill_list, tmp)

        return
    end subroutine swiftest_coarray_distribute_system


end submodule s_swiftest_coarray
