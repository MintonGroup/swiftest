! Copyright 2025 - David Minton
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

    module subroutine swiftest_coarray_balance_system(self, param)
        !! author: David A. Minton
        !!
        !! Checks whether or not the system needs to be rebalance. Rebalancing occurs when the sum of the absolute difference 
        !! between !! the number of test particles on each image and the average number distributed across all images is larger than
        !! the number of images. 
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: self
            !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param       
             !! Current run configuration parameters 
        ! Internals
        integer(I4B), codimension[*], save :: ntp, nsum, ndiff
        integer(I4B) :: img, nimg, navg
        character(len=NAMELEN) :: ntp_str, nsum_str, ndiff_str, navg_str, nimg_str

        if (.not.param%lcoarray) return
        img = this_image()
        nimg = num_images()
        ntp = self%tp%nbody
        sync all
        write(param%display_unit,*) "Checking whether test particles need to be reblanced."
        nsum = ntp
        ndiff = ntp
        call co_sum(nsum)
        navg = floor(real(nsum) / nimg)
        ndiff = abs(ntp - navg)
        call co_sum(ndiff)

        write(ntp_str,*) ntp
        write(nsum_str,*) nsum
        write(ndiff_str,*) ndiff
        write(navg_str,*) navg
        write(nimg_str,*) nimg

        if (ndiff >= nimg) then
            write(param%display_unit,*) trim(adjustl(ndiff_str)) // ">=" // trim(adjustl(nimg_str)) // ": Rebalancing"
            flush(param%display_unit)
            call self%coarray_collect(param)
            call self%coarray_distribute(param)
            sync all
            write(param%display_unit,*) "Rebalancing complete"
        else
            write(param%display_unit,*) trim(adjustl(ndiff_str)) // "<" // trim(adjustl(nimg_str)) // ": No rebalancing needed"
        end if
        flush(param%display_unit)
        return
    end subroutine swiftest_coarray_balance_system


    module subroutine swiftest_coarray_coclone_param(param)
        !! author: David A. Minton 
        !! 
        !! Broadcasts the image 1 parameter to all other images in a parameter coarray 
        implicit none
        ! Arguments
        type(swiftest_parameters),intent(inout),codimension[*]  :: param  
            !! Collection of parameters 


        call co_broadcast(param%integrator,1)
        call co_broadcast(param%param_file_name,1)
        call co_broadcast(param%t0,1)
        call co_broadcast(param%tstart,1)
        call co_broadcast(param%tstop,1)
        call co_broadcast(param%dt,1)
        call co_broadcast(param%iloop,1)
        call co_broadcast(param%nloops,1)
        call co_broadcast(param%incbfile,1)
        call co_broadcast(param%inplfile,1)
        call co_broadcast(param%intpfile,1)
        call co_broadcast(param%nc_in,1)
        call co_broadcast(param%in_type,1)
        call co_broadcast(param%in_form,1)
        call co_broadcast(param%istep_out,1)
        call co_broadcast(param%nstep_out,1)
        call co_broadcast(param%fstep_out,1)
        call co_broadcast(param%ltstretch,1)
        call co_broadcast(param%outfile,1)
        call co_broadcast(param%out_type,1)
        call co_broadcast(param%out_form,1)
        call co_broadcast(param%out_stat,1)
        call co_broadcast(param%dump_cadence,1)
        call co_broadcast(param%rmin,1)
        call co_broadcast(param%rmax,1)
        call co_broadcast(param%rmaxu,1)
        call co_broadcast(param%qmin,1)
        call co_broadcast(param%qmin_coord,1)
        call co_broadcast(param%qmin_alo,1)
        call co_broadcast(param%qmin_ahi,1)
        call co_broadcast(param%MU2KG,1)
        call co_broadcast(param%TU2S,1)
        call co_broadcast(param%DU2M,1)
        call co_broadcast(param%GU,1)
        call co_broadcast(param%inv_c2,1)
        call co_broadcast(param%GMTINY,1)
        call co_broadcast(param%min_GMfrag,1)
        call co_broadcast(param%nfrag_reduction,1)
        call co_broadcast(param%lmtiny_pl,1)
        call co_broadcast(param%collision_model,1)
        call co_broadcast(param%encounter_save,1)
        call co_broadcast(param%lenc_save_trajectory,1)
        call co_broadcast(param%lenc_save_closest,1)
        call co_broadcast(param%interaction_loops,1)
        call co_broadcast(param%encounter_check_plpl,1)
        call co_broadcast(param%encounter_check_pltp,1)
        call co_broadcast(param%lflatten_interactions,1)
        call co_broadcast(param%lencounter_sas_plpl,1)
        call co_broadcast(param%lencounter_sas_pltp,1)
        call co_broadcast(param%lrhill_present,1)
        call co_broadcast(param%lextra_force,1)
        call co_broadcast(param%lbig_discard,1)
        call co_broadcast(param%lclose,1)
        call co_broadcast(param%lenergy,1)
        call co_broadcast(param%lnon_spherical_cb,1)
        call co_broadcast(param%lrotation,1)
        call co_broadcast(param%ltides,1)
        call co_broadcast(param%E_orbit_orig,1)
        call co_broadcast(param%GMtot_orig,1)
        call co_broadcast(param%L_total_orig,1)
        call co_broadcast(param%L_orbit_orig,1)
        call co_broadcast(param%L_rot_orig,1)
        call co_broadcast(param%L_escape,1)
        call co_broadcast(param%GMescape,1)
        call co_broadcast(param%E_collisions,1)
        call co_broadcast(param%E_untracked,1)
        call co_broadcast(param%lfirstenergy,1)
        call co_broadcast(param%lfirstkick,1)
        call co_broadcast(param%lrestart,1)
        call co_broadcast(param%display_style,1)
        call co_broadcast(param%display_unit,1)
        call co_broadcast(param%log_output ,1)
        call co_broadcast(param%lgr,1)
        call co_broadcast(param%lyarkovsky,1)
        call co_broadcast(param%lyorp,1)
        call coclone(param%seed)
        call co_broadcast(param%lcoarray,1)

        return
    end subroutine swiftest_coarray_coclone_param 


    module subroutine swiftest_coarray_component_clone_info(var,src_img)
        !! author: David A. Minton
        !!
        !! Copies a component of a coarray derived type from the specified source image to the current local one. 
        !! The default source image is 1
        !! swiftest_particle_info scalar version
        implicit none
        ! Arguments
        type(swiftest_particle_info), intent(inout) :: var
            !! Variable to be cloned
        integer(I4B), intent(in),optional :: src_img
            !! Source image to clone from
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
        !! Copies a component of a coarray derived type from the specified source image to the current local one. 
        !! The default source image is 1
        !! swiftest_particle_info 1D allocatable array version
        implicit none
        ! Arguments
        type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: var
            !! Variable to be cloned
        integer(I4B), intent(in),optional :: src_img
            !! Source image to clone from
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


    module subroutine swiftest_coarray_component_collect_info_arr1D(var,dest_img)
        !! author: David A. Minton
        !!
        !! Collects components of a coarray derived type from all images and combines them into destination image component. 
        !! The default destination image is 1
        !! swiftest_particle_info 1D allocatable array version
        implicit none
        ! Arguments
        type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: var
            !! Variable to be collected
        integer(I4B), intent(in),optional :: dest_img
            !! Destination image to collect to
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
    end subroutine swiftest_coarray_component_collect_info_arr1D


    module subroutine swiftest_coarray_collect_system(self, param)
        !! author: David A. Minton
        !!
        !! Collects all the test particles from other images into the image #1 test particle system
        use whm, only: whm_tp, cocollect
        use rmvs, only: rmvs_tp, cocollect
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: self
            !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param        
            !! Current run configuration parameters 
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

        select type(tp => self%tp)
        class is (rmvs_tp)
            allocate(rmvs_cotp[*], source=tp)
            call cocollect(rmvs_cotp)
            if (this_image() == 1) then
                deallocate(self%tp)
                allocate(self%tp, source=rmvs_cotp)
            end if
            deallocate(rmvs_cotp)
        class is (whm_tp)
            allocate(whm_cotp[*], source=tp)
            call cocollect(whm_cotp)
            if (this_image() == 1) then
                deallocate(self%tp)
                allocate(self%tp, source=whm_cotp)
            end if
            deallocate(whm_cotp)
        end select

        return
    end subroutine swiftest_coarray_collect_system
 
 
    module subroutine swiftest_coarray_distribute_system(self, param)
        !! author: David A. Minton
        !!
        !! Distributes test particles from image #1 out to all images.
        use whm, only: whm_tp, coclone
        use rmvs, only: rmvs_tp, coclone
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: self
            !! Swiftest nbody system 
        class(swiftest_parameters),   intent(inout) :: param        
            !! Current run configuration parameters 
        ! Internals
        integer(I4B) :: istart, iend, ntot, num_per_image, nremaining, i
        logical, dimension(:), allocatable :: lspill_list
        integer(I4B), codimension[:], allocatable  :: ntp
        character(len=NAMELEN) :: image_num_char, ntp_num_char
        type(whm_tp), allocatable, codimension[:] :: whm_cotp
        type(rmvs_tp), allocatable, codimension[:] :: rmvs_cotp
        class(swiftest_tp), allocatable :: tmp

        if (.not.param%lcoarray) return
        associate(nbody_system => self)

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
            num_per_image = ntot / num_images()
            nremaining = mod(ntot, num_images())
            istart = 1
    
            ! Distribute the particles as evenly as possible
            lspill_list(:) = .true.
            if (this_image() <= nremaining) then
                istart = (this_image() - 1) * (num_per_image + 1) + 1
                iend = istart + num_per_image
            else
                istart = nremaining * (num_per_image + 1) + (this_image() - nremaining - 1) * num_per_image + 1
                iend = istart + num_per_image - 1
            end if
    
            if (istart <= iend) then
                lspill_list(istart:iend) = .false.
            end if

            select type(tp => nbody_system%tp)
            class is (rmvs_tp)
                allocate(rmvs_cotp[*], source=tp)
                call coclone(rmvs_cotp)
                if (this_image() /= 1) then
                    deallocate(nbody_system%tp)
                    allocate(nbody_system%tp, source=rmvs_cotp)
                end if
                deallocate(rmvs_cotp)
                allocate(tmp, mold=tp)
            class is (whm_tp) 
                allocate(whm_cotp[*], source=tp)
                call coclone(whm_cotp)
                if (this_image() /= 1) then
                    deallocate(nbody_system%tp)
                    allocate(nbody_system%tp, source=whm_cotp)
                end if
                deallocate(whm_cotp)
                allocate(tmp, mold=tp)
            end select

            call nbody_system%tp%spill(tmp, lspill_list(:), ldestructive=.true.)

            write(image_num_char,*) this_image()
            write(ntp_num_char,*) nbody_system%tp%nbody
            write(param%display_unit,*) "Image " // trim(adjustl(image_num_char)) // " ntp: " // trim(adjustl(ntp_num_char))
            if (param%log_output) flush(param%display_unit)

            deallocate(ntp, lspill_list, tmp)
            if (param%log_output) flush(param%display_unit)
            sync all
        end associate

        return
    end subroutine swiftest_coarray_distribute_system


end submodule s_swiftest_coarray
