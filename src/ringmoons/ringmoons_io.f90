! Copyright 2026 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(ringmoons) s_ringmoons_io
    use swiftest
    use netcdf
contains        

    module subroutine ringmoons_io_netcdf_flush(self, param)
        !! author: David A. Minton
        !!
        !! Flushes the current buffer to disk by closing and re-opening the file.
        !!    
        implicit none
        ! Arguments
        class(ringmoons_netcdf_parameters), intent(inout) :: self 
            !! Parameters used to identify a particular NetCDF dataset
        class(swiftest_parameters),         intent(inout) :: param 
            !! Current run configuration parameters 

        call self%close()
        call self%open(param,readonly=.false.)

        return
    end subroutine ringmoons_io_netcdf_flush

    module subroutine ringmoons_io_netcdf_open(self, param, readonly)
        !! author: David A. Minton
        !!
        !! Open up and existing Ringmoons file and create any missing variables if needed.
        use, intrinsic :: ieee_arithmetic
        implicit none
        class(ringmoons_netcdf_parameters), intent(inout) :: self     
            !! Parameters used to identify a particular NetCDF dataset
        class(swiftest_parameters),         intent(in)    :: param    
            !! Current run configuration parameters
        logical, optional,                  intent(in)    :: readonly 
            !! Logical flag indicating that this should be open read only
        ! Internals
        integer(I4B), parameter :: NO_FILL = 0
        integer(I4B) :: mode, status
        character(len=STRMAX) :: errmsg
        logical :: fileExists
        real(DP) :: dfill
        real(SP) :: sfill

        associate(nc => self)
            mode = NF90_WRITE
            if (present(readonly)) then
                if (readonly) mode = NF90_NOWRITE
            end if
    
            dfill = ieee_value(dfill, IEEE_QUIET_NAN)
            sfill = ieee_value(sfill, IEEE_QUIET_NAN)

            select case (param%out_type)
            case("NETCDF_FLOAT")
                nc%out_type = NF90_FLOAT
            case("NETCDF_DOUBLE")
                nc%out_type = NF90_DOUBLE
            case default
                write(*,*) trim(adjustl(param%out_type)), " is an invalid OUT_TYPE"
            end select


            inquire(file=nc%file_name, exist=fileExists)
            if (.not.fileExists) then
                write(*,*) "File not found: ", trim(adjustl(nc%file_name))
                call base_util_exit(FAILURE, param%display_unit) 
            end if

            write(errmsg,*) "ringmoons_io_netcdf_open nf90_open ",trim(adjustl(nc%file_name))
            call netcdf_io_check( nf90_open(nc%file_name, mode, nc%id), errmsg)
            self%lfile_is_open = .true.

            ! Dimensions
            call netcdf_io_check( nf90_inq_dimid(nc%id, nc%time_dimname, nc%time_dimid), &
                                    "ringmoons_io_netcdf_open nf90_inq_dimid time_dimid"  )
            call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%time_dimid, nc%time_dimname, len=nc%max_tslot), &
                                    "ringmoons_io_netcdf_open nf90_inquire_dimension max_tslot"  )
            call netcdf_io_check( nf90_inq_dimid(nc%id, nc%ringbin_dimname, nc%ringbin_dimid), &
                                    "ringmoons_io_netcdf_open nf90_inq_dimid ringbin_dimid"  )
            call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%ringbin_dimid, nc%ringbin_dimname, len=nc%nbins), &
                                    "ringmoons_io_netcdf_open nf90_inquire_dimension nbin"  )

            ! Dimension coordinates
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%time_dimname, nc%time_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid time_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%ringbin_dimname, nc%ringbin_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid ringbin_varid" )

            ! The following variables are required for initial conditions
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%r_varname, nc%r_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid r_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%sigma_varname, nc%sigma_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid sigma_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%r_p_varname, nc%r_p_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid r_p_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%m_p_varname, nc%m_p_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid m_p_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%r_inner_varname, nc%r_inner_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid r_inner_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%r_outer_varname, nc%r_outer_varid), &
                                    "ringmoons_io_netcdf_open nf90_inq_varid r_outer_varid" )

            ! The following variables are not required for initial conditions, so we must create them if they don't exist yet
            call nc%add_new_var(nc%tau_varname, nc%out_type, [nc%ringbin_dimid, nc%time_dimid], nc%tau_varid, &
                                 "ringmoons_io_netcdf_open add_new_var tau_varid")
            call nc%add_new_var(nc%nu_varname, nc%out_type, [nc%ringbin_dimid, nc%time_dimid], nc%nu_varid, &
                                 "ringmoons_io_netcdf_open add_new_var nu_varid")
            call nc%add_new_var(nc%toomre_varname, nc%out_type, [nc%ringbin_dimid, nc%time_dimid], nc%toomre_varid, &
                                 "ringmoons_io_netcdf_open add_new_var toomre_varid")
            call nc%add_new_var(nc%vrel_p_varname, nc%out_type, [nc%ringbin_dimid, nc%time_dimid], nc%vrel_p_varid, &
                                 "ringmoons_io_netcdf_open add_new_var vrel_p_varid")
        end associate
        return
    end subroutine ringmoons_io_netcdf_open

    module subroutine ringmoons_io_read_frame_ring(self, t, param)
        !! author: David A. Minton
        !!
        !! Read in ring data from a NetCDF file and initialize the ring from the input data
        implicit none
        class(ringmoons_ring), intent(inout) :: self
        real(DP), intent(in)                   :: t
        class(swiftest_parameters), intent(in) :: param
        ! Internals
        integer(I4B)                              :: i, nbin
        real(DP), dimension(:), allocatable       :: rtemp
        real(DP), dimension(:,:), allocatable     :: vectemp
        integer(I4B), dimension(:), allocatable   :: itemp
        real(DP), dimension(1)                    :: tmp_scalar

        if (param%in_type == "ASCII") then
            write(*,*) "ASCII not supported for RINGMOONS"
            call base_util_exit(FAILURE,param%display_unit)
        else
            associate(nc => self%nc, tslot => self%nc%tslot)
                call nc%open(param, readonly=.false.) ! Set to False so that we can add any missing variables
                call nc%find_tslot(t, tslot)
                self%t = t
                call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%time_dimid, len=nc%max_tslot), &
                                    "ringmoons_io_read_frame_ring nf90_inquire_dimension time_dimid"  )
                tslot = min(tslot, nc%max_tslot)

                call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%ringbin_dimid, len=nbin), &
                                    "ringmoons_io_read_frame_ring nf90_inquire_dimension ringbin_dimid"  )
                if (nbin == 0) then
                    write(*,*) "No ring bins found in NetCDF file"
                    call base_util_exit(FAILURE,param%display_unit)
                end if
                call self%setup(nbin, param)
                call netcdf_io_check( nf90_get_var(nc%id, nc%r_varid, self%r(1:nbin), start=[1, tslot], count=[nbin,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar r_varid"  )
                call netcdf_io_check( nf90_get_var(nc%id, nc%sigma_varid, self%sigma(1:nbin), start=[1, tslot], count=[nbin,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar sigma_varid"  )
                call netcdf_io_check( nf90_get_var(nc%id, nc%r_p_varid, self%r_p(1:nbin), start=[1, tslot], count=[nbin,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar r_p_varid"  )
                call netcdf_io_check( nf90_get_var(nc%id, nc%m_p_varid, self%m_p(1:nbin), start=[1, tslot], count=[nbin,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar m_p_varid"  )
                call netcdf_io_check( nf90_get_var(nc%id, nc%r_inner_varid, tmp_scalar, start=[1], count=[1]), &
                                  "netcdf_io_read_frame_system nf90_getvar r_inner_varid"  )
                self%r_inner = tmp_scalar(1)
                call netcdf_io_check( nf90_get_var(nc%id, nc%r_outer_varid, tmp_scalar, start=[1], count=[1]), &
                                  "netcdf_io_read_frame_system nf90_getvar r_outer_varid"  )
                self%r_outer = tmp_scalar(1)
                call nc%close()
            end associate
        end if

        return
    end subroutine ringmoons_io_read_frame_ring

    module subroutine ringmoons_io_read_frame_seed(self, t, nc, param) 
        implicit none
        class(ringmoons_seed),             intent(inout) :: self
        real(DP),                          intent(in)    :: t  
        class(swiftest_netcdf_parameters), intent(inout) :: nc
        class(swiftest_parameters),        intent(inout) :: param
        ! Internals
        integer(I4B) :: i, nseed, idmax, readstat
        character(len=NAMELEN), dimension(:), allocatable :: carrtemp
        logical, dimension(:), allocatable :: lvalid, seedmask
        integer(I4B), dimension(:), allocatable :: status
        character(len=NAMELEN) :: ctemp

        associate(tslot => nc%tslot)
            call nc%open(param, readonly=.false.) ! Set to False so that we can add any missing variables
            call nc%find_tslot(t, tslot)
            readstat = nf90_inq_varid(nc%id, nc%nseed_varname, nc%nseed_varid)
            if (readstat == NF90_NOERR) then
                call netcdf_io_check( nf90_get_var(nc%id, nc%nseed_varid,  nseed, start=[tslot]), &
                                    "ringmoons_io_read_frame_seed nf90_getvar nseed_varid"  )
            else
                nseed = 0 
            end if         
            call self%setup(nseed, param)
            if (nseed == 0) then
                call nc%close()
                return
            end if
            call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%name_dimid, len=idmax), &
                                  "ringmoons_io_read_frame_seed  nf90_inquire_dimension name_dimid"  )

            allocate(carrtemp(idmax))
            allocate(lvalid(idmax))
            allocate(seedmask(idmax))
            allocate(status(idmax))
            seedmask(:) = .false.
            lvalid(:) = .false.
            status(:) = INACTIVE
            call netcdf_io_check( nf90_get_var(nc%id, nc%ptype_varid, carrtemp, count=[NAMELEN, idmax]), &
                                    "ringmoons_io_read_frame_seed nf90_get_var ptype_varid"  )
            where(carrtemp(:) == SEED_TYPE_NAME) seedmask(:) = .true.

            call netcdf_io_check( nf90_get_var(nc%id, nc%status_varid, status, count=[tslot, idmax]), &
                                    "ringmoons_io_read_frame_seed nf90_get_var status_varid"  )
            where(seedmask(:)) lvalid(:) = status(:) == ACTIVE

            do i = 1, idmax
                if (lvalid(i)) then
                    call netcdf_io_check( nf90_get_var(nc%id, nc%id_varid, self%id(i), start=[i]), &
                                    "ringmoons_io_read_frame_seed nf90_get_var id_varid"  )
                    call netcdf_io_check( nf90_get_var(nc%id, nc%a_varid, self%a(i), start=[i, tslot]), &
                                    "ringmoons_io_read_frame_seed nf90_get_var a_varid"  )
                    call netcdf_io_check( nf90_get_var(nc%id, nc%name_varid, ctemp, start=[1, i], count=[NAMELEN, 1]), &
                                    "ringmoons_io_read_frame_seed nf90_get_var name_varid"  )
                    self%info(i)%name = trim(adjustl(ctemp))
                    call netcdf_io_check( nf90_get_var(nc%id, nc%mass_varid, self%mass(i), start=[i, tslot]), &
                                    "ringmoons_io_read_frame_seed nf90_get_var mass_varid"  )  
                    call netcdf_io_check( nf90_get_var(nc%id, nc%ptype_varid, ctemp, start=[1, i], count=[NAMELEN, 1]), &
                                    "ringmoons_io_read_frame_seed nf90_get_var particle_type_varid"  )
                    self%info(i)%particle_type = trim(adjustl(ctemp))
                end if
            end do
        end associate
        return
    end subroutine ringmoons_io_read_frame_seed

    module subroutine ringmoons_io_write_frame_ring(self, param) 
        !! author: David A. Minton
        !!
        !! Writes a frame of ring output to the binary output file.
        implicit none
        class(ringmoons_ring), intent(inout) :: self
        class(swiftest_parameters), intent(in) :: param
        ! Internals
        integer(I4B) :: old_mode, tmp

        associate(nc => self%nc, ring => self, tslot => self%nc%tslot, t => self%t, nbins => self%nbins)
            if (.not.nc%lfile_is_open) then
                call nc%open(param, readonly=.false.)
            end if
            call nc%find_tslot(t, tslot)
            call netcdf_io_check( nf90_set_fill(nc%id, NF90_NOFILL, old_mode), "ringmoons_io_write_frame_ring nf90_set_fill" )
            call netcdf_io_check( nf90_put_var(nc%id, nc%time_varid, t, start=[tslot]), &
                                    "ringmoons_io_write_frame_ring nf90_put_var time_varid" )
            call netcdf_io_check( nf90_put_var(nc%id, nc%r_varid, self%r(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var r_varid"  ) 
            call netcdf_io_check( nf90_put_var(nc%id, nc%sigma_varid, self%sigma(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var sigma_varid"  ) 
            call netcdf_io_check( nf90_put_var(nc%id, nc%r_p_varid, self%r_p(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var r_p_varid"  ) 
            call netcdf_io_check( nf90_put_var(nc%id, nc%m_p_varid, self%m_p(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var m_p_varid"  ) 
            call netcdf_io_check( nf90_put_var(nc%id, nc%tau_varid, self%tau(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var tau_varid"  ) 
            call netcdf_io_check( nf90_put_var(nc%id, nc%nu_varid, self%nu(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var nu_varid"  ) 
            call netcdf_io_check( nf90_put_var(nc%id, nc%toomre_varid, self%Q(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var toomre_varid"  ) 
            call netcdf_io_check( nf90_put_var(nc%id, nc%vrel_p_varid, self%vrel_p(1:nbins), start=[1,tslot], count=[nbins,1]), &
                                  "ringmoons_io_write_frame_ring nf90_put_var vrel_p_varid"  ) 

            call netcdf_io_check( nf90_set_fill(nc%id, old_mode, tmp), "ringmoons_io_write_frame_body nf90_set_fill old_mode" )
        end associate

        return
    end subroutine ringmoons_io_write_frame_ring

    module subroutine ringmoons_io_write_frame_seed(self, nc, param) 
        implicit none
        ! Arguments
        class(ringmoons_seed), intent(inout) :: self
            !! Ringmoons seed object
        class(swiftest_netcdf_parameters), intent(inout) :: nc
            !! Parameters used to write a NetCDF dataset to file
        class(swiftest_parameters), intent(in) :: param
            !! Current run configuration parameters
        ! Internals
        integer(I4B) :: i, j, idslot, old_mode, tmp
        integer(I4B), dimension(:), allocatable   :: ind
        real(DP), dimension(:), allocatable       :: mu,a, e, inc, capom, omega, capm
        real(DP), dimension(:,:), allocatable     :: rh, vh, pvh
        real(DP) :: rtmp
        character(len=NAMELEN) :: charstring

        call netcdf_io_check( nf90_set_fill(nc%id, NF90_NOFILL, old_mode), "ringmoons_io_write_frame_seed nf90_set_fill" )
        associate(n => self%nbody, tslot => nc%tslot)
            call nc%add_new_var(nc%nseed_varname, nc%out_type, [nc%time_dimid], nc%nseed_varid, &
                                 "ringmoons_io_write_frame_seed add_new_var nseed_varid")
            call netcdf_io_check( nf90_put_var(nc%id, nc%nseed_varid, n, start=[tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var nseed_varid"  )

            if (n == 0) return
            ! Seeds have no true orbital elements, but we can fill in the blanks so that we have something to use if we want to 
            ! convert between cartesian and orbital element
            allocate(a(n),e(n),inc(n),capom(n),omega(n),capm(n),mu(n))
            allocate(rh(NDIM,n),vh(NDIM,n))
            a(1:n) = self%a(1:n)
            mu(1:n) = self%mu(1:n)
            e(:) = 0.0_DP
            inc(:) = 0.0_DP
            call random_number(capom(:))
            call random_number(omega(:))
            call random_number(capm(:))
            capom(:) = capom(:) * 360.0_DP
            omega(:) = omega(:) * 360.0_DP
            capm(:) = capm(:) * 360.0_DP
#ifdef DOCONLOC
            do concurrent (i = 1:n) shared(mu, a, e, inc, capom, omega, capm, rh, vh)
#else
            do concurrent (i = 1:n)
#endif
                call swiftest_orbel_el2xv(mu(i), a(i), e(i), inc(i), capom(i), omega(i), capm(i), & 
                                            rh(1,i), rh(2,i), rh(3,i), vh(1,i), vh(2,i), vh(3,i)) 
            end do
            if (param%lgr) allocate(pvh, source=vh)
            ! Assign id values to the seeds
            do i =1, n
                if (self%id(i) < 0) then
                    self%maxid = self%maxid + 1
                    self%id(i) = self%maxid 
                end if
            end do

            call util_sort(self%id(1:n), ind)

            do i = 1, n
                j = ind(i)
                call nc%find_idslot(self%id(j), idslot) 
                call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid, self%id(j), start=[idslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var id_varid"  )
                call netcdf_io_check( nf90_put_var(nc%id, nc%status_varid, self%status(j), start=[idslot,tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var status_varid"  )
                call netcdf_io_check( nf90_put_var(nc%id, nc%mass_varid, self%mass(j), start=[idslot,tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var mass_varid"  )
                call netcdf_io_check( nf90_put_var(nc%id, nc%Gmass_varid, self%Gmass(j), start=[idslot,tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var Gmass_varid"  )
                call netcdf_io_check( nf90_put_var(nc%id, nc%radius_varid, self%radius(j), start=[idslot,tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var radius_varid"  )
                call netcdf_io_check( nf90_put_var(nc%id, nc%rhill_varid, self%rhill(j), start=[idslot,tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var rhill_varid"  )
                charstring = trim(adjustl(self%info(j)%name))
                call netcdf_io_check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var name_varid"  )

                charstring = trim(adjustl(self%info(j)%particle_type))
                call netcdf_io_check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var particle_type_varid"  )

                if (param%lgr) call swiftest_gr_pseudovel2vel(param, mu(j), rh(:, j), pvh(:, j), vh(:,j))

                if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
                    call netcdf_io_check( nf90_put_var(nc%id, nc%rh_varid, rh(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]),&
                                    "ringmoons_io_write_frame_seed nf90_put_var rh_varid"  )
                    if (param%lgr) then ! Convert from pseudovelocity to heliocentric without replacing the current value of 
                                        !  pseudovelocity
                        call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid, vh(:,j), start=[1,idslot, tslot], count=[NDIM,1,1]),&
                                    "ringmoons_io_write_frame_seed nf90_put_var vh_varid"  )
                        call netcdf_io_check( nf90_put_var(nc%id, nc%gr_pseudo_vh_varid, pvh(:, j), start=[1,idslot, tslot], &
                                                            count=[NDIM,1,1]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var gr_pseudo_vhx_varid"  )

                    else
                        call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid, vh(:, j), start=[1,idslot, tslot], &
                                                            count=[NDIM,1,1]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var vh_varid"  )
                    end if
                end if

                if ((param%out_form == "EL") .or. (param%out_form == "XVEL")) then
                    call netcdf_io_check( nf90_put_var(nc%id, nc%a_varid, a(j), start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body a_varid"  )
                    call netcdf_io_check( nf90_put_var(nc%id, nc%e_varid, e(j), start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body e_varid"  )
                    call netcdf_io_check( nf90_put_var(nc%id, nc%inc_varid, inc(j), start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body inc_varid"  )
                    call netcdf_io_check( nf90_put_var(nc%id, nc%capom_varid, capom(j), start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body capom_varid"  )
                    call netcdf_io_check( nf90_put_var(nc%id, nc%omega_varid, omega(j), start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body omega_varid"  )
                    call netcdf_io_check( nf90_put_var(nc%id, nc%capm_varid, capm(j), start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body capm_varid"  ) 
                    rtmp = mod(omega(j) + capom(j), 360.0_DP)
                    call netcdf_io_check( nf90_put_var(nc%id, nc%varpi_varid,rtmp, start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body varpi_varid"  ) 
                    rtmp = mod(capm(j)+rtmp, 360.0_DP)
                    call netcdf_io_check( nf90_put_var(nc%id, nc%lam_varid, rtmp, start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body lam_varid"  ) 
                    call netcdf_io_check( nf90_put_var(nc%id, nc%f_varid, capm, start=[idslot, tslot]), &
                                    "ringmoons_io_write_frame_seed nf90_put_var body f_varid"  ) 
                end if

            end do
            deallocate(mu,a,e,inc,capom,omega,capm)
            deallocate(rh,vh)
            if (param%lgr) deallocate(pvh)

        end associate

        return
    end subroutine ringmoons_io_write_frame_seed

end submodule s_ringmoons_io
