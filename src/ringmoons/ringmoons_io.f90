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
        class(ringmoons_seed), intent(inout) :: self
        class(swiftest_netcdf_parameters), intent(inout) :: nc
        class(swiftest_parameters), intent(in) :: param
        return
    end subroutine ringmoons_io_write_frame_seed

end submodule s_ringmoons_io
