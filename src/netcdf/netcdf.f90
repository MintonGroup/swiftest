!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_netcdf
   use swiftest
   use netcdf
contains

   module subroutine check(status, call_identifier)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Checks the status of all NetCDF operations to catch errors
      implicit none
      ! Arguments
      integer, intent (in) :: status !! The status code returned by a NetCDF function
      character(len=*), intent(in), optional :: call_identifier !! String that indicates which calling function caused the error for diagnostic purposes

      if(status /= nf90_noerr) then
         if (present(call_identifier)) write(*,*) "NetCDF error in ",trim(call_identifier)
         write(*,*) trim(nf90_strerror(status))
         call util_exit(FAILURE)
      end if

      return
   end subroutine check


   module subroutine netcdf_close(self)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Closes a NetCDF file
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset

      call check( nf90_close(self%id), "netcdf_close" )

      return
   end subroutine netcdf_close


   module subroutine netcdf_flush(self, param)
      !! author: David A. Minton
      !!
      !! Flushes the current buffer to disk by closing and re-opening the file.
      !!    
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 

      call self%close()
      call self%open(param)

      return
   end subroutine netcdf_flush


   module function netcdf_get_old_t_final_system(self, param) result(old_t_final)
      !! author: David A. Minton
      !!
      !! Validates the dump file to check whether the dump file initial conditions duplicate the last frame of the netcdf output.
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self
      class(swiftest_parameters),   intent(inout) :: param
      ! Result
      real(DP)                                    :: old_t_final
      ! Internals
      integer(I4B)                              :: itmax, idmax
      real(DP), dimension(:), allocatable       :: vals
      real(DP), dimension(1)                    :: rtemp
      real(DP), dimension(NDIM)                 :: rot0, Ip0, Lnow
      real(DP) :: KE_orb_orig, KE_spin_orig, PE_orig

      associate (nc => param%system_history%nc)
         call nc%open(param)
         call check( nf90_inquire_dimension(nc%id, nc%time_dimid, len=itmax), "netcdf_get_old_t_final_system time_dimid" )
         call check( nf90_inquire_dimension(nc%id, nc%name_dimid, len=idmax), "netcdf_get_old_t_final_system name_dimid" )
         allocate(vals(idmax))
         call check( nf90_get_var(nc%id, nc%time_varid, rtemp, start=[1], count=[1]), "netcdf_get_old_t_final_system time_varid" )

         !old_t_final = rtemp(1)
         old_t_final = param%t0 ! For NetCDF it is safe to overwrite the final t value on a restart

         if (param%lenergy) then
            call check( nf90_get_var(nc%id, nc%KE_orb_varid, rtemp, start=[1], count=[1]), "netcdf_get_old_t_final_system KE_orb_varid" )
            KE_orb_orig = rtemp(1)

            call check( nf90_get_var(nc%id, nc%KE_spin_varid, rtemp, start=[1], count=[1]), "netcdf_get_old_t_final_system KE_spin_varid" )
            KE_spin_orig = rtemp(1)

            call check( nf90_get_var(nc%id, nc%PE_varid, rtemp, start=[1], count=[1]), "netcdf_get_old_t_final_system PE_varid" )
            PE_orig = rtemp(1)

            call check( nf90_get_var(nc%id, nc%Ecollisions_varid, self%Ecollisions, start=[1]), "netcdf_get_old_t_final_system Ecollisions_varid" )
            call check( nf90_get_var(nc%id, nc%Euntracked_varid,  self%Euntracked,  start=[1]), "netcdf_get_old_t_final_system Euntracked_varid" )

            self%Eorbit_orig = KE_orb_orig + KE_spin_orig + PE_orig + self%Ecollisions + self%Euntracked

            call check( nf90_get_var(nc%id, nc%L_orb_varid, self%Lorbit_orig(:), start=[1,1], count=[NDIM,1]), "netcdf_get_old_t_final_system L_orb_varid" )
            call check( nf90_get_var(nc%id, nc%L_spin_varid, self%Lspin_orig(:), start=[1,1], count=[NDIM,1]), "netcdf_get_old_t_final_system L_spin_varid" )
            call check( nf90_get_var(nc%id, nc%L_escape_varid, self%Lescape(:),  start=[1,1], count=[NDIM,1]), "netcdf_get_old_t_final_system L_escape_varid" )

            self%Ltot_orig(:) = self%Lorbit_orig(:) + self%Lspin_orig(:) + self%Lescape(:)

            call check( nf90_get_var(nc%id, nc%Gmass_varid, vals, start=[1,1], count=[idmax,1]), "netcdf_get_old_t_final_system Gmass_varid" )
            call check( nf90_get_var(nc%id, nc%GMescape_varid,    self%GMescape,    start=[1]), "netcdf_get_old_t_final_system GMescape_varid" )
            self%GMtot_orig = vals(1) + sum(vals(2:idmax), vals(2:idmax) == vals(2:idmax)) + self%GMescape

            select type(cb => self%cb)
            class is (symba_cb)
               cb%GM0 = vals(1)
               cb%dGM = cb%Gmass - cb%GM0

               call check( nf90_get_var(nc%id, nc%radius_varid, rtemp, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system radius_varid" )
               cb%R0 = rtemp(1) 

               if (param%lrotation) then

                  call check( nf90_get_var(nc%id, nc%rot_varid, rot0, start=[1,1,1], count=[NDIM,1,1]), "netcdf_get_old_t_final_system rot_varid" )
                  call check( nf90_get_var(nc%id, nc%Ip_varid, Ip0, start=[1,1,1], count=[NDIM,1,1]), "netcdf_get_old_t_final_system Ip_varid" )

                  cb%L0(:) = Ip0(3) * cb%GM0 * cb%R0**2 * rot0(:)

                  Lnow(:) = cb%Ip(3) * cb%Gmass * cb%radius**2 * cb%rot(:)
                  cb%dL(:) = Lnow(:) - cb%L0(:)
               end if
            end select

         end if

         deallocate(vals)
      end associate
      
      return
   end function netcdf_get_old_t_final_system


   module subroutine netcdf_initialize_output(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Initialize a NetCDF file system and defines all variables.
      use, intrinsic :: ieee_arithmetic
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param   !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: nvar, varid, vartype
      real(DP) :: dfill
      real(SP) :: sfill
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B) :: ndims

      associate(nc => self)

         dfill = ieee_value(dfill, IEEE_QUIET_NAN)
         sfill = ieee_value(sfill, IEEE_QUIET_NAN)

         select case (param%out_type)
         case("NETCDF_FLOAT")
            nc%out_type = NF90_FLOAT
         case("NETCDF_DOUBLE")
            nc%out_type = NF90_DOUBLE
         end select

         ! Check if the file exists, and if it does, delete it
         inquire(file=nc%file_name, exist=fileExists)
         if (fileExists) then
            open(unit=LUN, file=nc%file_name, status="old", err=667, iomsg=errmsg)
            close(unit=LUN, status="delete")
         end if

         ! Create the file
         call check( nf90_create(nc%file_name, NF90_NETCDF4, nc%id), "netcdf_initialize_output nf90_create" )

         ! Dimensions
         call check( nf90_def_dim(nc%id, nc%time_dimname, NF90_UNLIMITED, nc%time_dimid), "netcdf_initialize_output nf90_def_dim time_dimid" ) ! Simulation time dimension
         call check( nf90_def_dim(nc%id, nc%space_dimname, NDIM, nc%space_dimid), "netcdf_initialize_output nf90_def_dim space_dimid" )           ! 3D space dimension
         call check( nf90_def_dim(nc%id, nc%name_dimname, NF90_UNLIMITED, nc%name_dimid), "netcdf_initialize_output nf90_def_dim name_dimid" )       ! dimension to store particle id numbers
         call check( nf90_def_dim(nc%id, nc%str_dimname, NAMELEN, nc%str_dimid), "netcdf_initialize_output nf90_def_dim str_dimid"  )          ! Dimension for string variables (aka character arrays)

         ! Dimension coordinates
         call check( nf90_def_var(nc%id, nc%time_dimname, nc%out_type, nc%time_dimid, nc%time_varid), "netcdf_initialize_output nf90_def_var time_varid"  )
         call check( nf90_def_var(nc%id, nc%space_dimname, NF90_CHAR, nc%space_dimid, nc%space_varid), "netcdf_initialize_output nf90_def_var space_varid"  )
         call check( nf90_def_var(nc%id, nc%name_dimname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], nc%name_varid), "netcdf_initialize_output nf90_def_var name_varid"  )

         ! Variables
         call check( nf90_def_var(nc%id, nc%id_varname, NF90_INT, nc%name_dimid, nc%id_varid), "netcdf_initialize_output nf90_def_var id_varid"  )
         call check( nf90_def_var(nc%id, nc%npl_varname, NF90_INT, nc%time_dimid, nc%npl_varid), "netcdf_initialize_output nf90_def_var npl_varid"  )
         call check( nf90_def_var(nc%id, nc%ntp_varname, NF90_INT, nc%time_dimid, nc%ntp_varid), "netcdf_initialize_output nf90_def_var ntp_varid"  )
         if (param%integrator == SYMBA) call check( nf90_def_var(nc%id, nc%nplm_varname, NF90_INT, nc%time_dimid, nc%nplm_varid), "netcdf_initialize_output nf90_def_var nplm_varid"  )
         call check( nf90_def_var(nc%id, nc%ptype_varname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], nc%ptype_varid), "netcdf_initialize_output nf90_def_var ptype_varid"  )

         if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
            call check( nf90_def_var(nc%id, nc%rh_varname,  nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%rh_varid), "netcdf_initialize_output nf90_def_var rh_varid"  )
            call check( nf90_def_var(nc%id, nc%vh_varname,  nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%vh_varid), "netcdf_initialize_output nf90_def_var vh_varid"  )

            !! When GR is enabled, we need to save the pseudovelocity vectors in addition to the true heliocentric velocity vectors, otherwise
            !! we cannnot expect bit-identical runs from restarted runs with GR enabled due to floating point errors during the conversion.
            if (param%lgr) then
               call check( nf90_def_var(nc%id, nc%gr_pseudo_vh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%gr_pseudo_vh_varid), "netcdf_initialize_output nf90_def_var gr_psuedo_vh_varid"  )
               nc%lpseudo_vel_exists = .true.
            end if

         end if
      
         if ((param%out_form == "EL") .or. (param%out_form == "XVEL")) then
            call check( nf90_def_var(nc%id, nc%a_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%a_varid), "netcdf_initialize_output nf90_def_var a_varid"  )
            call check( nf90_def_var(nc%id, nc%e_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%e_varid), "netcdf_initialize_output nf90_def_var e_varid"  )
            call check( nf90_def_var(nc%id, nc%inc_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%inc_varid), "netcdf_initialize_output nf90_def_var inc_varid"  )
            call check( nf90_def_var(nc%id, nc%capom_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%capom_varid), "netcdf_initialize_output nf90_def_var capom_varid"  )
            call check( nf90_def_var(nc%id, nc%omega_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%omega_varid), "netcdf_initialize_output nf90_def_var omega_varid"  )
            call check( nf90_def_var(nc%id, nc%capm_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%capm_varid), "netcdf_initialize_output nf90_def_var capm_varid"  )
            call check( nf90_def_var(nc%id, nc%varpi_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%varpi_varid), "netcdf_initialize_output nf90_def_var varpi_varid"  )
            call check( nf90_def_var(nc%id, nc%lam_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%lam_varid), "netcdf_initialize_output nf90_def_var lam_varid"  )
            call check( nf90_def_var(nc%id, nc%f_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%f_varid), "netcdf_initialize_output nf90_def_var f_varid"  )
            call check( nf90_def_var(nc%id, nc%cape_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%cape_varid), "netcdf_initialize_output nf90_def_var cape_varid"  )
         end if

         call check( nf90_def_var(nc%id, nc%gmass_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%Gmass_varid), "netcdf_initialize_output nf90_def_var Gmass_varid"  )

         if (param%lrhill_present) then
            call check( nf90_def_var(nc%id, nc%rhill_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%rhill_varid), "netcdf_initialize_output nf90_def_var rhill_varid"  )
         end if

         if (param%lclose) then
            call check( nf90_def_var(nc%id, nc%radius_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%radius_varid), "netcdf_initialize_output nf90_def_var radius_varid"  )

            call check( nf90_def_var(nc%id, nc%origin_time_varname, nc%out_type, nc%name_dimid, nc%origin_time_varid), "netcdf_initialize_output nf90_def_var origin_time_varid"  )
            call check( nf90_def_var(nc%id, nc%origin_type_varname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], &
                                    nc%origin_type_varid), "netcdf_initialize_output nf90_create"  )
            call check( nf90_def_var(nc%id, nc%origin_rh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], nc%origin_rh_varid), "netcdf_initialize_output nf90_def_var origin_rh_varid"  )
            call check( nf90_def_var(nc%id, nc%origin_vh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], nc%origin_vh_varid), "netcdf_initialize_output nf90_def_var origin_vh_varid"  )

            call check( nf90_def_var(nc%id, nc%collision_id_varname, NF90_INT, nc%name_dimid, nc%collision_id_varid), "netcdf_initialize_output nf90_def_var collision_id_varid"  )
            call check( nf90_def_var(nc%id, nc%discard_time_varname, nc%out_type, nc%name_dimid, nc%discard_time_varid), "netcdf_initialize_output nf90_def_var discard_time_varid"  )
            call check( nf90_def_var(nc%id, nc%discard_rh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], nc%discard_rh_varid), "netcdf_initialize_output nf90_def_var discard_rh_varid"  )
            call check( nf90_def_var(nc%id, nc%discard_vh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], nc%discard_vh_varid), "netcdf_initialize_output nf90_def_var discard_vh_varid"  )
            call check( nf90_def_var(nc%id, nc%discard_body_id_varname, NF90_INT, nc%name_dimid, nc%discard_body_id_varid), "netcdf_initialize_output nf90_def_var discard_body_id_varid"  )
         end if

         if (param%lrotation) then
            call check( nf90_def_var(nc%id, nc%Ip_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%Ip_varid), "netcdf_initialize_output nf90_def_var Ip_varid"  )
            call check( nf90_def_var(nc%id, nc%rot_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%rot_varid), "netcdf_initialize_output nf90_def_var rot_varid"  )
         end if

         ! if (param%ltides) then
         !    call check( nf90_def_var(nc%id, nc%k2_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%k2_varid), "netcdf_initialize_output nf90_def_var k2_varid"  )
         !    call check( nf90_def_var(nc%id, nc%q_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%Q_varid), "netcdf_initialize_output nf90_def_var Q_varid"  )
         ! end if

         if (param%lenergy) then
            call check( nf90_def_var(nc%id, nc%ke_orb_varname, nc%out_type, nc%time_dimid, nc%KE_orb_varid), "netcdf_initialize_output nf90_def_var KE_orb_varid"  )
            call check( nf90_def_var(nc%id, nc%ke_spin_varname, nc%out_type, nc%time_dimid, nc%KE_spin_varid), "netcdf_initialize_output nf90_def_var KE_spin_varid"  )
            call check( nf90_def_var(nc%id, nc%pe_varname, nc%out_type, nc%time_dimid, nc%PE_varid), "netcdf_initialize_output nf90_def_var PE_varid"  )
            call check( nf90_def_var(nc%id, nc%L_orb_varname, nc%out_type, [nc%space_dimid, nc%time_dimid], nc%L_orb_varid), "netcdf_initialize_output nf90_def_var L_orb_varid"  )
            call check( nf90_def_var(nc%id, nc%L_spin_varname, nc%out_type, [nc%space_dimid, nc%time_dimid], nc%L_spin_varid), "netcdf_initialize_output nf90_def_var L_spin_varid"  )
            call check( nf90_def_var(nc%id, nc%L_escape_varname, nc%out_type, [nc%space_dimid, nc%time_dimid], nc%L_escape_varid), "netcdf_initialize_output nf90_def_var L_escape_varid"  )
            call check( nf90_def_var(nc%id, nc%Ecollisions_varname, nc%out_type, nc%time_dimid, nc%Ecollisions_varid), "netcdf_initialize_output nf90_def_var Ecollisions_varid"  )
            call check( nf90_def_var(nc%id, nc%Euntracked_varname, nc%out_type, nc%time_dimid, nc%Euntracked_varid), "netcdf_initialize_output nf90_def_var Euntracked_varid"  )
            call check( nf90_def_var(nc%id, nc%GMescape_varname, nc%out_type, nc%time_dimid, nc%GMescape_varid), "netcdf_initialize_output nf90_def_var GMescape_varid"  )
         end if

         call check( nf90_def_var(nc%id, nc%j2rp2_varname, nc%out_type, nc%time_dimid, nc%j2rp2_varid), "netcdf_initialize_output nf90_def_var j2rp2_varid"  )
         call check( nf90_def_var(nc%id, nc%j4rp4_varname, nc%out_type, nc%time_dimid, nc%j4rp4_varid), "netcdf_initialize_output nf90_def_var j4rp4_varid"  )


         ! Set fill mode to NaN for all variables
         call check( nf90_inquire(nc%id, nVariables=nvar), "netcdf_initialize_output nf90_inquire nVariables"  )
         do varid = 1, nvar
            call check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), "netcdf_initialize_output nf90_inquire_variable"  )
            select case(vartype)
            case(NF90_INT)
               call check( nf90_def_var_fill(nc%id, varid, 0, NF90_FILL_INT), "netcdf_initialize_output nf90_def_var_fill NF90_INT"  )
            case(NF90_FLOAT)
               call check( nf90_def_var_fill(nc%id, varid, 0, sfill), "netcdf_initialize_output nf90_def_var_fill NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call check( nf90_def_var_fill(nc%id, varid, 0, dfill), "netcdf_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
            case(NF90_CHAR)
               call check( nf90_def_var_fill(nc%id, varid, 0, 0), "netcdf_initialize_output nf90_def_var_fill NF90_CHAR"  )
            end select
         end do

         ! Take the file out of define mode
         call check( nf90_enddef(nc%id), "netcdf_initialize_output nf90_enddef"  )

         ! Add in the space dimension coordinates
         call check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), "netcdf_initialize_output nf90_put_var space"  )

      end associate
      return

      667 continue
      write(*,*) "Error creating NetCDF output file. " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine netcdf_initialize_output


   module subroutine netcdf_open(self, param, readonly)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Opens a NetCDF file and does the variable inquiries to activate variable ids
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param    !! Current run configuration parameters
      logical, optional,          intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      ! Internals
      integer(I4B) :: mode, status
      character(len=STRMAX) :: errmsg

      mode = NF90_WRITE
      if (present(readonly)) then
         if (readonly) mode = NF90_NOWRITE
      end if

      associate(nc => self)

         write(errmsg,*) "netcdf_open nf90_open ",trim(adjustl(nc%file_name))
         call check( nf90_open(nc%file_name, mode, nc%id), errmsg)

         ! Dimensions
         call check( nf90_inq_dimid(nc%id, nc%time_dimname, nc%time_dimid), "netcdf_open nf90_inq_dimid time_dimid"  )
         call check( nf90_inq_dimid(nc%id, nc%space_dimname, nc%space_dimid), "netcdf_open nf90_inq_dimid space_dimid"  )
         call check( nf90_inq_dimid(nc%id, nc%name_dimname, nc%name_dimid), "netcdf_open nf90_inq_dimid name_dimid"  )
         call check( nf90_inq_dimid(nc%id, nc%str_dimname, nc%str_dimid), "netcdf_open nf90_inq_dimid str_dimid"  )

         ! Dimension coordinates
         call check( nf90_inq_varid(nc%id, nc%time_dimname, nc%time_varid), "netcdf_open nf90_inq_varid time_varid" )
         call check( nf90_inq_varid(nc%id, nc%space_dimname, nc%space_varid), "netcdf_open nf90_inq_varid space_varid" )
         call check( nf90_inq_varid(nc%id, nc%name_dimname, nc%name_varid), "netcdf_open nf90_inq_varid name_varid" )

         ! Required Variables
         call check( nf90_inq_varid(nc%id, nc%id_varname, nc%id_varid), "netcdf_open nf90_inq_varid name_varid" )
         call check( nf90_inq_varid(nc%id, nc%gmass_varname, nc%Gmass_varid), "netcdf_open nf90_inq_varid Gmass_varid" )

         if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
            call check( nf90_inq_varid(nc%id, nc%rh_varname, nc%rh_varid), "netcdf_open nf90_inq_varid rh_varid" )
            call check( nf90_inq_varid(nc%id, nc%vh_varname, nc%vh_varid), "netcdf_open nf90_inq_varid vh_varid" )

            if (param%lgr) then
               !! check if pseudovelocity vectors exist in this file. If they are, set the correct flag so we know whe should not do the conversion.
               status = nf90_inq_varid(nc%id, nc%gr_pseudo_vh_varname, nc%gr_pseudo_vh_varid)
               nc%lpseudo_vel_exists = (status == nf90_noerr)
               if (param%lrestart .and. .not.nc%lpseudo_vel_exists) then
                  write(*,*) "Warning! Pseudovelocity not found in input file for GR enabled run. If this is a restarted run, bit-identical trajectories are not guarunteed!"
               end if

            end if
         end if

         if ((param%out_form == "EL") .or. (param%out_form == "XVEL")) then
            call check( nf90_inq_varid(nc%id, nc%a_varname, nc%a_varid), "netcdf_open nf90_inq_varid a_varid" )
            call check( nf90_inq_varid(nc%id, nc%e_varname, nc%e_varid), "netcdf_open nf90_inq_varid e_varid" )
            call check( nf90_inq_varid(nc%id, nc%inc_varname, nc%inc_varid), "netcdf_open nf90_inq_varid inc_varid" )
            call check( nf90_inq_varid(nc%id, nc%capom_varname, nc%capom_varid), "netcdf_open nf90_inq_varid capom_varid" )
            call check( nf90_inq_varid(nc%id, nc%omega_varname, nc%omega_varid), "netcdf_open nf90_inq_varid omega_varid" )
            call check( nf90_inq_varid(nc%id, nc%capm_varname, nc%capm_varid), "netcdf_open nf90_inq_varid capm_varid" )
         end if

         if (param%lclose) then
            call check( nf90_inq_varid(nc%id, nc%radius_varname, nc%radius_varid), "netcdf_open nf90_inq_varid radius_varid" )
         end if 

         if (param%lrotation) then
            call check( nf90_inq_varid(nc%id, nc%Ip_varname, nc%Ip_varid), "netcdf_open nf90_inq_varid Ip_varid" )
            call check( nf90_inq_varid(nc%id, nc%rot_varname, nc%rot_varid), "netcdf_open nf90_inq_varid rot_varid" )
         end if

         ! if (param%ltides) then
         !    call check( nf90_inq_varid(nc%id, nc%k2_varname, nc%k2_varid), "netcdf_open nf90_inq_varid k2_varid" )
         !    call check( nf90_inq_varid(nc%id, nc%q_varname, nc%Q_varid), "netcdf_open nf90_inq_varid Q_varid" )
         ! end if

         ! Optional Variables
         if (param%lrhill_present) then
            status = nf90_inq_varid(nc%id, nc%rhill_varname, nc%rhill_varid)
            if (status /= nf90_noerr) write(*,*) "Warning! RHILL variable not set in input file. Calculating."
         end if

         ! Optional variables The User Doesn't Need to Know About
         status = nf90_inq_varid(nc%id, nc%npl_varname, nc%npl_varid)
         status = nf90_inq_varid(nc%id, nc%ntp_varname, nc%ntp_varid)
         status = nf90_inq_varid(nc%id, nc%j2rp2_varname, nc%j2rp2_varid)
         status = nf90_inq_varid(nc%id, nc%j4rp4_varname, nc%j4rp4_varid)
         status = nf90_inq_varid(nc%id, nc%ptype_varname, nc%ptype_varid)
         status = nf90_inq_varid(nc%id, nc%varpi_varname, nc%varpi_varid)
         status = nf90_inq_varid(nc%id, nc%lam_varname, nc%lam_varid)
         status = nf90_inq_varid(nc%id, nc%f_varname, nc%f_varid)
         status = nf90_inq_varid(nc%id, nc%cape_varname, nc%cape_varid)

         if (param%integrator == SYMBA) then
            status = nf90_inq_varid(nc%id, nc%nplm_varname, nc%nplm_varid)
         end if

         if (param%lclose) then
            status = nf90_inq_varid(nc%id, nc%origin_type_varname, nc%origin_type_varid)
            status = nf90_inq_varid(nc%id, nc%origin_time_varname, nc%origin_time_varid)
            status = nf90_inq_varid(nc%id, nc%origin_rh_varname, nc%origin_rh_varid)
            status = nf90_inq_varid(nc%id, nc%origin_vh_varname, nc%origin_vh_varid)
            status = nf90_inq_varid(nc%id, nc%collision_id_varname, nc%collision_id_varid)
            status = nf90_inq_varid(nc%id, nc%discard_time_varname, nc%discard_time_varid)
            status = nf90_inq_varid(nc%id, nc%discard_rh_varname, nc%discard_rh_varid)
            status = nf90_inq_varid(nc%id, nc%discard_vh_varname, nc%discard_vh_varid)
            status = nf90_inq_varid(nc%id, nc%discard_body_id_varname, nc%discard_body_id_varid)
         end if

         if (param%lenergy) then
            status = nf90_inq_varid(nc%id, nc%ke_orb_varname, nc%KE_orb_varid)
            status = nf90_inq_varid(nc%id, nc%ke_spin_varname, nc%KE_spin_varid)
            status = nf90_inq_varid(nc%id, nc%pe_varname, nc%PE_varid)
            status = nf90_inq_varid(nc%id, nc%L_orb_varname, nc%L_orb_varid)
            status = nf90_inq_varid(nc%id, nc%L_spin_varname, nc%L_spin_varid)
            status = nf90_inq_varid(nc%id, nc%L_escape_varname, nc%L_escape_varid)
            status = nf90_inq_varid(nc%id, nc%Ecollisions_varname, nc%Ecollisions_varid)
            status = nf90_inq_varid(nc%id, nc%Euntracked_varname, nc%Euntracked_varid)
            status = nf90_inq_varid(nc%id, nc%GMescape_varname, nc%GMescape_varid)
         end if

      end associate

      return
   end subroutine netcdf_open


   module function netcdf_read_frame_system(self, nc, param) result(ierr)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read a frame (header plus records for each massive body and active test particle) from an output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
      class(netcdf_parameters),     intent(inout) :: nc    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      ! Return
      integer(I4B)                                :: ierr  !! Error code: returns 0 if the read is successful
      ! Internals
      integer(I4B)                              :: i, tslot, idmax, npl_check, ntp_check, nplm_check, t_max, str_max, status
      real(DP), dimension(:), allocatable       :: rtemp
      real(DP), dimension(:,:), allocatable     :: vectemp
      integer(I4B), dimension(:), allocatable   :: itemp
      logical, dimension(:), allocatable        :: validmask, tpmask, plmask

      tslot = param%ioutput

      call nc%open(param, readonly=.true.)
      call self%read_hdr(nc, param)

      associate(cb => self%cb, pl => self%pl, tp => self%tp, npl => self%pl%nbody, ntp => self%tp%nbody)

         call pl%setup(npl, param)
         call tp%setup(ntp, param)

         call check( nf90_inquire_dimension(nc%id, nc%name_dimid, len=idmax), "netcdf_read_frame_system nf90_inquire_dimension name_dimid"  )
         allocate(rtemp(idmax))
         allocate(vectemp(NDIM,idmax))
         allocate(itemp(idmax))
         allocate(validmask(idmax))
         allocate(tpmask(idmax))
         allocate(plmask(idmax))
         call check( nf90_inquire_dimension(nc%id, nc%time_dimid, len=t_max), "netcdf_read_frame_system nf90_inquire_dimension time_dimid"  )
         call check( nf90_inquire_dimension(nc%id, nc%str_dimid, len=str_max), "netcdf_read_frame_system nf90_inquire_dimension str_dimid"  )

         ! First filter out only the id slots that contain valid bodies
         if (param%in_form == "XV") then
            call check( nf90_get_var(nc%id, nc%rh_varid, vectemp(:,:), start=[1, 1, tslot]), "netcdf_read_frame_system filter pass nf90_getvar rh_varid"  )
            validmask(:) = vectemp(1,:) == vectemp(1,:)
         else
            call check( nf90_get_var(nc%id, nc%a_varid, rtemp(:), start=[1, tslot]), "netcdf_read_frame_system filter pass nf90_getvar a_varid"  )
            validmask(:) = rtemp(:) == rtemp(:)
         end if

         ! Next, filter only bodies that don't have mass (test particles)
         call check( nf90_get_var(nc%id, nc%Gmass_varid, rtemp(:), start=[1, tslot]), "netcdf_read_frame_system nf90_getvar tp finder Gmass_varid"  )
         plmask(:) = rtemp(:) == rtemp(:) .and. validmask(:)
         tpmask(:) = .not. plmask(:) .and. validmask(:)
         plmask(1) = .false. ! This is the central body

         ! Check to make sure the number of bodies is correct
         npl_check = count(plmask(:))
         ntp_check = count(tpmask(:))

         if (npl_check /= npl) then
            write(*,*) "Error reading in NetCDF file: The recorded value of npl does not match the number of active massive bodies"
            call util_exit(failure)
         end if

         if (ntp_check /= ntp) then
            write(*,*) "Error reading in NetCDF file: The recorded value of ntp does not match the number of active test particles"
            call util_exit(failure)
         end if

         select type (pl)
         class is (symba_pl)
            select type (param)
            class is (symba_parameters)
               nplm_check = count(pack(rtemp,plmask) > param%GMTINY )
               if (nplm_check /= pl%nplm) then
                  write(*,*) "Error reading in NetCDF file: The recorded value of nplm does not match the number of active fully interacting massive bodies"
                  call util_exit(failure)
               end if
            end select
         end select

         ! Now read in each variable and split the outputs by body type
         if ((param%in_form == "XV") .or. (param%in_form == "XVEL")) then
            call check( nf90_get_var(nc%id, nc%rh_varid, vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), "netcdf_read_frame_system nf90_getvar rh_varid"  )
            do i = 1, NDIM
               if (npl > 0) pl%rh(i,:) = pack(vectemp(i,:), plmask(:))
               if (ntp > 0) tp%rh(i,:) = pack(vectemp(i,:), tpmask(:))
            end do

            if (param%lgr .and. nc%lpseudo_vel_exists) then
               call check( nf90_get_var(nc%id, nc%gr_pseudo_vh_varid, vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), "netcdf_read_frame_system nf90_getvar gr_pseudo_vh_varid"  )
               do i = 1, NDIM
                  if (npl > 0) pl%vh(i,:) = pack(vectemp(i,:), plmask(:))
                  if (ntp > 0) tp%vh(i,:) = pack(vectemp(i,:), tpmask(:))
               end do
            else
               call check( nf90_get_var(nc%id, nc%vh_varid, vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), "netcdf_read_frame_system nf90_getvar vh_varid"  )
               do i = 1, NDIM
                  if (npl > 0) pl%vh(i,:) = pack(vectemp(i,:), plmask(:))
                  if (ntp > 0) tp%vh(i,:) = pack(vectemp(i,:), tpmask(:))
               end do
            end if
         end if

         if ((param%in_form == "EL")  .or. (param%in_form == "XVEL")) then
            call check( nf90_get_var(nc%id, nc%a_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar a_varid"  )
            if (.not.allocated(pl%a)) allocate(pl%a(npl))
            if (.not.allocated(tp%a)) allocate(tp%a(ntp))
            if (npl > 0) pl%a(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%a(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(nc%id, nc%e_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar e_varid"  )
            if (.not.allocated(pl%e)) allocate(pl%e(npl))
            if (.not.allocated(tp%e)) allocate(tp%e(ntp))
            if (npl > 0) pl%e(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%e(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(nc%id, nc%inc_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar inc_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%inc)) allocate(pl%inc(npl))
            if (.not.allocated(tp%inc)) allocate(tp%inc(ntp))
            if (npl > 0) pl%inc(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%inc(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(nc%id, nc%capom_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar capom_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%capom)) allocate(pl%capom(npl))
            if (.not.allocated(tp%capom)) allocate(tp%capom(ntp))
            if (npl > 0) pl%capom(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%capom(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(nc%id, nc%omega_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar omega_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%omega)) allocate(pl%omega(npl))
            if (.not.allocated(tp%omega)) allocate(tp%omega(ntp))
            if (npl > 0) pl%omega(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%omega(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(nc%id, nc%capm_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar capm_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%capm)) allocate(pl%capm(npl))
            if (.not.allocated(tp%capm)) allocate(tp%capm(ntp))
            if (npl > 0) pl%capm(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%capm(:) = pack(rtemp, tpmask)

         end if
      
         call check( nf90_get_var(nc%id, nc%Gmass_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar Gmass_varid"  )
         cb%Gmass = rtemp(1)
         cb%mass = cb%Gmass / param%GU

         ! Set initial central body mass for Helio bookkeeping
         select type(cb)
            class is (symba_cb)
               cb%GM0 = cb%Gmass
         end select
            

         if (npl > 0) then
            pl%Gmass(:) = pack(rtemp, plmask)
            pl%mass(:) = pl%Gmass(:) / param%GU

            if (param%lrhill_present) then 
               call check( nf90_get_var(nc%id, nc%rhill_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar rhill_varid"  )
               pl%rhill(:) = pack(rtemp, plmask)
            end if
         end if

         if (param%lclose) then
            call check( nf90_get_var(nc%id, nc%radius_varid, rtemp, start=[1, tslot], count=[idmax,1]), "netcdf_read_frame_system nf90_getvar radius_varid"  )
            cb%radius = rtemp(1)

            ! Set initial central body radius for SyMBA bookkeeping
            select type(cb)
               class is (symba_cb)
                  cb%R0 = cb%radius
            end select
            if (npl > 0) pl%radius(:) = pack(rtemp, plmask)
         else
            cb%radius = param%rmin
            if (npl > 0) pl%radius(:) = 0.0_DP
         end if

         if (param%lrotation) then
            call check( nf90_get_var(nc%id, nc%Ip_varid,  vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), "netcdf_read_frame_system nf90_getvar Ip_varid"  )
            cb%Ip(:) = vectemp(:,1)
            do i = 1, NDIM
               if (npl > 0) pl%Ip(i,:) = pack(vectemp(i,:), plmask(:))
            end do

            call check( nf90_get_var(nc%id, nc%rot_varid, vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), "netcdf_read_frame_system nf90_getvar rot_varid"  )
            cb%rot(:) = vectemp(:,1)
            do i = 1, NDIM
               if (npl > 0) pl%rot(i,:) = pack(vectemp(i,:), plmask(:))
            end do

            ! Set initial central body angular momentum for Helio bookkeeping
            select type(cb)
               class is (symba_cb)
                  cb%L0(:) = cb%Ip(3) * cb%GM0 * cb%R0**2 * cb%rot(:)         
            end select
         end if

         ! if (param%ltides) then
         !    call check( nf90_get_var(nc%id, nc%k2_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar k2_varid"  )
         !    cb%k2 = rtemp(1)
         !    if (npl > 0) pl%k2(:) = pack(rtemp, plmask)

         !    call check( nf90_get_var(nc%id, nc%Q_varid,  rtemp,  start=[1, tslot]), "netcdf_read_frame_system nf90_getvar Q_varid"  )
         !    cb%Q = rtemp(1)
         !    if (npl > 0) pl%Q(:) = pack(rtemp, plmask)
         ! end if

         status = nf90_inq_varid(nc%id, nc%j2rp2_varname, nc%j2rp2_varid)
         if (status == nf90_noerr) then
            call check( nf90_get_var(nc%id, nc%j2rp2_varid, cb%j2rp2, start=[tslot]), "netcdf_read_frame_system nf90_getvar j2rp2_varid"  )
         else 
            cb%j2rp2 = 0.0_DP
         end if

         status = nf90_inq_varid(nc%id, nc%j4rp4_varname, nc%j4rp4_varid)   
         if (status == nf90_noerr) then      
            call check( nf90_get_var(nc%id, nc%j4rp4_varid, cb%j4rp4, start=[tslot]), "netcdf_read_frame_system nf90_getvar j4rp4_varid"  )
         else 
            cb%j4rp4 = 0.0_DP
         end if

         call self%read_particle_info(nc, param, plmask, tpmask) 

         if (param%in_form == "EL") then
            call pl%el2xv(cb)
            call tp%el2xv(cb)
         end if
         ! if this is a GR-enabled run, check to see if we got the pseudovelocities in. Otherwise, we'll need to generate them.
         if (param%lgr .and. .not.(nc%lpseudo_vel_exists)) then
            call pl%set_mu(cb)
            call tp%set_mu(cb)
            call pl%v2pv(param)
            call tp%v2pv(param)
         end if
         
      end associate

      call nc%close()

      ierr = 0
      return

      667 continue
      write(*,*) "Error reading system frame in netcdf_read_frame_system"

   end function netcdf_read_frame_system


   module subroutine netcdf_read_hdr_system(self, nc, param) 
      !! author: David A. Minton
      !!
      !! Reads header information (variables that change with time, but not particle id). 
      !! This subroutine significantly improves the output over the original binary file, allowing us to track energy, momentum, and other quantities that 
      !! previously were handled as separate output files.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
      class(netcdf_parameters),     intent(inout) :: nc      !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: tslot, status, idmax
      real(DP), dimension(:), allocatable       :: gmtemp
      logical, dimension(:), allocatable        :: plmask, tpmask, plmmask


      tslot = param%ioutput
      call check( nf90_inquire_dimension(nc%id, nc%name_dimid, len=idmax), "netcdf_read_hdr_system nf90_inquire_dimension name_dimid"  )
      call check( nf90_get_var(nc%id, nc%time_varid, self%t, start=[tslot]), "netcdf_read_hdr_system nf90_getvar time_varid"  )

      allocate(gmtemp(idmax))
      allocate(tpmask(idmax))
      allocate(plmask(idmax))
      allocate(plmmask(idmax))

      call check( nf90_get_var(nc%id, nc%Gmass_varid, gmtemp, start=[1,1], count=[idmax,1]), "netcdf_read_hdr_system nf90_getvar Gmass_varid"  )

      plmask(:) = gmtemp(:) == gmtemp(:)
      tpmask(:) = .not. plmask(:)
      plmask(1) = .false. ! This is the central body
      select type (param)
      class is (symba_parameters)
         plmmask(:) = plmask(:)
         where(plmask(:))
            plmmask(:) = gmtemp(:) > param%GMTINY
         endwhere
      end select

      status = nf90_inq_varid(nc%id, nc%npl_varname, nc%npl_varid)
      if (status == nf90_noerr) then
         call check( nf90_get_var(nc%id, nc%npl_varid,  self%pl%nbody, start=[tslot]), "netcdf_read_hdr_system nf90_getvar npl_varid"  )
      else
         self%pl%nbody = count(plmask(:))
      end if

      status = nf90_inq_varid(nc%id, nc%ntp_varname, nc%ntp_varid)
      if (status == nf90_noerr) then
         call check( nf90_get_var(nc%id, nc%ntp_varid,  self%tp%nbody, start=[tslot]), "netcdf_read_hdr_system nf90_getvar ntp_varid"  )
      else
         self%tp%nbody = count(tpmask(:))
      end if

      if (param%integrator == SYMBA) then
         status = nf90_inq_varid(nc%id, nc%nplm_varname, nc%nplm_varid)
         select type(pl => self%pl)
         class is (symba_pl)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%nplm_varid,  pl%nplm, start=[tslot]), "netcdf_read_hdr_system nf90_getvar nplm_varid"  )
            else
               pl%nplm = count(plmmask(:))
            end if
         end select
      end if

      if (param%lenergy) then
         status = nf90_inq_varid(nc%id, nc%ke_orb_varname, nc%KE_orb_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%KE_orb_varid,      self%ke_orbit,    start=[tslot]), "netcdf_read_hdr_system nf90_getvar KE_orb_varid"  )
         status = nf90_inq_varid(nc%id, nc%ke_spin_varname, nc%KE_spin_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%KE_spin_varid,     self%ke_spin,     start=[tslot]), "netcdf_read_hdr_system nf90_getvar KE_spin_varid"  )
         status = nf90_inq_varid(nc%id, nc%pe_varname, nc%PE_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%PE_varid,          self%pe,          start=[tslot]), "netcdf_read_hdr_system nf90_getvar PE_varid"  )
         status = nf90_inq_varid(nc%id, nc%L_orb_varname, nc%L_orb_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%L_orb_varid,      self%Lorbit(:),   start=[1,tslot], count=[NDIM,1]), "netcdf_read_hdr_system nf90_getvar L_orb_varid"  )
         status = nf90_inq_varid(nc%id, nc%L_spin_varname, nc%L_spin_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%L_spin_varid,     self%Lspin(:),    start=[1,tslot], count=[NDIM,1]), "netcdf_read_hdr_system nf90_getvar L_spin_varid"  )
         status = nf90_inq_varid(nc%id, nc%L_escape_varname, nc%L_escape_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%L_escape_varid,   self%Lescape(:),  start=[1, tslot], count=[NDIM,1]), "netcdf_read_hdr_system nf90_getvar L_escape_varid"  )
         status = nf90_inq_varid(nc%id, nc%Ecollisions_varname, nc%Ecollisions_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%Ecollisions_varid, self%Ecollisions, start=[tslot]), "netcdf_read_hdr_system nf90_getvar Ecollisions_varid"  )
         status = nf90_inq_varid(nc%id, nc%Euntracked_varname, nc%Euntracked_varid)
         if (status == nf90_noerr) call check( nf90_get_var(nc%id, nc%Euntracked_varid,  self%Euntracked,  start=[tslot]), "netcdf_read_hdr_system nf90_getvar Euntracked_varid"  )
         status = nf90_inq_varid(nc%id, nc%GMescape_varname, nc%GMescape_varid)
         if (status == nf90_noerr)  call check( nf90_get_var(nc%id, nc%GMescape_varid,    self%GMescape,    start=[tslot]), "netcdf_read_hdr_system nf90_getvar GMescape_varid"  )
      end if

      return
   end subroutine netcdf_read_hdr_system


   module subroutine netcdf_read_particle_info_system(self, nc, param, plmask, tpmask)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Reads particle information metadata from file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(netcdf_parameters),     intent(inout) :: nc       !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      logical, dimension(:),        intent(in)    :: plmask !! Logical array indicating which index values belong to massive bodies
      logical, dimension(:),        intent(in)    :: tpmask !! Logical array indicating which index values belong to test particles
      ! Internals
      integer(I4B)                                :: i, idmax, status
      real(DP), dimension(:), allocatable         :: rtemp
      real(DP), dimension(:,:), allocatable       :: vectemp
      integer(I4B), dimension(:), allocatable     :: itemp
      character(len=NAMELEN), dimension(:), allocatable :: ctemp
      integer(I4B), dimension(:), allocatable     :: plind, tpind

      ! This string of spaces of length NAMELEN is used to clear out any old data left behind inside the string variables
      idmax = size(plmask)
      allocate(rtemp(idmax))
      allocate(vectemp(NDIM,idmax))
      allocate(itemp(idmax))
      allocate(ctemp(idmax))

      associate(cb => self%cb, pl => self%pl, tp => self%tp, npl => self%pl%nbody, ntp => self%tp%nbody)

         if (npl > 0) then
            pl%status(:) = ACTIVE
            pl%lmask(:) = .true.
            do i = 1, npl
               call pl%info(i)%set_value(status="ACTIVE")
            end do
            allocate(plind(npl))
            plind(:) = pack([(i, i = 1, idmax)], plmask(:))
         end if
         if (ntp > 0) then
            tp%status(:) = ACTIVE
            tp%lmask(:) = .true.
            do i = 1, ntp
               call tp%info(i)%set_value(status="ACTIVE")
            end do
            allocate(tpind(ntp))
            tpind(:) = pack([(i, i = 1, idmax)], tpmask(:))
         end if

         call check( nf90_get_var(nc%id, nc%id_varid, itemp), "netcdf_read_particle_info_system nf90_getvar id_varid"  )
         cb%id = itemp(1)
         pl%id(:) = pack(itemp, plmask)
         tp%id(:) = pack(itemp, tpmask)
         cb%id = 0
         pl%id(:) = pack([(i,i=0,idmax-1)],plmask)
         tp%id(:) = pack([(i,i=0,idmax-1)],tpmask)

         call check( nf90_get_var(nc%id, nc%name_varid, ctemp, count=[NAMELEN, idmax]), "netcdf_read_particle_info_system nf90_getvar name_varid"  )
         call cb%info%set_value(name=ctemp(1))
         do i = 1, npl
            call pl%info(i)%set_value(name=ctemp(plind(i)))
         end do
         do i = 1, ntp
            call tp%info(i)%set_value(name=ctemp(tpind(i)))
         end do

         status = nf90_get_var(nc%id, nc%ptype_varid, ctemp, count=[NAMELEN, idmax])
         if (status /= nf90_noerr) then ! Set default particle types
            call cb%info%set_value(particle_type=CB_TYPE_NAME)

            ! Handle semi-interacting bodies in SyMBA
            select type(pl)
            class is (symba_pl)
               select type (param)
               class is (symba_parameters)
                  do i = 1, npl
                     if (pl%Gmass(i) < param%GMTINY) then
                        call pl%info(i)%set_value(particle_type=PL_TINY_TYPE_NAME)
                     else
                        call pl%info(i)%set_value(particle_type=PL_TYPE_NAME)
                     end if
                  end do
               end select
            class default ! Non-SyMBA massive bodies
               do i = 1, npl
                  call pl%info(i)%set_value(particle_type=PL_TYPE_NAME)
               end do
            end select
            do i = 1, ntp
               call tp%info(i)%set_value(particle_type=TP_TYPE_NAME)
            end do
         else ! Use particle types defined in input file
            call cb%info%set_value(particle_type=ctemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(particle_type=ctemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(particle_type=ctemp(tpind(i)))
            end do
         end if

         call cb%info%set_value(status="ACTIVE")

         if (param%lclose) then

            status = nf90_inq_varid(nc%id, nc%origin_type_varname, nc%origin_type_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%origin_type_varid, ctemp, count=[NAMELEN, idmax]), "netcdf_read_particle_info_system nf90_getvar origin_type_varid"  )
            else
               ctemp = "Initial Conditions"
            end if

            call cb%info%set_value(origin_type=ctemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(origin_type=ctemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_type=ctemp(tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%origin_time_varname, nc%origin_time_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%origin_time_varid, rtemp), "netcdf_read_particle_info_system nf90_getvar origin_time_varid"  )
            else
               rtemp = param%t0
            end if

            call cb%info%set_value(origin_time=rtemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(origin_time=rtemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_time=rtemp(tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%origin_rh_varname, nc%origin_rh_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%origin_rh_varid, vectemp(:,:)), "netcdf_read_particle_info_system nf90_getvar origin_rh_varid"  )
            else if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
               call check( nf90_get_var(nc%id, nc%rh_varid, vectemp(:,:)), "netcdf_read_particle_info_system nf90_getvar rh_varid"  )
            else 
               vectemp(:,:) = 0._DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(origin_rh=vectemp(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_rh=vectemp(:,tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%origin_vh_varname, nc%origin_vh_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%origin_vh_varid, vectemp(:,:)), "netcdf_read_particle_info_system nf90_getvar origin_vh_varid"  )
            else if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
               call check( nf90_get_var(nc%id, nc%vh_varid, vectemp(:,:)), "netcdf_read_particle_info_system nf90_getvar vh_varid"  )
            else
               vectemp(:,:) = 0._DP
            end if 
            
            do i = 1, npl
               call pl%info(i)%set_value(origin_vh=vectemp(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_vh=vectemp(:,tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%collision_id_varname, nc%collision_id_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%collision_id_varid, itemp), "netcdf_read_particle_info_system nf90_getvar collision_id_varid"  )
            else
               itemp = 0.0_DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(collision_id=itemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(collision_id=itemp(tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%discard_time_varname, nc%discard_time_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%discard_time_varid, rtemp), "netcdf_read_particle_info_system nf90_getvar discard_time_varid"  )
            else
               rtemp = 0.0_DP
            end if 

            call cb%info%set_value(discard_time=rtemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(discard_time=rtemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(discard_time=rtemp(tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%discard_rh_varname, nc%discard_rh_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%discard_rh_varid, vectemp(:,:)), "netcdf_read_particle_info_system nf90_getvar discard_rh_varid"  )
            else
               vectemp(:,:) = 0.0_DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(discard_rh=vectemp(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(discard_rh=vectemp(:,tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%discard_vh_varname, nc%discard_vh_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(nc%id, nc%discard_vh_varid, vectemp(:,:)), "netcdf_read_particle_info_system nf90_getvar discard_vh_varid"  )
            else
               vectemp(:,:) = 0.0_DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(discard_vh=vectemp(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(discard_vh=vectemp(:,tpind(i)))
            end do
         end if

      end associate

      return
   end subroutine netcdf_read_particle_info_system


   module subroutine netcdf_sync(self)
      !! author: David A. Minton
      !!
      !! Syncrhonize the disk and memory buffer of the NetCDF file (e.g. commit the frame files stored in memory to disk) 
      !!    
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset

      call check( nf90_sync(self%id), "netcdf_sync nf90_sync"  )

      return
   end subroutine netcdf_sync


   module subroutine netcdf_write_frame_base(self, nc, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      implicit none
      ! Arguments
      class(swiftest_base),       intent(in)    :: self  !! Swiftest particle object
      class(netcdf_parameters),   intent(inout) :: nc    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, tslot, idslot, old_mode
      integer(I4B), dimension(:), allocatable   :: ind
      real(DP), dimension(NDIM)                 :: vh !! Temporary variable to store heliocentric velocity values when converting from pseudovelocity in GR-enabled runs
      real(DP)                                  :: a, e, inc, omega, capom, capm, varpi, lam, f, cape, capf

      tslot = param%ioutput

      call self%write_info(nc, param)

      call check( nf90_set_fill(nc%id, nf90_nofill, old_mode), "netcdf_write_frame_base nf90_set_fill"  )
      select type(self)
         class is (swiftest_body)
         associate(n => self%nbody)
            if (n == 0) return

            call util_sort(self%id(1:n), ind)

            do i = 1, n
               j = ind(i)
               idslot = self%id(j) + 1

               !! Convert from pseudovelocity to heliocentric without replacing the current value of pseudovelocity 
               if (param%lgr) call gr_pseudovel2vel(param, self%mu(j), self%rh(:, j), self%vh(:, j), vh(:))

               if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
                  call check( nf90_put_var(nc%id, nc%rh_varid, self%rh(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var rh_varid"  )
                  if (param%lgr) then !! Convert from pseudovelocity to heliocentric without replacing the current value of pseudovelocity
                     call check( nf90_put_var(nc%id, nc%vh_varid, vh(:), start=[1,idslot, tslot], count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var vh_varid"  )
                     call check( nf90_put_var(nc%id, nc%gr_pseudo_vh_varid, self%vh(:, j), start=[1,idslot, tslot],count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var gr_pseudo_vhx_varid"  )

                  else
                     call check( nf90_put_var(nc%id, nc%vh_varid, self%vh(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var vh_varid"  )
                  end if
               end if

               if ((param%out_form == "EL") .or. (param%out_form == "XVEL")) then
                  if (param%lgr) then !! For GR-enabled runs, use the true value of velocity computed above
                     call orbel_xv2el(self%mu(j), self%rh(1,j), self%rh(2,j), self%rh(3,j), &
                                       vh(1), vh(2), vh(3), &
                                       a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf)
                  else !! For non-GR runs just convert from the velocity we have
                     call orbel_xv2el(self%mu(j), self%rh(1,j), self%rh(2,j), self%rh(3,j), &
                                       self%vh(1,j), self%vh(2,j), self%vh(3,j), &
                                       a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf)
                  end if
                  call check( nf90_put_var(nc%id, nc%a_varid, a, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body a_varid"  )
                  call check( nf90_put_var(nc%id, nc%e_varid, e, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body e_varid"  )
                  call check( nf90_put_var(nc%id, nc%inc_varid, inc * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body inc_varid"  )
                  call check( nf90_put_var(nc%id, nc%capom_varid, capom * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body capom_varid"  )
                  call check( nf90_put_var(nc%id, nc%omega_varid, omega * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body omega_varid"  )
                  call check( nf90_put_var(nc%id, nc%capm_varid, capm * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body capm_varid"  ) 
                  call check( nf90_put_var(nc%id, nc%varpi_varid, varpi * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body varpi_varid"  ) 
                  call check( nf90_put_var(nc%id, nc%lam_varid, lam * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body lam_varid"  ) 
                  call check( nf90_put_var(nc%id, nc%f_varid, f * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body f_varid"  ) 
                  if (e < 1.0_DP) then
                     call check( nf90_put_var(nc%id, nc%cape_varid, cape * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body cape_varid"  ) 
                  else if (e > 1.0_DP) then
                     call check( nf90_put_var(nc%id, nc%cape_varid, capf * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body (capf) cape_varid"  ) 
                  end if
               end if

               select type(self)  
               class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
                  call check( nf90_put_var(nc%id, nc%Gmass_varid, self%Gmass(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body Gmass_varid"  )
                  if (param%lrhill_present) then
                     call check( nf90_put_var(nc%id, nc%rhill_varid, self%rhill(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body rhill_varid"  )
                  end if
                  if (param%lclose) call check( nf90_put_var(nc%id, nc%radius_varid, self%radius(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body radius_varid"  )
                  if (param%lrotation) then
                     call check( nf90_put_var(nc%id, nc%Ip_varid, self%Ip(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var body Ip_varid"  )
                     call check( nf90_put_var(nc%id, nc%rot_varid, self%rot(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var body rotx_varid"  )
                  end if
                  ! if (param%ltides) then
                  !    call check( nf90_put_var(nc%id, nc%k2_varid, self%k2(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body k2_varid"  )
                  !    call check( nf90_put_var(nc%id, nc%Q_varid, self%Q(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var body Q_varid"  )
                  ! end if

               end select
            end do
         end associate
      class is (swiftest_cb)
         idslot = self%id + 1
         call check( nf90_put_var(nc%id, nc%id_varid, self%id, start=[idslot]), "netcdf_write_frame_base nf90_put_var cb id_varid"  )

         call check( nf90_put_var(nc%id, nc%Gmass_varid, self%Gmass, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb Gmass_varid"  )
         if (param%lclose) call check( nf90_put_var(nc%id, nc%radius_varid, self%radius, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb radius_varid"  )
         call check( nf90_put_var(nc%id, nc%j2rp2_varid, self%j2rp2, start=[tslot]), "netcdf_write_frame_base nf90_put_var cb j2rp2_varid" )
         call check( nf90_put_var(nc%id, nc%j4rp4_varid, self%j4rp4, start=[tslot]), "netcdf_write_frame_base nf90_put_var cb j4rp4_varid" )
         if (param%lrotation) then
            call check( nf90_put_var(nc%id, nc%Ip_varid, self%Ip(:), start=[1, idslot, tslot], count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var cb Ip_varid"  )
            call check( nf90_put_var(nc%id, nc%rot_varid, self%rot(:), start=[1, idslot, tslot], count=[NDIM,1,1]), "netcdf_write_frame_base nf90_put_var cb rot_varid"  )
         end if
         ! if (param%ltides) then
         !    call check( nf90_put_var(nc%id, nc%k2_varid, self%k2, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb k2_varid"  )
         !    call check( nf90_put_var(nc%id, nc%Q_varid, self%Q, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb Q_varid"  )
         ! end if

      end select
      call check( nf90_set_fill(nc%id, old_mode, old_mode), "netcdf_write_frame_base nf90_set_fill old_mode"  )

      return
   end subroutine netcdf_write_frame_base


   module subroutine netcdf_write_frame_system(self, nc, param)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to a output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
      class(netcdf_parameters),     intent(inout) :: nc      !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 

      call self%write_hdr(nc, param)
      call self%cb%write_frame(nc, param)
      call self%pl%write_frame(nc, param)
      call self%tp%write_frame(nc, param)

      return
   end subroutine netcdf_write_frame_system


   module subroutine netcdf_write_info_base(self, nc, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write all current particle to file
      implicit none
      ! Arguments
      class(swiftest_base),       intent(in)    :: self   !! Swiftest particle object
      class(netcdf_parameters),   intent(inout) :: nc       !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, idslot, old_mode
      integer(I4B), dimension(:), allocatable   :: ind
      character(len=:), allocatable             :: charstring

      ! This string of spaces of length NAMELEN is used to clear out any old data left behind inside the string variables
      call check( nf90_set_fill(nc%id, nf90_nofill, old_mode), "netcdf_write_info_base nf90_set_fill nf90_nofill"  )

      select type(self)
         class is (swiftest_body)
         associate(n => self%nbody)
            if (n == 0) return
            call util_sort(self%id(1:n), ind)

            do i = 1, n
               j = ind(i)
               idslot = self%id(j) + 1
               call check( nf90_put_var(nc%id, nc%id_varid, self%id(j), start=[idslot]), "netcdf_write_info_base nf90_put_var id_varid"  )

               charstring = trim(adjustl(self%info(j)%name))
               call check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "netcdf_write_info_base nf90_put_var name_varid"  )

               charstring = trim(adjustl(self%info(j)%particle_type))
               call check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "netcdf_write_info_base nf90_put_var particle_type_varid"  )

               if (param%lclose) then
                  charstring = trim(adjustl(self%info(j)%origin_type))
                  call check( nf90_put_var(nc%id, nc%origin_type_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "netcdf_write_info_base nf90_put_var origin_type_varid"  )
                  call check( nf90_put_var(nc%id, nc%origin_time_varid,  self%info(j)%origin_time,  start=[idslot]), "netcdf_write_info_base nf90_put_var origin_time_varid"  )
                  call check( nf90_put_var(nc%id, nc%origin_rh_varid,    self%info(j)%origin_rh(:), start=[1,idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var origin_rh_varid"  )
                  call check( nf90_put_var(nc%id, nc%origin_vh_varid,    self%info(j)%origin_vh(:), start=[1,idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var origin_vh_varid"  )
   
                  call check( nf90_put_var(nc%id, nc%collision_id_varid, self%info(j)%collision_id, start=[idslot]), "netcdf_write_info_base nf90_put_var collision_id_varid"  )
                  call check( nf90_put_var(nc%id, nc%discard_time_varid, self%info(j)%discard_time, start=[idslot]), "netcdf_write_info_base nf90_put_var discard_time_varid"  )
                  call check( nf90_put_var(nc%id, nc%discard_rh_varid,   self%info(j)%discard_rh(:), start=[1,idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var discard_rh_varid"  )
                  call check( nf90_put_var(nc%id, nc%discard_vh_varid,   self%info(j)%discard_vh(:), start=[1,idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var discard_vh_varid"  )
               end if

            end do
         end associate

      class is (swiftest_cb)
         idslot = self%id + 1
         call check( nf90_put_var(nc%id, nc%id_varid, self%id, start=[idslot]), "netcdf_write_info_base nf90_put_var cb id_varid"  )

         charstring = trim(adjustl(self%info%name))
         call check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "netcdf_write_info_base nf90_put_var cb name_varid"  )

         charstring = trim(adjustl(self%info%particle_type))
         call check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "netcdf_write_info_base nf90_put_var cb ptype_varid"  )

         if (param%lclose) then
            charstring = trim(adjustl(self%info%origin_type))
            call check( nf90_put_var(nc%id, nc%origin_type_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "netcdf_write_info_base nf90_put_var cb origin_type_varid"  )

            call check( nf90_put_var(nc%id, nc%origin_time_varid, self%info%origin_time, start=[idslot]), "netcdf_write_info_base nf90_put_var cb origin_time_varid"  )
            call check( nf90_put_var(nc%id, nc%origin_rh_varid, self%info%origin_rh(:), start=[1, idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var cb origin_rh_varid"  )
            call check( nf90_put_var(nc%id, nc%origin_vh_varid, self%info%origin_vh(:), start=[1, idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var cb origin_vh_varid"  )
   
            call check( nf90_put_var(nc%id, nc%collision_id_varid, self%info%collision_id, start=[idslot]), "netcdf_write_info_base nf90_put_var cb collision_id_varid"  )
            call check( nf90_put_var(nc%id, nc%discard_time_varid, self%info%discard_time, start=[idslot]), "netcdf_write_info_base nf90_put_var cb discard_time_varid"  )
            call check( nf90_put_var(nc%id, nc%discard_rh_varid, self%info%discard_rh(:), start=[1, idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var cb discard_rh_varid"  )
            call check( nf90_put_var(nc%id, nc%discard_vh_varid, self%info%discard_vh(:), start=[1, idslot], count=[NDIM,1]), "netcdf_write_info_base nf90_put_var cb discard_vh_varid"  )
         end if

      end select

      call check( nf90_set_fill(nc%id, old_mode, old_mode) )
      return
   end subroutine netcdf_write_info_base


   module subroutine netcdf_write_hdr_system(self, nc, param) 
      !! author: David A. Minton
      !!
      !! Writes header information (variables that change with time, but not particle id). 
      !! This subroutine significantly improves the output over the original binary file, allowing us to track energy, momentum, and other quantities that 
      !! previously were handled as separate output files.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(in)    :: self  !! Swiftest nbody system object
      class(netcdf_parameters),     intent(inout) :: nc      !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: tslot

      tslot = param%ioutput

      call check( nf90_put_var(nc%id, nc%time_varid, self%t, start=[tslot]), "netcdf_write_hdr_system nf90_put_var time_varid"  )
      call check( nf90_put_var(nc%id, nc%npl_varid, self%pl%nbody, start=[tslot]), "netcdf_write_hdr_system nf90_put_var npl_varid"  )
      call check( nf90_put_var(nc%id, nc%ntp_varid, self%tp%nbody, start=[tslot]), "netcdf_write_hdr_system nf90_put_var ntp_varid"  )
      select type(pl => self%pl)
      class is (symba_pl)
         call check( nf90_put_var(nc%id, nc%nplm_varid, pl%nplm, start=[tslot]), "netcdf_write_hdr_system nf90_put_var nplm_varid"  )
      end select

      if (param%lenergy) then
         call check( nf90_put_var(nc%id, nc%KE_orb_varid, self%ke_orbit, start=[tslot]), "netcdf_write_hdr_system nf90_put_var KE_orb_varid"  )
         call check( nf90_put_var(nc%id, nc%KE_spin_varid, self%ke_spin, start=[tslot]), "netcdf_write_hdr_system nf90_put_var KE_spin_varid"  )
         call check( nf90_put_var(nc%id, nc%PE_varid, self%pe, start=[tslot]), "netcdf_write_hdr_system nf90_put_var PE_varid"  )
         call check( nf90_put_var(nc%id, nc%L_orb_varid, self%Lorbit(:), start=[1,tslot], count=[NDIM,1]), "netcdf_write_hdr_system nf90_put_var L_orb_varid"  )
         call check( nf90_put_var(nc%id, nc%L_spin_varid, self%Lspin(:), start=[1,tslot], count=[NDIM,1]), "netcdf_write_hdr_system nf90_put_var L_spin_varid"  )
         call check( nf90_put_var(nc%id, nc%L_escape_varid, self%Lescape(:), start=[1,tslot], count=[NDIM,1]), "netcdf_write_hdr_system nf90_put_var L_escape_varid"  )
         call check( nf90_put_var(nc%id, nc%Ecollisions_varid, self%Ecollisions, start=[tslot]), "netcdf_write_hdr_system nf90_put_var Ecollisions_varid"  )
         call check( nf90_put_var(nc%id, nc%Euntracked_varid, self%Euntracked, start=[tslot]), "netcdf_write_hdr_system nf90_put_var Euntracked_varid"  )
         call check( nf90_put_var(nc%id, nc%GMescape_varid, self%GMescape, start=[tslot]), "netcdf_write_hdr_system nf90_put_var GMescape_varid"  )
      end if

      return
   end subroutine netcdf_write_hdr_system

end submodule s_netcdf
