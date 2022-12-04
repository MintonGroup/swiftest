!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (encounter_classes) s_encounter_io
   use swiftest
   use netcdf
contains

   module subroutine encounter_io_dump_storage_list(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the time history of an encounter to file.
      implicit none
      ! Arguments
      class(encounter_storage(*)),  intent(inout)        :: self   !! Encounter storage object
      class(swiftest_parameters),   intent(inout)        :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i

      ! Most of this is just temporary test code just to get something working. Eventually this should get cleaned up.
      call self%nciu%initialize(param)
      do i = 1, self%nframes
         if (allocated(self%frame(i)%item)) then
            select type(plplenc_list => self%frame(i)%item)
            class is (symba_plplenc)
               self%nciu%ienc_frame = i
               call plplenc_list%write_frame(self%nciu,param)
            end select
         end if
      end do
      call self%nciu%close()


      return
   end subroutine encounter_io_dump_storage_list


   module subroutine encounter_io_initialize_output(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF encounter file system. This is a simplified version of the main simulation output NetCDF file, but with fewer variables.
      use, intrinsic :: ieee_arithmetic
      implicit none
      ! Arguments
      class(encounter_io_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),     intent(in)    :: param   !! Current run configuration parameters
      ! Internals
      integer(I4B) :: nvar, varid, vartype
      real(DP) :: dfill
      real(SP) :: sfill
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B) :: ndims


      associate(nciu => self)
         dfill = ieee_value(dfill, IEEE_QUIET_NAN)
         sfill = ieee_value(sfill, IEEE_QUIET_NAN)

         select case (param%out_type)
         case("NETCDF_FLOAT")
            self%out_type = NF90_FLOAT
         case("NETCDF_DOUBLE")
            self%out_type = NF90_DOUBLE
         end select


         ! Check if the file exists, and if it does, delete it
         inquire(file=nciu%enc_file, exist=fileExists)
         if (fileExists) then
            open(unit=LUN, file=nciu%enc_file, status="old", err=667, iomsg=errmsg)
            close(unit=LUN, status="delete")
         end if

         call check( nf90_create(nciu%enc_file, NF90_NETCDF4, nciu%id), "encounter_io_initialize_output nf90_create" )

         ! Dimensions
         call check( nf90_def_dim(nciu%id, nciu%time_dimname, NF90_UNLIMITED, nciu%time_dimid), "encounter_io_initialize_output nf90_def_dim time_dimid" ) ! Simulation time dimension
         call check( nf90_def_dim(nciu%id, nciu%space_dimname, NDIM, nciu%space_dimid), "encounter_io_initialize_output nf90_def_dim space_dimid" )           ! 3D space dimension
         call check( nf90_def_dim(nciu%id, nciu%id_dimname, NF90_UNLIMITED, nciu%id_dimid), "encounter_io_initialize_output nf90_def_dim id_dimid" )       ! dimension to store particle id numbers
         call check( nf90_def_dim(nciu%id, nciu%str_dimname, NAMELEN, nciu%str_dimid), "encounter_io_initialize_output nf90_def_dim str_dimid"  )          ! Dimension for string variables (aka character arrays)

         ! Dimension coordinates
         call check( nf90_def_var(nciu%id, nciu%time_dimname, nciu%out_type, nciu%time_dimid, nciu%time_varid), "encounter_io_initialize_output nf90_def_var time_varid"  )
         call check( nf90_def_var(nciu%id, nciu%space_dimname, NF90_CHAR, nciu%space_dimid, nciu%space_varid), "encounter_io_initialize_output nf90_def_var space_varid"  )
         call check( nf90_def_var(nciu%id, nciu%id_dimname, NF90_INT, nciu%id_dimid, nciu%id_varid), "encounter_io_initialize_output nf90_def_var id_varid"  )
         call check( nf90_def_var(nciu%id, nciu%name_varname, NF90_CHAR, [nciu%str_dimid, nciu%id_dimid], nciu%name_varid), "encounter_io_initialize_output nf90_def_var name_varid"  )
      
         ! Variables
         call check( nf90_def_var(nciu%id, nciu%name_varname, NF90_CHAR, [nciu%str_dimid, nciu%id_dimid], nciu%name_varid), "encounter_io_initialize_output nf90_def_var name_varid"  )
         call check( nf90_def_var(nciu%id, nciu%ptype_varname, NF90_CHAR, [nciu%str_dimid, nciu%id_dimid], nciu%ptype_varid), "encounter_io_initialize_output nf90_def_var ptype_varid"  )
         call check( nf90_def_var(nciu%id, nciu%rh_varname,  nciu%out_type, [nciu%space_dimid, nciu%id_dimid, nciu%time_dimid], nciu%rh_varid), "encounter_io_initialize_output nf90_def_var rh_varid"  )
         call check( nf90_def_var(nciu%id, nciu%vh_varname,  nciu%out_type, [nciu%space_dimid, nciu%id_dimid, nciu%time_dimid], nciu%vh_varid), "encounter_io_initialize_output nf90_def_var vh_varid"  )
         call check( nf90_def_var(nciu%id, nciu%gmass_varname, nciu%out_type, [nciu%id_dimid, nciu%time_dimid], nciu%Gmass_varid), "encounter_io_initialize_output nf90_def_var Gmass_varid"  )
         if (param%lclose) then
            call check( nf90_def_var(nciu%id, nciu%radius_varname, nciu%out_type, [nciu%id_dimid, nciu%time_dimid], nciu%radius_varid), "encounter_io_initialize_output nf90_def_var radius_varid"  )
         end if
         if (param%lrotation) then
            call check( nf90_def_var(nciu%id, nciu%Ip_varname, nciu%out_type, [nciu%space_dimid, nciu%id_dimid, nciu%time_dimid], nciu%Ip_varid), "encounter_io_initialize_output nf90_def_var Ip_varid"  )
            call check( nf90_def_var(nciu%id, nciu%rot_varname, nciu%out_type, [nciu%space_dimid, nciu%id_dimid, nciu%time_dimid], nciu%rot_varid), "encounter_io_initialize_output nf90_def_var rot_varid"  )
         end if
         call check( nf90_inquire(nciu%id, nVariables=nvar), "encounter_io_initialize_output nf90_inquire nVariables"  )
         do varid = 1, nvar
            call check( nf90_inquire_variable(nciu%id, varid, xtype=vartype, ndims=ndims), "encounter_io_initialize_output nf90_inquire_variable"  )
            select case(vartype)
            case(NF90_INT)
               call check( nf90_def_var_fill(nciu%id, varid, 0, NF90_FILL_INT), "encounter_io_initialize_output nf90_def_var_fill NF90_INT"  )
            case(NF90_FLOAT)
               call check( nf90_def_var_fill(nciu%id, varid, 0, sfill), "encounter_io_initialize_output nf90_def_var_fill NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call check( nf90_def_var_fill(nciu%id, varid, 0, dfill), "encounter_io_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
            case(NF90_CHAR)
               call check( nf90_def_var_fill(nciu%id, varid, 0, 0), "encounter_io_initialize_output nf90_def_var_fill NF90_CHAR"  )
            end select
         end do

         ! Take the file out of define mode
         call check( nf90_enddef(nciu%id), "encounter_io_initialize_output nf90_enddef"  )
      end associate

      return

      667 continue
      write(*,*) "Error creating encounter output file. " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine encounter_io_initialize_output


   module subroutine encounter_io_open_file(self, param, readonly)
      !! author: David A. Minton
      !!
      !! Opens a NetCDF encounter file and does the variable inquiries to activate variable ids
      implicit none
      ! Arguments
      class(encounter_io_parameters), intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),     intent(in)    :: param    !! Current run configuration parameters
      logical, optional,              intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      ! Internals
      integer(I4B) :: mode
      character(len=STRMAX) :: errmsg

      mode = NF90_WRITE
      if (present(readonly)) then
         if (readonly) mode = NF90_NOWRITE
      end if

      write(errmsg,*) "encounter_io_open_file nf90_open ",trim(adjustl(param%outfile))
      call check( nf90_open(self%enc_file, mode, self%id), errmsg)

      return
   end subroutine encounter_io_open_file

   module subroutine encounter_io_write_frame(self, iu, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of an encounter list structure.
      implicit none
      ! Arguments
      class(encounter_list),          intent(in)    :: self   !! Swiftest encounter structure
      class(encounter_io_parameters), intent(inout) :: iu     !! Parameters used to identify a particular encounter io NetCDF dataset
      class(swiftest_parameters),     intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                                  :: i,old_mode, n

      i = iu%ienc_frame
      n = int(self%nenc, kind=I4B)
      call check( nf90_set_fill(iu%id, nf90_nofill, old_mode), "encounter_io_write_frame_base nf90_set_fill"  )
      call check( nf90_put_var(iu%id, iu%time_varid, self%t, start=[i]), "encounter_io_write_frame nf90_put_var time_varid"  )

      ! call check( nf90_put_var(iu%id, iu%nenc_varid, self%nenc, start=[i]), "encounter_io_frame nf90_put_var nenc_varid"  )
      ! call check( nf90_put_var(iu%id, iu%name_varid, self%name1(1:n), start=[1, 1, i], count=[NAMELEN,1,1]), "netcdf_write_frame nf90_put_var name 1"  )
      ! call check( nf90_put_var(iu%id, iu%name_varid, self%name2(1:n), start=[1, 2, i], count=[NAMELEN,1,1]), "netcdf_write_frame nf90_put_var name 2"  )
      ! call check( nf90_put_var(iu%id, iu%xhx_varid, self%x1(1, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhx_varid 1"  )
      ! call check( nf90_put_var(iu%id, iu%xhy_varid, self%x1(2, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhy_varid 1"  )
      ! call check( nf90_put_var(iu%id, iu%xhz_varid, self%x1(3, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhz_varid 1"  )
      ! call check( nf90_put_var(iu%id, iu%xhx_varid, self%x2(1, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhx_varid 2"  )
      ! call check( nf90_put_var(iu%id, iu%xhy_varid, self%x2(2, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhy_varid 2"  )
      ! call check( nf90_put_var(iu%id, iu%xhz_varid, self%x2(3, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhz_varid 2"  )
      ! call check( nf90_put_var(iu%id, iu%vhx_varid, self%v1(1, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhx_varid 1"  )
      ! call check( nf90_put_var(iu%id, iu%vhy_varid, self%v1(2, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhy_varid 1"  )
      ! call check( nf90_put_var(iu%id, iu%vhz_varid, self%v1(3, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhz_varid 1"  )
      ! call check( nf90_put_var(iu%id, iu%vhx_varid, self%v2(1, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhx_varid 2"  )
      ! call check( nf90_put_var(iu%id, iu%vhy_varid, self%v2(2, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhy_varid 2"  )
      ! call check( nf90_put_var(iu%id, iu%vhz_varid, self%v2(3, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhz_varid 2"  )
      ! call check( nf90_put_var(iu%id, iu%Gmass_varid, self%Gmass1(1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var Gmass 1"  )
      ! call check( nf90_put_var(iu%id, iu%Gmass_varid, self%Gmass2(1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var Gmass 2"  )
      ! call check( nf90_put_var(iu%id, iu%radius_varid, self%radius1(1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var radius 1"  )
      ! call check( nf90_put_var(iu%id, iu%radius_varid, self%radius2(1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var radius 2"  )
      ! select type(self)
      ! class is (symba_encounter)
      !    call check( nf90_put_var(iu%id, iu%level_varid, self%level(1:n), start=[1, i], count=[n,1]), "netcdf_write_frame nf90_put_var level"  )
      ! end select

      return
   end subroutine encounter_io_write_frame

end submodule s_encounter_io