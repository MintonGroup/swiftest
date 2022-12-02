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
      class(encounter_storage(*)), intent(inout) :: self   !! Encounter storage object
      class(swiftest_parameters),  intent(inout) :: param  !! Current run configuration parameters 
   end subroutine encounter_io_dump_storage_list


   module subroutine encounter_io_initialize_output(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF encounter file system and defines all variables.
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

      dfill = ieee_value(dfill, IEEE_QUIET_NAN)
      sfill = ieee_value(sfill, IEEE_QUIET_NAN)

      ! Check if the file exists, and if it does, delete it
      inquire(file=param%outfile, exist=fileExists)
      if (fileExists) then
         open(unit=LUN, file=self%outfile, status="old", err=667, iomsg=errmsg)
         close(unit=LUN, status="delete")
      end if

      call check( nf90_create(self%outfile, NF90_NETCDF4, self%ncid), "encounter_io_initialize_output nf90_create" )

      call check( nf90_def_dim(self%ncid, ENCID_DIMNAME, NF90_UNLIMITED, self%encid_dimid), "encounter_io_initialize_output nf90_def_dim encid_dimid" )    
      call check( nf90_def_dim(self%ncid, STR_DIMNAME, NAMELEN, self%str_dimid), "encounter_io_initialize_output nf90_def_dim str_dimid"  ) ! Dimension for string variables (aka character arrays)
      call check( nf90_def_dim(self%ncid, TIME_DIMNAME, NF90_UNLIMITED, self%time_dimid), "encounter_io_initialize_output nf90_def_dim time_dimid" ) ! 'y' dimension

      select case (param%out_type)
      case("NETCDF_FLOAT")
         self%out_type = NF90_FLOAT
      case("NETCDF_DOUBLE")
         self%out_type = NF90_DOUBLE
      end select

      call check( nf90_def_var(self%ncid, TIME_DIMNAME, self%out_type, self%time_dimid, self%time_varid), "encounter_io_initialize_output nf90_def_var time_varid"  )
      call check( nf90_def_var(self%ncid, ENCID_DIMNAME, NF90_INT, self%encid_dimid, self%encid_varid), "encounter_io_initialize_output nf90_def_var encid_varid"  )
      call check( nf90_def_var(self%ncid, NENC_VARNAME, NF90_INT, self%time_dimid, self%nenc_varid), "encounter_io_initialize_output nf90_def_var nenc_varid"  )
      call check( nf90_def_var(self%ncid, ID1_VARNAME, NF90_INT, [self%encid_dimid, self%time_dimid], self%id1_varid), "encounter_io_initialize_output nf90_def_var id1_varid"  )
      call check( nf90_def_var(self%ncid, ID2_VARNAME, NF90_INT, [self%encid_dimid, self%time_dimid], self%id2_varid), "encounter_io_initialize_output nf90_def_var id2_varid"  )
      call check( nf90_def_var(self%ncid, X1X_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%x1x_varid), "encounter_io_initialize_output nf90_def_var x1x_varid"  )
      call check( nf90_def_var(self%ncid, X1Y_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%x1y_varid), "encounter_io_initialize_output nf90_def_var x1y_varid"  )
      call check( nf90_def_var(self%ncid, X1Z_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%x1z_varid), "encounter_io_initialize_output nf90_def_var x1z_varid"  )
      call check( nf90_def_var(self%ncid, X2X_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%x2x_varid), "encounter_io_initialize_output nf90_def_var x2x_varid"  )
      call check( nf90_def_var(self%ncid, X2Y_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%x2y_varid), "encounter_io_initialize_output nf90_def_var x2y_varid"  )
      call check( nf90_def_var(self%ncid, X2Z_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%x2z_varid), "encounter_io_initialize_output nf90_def_var x2z_varid"  )
      call check( nf90_def_var(self%ncid, V1X_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%v1x_varid), "encounter_io_initialize_output nf90_def_var v1x_varid"  )
      call check( nf90_def_var(self%ncid, V1Y_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%v1y_varid), "encounter_io_initialize_output nf90_def_var v1y_varid"  )
      call check( nf90_def_var(self%ncid, V1Z_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%v1z_varid), "encounter_io_initialize_output nf90_def_var v1z_varid"  )
      call check( nf90_def_var(self%ncid, V2X_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%v2x_varid), "encounter_io_initialize_output nf90_def_var v2x_varid"  )
      call check( nf90_def_var(self%ncid, V2Y_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%v2y_varid), "encounter_io_initialize_output nf90_def_var v2y_varid"  )
      call check( nf90_def_var(self%ncid, V2Z_VARNAME, self%out_type, [self%encid_dimid, self%time_dimid], self%v2z_varid), "encounter_io_initialize_output nf90_def_var v2z_varid"  )
      call check( nf90_def_var(self%ncid, LEVEL_VARNAME, NF90_INT, [self%encid_dimid, self%time_dimid], self%level_varid), "encounter_io_initialize_output nf90_def_var level_varid"  )


      ! Take the file out of define mode
      call check( nf90_enddef(self%ncid), "encounter_io_initialize_output nf90_enddef"  )

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
      integer(I4B) :: mode, status
      character(len=NF90_MAX_NAME) :: str_dim_name
      character(len=STRMAX) :: errmsg

      mode = NF90_WRITE
      if (present(readonly)) then
         if (readonly) mode = NF90_NOWRITE
      end if

      write(errmsg,*) "netcdf_open nf90_open ",trim(adjustl(param%outfile))
      call check( nf90_open(self%outfile, mode, self%ncid), errmsg)

      return
   end subroutine encounter_io_open_file

end submodule s_encounter_io