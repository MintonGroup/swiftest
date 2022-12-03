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
      !! Initialize a NetCDF encounter file system and defines all variables.
      use, intrinsic :: ieee_arithmetic
      implicit none
      ! Arguments
      class(encounter_io_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),     intent(in)    :: param   !! Current run configuration parameters
      ! Internals
      real(DP) :: dfill
      real(SP) :: sfill
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B), dimension(2), parameter :: collider_dimension = [1,2]

      dfill = ieee_value(dfill, IEEE_QUIET_NAN)
      sfill = ieee_value(sfill, IEEE_QUIET_NAN)

      ! Check if the file exists, and if it does, delete it
      inquire(file=self%enc_file, exist=fileExists)
      if (fileExists) then
         open(unit=LUN, file=self%enc_file, status="old", err=667, iomsg=errmsg)
         close(unit=LUN, status="delete")
      end if

      call check( nf90_create(self%enc_file, NF90_NETCDF4, self%id), "encounter_io_initialize_output nf90_create" )

      call check( nf90_def_dim(self%id, self%eid_dimname, NF90_UNLIMITED, self%eid_dimid), "encounter_io_initialize_output nf90_def_dim eid_dimid" )    
      call check( nf90_def_dim(self%id, self%str_dimname, NAMELEN, self%str_dimid), "encounter_io_initialize_output nf90_def_dim str_dimid"  ) ! Dimension for string variables (aka character arrays)
      call check( nf90_def_dim(self%id, self%time_dimname, NF90_UNLIMITED, self%time_dimid), "encounter_io_initialize_output nf90_def_dim time_dimid" ) ! 'y' dimension
      call check( nf90_def_dim(self%id, self%collider_dimname, self%collider_dim_size, self%collider_dimid), "encounter_io_initialize_output nf90_def_dim time_dimid" ) ! 'y' dimension

      select case (param%out_type)
      case("NETCDF_FLOAT")
         self%out_type = NF90_FLOAT
      case("NETCDF_DOUBLE")
         self%out_type = NF90_DOUBLE
      end select

      call check( nf90_def_var(self%id, self%time_dimname, self%out_type, self%time_dimid, self%time_varid), "encounter_io_initialize_output nf90_def_var time_varid"  )
      call check( nf90_def_var(self%id, self%nenc_varname, NF90_INT, self%time_dimid, self%nenc_varid), "encounter_io_initialize_output nf90_def_var nenc_varid"  )
      call check( nf90_def_var(self%id, self%name_varname, NF90_CHAR, [self%str_dimid, self%collider_dimid, self%eid_dimid], self%name_varid), "encounter_io_initialize_output nf90_def_var name_varid"  )
      call check( nf90_def_var(self%id, self%id_dimname, NF90_INT, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%id_varid), "encounter_io_initialize_output nf90_def_var id_varid"  )
      call check( nf90_def_var(self%id, self%xhx_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%xhx_varid), "encounter_io_initialize_output nf90_def_var xhx_varid"  )
      call check( nf90_def_var(self%id, self%xhy_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%xhy_varid), "encounter_io_initialize_output nf90_def_var xhy_varid"  )
      call check( nf90_def_var(self%id, self%xhz_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%xhz_varid), "encounter_io_initialize_output nf90_def_var xhz_varid"  )
      call check( nf90_def_var(self%id, self%vhx_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%vhx_varid), "encounter_io_initialize_output nf90_def_var vhx_varid"  )
      call check( nf90_def_var(self%id, self%vhy_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%vhy_varid), "encounter_io_initialize_output nf90_def_var vhy_varid"  )
      call check( nf90_def_var(self%id, self%vhz_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%vhz_varid), "encounter_io_initialize_output nf90_def_var vhz_varid"  )
      call check( nf90_def_var(self%id, self%level_varname, NF90_INT, [self%eid_dimid, self%time_dimid], self%level_varid), "encounter_io_initialize_output nf90_def_var level_varid"  )
      call check( nf90_def_var(self%id, self%gmass_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%Gmass_varid), "encounter_io_initialize_output nf90_def_var Gmass_varid"  )
      call check( nf90_def_var(self%id, self%radius_varname, self%out_type, [self%collider_dimid, self%eid_dimid, self%time_dimid], self%radius_varid), "encounter_io_initialize_output nf90_def_var radius_varid"  )


      ! Take the file out of define mode
      call check( nf90_enddef(self%id), "encounter_io_initialize_output nf90_enddef"  )

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

      call check( nf90_inq_dimid(self%id, self%time_dimname, self%time_dimid), "encounter_io_open_file nf90_inq_dimid time_dimid"  )
      call check( nf90_inq_dimid(self%id, self%eid_dimname, self%eid_dimid), "encounter_io_open_file nf90_inq_dimid eid_dimid"  )
      call check( nf90_inq_dimid(self%id, self%collider_dimname, self%collider_dimid), "encounter_io_open_file nf90_inq_dimid collider_dimid"  )
      call check( nf90_inq_dimid(self%id, self%str_dimname, self%str_dimid), "encounter_io_open_file nf90_inq_dimid collider_str"  )

      call check( nf90_inq_varid(self%id, self%time_dimname, self%time_varid), "encounter_io_open_file nf90_inq_varid time_varid" )
      call check( nf90_inq_varid(self%id, self%name_varname, self%name_varid), "encounter_io_open_file nf90_inq_varid name_varid" )
      call check( nf90_inq_varid(self%id, self%nenc_varname, self%nenc_varid), "encounter_io_open_file nf90_inq_varid nenc_varid" )

      call check( nf90_inq_varid(self%id, self%xhx_varname, self%xhx_varid), "encounter_io_open_file nf90_inq_varid xhx_varid" )
      call check( nf90_inq_varid(self%id, self%xhy_varname, self%xhy_varid), "encounter_io_open_file nf90_inq_varid xhy_varid" )
      call check( nf90_inq_varid(self%id, self%xhz_varname, self%xhz_varid), "encounter_io_open_file nf90_inq_varid xhz_varid" )
      call check( nf90_inq_varid(self%id, self%vhx_varname, self%vhx_varid), "encounter_io_open_file nf90_inq_varid vhx_varid" )
      call check( nf90_inq_varid(self%id, self%vhy_varname, self%vhy_varid), "encounter_io_open_file nf90_inq_varid vhy_varid" )
      call check( nf90_inq_varid(self%id, self%vhz_varname, self%vhz_varid), "encounter_io_open_file nf90_inq_varid vhz_varid" )
      call check( nf90_inq_varid(self%id, self%level_varname, self%level_varid), "encounter_io_open_file nf90_inq_varid level_varid" )
      call check( nf90_inq_varid(self%id, self%gmass_varname, self%Gmass_varid), "encounter_io_open_file nf90_inq_varid Gmass_varid" )
      call check( nf90_inq_varid(self%id, self%radius_varname, self%radius_varid), "encounter_io_open_file nf90_inq_varid radius_varid" )

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

      call check( nf90_put_var(iu%id, iu%nenc_varid, self%nenc, start=[i]), "encounter_io_frame nf90_put_var nenc_varid"  )
      call check( nf90_put_var(iu%id, iu%name_varid, self%name1(1:n), start=[1, 1, i], count=[NAMELEN,1,1]), "netcdf_write_frame nf90_put_var name 1"  )
      call check( nf90_put_var(iu%id, iu%name_varid, self%name2(1:n), start=[1, 2, i], count=[NAMELEN,1,1]), "netcdf_write_frame nf90_put_var name 2"  )
      call check( nf90_put_var(iu%id, iu%xhx_varid, self%x1(1, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhx_varid 1"  )
      call check( nf90_put_var(iu%id, iu%xhy_varid, self%x1(2, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhy_varid 1"  )
      call check( nf90_put_var(iu%id, iu%xhz_varid, self%x1(3, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhz_varid 1"  )
      call check( nf90_put_var(iu%id, iu%xhx_varid, self%x2(1, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhx_varid 2"  )
      call check( nf90_put_var(iu%id, iu%xhy_varid, self%x2(2, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhy_varid 2"  )
      call check( nf90_put_var(iu%id, iu%xhz_varid, self%x2(3, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var xhz_varid 2"  )
      call check( nf90_put_var(iu%id, iu%vhx_varid, self%v1(1, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhx_varid 1"  )
      call check( nf90_put_var(iu%id, iu%vhy_varid, self%v1(2, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhy_varid 1"  )
      call check( nf90_put_var(iu%id, iu%vhz_varid, self%v1(3, 1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhz_varid 1"  )
      call check( nf90_put_var(iu%id, iu%vhx_varid, self%v2(1, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhx_varid 2"  )
      call check( nf90_put_var(iu%id, iu%vhy_varid, self%v2(2, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhy_varid 2"  )
      call check( nf90_put_var(iu%id, iu%vhz_varid, self%v2(3, 1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var vhz_varid 2"  )
      call check( nf90_put_var(iu%id, iu%Gmass_varid, self%Gmass1(1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var Gmass 1"  )
      call check( nf90_put_var(iu%id, iu%Gmass_varid, self%Gmass2(1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var Gmass 2"  )
      call check( nf90_put_var(iu%id, iu%radius_varid, self%radius1(1:n), start=[1, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var radius 1"  )
      call check( nf90_put_var(iu%id, iu%radius_varid, self%radius2(1:n), start=[2, 1, i], count=[1,n,1]), "netcdf_write_frame nf90_put_var radius 2"  )
      select type(self)
      class is (symba_encounter)
         call check( nf90_put_var(iu%id, iu%level_varid, self%level(1:n), start=[1, i], count=[n,1]), "netcdf_write_frame nf90_put_var level"  )
      end select

      return
   end subroutine encounter_io_write_frame

end submodule s_encounter_io