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


   module subroutine encounter_io_write_frame(self, nciu, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of an encounter list structure.
      implicit none
      ! Arguments
      class(encounter_list),          intent(in)    :: self   !! Swiftest encounter structure
      class(encounter_io_parameters), intent(inout) :: nciu     !! Parameters used to identify a particular encounter io NetCDF dataset
      class(swiftest_parameters),     intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                             :: tslot,i,old_mode, n
      character(len=NAMELEN)                   :: charstring

      tslot = nciu%ienc_frame
      call check( nf90_set_fill(nciu%id, nf90_nofill, old_mode), "encounter_io_write_frame_base nf90_set_fill"  )
      call check( nf90_put_var(nciu%id, nciu%time_varid, self%t, start=[tslot]), "encounter_io_write_frame nf90_put_var time_varid"  )

      ! charstring = trim(adjustl(self%info(j)%name))
      ! call check( nf90_put_var(nciu%id, nciu%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "encounter_io_write_info_base nf90_put_var name_varid"  )

      ! charstring = trim(adjustl(self%info(j)%particle_type))
      ! call check( nf90_put_var(nciu%id, nciu%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "encounter_io_write_info_base nf90_put_var particle_type_varid"  )


      ! call check( nf90_put_var(nciu%id, nciu%rh_varid, self%rh(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame_base nf90_put_var rh_varid"  )
      ! call check( nf90_put_var(nciu%id, nciu%vh_varid, self%vh(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame_base nf90_put_var vh_varid"  )
      ! call check( nf90_put_var(nciu%id, nciu%Gmass_varid, self%Gmass(j), start=[idslot, tslot]), "encounter_io_write_frame_base nf90_put_var body Gmass_varid"  )
      ! if (param%lclose) call check( nf90_put_var(nciu%id, nciu%radius_varid, self%radius(j), start=[idslot, tslot]), "encounter_io_write_frame_base nf90_put_var body radius_varid"  )
      ! if (param%lrotation) then
      !    call check( nf90_put_var(nciu%id, nciu%Ip_varid, self%Ip(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame_base nf90_put_var body Ip_varid"  )
      !    call check( nf90_put_var(nciu%id, nciu%rot_varid, self%rot(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame_base nf90_put_var body rotx_varid"  )
      ! end if

      return
   end subroutine encounter_io_write_frame

end submodule s_encounter_io