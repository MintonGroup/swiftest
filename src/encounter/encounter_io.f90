!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (encounter) s_encounter_io
   use swiftest
contains

   module subroutine encounter_io_dump(self, param)
      ! author: David A. Minton
      !!
      !! Dumps the time history of an encounter to file.
      implicit none
      ! Arguments
      class(encounter_storage(*)),  intent(inout)        :: self   !! Encounter storage object
      class(base_parameters),   intent(inout)        :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i

      select type(nc => self%nc)
      class is (encounter_netcdf_parameters)
         if (self%iframe > 0) then
            ! Create and save the output files for this encounter and fragmentation
            nc%file_number = nc%file_number + 1 
            call self%make_index_map()
            nc%time_dimsize = self%nt
            nc%name_dimsize = self%nid
            write(nc%file_name, '("encounter_",I0.6,".nc")') nc%file_number
            call nc%initialize(param)

            do i = 1, self%nframes
               if (allocated(self%frame(i)%item)) then
                  select type(snapshot => self%frame(i)%item)
                  class is (encounter_snapshot)
                     param%ioutput = self%tmap(i)
                     call snapshot%write_frame(self,param)
                  end select
               else
                  exit
               end if
            end do

            call nc%close()
            call self%reset()
         end if
      end select

      return
   end subroutine encounter_io_dump


   module subroutine encounter_io_initialize_output(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF encounter file system. This is a simplified version of the main simulation output NetCDF file, but with fewer variables.
      use, intrinsic :: ieee_arithmetic
      use netcdf
      implicit none
      ! Arguments
      class(encounter_netcdf_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(base_parameters),     intent(in)    :: param   !! Current run configuration parameters
      ! Internals
      integer(I4B) :: nvar, varid, vartype
      real(DP) :: dfill
      real(SP) :: sfill
      integer(I4B), parameter :: NO_FILL = 0
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B) :: ndims

      associate(nc => self)
         dfill = ieee_value(dfill, IEEE_QUIET_NAN)
         sfill = ieee_value(sfill, IEEE_QUIET_NAN)

         select case (param%out_type)
         case("NETCDF_FLOAT")
            self%out_type = NF90_FLOAT
         case("NETCDF_DOUBLE")
            self%out_type = NF90_DOUBLE
         end select

         ! Check if the file exists, and if it does, delete it
         inquire(file=nc%file_name, exist=fileExists)
         if (fileExists) then
            open(unit=LUN, file=nc%file_name, status="old", err=667, iomsg=errmsg)
            close(unit=LUN, status="delete")
         end if

         call netcdf_io_check( nf90_create(nc%file_name, NF90_NETCDF4, nc%id), "encounter_io_initialize_output nf90_create" )

         ! Dimensions
         call netcdf_io_check( nf90_def_dim(nc%id, nc%time_dimname, nc%time_dimsize, nc%time_dimid), "encounter_io_initialize_output nf90_def_dim time_dimid" ) ! Simulation time dimension
         call netcdf_io_check( nf90_def_dim(nc%id, nc%space_dimname, NDIM, nc%space_dimid), "encounter_io_initialize_output nf90_def_dim space_dimid" )           ! 3D space dimension
         call netcdf_io_check( nf90_def_dim(nc%id, nc%name_dimname, nc%name_dimsize, nc%name_dimid), "encounter_io_initialize_output nf90_def_dim name_dimid" )       ! dimension to store particle id numbers
         call netcdf_io_check( nf90_def_dim(nc%id, nc%str_dimname, NAMELEN, nc%str_dimid), "encounter_io_initialize_output nf90_def_dim str_dimid"  )          ! Dimension for string variables (aka character arrays)

         ! Dimension coordinates
         call netcdf_io_check( nf90_def_var(nc%id, nc%time_dimname, nc%out_type, nc%time_dimid, nc%time_varid), "encounter_io_initialize_output nf90_def_var time_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%space_dimname, NF90_CHAR, nc%space_dimid, nc%space_varid), "encounter_io_initialize_output nf90_def_var space_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%name_dimname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], nc%name_varid), "encounter_io_initialize_output nf90_def_var id_varid"  )
      
         ! Variables
         call netcdf_io_check( nf90_def_var(nc%id, nc%id_varname, NF90_INT, nc%name_dimid, nc%id_varid), "encounter_io_initialize_output nf90_def_var id_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%ptype_varname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], nc%ptype_varid), "encounter_io_initialize_output nf90_def_var ptype_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%rh_varname,  nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%rh_varid), "encounter_io_initialize_output nf90_def_var rh_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%vh_varname,  nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%vh_varid), "encounter_io_initialize_output nf90_def_var vh_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%Gmass_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%Gmass_varid), "encounter_io_initialize_output nf90_def_var Gmass_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%loop_varname, NF90_INT, [nc%time_dimid], nc%loop_varid), "encounter_io_initialize_output nf90_def_var loop_varid"  )
         if (param%lclose) then
            call netcdf_io_check( nf90_def_var(nc%id, nc%radius_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%radius_varid), "encounter_io_initialize_output nf90_def_var radius_varid"  )
         end if
         if (param%lrotation) then
            call netcdf_io_check( nf90_def_var(nc%id, nc%Ip_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%Ip_varid), "encounter_io_initialize_output nf90_def_var Ip_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%rot_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%rot_varid), "encounter_io_initialize_output nf90_def_var rot_varid"  )
         end if

         call netcdf_io_check( nf90_inquire(nc%id, nVariables=nvar), "encounter_io_initialize_output nf90_inquire nVariables"  )
         do varid = 1, nvar
            call netcdf_io_check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), "encounter_io_initialize_output nf90_inquire_variable"  )
            select case(vartype)
            case(NF90_INT)
               call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, NF90_FILL_INT), "encounter_io_initialize_output nf90_def_var_fill NF90_INT"  )
            case(NF90_FLOAT)
               call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, sfill), "encounter_io_initialize_output nf90_def_var_fill NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, dfill), "encounter_io_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
            case(NF90_CHAR)
               call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, 0), "encounter_io_initialize_output nf90_def_var_fill NF90_CHAR"  )
            end select
         end do

         ! Take the file out of define mode
         call netcdf_io_check( nf90_enddef(nc%id), "encounter_io_initialize_output nf90_enddef"  )

         ! Add in the space dimension coordinates
         call netcdf_io_check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), "encounter_io_initialize_output nf90_put_var space"  )

      end associate

      return

      667 continue
      write(*,*) "Error creating encounter output file. " // trim(adjustl(errmsg))
      call swiftest_util_exit(FAILURE)
   end subroutine encounter_io_initialize_output


   module subroutine encounter_io_write_frame_snapshot(self, history, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of an encounter trajectory.
      use netcdf
      implicit none
      ! Arguments
      class(encounter_snapshot),   intent(in)    :: self              !! Swiftest encounter structure
      class(encounter_storage(*)), intent(inout) :: history !! Encounter storage object
      class(base_parameters),      intent(inout) :: param             !! Current run configuration parameters
 
      ! Internals
      integer(I4B)           :: i, idslot, old_mode, npl, ntp
      character(len=:), allocatable :: charstring

      select type(param)
      class is (swiftest_parameters)
      select type(pl => self%pl)
      class is (swiftest_pl)
      select type(tp => self%tp)
      class is (swiftest_pl)
      select type (nc => history%nc)
      class is (encounter_netcdf_parameters)
         associate(tslot => param%ioutput)
            call netcdf_io_check( nf90_set_fill(nc%id, nf90_nofill, old_mode), "encounter_io_write_frame_snapshot nf90_set_fill"  )
      
            call netcdf_io_check( nf90_put_var(nc%id, nc%time_varid, self%t, start=[tslot]), "encounter_io_write_frame_snapshot nf90_put_var time_varid"  )
            call netcdf_io_check( nf90_put_var(nc%id, nc%loop_varid, int(self%iloop,kind=I4B), start=[tslot]), "encounter_io_write_frame_snapshot nf90_put_var pl loop_varid"  )

            npl = pl%nbody
            do i = 1, npl
               idslot = findloc(history%idvals,pl%id(i),dim=1)
               call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid, pl%id(i),   start=[idslot]), "encounter_io_write_frame_snapshot nf90_put_var pl id_varid"  )
               call netcdf_io_check( nf90_put_var(nc%id, nc%rh_varid, pl%rh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame_snapshot nf90_put_var pl rh_varid"  )
               call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid, pl%vh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame_snapshot nf90_put_var pl vh_varid"  )
               call netcdf_io_check( nf90_put_var(nc%id, nc%Gmass_varid, pl%Gmass(i), start=[idslot, tslot]), "encounter_io_write_frame_snapshot nf90_put_var pl Gmass_varid"  )

               if (param%lclose) call netcdf_io_check( nf90_put_var(nc%id, nc%radius_varid, pl%radius(i), start=[idslot, tslot]), "encounter_io_write_frame_snapshot nf90_put_var pl radius_varid"  )

               if (param%lrotation) then
                  call netcdf_io_check( nf90_put_var(nc%id, nc%Ip_varid, pl%Ip(:,i), start=[1, idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame_snapshot nf90_put_var pl Ip_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%rot_varid, pl%rot(:,i), start=[1,idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame_snapshot nf90_put_var pl rotx_varid"  )
               end if

               charstring = trim(adjustl(pl%info(i)%name))
               call netcdf_io_check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame_snapshot nf90_put_var pl name_varid"  )
               charstring = trim(adjustl(pl%info(i)%particle_type))
               call netcdf_io_check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame_snapshot nf90_put_var pl particle_type_varid"  )
            end do

            ntp = tp%nbody
            do i = 1, ntp
               idslot = findloc(history%idvals,tp%id(i),dim=1)
               call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid, tp%id(i), start=[idslot]), "encounter_io_write_frame_snapshot nf90_put_var tp id_varid"  )
               call netcdf_io_check( nf90_put_var(nc%id, nc%rh_varid, tp%rh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame_snapshot nf90_put_var tp rh_varid"  )
               call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid, tp%vh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame_snapshot nf90_put_var tp vh_varid"  )

               charstring = trim(adjustl(tp%info(i)%name))
               call netcdf_io_check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame_snapshot nf90_put_var tp name_varid"  )
               charstring = trim(adjustl(tp%info(i)%particle_type))
               call netcdf_io_check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame_snapshot nf90_put_var tp particle_type_varid"  )
            end do

            call netcdf_io_check( nf90_set_fill(nc%id, old_mode, old_mode) )
         end associate
      end select
      end select
      end select
      end select

      return
   end subroutine encounter_io_write_frame_snapshot




end submodule s_encounter_io