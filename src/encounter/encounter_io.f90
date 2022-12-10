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
contains


   module subroutine encounter_io_dump(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the time history of an encounter to file.
      implicit none
      ! Arguments
      class(encounter_storage(*)),  intent(inout)        :: self   !! Encounter storage object
      class(swiftest_parameters),   intent(inout)        :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i

      do i = 1, self%nframes
         if (allocated(self%frame(i)%item)) then
            select type(snapshot => self%frame(i)%item)
            class is (fraggle_collision_snapshot)
               param%ioutput = i
               call snapshot%write_frame(self%nc,param)
            class is (encounter_snapshot)
               param%ioutput = self%tslot(i)
               call snapshot%write_frame(self%nc,param)
            end select
         else
            exit
         end if
      end do


      return
   end subroutine encounter_io_dump


   module subroutine encounter_io_initialize(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF encounter file system. This is a simplified version of the main simulation output NetCDF file, but with fewer variables.
      use, intrinsic :: ieee_arithmetic
      use netcdf
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
      integer(I4B) :: ndims, i

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

         call check( nf90_create(nc%file_name, NF90_NETCDF4, nc%id), "encounter_io_initialize nf90_create" )

         ! Dimensions
         call check( nf90_def_dim(nc%id, nc%time_dimname, nc%time_dimsize, nc%time_dimid), "encounter_io_initialize nf90_def_dim time_dimid" ) ! Simulation time dimension
         call check( nf90_def_dim(nc%id, nc%space_dimname, NDIM , nc%space_dimid), "encounter_io_initialize nf90_def_dim space_dimid" )           ! 3D space dimension
         call check( nf90_def_dim(nc%id, nc%id_dimname, param%maxid+1, nc%id_dimid), "encounter_io_initialize nf90_def_dim id_dimid" )       ! dimension to store particle id numbers
         call check( nf90_def_dim(nc%id, nc%str_dimname, NAMELEN, nc%str_dimid), "encounter_io_initialize nf90_def_dim str_dimid"  )          ! Dimension for string variables (aka character arrays)

         ! Dimension coordinates
         call check( nf90_def_var(nc%id, nc%time_dimname, nc%out_type, nc%time_dimid, nc%time_varid), "encounter_io_initialize nf90_def_var time_varid"  )
         call check( nf90_def_var(nc%id, nc%space_dimname, NF90_CHAR, nc%space_dimid, nc%space_varid), "encounter_io_initialize nf90_def_var space_varid"  )
         call check( nf90_def_var(nc%id, nc%id_dimname, NF90_INT, nc%id_dimid, nc%id_varid), "encounter_io_initialize nf90_def_var id_varid"  )
      
         ! Variables
         call check( nf90_def_var(nc%id, nc%name_varname, NF90_CHAR, [nc%str_dimid, nc%id_dimid], nc%name_varid), "encounter_io_initialize nf90_def_var name_varid"  )
         call check( nf90_def_var(nc%id, nc%ptype_varname, NF90_CHAR, [nc%str_dimid, nc%id_dimid], nc%ptype_varid), "encounter_io_initialize nf90_def_var ptype_varid"  )
         call check( nf90_def_var(nc%id, nc%rh_varname,  nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%rh_varid), "encounter_io_initialize nf90_def_var rh_varid"  )
         call check( nf90_def_var(nc%id, nc%vh_varname,  nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%vh_varid), "encounter_io_initialize nf90_def_var vh_varid"  )
         call check( nf90_def_var(nc%id, nc%Gmass_varname, nc%out_type, [nc%id_dimid, nc%time_dimid], nc%Gmass_varid), "encounter_io_initialize nf90_def_var Gmass_varid"  )
         call check( nf90_def_var(nc%id, nc%loop_varname, NF90_INT, [nc%time_dimid], nc%loop_varid), "encounter_io_initialize nf90_def_var loop_varid"  )
         if (param%lclose) then
            call check( nf90_def_var(nc%id, nc%radius_varname, nc%out_type, [nc%id_dimid, nc%time_dimid], nc%radius_varid), "encounter_io_initialize nf90_def_var radius_varid"  )
         end if
         if (param%lrotation) then
            call check( nf90_def_var(nc%id, nc%Ip_varname, nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%Ip_varid), "encounter_io_initialize nf90_def_var Ip_varid"  )
            call check( nf90_def_var(nc%id, nc%rot_varname, nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%rot_varid), "encounter_io_initialize nf90_def_var rot_varid"  )
         end if

         call check( nf90_inquire(nc%id, nVariables=nvar), "encounter_io_initialize nf90_inquire nVariables"  )
         do varid = 1, nvar
            call check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), "encounter_io_initialize nf90_inquire_variable"  )
            select case(vartype)
            case(NF90_INT)
               call check( nf90_def_var_fill(nc%id, varid, 0, NF90_FILL_INT), "encounter_io_initialize nf90_def_var_fill NF90_INT"  )
            case(NF90_FLOAT)
               call check( nf90_def_var_fill(nc%id, varid, 0, sfill), "encounter_io_initialize nf90_def_var_fill NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call check( nf90_def_var_fill(nc%id, varid, 0, dfill), "encounter_io_initialize nf90_def_var_fill NF90_DOUBLE"  )
            case(NF90_CHAR)
               call check( nf90_def_var_fill(nc%id, varid, 0, 0), "encounter_io_initialize nf90_def_var_fill NF90_CHAR"  )
            end select
         end do

         ! Take the file out of define mode
         call check( nf90_enddef(nc%id), "encounter_io_initialize nf90_enddef"  )

         ! Add in the space dimension coordinates
         call check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), "encounter_io_initialize nf90_put_var space"  )

         ! Pre-fill name slots with ids 
         call check( nf90_put_var(nc%id, nc%id_varid, [(-1,i=1,param%maxid+1)], start=[1], count=[param%maxid+1]), "encounter_io_initialize nf90_put_var pl id_varid"  )
      end associate

      return

      667 continue
      write(*,*) "Error creating encounter output file. " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine encounter_io_initialize


   module subroutine encounter_io_write_frame(self, nc, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of an encounter trajectory.
      use netcdf
      implicit none
      ! Arguments
      class(encounter_snapshot),      intent(in)    :: self  !! Swiftest encounter structure
      class(encounter_io_parameters), intent(inout) :: nc   !! Parameters used to identify a particular encounter io NetCDF dataset
      class(swiftest_parameters),     intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B)           :: i, tslot, idslot, old_mode, npl, ntp
      character(len=:), allocatable :: charstring

      tslot = param%ioutput
      associate(pl => self%pl, tp => self%tp)

         call check( nf90_set_fill(nc%id, nf90_nofill, old_mode), "encounter_io_write_frame nf90_set_fill"  )
   
         call check( nf90_put_var(nc%id, nc%time_varid, self%t, start=[tslot]), "encounter_io_write_frame nf90_put_var time_varid"  )
         call check( nf90_put_var(nc%id, nc%loop_varid, int(self%iloop,kind=I4B), start=[tslot]), "encounter_io_write_frame nf90_put_var pl loop_varid"  )

         npl = pl%nbody
         do i = 1, npl
            idslot = pl%id(i) + 1
            call check( nf90_put_var(nc%id, nc%id_varid, pl%id(i),   start=[idslot]), "encounter_io_write_frame nf90_put_var pl id_varid"  )
            call check( nf90_put_var(nc%id, nc%rh_varid, pl%rh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame nf90_put_var pl rh_varid"  )
            call check( nf90_put_var(nc%id, nc%vh_varid, pl%vh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame nf90_put_var pl vh_varid"  )
            call check( nf90_put_var(nc%id, nc%Gmass_varid, pl%Gmass(i), start=[idslot, tslot]), "encounter_io_write_frame nf90_put_var pl Gmass_varid"  )

            if (param%lclose) call check( nf90_put_var(nc%id, nc%radius_varid, pl%radius(i), start=[idslot, tslot]), "encounter_io_write_frame nf90_put_var pl radius_varid"  )

            if (param%lrotation) then
               call check( nf90_put_var(nc%id, nc%Ip_varid, pl%Ip(:,i), start=[1, idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame nf90_put_var pl Ip_varid"  )
               call check( nf90_put_var(nc%id, nc%rot_varid, pl%rot(:,i), start=[1,idslot, tslot], count=[NDIM,1,1]), "encounter_io_write_frame nf90_put_var pl rotx_varid"  )
            end if

            charstring = trim(adjustl(pl%info(i)%name))
            call check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame nf90_put_var pl name_varid"  )
            charstring = trim(adjustl(pl%info(i)%particle_type))
            call check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame nf90_put_var pl particle_type_varid"  )
         end do

         ntp = tp%nbody
         do i = 1, ntp
            idslot = tp%id(i) + 1
            call check( nf90_put_var(nc%id, nc%id_varid, tp%id(i), start=[idslot]), "encounter_io_write_frame nf90_put_var tp id_varid"  )
            call check( nf90_put_var(nc%id, nc%rh_varid, tp%rh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame nf90_put_var tp rh_varid"  )
            call check( nf90_put_var(nc%id, nc%vh_varid, tp%vh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "encounter_io_write_frame nf90_put_var tp vh_varid"  )

            charstring = trim(adjustl(tp%info(i)%name))
            call check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame nf90_put_var tp name_varid"  )
            charstring = trim(adjustl(tp%info(i)%particle_type))
            call check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[len(charstring), 1]), "encounter_io_write_frame nf90_put_var tp particle_type_varid"  )
         end do

         call check( nf90_set_fill(nc%id, old_mode, old_mode) )
      end associate

      return
   end subroutine encounter_io_write_frame




end submodule s_encounter_io