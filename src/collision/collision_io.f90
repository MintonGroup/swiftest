!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(collision) s_collision_netcdf_io
   use swiftest

contains

   module subroutine collision_netcdf_io_dump(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the time history of an encounter to file.
      implicit none
      ! Arguments
      class(collision_storage(*)),  intent(inout)  :: self   !! Encounter storage object
      class(base_parameters), intent(inout)    :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i

      select type(nc => self%nc)
      class is (collision_netcdf_parameters)
      select type(param)
      class is (swiftest_parameters)
         if (self%iframe > 0) then
            nc%file_number = nc%file_number + 1 
            call self%make_index_map()
            nc%event_dimsize = self%nt
            nc%name_dimsize = self%nid

            write(nc%file_name, '("collision_",I0.6,".nc")') nc%file_number
            call nc%initialize(param)

            do i = 1, self%nframes
               if (allocated(self%frame(i)%item)) then
                  select type(snapshot => self%frame(i)%item)
                  class is (collision_snapshot)
                     param%ioutput = i
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
      end select

      return
   end subroutine collision_netcdf_io_dump

   module subroutine collision_netcdf_io_initialize_output(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF fragment history file system. This is a simplified version of the main simulation output NetCDF file, but with fewer variables.
      use, intrinsic :: ieee_arithmetic
      use netcdf
      implicit none
      ! Arguments
      class(collision_netcdf_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(base_parameters),     intent(in)    :: param   
      ! Internals
      integer(I4B) :: nvar, varid, vartype
      real(DP) :: dfill
      real(SP) :: sfill
      integer(I4B), parameter :: NO_FILL = 0
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B) :: ndims

      select type(param)
      class is (base_parameters)
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

            call netcdf_io_check( nf90_create(nc%file_name, NF90_NETCDF4, nc%id), "collision_netcdf_io_initialize_output nf90_create" )

            ! Dimensions
            call netcdf_io_check( nf90_def_dim(nc%id, nc%event_dimname, nc%event_dimsize, nc%event_dimid), "collision_netcdf_io_initialize_output nf90_def_dim event_dimid"  ) ! Dimension to store individual collision events
            call netcdf_io_check( nf90_def_dim(nc%id, nc%space_dimname, NDIM,             nc%space_dimid), "collision_netcdf_io_initialize_output nf90_def_dim space_dimid" )  ! 3D space dimension
            call netcdf_io_check( nf90_def_dim(nc%id, nc%name_dimname,  nc%name_dimsize,  nc%name_dimid),    "collision_netcdf_io_initialize_output nf90_def_dim name_dimid" )     ! Dimension to store particle id numbers
            call netcdf_io_check( nf90_def_dim(nc%id, nc%str_dimname,   NAMELEN,          nc%str_dimid),   "collision_netcdf_io_initialize_output nf90_def_dim str_dimid"  )   ! Dimension for string variables (aka character arrays)
            call netcdf_io_check( nf90_def_dim(nc%id, nc%stage_dimname, 2,                nc%stage_dimid), "collision_netcdf_io_initialize_output nf90_def_dim stage_dimid"  ) ! Dimension for stage variables (aka "before" vs. "after"

            ! Dimension coordinates
            call netcdf_io_check( nf90_def_var(nc%id, nc%space_dimname, NF90_CHAR,  nc%space_dimid, nc%space_varid), "collision_netcdf_io_initialize_output nf90_def_var space_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%name_dimname,  NF90_CHAR, [nc%str_dimid, nc%name_dimid], nc%name_varid),   "collision_netcdf_io_initialize_output nf90_def_var name_varid")
            call netcdf_io_check( nf90_def_var(nc%id, nc%stage_dimname, NF90_CHAR, [nc%str_dimid, nc%stage_dimid], nc%stage_varid), "collision_netcdf_io_initialize_output nf90_def_var stage_varid"  )
         
            ! Variables
            call netcdf_io_check( nf90_def_var(nc%id, nc%id_varname,    NF90_INT,   nc%name_dimid,    nc%id_varid),    "collision_netcdf_io_initialize_output nf90_def_var id_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%time_dimname,    nc%out_type, &
                                                                                 nc%event_dimid, nc%time_varid),    "collision_netcdf_io_initialize_output nf90_def_var time_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%regime_varname,  NF90_CHAR,  &
               [nc%str_dimid,                                                    nc%event_dimid], nc%regime_varid), "collision_netcdf_io_initialize_output nf90_def_var regime_varid")
            call netcdf_io_check( nf90_def_var(nc%id, nc%Qloss_varname,  nc%out_type,  &
               [                                                                 nc%event_dimid], nc%Qloss_varid),  "collision_netcdf_io_initialize_output nf90_def_var Qloss_varid")
            
            call netcdf_io_check( nf90_def_var(nc%id, nc%ptype_varname,   NF90_CHAR,  &
               [nc%str_dimid,                   nc%name_dimid,    nc%stage_dimid,  nc%event_dimid], nc%ptype_varid),  "collision_netcdf_io_initialize_output nf90_def_var ptype_varid")

            call netcdf_io_check( nf90_def_var(nc%id, nc%loop_varname,    NF90_INT, &
               [                                                                 nc%event_dimid], nc%loop_varid),   "collision_netcdf_io_initialize_output nf90_def_var loop_varid")   

            call netcdf_io_check( nf90_def_var(nc%id, nc%rh_varname,      nc%out_type,&
               [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%event_dimid], nc%rh_varid),     "collision_netcdf_io_initialize_output nf90_def_var rh_varid")

            call netcdf_io_check( nf90_def_var(nc%id, nc%vh_varname,      nc%out_type,&
               [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%event_dimid], nc%vh_varid),     "collision_netcdf_io_initialize_output nf90_def_var vh_varid")

            call netcdf_io_check( nf90_def_var(nc%id, nc%Gmass_varname,   nc%out_type,&
               [                                nc%name_dimid,    nc%stage_dimid,  nc%event_dimid], nc%Gmass_varid),  "collision_netcdf_io_initialize_output nf90_def_var Gmass_varid")


            call netcdf_io_check( nf90_def_var(nc%id, nc%radius_varname,  nc%out_type,&
               [                                nc%name_dimid,    nc%stage_dimid,  nc%event_dimid], nc%radius_varid), "collision_netcdf_io_initialize_output nf90_def_var radius_varid")

            call netcdf_io_check( nf90_def_var(nc%id, nc%Ip_varname,      nc%out_type,&
               [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%event_dimid], nc%Ip_varid),     "collision_netcdf_io_initialize_output nf90_def_var Ip_varid")

            call netcdf_io_check( nf90_def_var(nc%id, nc%rot_varname,     nc%out_type,&
               [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%event_dimid], nc%rot_varid),    "collision_netcdf_io_initialize_output nf90_def_var rot_varid")
            
            if (param%lenergy) then

               call netcdf_io_check( nf90_def_var(nc%id, nc%ke_orb_varname,  nc%out_type,&
                  [                                                nc%stage_dimid,  nc%event_dimid], nc%KE_orb_varid), "collision_netcdf_io_initialize_output nf90_def_var KE_orb_varid")

               call netcdf_io_check( nf90_def_var(nc%id, nc%ke_spin_varname, nc%out_type,&
                  [                                                nc%stage_dimid,  nc%event_dimid], nc%KE_spin_varid), "collision_netcdf_io_initialize_output nf90_def_var KE_spin_varid"  )

               call netcdf_io_check( nf90_def_var(nc%id, nc%pe_varname,      nc%out_type,&
                  [                                                nc%stage_dimid,  nc%event_dimid], nc%PE_varid),      "collision_netcdf_io_initialize_output nf90_def_var PE_varid"  )

               call netcdf_io_check( nf90_def_var(nc%id, nc%L_orb_varname,   nc%out_type, &
                  [              nc%space_dimid,                   nc%stage_dimid,  nc%event_dimid], nc%L_orb_varid),   "collision_netcdf_io_initialize_output nf90_def_var L_orb_varid"  )

               call netcdf_io_check( nf90_def_var(nc%id, nc%Lspin_varname,  nc%out_type,&
                  [              nc%space_dimid,                   nc%stage_dimid,  nc%event_dimid], nc%Lspin_varid),  "collision_netcdf_io_initialize_output nf90_def_var Lspin_varid"  )
            end if 

            call netcdf_io_check( nf90_inquire(nc%id, nVariables=nvar), "collision_netcdf_io_initialize_output nf90_inquire nVariables"  )
            do varid = 1, nvar
               call netcdf_io_check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), "collision_netcdf_io_initialize_output nf90_inquire_variable"  )
               select case(vartype)
               case(NF90_INT)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, NF90_FILL_INT), "collision_netcdf_io_initialize_output nf90_def_var_fill NF90_INT"  )
               case(NF90_FLOAT)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, sfill), "collision_netcdf_io_initialize_output nf90_def_var_fill NF90_FLOAT"  )
               case(NF90_DOUBLE)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, dfill), "collision_netcdf_io_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
               case(NF90_CHAR)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, 0), "collision_netcdf_io_initialize_output nf90_def_var_fill NF90_CHAR"  )
               end select
            end do
            ! Take the file out of define mode
            call netcdf_io_check( nf90_enddef(nc%id), "collision_netcdf_io_initialize_output nf90_enddef"  )

            ! Add in the space and stage dimension coordinates
            call netcdf_io_check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), "collision_netcdf_io_initialize_output nf90_put_var space"  )
            call netcdf_io_check( nf90_put_var(nc%id, nc%stage_varid, nc%stage_coords(1), start=[1,1], count=[len(nc%stage_coords(1)),1]), "collision_netcdf_io_initialize_output nf90_put_var stage 1"  )
            call netcdf_io_check( nf90_put_var(nc%id, nc%stage_varid, nc%stage_coords(2), start=[1,2], count=[len(nc%stage_coords(2)),1]), "collision_netcdf_io_initialize_output nf90_put_var stage 2"  )

         end associate
      end select

      return

      667 continue
      write(*,*) "Error creating fragmentation output file. " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine collision_netcdf_io_initialize_output


   module subroutine collision_netcdf_io_write_frame_snapshot(self, history, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of a collision result
      use netcdf
      implicit none
      ! Arguments
      class(collision_snapshot),   intent(in)    :: self    !! Swiftest encounter structure
      class(encounter_storage(*)), intent(inout) :: history !! Collision history object
      class(base_parameters),      intent(inout) :: param   !! Current run configuration parameters
      ! Internals
      integer(I4B)           :: i, idslot, old_mode, npl, stage
      character(len=:), allocatable :: charstring
      class(swiftest_pl), allocatable :: pl

      select type(nc => history%nc)
      class is (collision_netcdf_parameters)
         associate(system => self%collision_system, impactors => self%collision_system%impactors, fragments => self%collision_system%fragments, eslot => param%ioutput)
            call netcdf_io_check( nf90_set_fill(nc%id, nf90_nofill, old_mode), "collision_netcdf_io_write_frame_snapshot nf90_set_fill" )

            call netcdf_io_check( nf90_put_var(nc%id, nc%time_varid, self%t,                   start=[eslot]), "collision_netcdf_io_write_frame_snapshot nf90_put_var time_varid" )
            call netcdf_io_check( nf90_put_var(nc%id, nc%loop_varid, int(self%iloop,kind=I4B), start=[eslot]), "collision_netcdf_io_write_frame_snapshot nf90_put_varloop_varid" )

            charstring = trim(adjustl(REGIME_NAMES(impactors%regime)))
            call netcdf_io_check( nf90_put_var(nc%id, nc%regime_varid, charstring,             start=[1, eslot], count=[len(charstring), 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var regime_varid" )
            call netcdf_io_check( nf90_put_var(nc%id, nc%Qloss_varid, impactors%Qloss,         start=[eslot] ), "collision_netcdf_io_write_frame_snapshot nf90_put_var Qloss_varid" )

            select type(before =>self%collision_system%before)
            class is (swiftest_nbody_system)
            select type(after =>self%collision_system%before)
            class is (swiftest_nbody_system)
               do stage = 1,2
                  if (allocated(pl)) deallocate(pl)
                  select case(stage)
                  case(1)
                     allocate(pl, source=before%pl)
                  case(2)
                     allocate(pl, source=after%pl)
                  end select
                  npl = pl%nbody
                  do i = 1, npl
                     idslot = findloc(history%idvals,pl%id(i),dim=1)
                     call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid,     pl%id(i),     start=[   idslot              ]), "collision_netcdf_io_write_frame_snapshot nf90_put_var id_varid"  )
                     charstring = trim(adjustl(pl%info(i)%name))
                     call netcdf_io_check( nf90_put_var(nc%id, nc%name_varid,   charstring,   start=[1, idslot              ], count=[len(charstring), 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var name_varid"  )
                     charstring = trim(adjustl(pl%info(i)%particle_type))
                     call netcdf_io_check( nf90_put_var(nc%id, nc%ptype_varid,  charstring,   start=[1, idslot, stage, eslot], count=[len(charstring), 1, 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var particle_type_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%rh_varid,     pl%rh(:,i),   start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var rh_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid,     pl%vh(:,i),   start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var vh_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%Gmass_varid,  pl%Gmass(i),  start=[   idslot, stage, eslot]), "collision_netcdf_io_write_frame_snapshot nf90_put_var Gmass_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%radius_varid, pl%radius(i), start=[   idslot, stage, eslot]), "collision_netcdf_io_write_frame_snapshot nf90_put_var radius_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%Ip_varid,     pl%Ip(:,i),   start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var Ip_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%rot_varid,    pl%rot(:,i),  start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var rotx_varid"  )
                  end do
               end do
            end select
            end select
            if (param%lenergy) then
               call netcdf_io_check( nf90_put_var(nc%id, nc%ke_orb_varid,  system%ke_orbit(:), start=[   1, eslot], count=[      2, 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var ke_orb_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%ke_spin_varid, system%ke_spin(:),  start=[   1, eslot], count=[      2, 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var ke_spin_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%pe_varid,      system%pe(:),       start=[   1, eslot], count=[      2, 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var pe_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%L_orb_varid,   system%Lorbit(:,:), start=[1, 1, eslot], count=[NDIM, 2, 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var L_orb_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%Lspin_varid,  system%Lspin(:,:),  start=[1, 1, eslot], count=[NDIM, 2, 1]), "collision_netcdf_io_write_frame_snapshot nf90_put_var Lspin_varid before" )
            end if
      
            call netcdf_io_check( nf90_set_fill(nc%id, old_mode, old_mode) )
         end associate
      end select
      return
   end subroutine collision_netcdf_io_write_frame_snapshot


end submodule s_collision_netcdf_io