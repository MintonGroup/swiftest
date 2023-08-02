!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(collision) s_collision_io
   use swiftest

contains

   module subroutine collision_io_collider_message(pl, collidx, collider_message)
      !! author: David A. Minton
      !!
      !! Prints a nicely formatted message about which bodies collided, including their names
      !! This subroutine appends the body names and ids to an input message.
      implicit none
      ! Arguments
      class(base_object),            intent(in)    :: pl            !! Swiftest massive body object
      integer(I4B),    dimension(:), intent(in)    :: collidx           !! Index of collisional colliders%idx members
      character(*),                  intent(inout) :: collider_message !! The message to print to the screen.
      ! Internals
      integer(I4B) :: i, n
      character(len=STRMAX) :: idstr
      
      n = size(collidx)
      if (n == 0) return

      select type(pl)
      class is (swiftest_pl)
         do i = 1, n
            if (i > 1) collider_message = trim(adjustl(collider_message)) // " and "
            collider_message = " " // trim(adjustl(collider_message)) // " " // trim(adjustl(pl%info(collidx(i))%name))
            write(idstr, '(I10)') pl%id(collidx(i))
            collider_message = trim(adjustl(collider_message)) // " (" // trim(adjustl(idstr)) // ") "
         end do

      end select

      return
   end subroutine collision_io_collider_message


   module subroutine collision_io_log_regime(impactors)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the collisional regime determination
      implicit none
      ! Arguments
      class(collision_impactors), intent(inout) :: impactors  !! Collision system object
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=LUN, file=COLLISION_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Collisional regime determination results"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "True number of impactors : ",impactors%ncoll
      write(LUN, *) "Index list of true impactors  : ",impactors%id(1:impactors%ncoll)
      select case(impactors%regime) 
      case(COLLRESOLVE_REGIME_MERGE)
         write(LUN, *) "Regime: Merge"
      case(COLLRESOLVE_REGIME_DISRUPTION)
         write(LUN, *) "Regime: Disruption"
      case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
         write(LUN, *) "Regime: Supercatastrophic disruption"
      case(COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
         write(LUN, *) "Regime: Graze and merge"
      case(COLLRESOLVE_REGIME_HIT_AND_RUN)
         write(LUN, *) "Regime: Hit and run"
      end select
      write(LUN, *) "Expected energy change       : ", -impactors%Qloss
      write(LUN, *) "--------------------------------------------------------------------"
      close(LUN)

      return
      667 continue
      write(*,*) "Error writing collision regime information to log file: " // trim(adjustl(errmsg))
   end subroutine collision_io_log_regime


   module subroutine collision_io_netcdf_dump(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the time history of an encounter to file.
      implicit none
      ! Arguments
      class(collision_storage),  intent(inout)  :: self   !! Encounter storage object
      class(base_parameters), intent(inout)    :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i

      select type(nc => self%nc)
      class is (collision_netcdf_parameters)
      select type(param)
      class is (base_parameters)
         if (self%iframe > 0) then
            call self%make_index_map()
            call nc%open(param)
            do i = 1, self%iframe
               if (allocated(self%frame(i)%item)) then
                  select type(snapshot => self%frame(i)%item)
                  class is (collision_snapshot)
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
   end subroutine collision_io_netcdf_dump


   module subroutine collision_io_netcdf_initialize_output(self, param)
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
            case default
               write(*,*) trim(adjustl(param%out_type)), " is an invalid OUT_TYPE"
            end select

            ! Check if the file exists, and if it does, delete it
            inquire(file=nc%file_name, exist=fileExists)
            if (fileExists) then
               open(unit=LUN, file=nc%file_name, status="old", err=667, iomsg=errmsg)
               close(unit=LUN, status="delete")
            end if

            call netcdf_io_check( nf90_create(nc%file_name, NF90_NETCDF4, nc%id), "collision_io_netcdf_initialize_output nf90_create" )
            nc%lfile_is_open = .true.

            ! Dimensions
            call netcdf_io_check( nf90_def_dim(nc%id, nc%collision_id_varname, NF90_UNLIMITED, nc%collision_id_dimid), "collision_io_netcdf_initialize_output nf90_def_dim collision_id_dimid"  ) ! Dimension to store individual collision events
            call netcdf_io_check( nf90_def_dim(nc%id, nc%space_dimname,        NDIM,           nc%space_dimid), "collision_io_netcdf_initialize_output nf90_def_dim space_dimid" )  ! 3D space dimension
            call netcdf_io_check( nf90_def_dim(nc%id, nc%name_dimname,         NF90_UNLIMITED, nc%name_dimid),    "collision_io_netcdf_initialize_output nf90_def_dim name_dimid" )     ! Dimension to store particle id numbers
            call netcdf_io_check( nf90_def_dim(nc%id, nc%str_dimname,          NAMELEN,        nc%str_dimid),   "collision_io_netcdf_initialize_output nf90_def_dim str_dimid"  )   ! Dimension for string variables (aka character arrays)
            call netcdf_io_check( nf90_def_dim(nc%id, nc%stage_dimname,        2,              nc%stage_dimid), "collision_io_netcdf_initialize_output nf90_def_dim stage_dimid"  ) ! Dimension for stage variables (aka "before" vs. "after"

            ! Dimension coordinates
            call netcdf_io_check( nf90_def_var(nc%id, nc%collision_id_varname,  NF90_INT,  [nc%collision_id_dimid], nc%collision_id_varid),  "collision_io_netcdf_initialize_output nf90_def_var collision_id_varid")
            call netcdf_io_check( nf90_def_var(nc%id, nc%space_dimname,         NF90_CHAR,  nc%space_dimid, nc%space_varid), "collision_io_netcdf_initialize_output nf90_def_var space_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%name_dimname,          NF90_CHAR, [nc%str_dimid,   nc%name_dimid], nc%name_varid),   "collision_io_netcdf_initialize_output nf90_def_var name_varid")
            call netcdf_io_check( nf90_def_var(nc%id, nc%stage_dimname,         NF90_CHAR, [nc%str_dimid,   nc%stage_dimid], nc%stage_varid), "collision_io_netcdf_initialize_output nf90_def_var stage_varid"  )
         
            ! Variables
            call netcdf_io_check( nf90_def_var(nc%id, nc%id_varname,    NF90_INT,   nc%name_dimid,    nc%id_varid),    "collision_io_netcdf_initialize_output nf90_def_var id_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%time_dimname,    nc%out_type, &
                                                                                 nc%collision_id_dimid, nc%time_varid),    "collision_io_netcdf_initialize_output nf90_def_var time_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%regime_varname,  NF90_CHAR,  &
               [nc%str_dimid,                                                    nc%collision_id_dimid], nc%regime_varid), "collision_io_netcdf_initialize_output nf90_def_var regime_varid")
            call netcdf_io_check( nf90_def_var(nc%id, nc%Qloss_varname,  nc%out_type,  &
               [                                                                 nc%collision_id_dimid], nc%Qloss_varid),  "collision_io_netcdf_initialize_output nf90_def_var Qloss_varid")
            
            call netcdf_io_check( nf90_def_var(nc%id, nc%ptype_varname,   NF90_CHAR,  &
               [nc%str_dimid,                   nc%name_dimid,    nc%stage_dimid,  nc%collision_id_dimid], nc%ptype_varid),  "collision_io_netcdf_initialize_output nf90_def_var ptype_varid")


            call netcdf_io_check( nf90_def_var(nc%id, nc%rh_varname,      nc%out_type,&
               [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%collision_id_dimid], nc%rh_varid),     "collision_io_netcdf_initialize_output nf90_def_var rh_varid")

            call netcdf_io_check( nf90_def_var(nc%id, nc%vh_varname,      nc%out_type,&
               [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%collision_id_dimid], nc%vh_varid),     "collision_io_netcdf_initialize_output nf90_def_var vh_varid")

            call netcdf_io_check( nf90_def_var(nc%id, nc%Gmass_varname,   nc%out_type,&
               [                                nc%name_dimid,    nc%stage_dimid,  nc%collision_id_dimid], nc%Gmass_varid),  "collision_io_netcdf_initialize_output nf90_def_var Gmass_varid")


            call netcdf_io_check( nf90_def_var(nc%id, nc%radius_varname,  nc%out_type,&
               [                                nc%name_dimid,    nc%stage_dimid,  nc%collision_id_dimid], nc%radius_varid), "collision_io_netcdf_initialize_output nf90_def_var radius_varid")

            if (param%lrotation) then
               call netcdf_io_check( nf90_def_var(nc%id, nc%Ip_varname,      nc%out_type,&
                  [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%collision_id_dimid], nc%Ip_varid),     "collision_io_netcdf_initialize_output nf90_def_var Ip_varid")

               call netcdf_io_check( nf90_def_var(nc%id, nc%rot_varname,     nc%out_type,&
                  [              nc%space_dimid,   nc%name_dimid,    nc%stage_dimid,  nc%collision_id_dimid], nc%rot_varid),    "collision_io_netcdf_initialize_output nf90_def_var rot_varid")
            end if

            if (param%lenergy) then
               call netcdf_io_check( nf90_def_var(nc%id, nc%ke_orb_varname,  nc%out_type,&
                  [                                                nc%stage_dimid,  nc%collision_id_dimid], nc%KE_orb_varid), "collision_io_netcdf_initialize_output nf90_def_var KE_orb_varid")

               call netcdf_io_check( nf90_def_var(nc%id, nc%ke_spin_varname, nc%out_type,&
                  [                                                nc%stage_dimid,  nc%collision_id_dimid], nc%KE_spin_varid), "collision_io_netcdf_initialize_output nf90_def_var KE_spin_varid"  )

               call netcdf_io_check( nf90_def_var(nc%id, nc%pe_varname,      nc%out_type,&
                  [                                                nc%stage_dimid,  nc%collision_id_dimid], nc%PE_varid),      "collision_io_netcdf_initialize_output nf90_def_var PE_varid"  )

               call netcdf_io_check( nf90_def_var(nc%id, nc%be_varname,      nc%out_type,&
                  [                                                nc%stage_dimid,  nc%collision_id_dimid], nc%BE_varid),      "collision_io_netcdf_initialize_output nf90_def_var BE_varid"  )
               call netcdf_io_check( nf90_def_var(nc%id, nc%te_varname,      nc%out_type,&
                  [                                                nc%stage_dimid,  nc%collision_id_dimid], nc%TE_varid),      "collision_io_netcdf_initialize_output nf90_def_var TE_varid"  )

               if (param%lrotation) then
                  call netcdf_io_check( nf90_def_var(nc%id, nc%L_orbit_varname,   nc%out_type, &
                     [              nc%space_dimid,                   nc%stage_dimid,  nc%collision_id_dimid], nc%L_orbit_varid),   "collision_io_netcdf_initialize_output nf90_def_var L_orbit_varid"  )

                  call netcdf_io_check( nf90_def_var(nc%id, nc%L_spin_varname,  nc%out_type,&
                     [              nc%space_dimid,                   nc%stage_dimid,  nc%collision_id_dimid], nc%L_spin_varid),  "collision_io_netcdf_initialize_output nf90_def_var L_spin_varid"  )
               end if
            end if 

            call netcdf_io_check( nf90_inquire(nc%id, nVariables=nvar), "collision_io_netcdf_initialize_output nf90_inquire nVariables"  )
            do varid = 1, nvar
               call netcdf_io_check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), "collision_io_netcdf_initialize_output nf90_inquire_variable"  )
               select case(vartype)
               case(NF90_INT)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, NF90_FILL_INT), "collision_io_netcdf_initialize_output nf90_def_var_fill NF90_INT"  )
               case(NF90_FLOAT)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, sfill), "collision_io_netcdf_initialize_output nf90_def_var_fill NF90_FLOAT"  )
               case(NF90_DOUBLE)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, dfill), "collision_io_netcdf_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
               case(NF90_CHAR)
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, 0), "collision_io_netcdf_initialize_output nf90_def_var_fill NF90_CHAR"  )
               end select
            end do

            ! Take the file out of define mode
            call netcdf_io_check( nf90_enddef(nc%id), "collision_io_netcdf_initialize_output nf90_enddef"  )

            ! Add in the space and stage dimension coordinates
            call netcdf_io_check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), "collision_io_netcdf_initialize_output nf90_put_var space"  )
            call netcdf_io_check( nf90_put_var(nc%id, nc%stage_varid, nc%stage_coords(1), start=[1,1], count=[len(nc%stage_coords(1)),1]), "collision_io_netcdf_initialize_output nf90_put_var stage 1"  )
            call netcdf_io_check( nf90_put_var(nc%id, nc%stage_varid, nc%stage_coords(2), start=[1,2], count=[len(nc%stage_coords(2)),1]), "collision_io_netcdf_initialize_output nf90_put_var stage 2"  )

         end associate
      end select

      return

      667 continue
      write(*,*) "Error creating fragmentation output file. " // trim(adjustl(errmsg))
      call base_util_exit(FAILURE,unit=param%display_unit)
   end subroutine collision_io_netcdf_initialize_output


   module subroutine collision_io_netcdf_open(self, param, readonly)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Opens a NetCDF file and does the variable inquiries to activate variable ids
      use netcdf
      implicit none
      ! Arguments
      class(collision_netcdf_parameters), intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
      class(base_parameters),             intent(in)    :: param    !! Current run configuration parameters
      logical, optional,                  intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      ! Internals
      integer(I4B) :: mode
      character(len=STRMAX) :: errmsg
      logical fileExists

      mode = NF90_WRITE
      if (present(readonly)) then
         if (readonly) mode = NF90_NOWRITE
      end if

      select type(param)
      class is (base_parameters)
         associate(nc => self)

            inquire(file=nc%file_name, exist=fileExists)
            if (.not.fileExists) then
               call nc%initialize(param)
               return
            end if

            write(errmsg,*) "collision_io_netcdf_open nf90_open ",trim(adjustl(nc%file_name))
            call netcdf_io_check( nf90_open(nc%file_name, mode, nc%id), errmsg)
            self%lfile_is_open = .true.

            ! Dimensions
            call netcdf_io_check( nf90_inq_dimid(nc%id, nc%collision_id_varname, nc%collision_id_dimid), "collision_io_netcdf_open nf90_inq_dimid collision_id_dimid"  )
            call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%collision_id_dimid, nc%collision_id_varname, len=nc%max_idslot), "collision_io_find_idslot nf90_inquire_dimension max_idslot"  )
            call netcdf_io_check( nf90_inq_dimid(nc%id, nc%space_dimname, nc%space_dimid), "collision_io_netcdf_open nf90_inq_dimid space_dimid"  )
            call netcdf_io_check( nf90_inq_dimid(nc%id, nc%name_dimname, nc%name_dimid), "collision_io_netcdf_open nf90_inq_dimid name_dimid"  )
            call netcdf_io_check( nf90_inq_dimid(nc%id, nc%str_dimname, nc%str_dimid), "collision_io_netcdf_open nf90_inq_dimid str_dimid"  )
            call netcdf_io_check( nf90_inq_dimid(nc%id, nc%stage_dimname, nc%stage_dimid), "collision_io_netcdf_open nf90_inq_dimid stage_dimid"  )

            ! Dimension coordinates
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%collision_id_varname, nc%collision_id_varid), "collision_io_netcdf_open nf90_inq_varid collision_id_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%space_dimname, nc%space_varid), "collision_io_netcdf_open nf90_inq_varid space_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%name_dimname, nc%name_varid), "collision_io_netcdf_open nf90_inq_varid name_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%stage_dimname, nc%stage_varid), "collision_io_netcdf_open nf90_inq_varid stage_varid" )

            ! Required Variables
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%id_varname, nc%id_varid), "collision_io_netcdf_open nf90_inq_varid name_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%time_dimname, nc%time_varid), "collision_io_netcdf_open nf90_inq_varid time_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%regime_varname, nc%regime_varid), "collision_io_netcdf_open nf90_inq_varid regime_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%Qloss_varname, nc%Qloss_varid), "collision_io_netcdf_open nf90_inq_varid Qloss_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%ptype_varname, nc%ptype_varid), "collision_io_netcdf_open nf90_inq_varid ptype_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%rh_varname, nc%rh_varid), "collision_io_netcdf_open nf90_inq_varid rh_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%vh_varname, nc%vh_varid), "collision_io_netcdf_open nf90_inq_varid vh_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%Gmass_varname, nc%Gmass_varid), "collision_io_netcdf_open nf90_inq_varid Gmass_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%radius_varname, nc%radius_varid), "collision_io_netcdf_open nf90_inq_varid radius_varid" )
            
            if (param%lrotation) then
               call netcdf_io_check( nf90_inq_varid(nc%id, nc%Ip_varname, nc%Ip_varid), "collision_io_netcdf_open nf90_inq_varid Ip_varid" )
               call netcdf_io_check( nf90_inq_varid(nc%id, nc%rot_varname, nc%rot_varid), "collision_io_netcdf_open nf90_inq_varid rot_varid" )
            end if

            if (param%lenergy) then
               call netcdf_io_check( nf90_inq_varid(nc%id, nc%ke_orb_varname, nc%ke_orb_varid), "collision_io_netcdf_open nf90_inq_varid ke_orb_varid" )
               call netcdf_io_check( nf90_inq_varid(nc%id, nc%pe_varname, nc%pe_varid), "collision_io_netcdf_open nf90_inq_varid pe_varid" )
               call netcdf_io_check( nf90_inq_varid(nc%id, nc%be_varname, nc%be_varid), "collision_io_netcdf_open nf90_inq_varid be_varid" )
               call netcdf_io_check( nf90_inq_varid(nc%id, nc%te_varname, nc%te_varid), "collision_io_netcdf_open nf90_inq_varid te_varid" )
               call netcdf_io_check( nf90_inq_varid(nc%id, nc%L_orbit_varname, nc%L_orbit_varid), "collision_io_netcdf_open nf90_inq_varid L_orbit_varid" )
               if (param%lrotation) then
                  call netcdf_io_check( nf90_inq_varid(nc%id, nc%ke_spin_varname, nc%ke_spin_varid), "collision_io_netcdf_open nf90_inq_varid ke_spin_varid" )
                  call netcdf_io_check( nf90_inq_varid(nc%id, nc%L_spin_varname, nc%L_spin_varid), "collision_io_netcdf_open nf90_inq_varid L_spin_varid" )
               end if
            end if

         end associate
      end select 

      return
   end subroutine collision_io_netcdf_open


   module subroutine collision_io_netcdf_write_frame_snapshot(self, history, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of a collision result
      use netcdf
      implicit none
      ! Arguments
      class(collision_snapshot),   intent(in)    :: self    !! Swiftest encounter structure
      class(encounter_storage), intent(inout) :: history !! Collision history object
      class(base_parameters),      intent(inout) :: param   !! Current run configuration parameters
      ! Internals
      integer(I4B)           :: i, idslot, old_mode, npl, stage, tmp
      character(len=NAMELEN) :: charstring
      class(swiftest_pl), allocatable :: pl

      select type(nc => history%nc)
      class is (collision_netcdf_parameters)
         associate(collider => self%collider, impactors => self%collider%impactors, fragments => self%collider%fragments, eslot => self%collider%collision_id)
            call netcdf_io_check( nf90_set_fill(nc%id, NF90_NOFILL, old_mode), "collision_io_netcdf_write_frame_snapshot nf90_set_fill" )

            call netcdf_io_check( nf90_put_var(nc%id, nc%collision_id_varid, eslot,            start=[eslot]), "collision_io_netcdf_write_frame_snapshot nf90_put_var collision_id_varid" )
            call netcdf_io_check( nf90_put_var(nc%id, nc%time_varid, self%t,                   start=[eslot]), "collision_io_netcdf_write_frame_snapshot nf90_put_var time_varid" )

            charstring = trim(adjustl(REGIME_NAMES(impactors%regime)))
            call netcdf_io_check( nf90_put_var(nc%id, nc%regime_varid, charstring,             start=[1, eslot], count=[NAMELEN, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var regime_varid" )
            call netcdf_io_check( nf90_put_var(nc%id, nc%Qloss_varid, impactors%Qloss,         start=[eslot] ), "collision_io_netcdf_write_frame_snapshot nf90_put_var Qloss_varid" )

            select type(before =>self%collider%before)
            class is (swiftest_nbody_system)
            select type(after =>self%collider%after)
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

                  ! This ensures that there first idslot will have the first body in it, not id 0 which is the default for a new idvals array
                  if (.not.allocated(nc%idvals)) allocate(nc%idvals, source=pl%id)
                  do i = 1, npl
                     call nc%find_idslot(pl%id(i), idslot) 
                     call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid,     pl%id(i),     start=[   idslot              ]), "collision_io_netcdf_write_frame_snapshot nf90_put_var id_varid"  )
                     charstring = trim(adjustl(pl%info(i)%name))
                     call netcdf_io_check( nf90_put_var(nc%id, nc%name_varid,   charstring,   start=[1, idslot              ], count=[NAMELEN, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var name_varid"  )
                     charstring = trim(adjustl(pl%info(i)%particle_type))
                     call netcdf_io_check( nf90_put_var(nc%id, nc%ptype_varid,  charstring,   start=[1, idslot, stage, eslot], count=[NAMELEN, 1, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var particle_type_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%rh_varid,     pl%rh(:,i),   start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var rh_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid,     pl%vh(:,i),   start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var vh_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%Gmass_varid,  pl%Gmass(i),  start=[   idslot, stage, eslot]), "collision_io_netcdf_write_frame_snapshot nf90_put_var Gmass_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%radius_varid, pl%radius(i), start=[   idslot, stage, eslot]), "collision_io_netcdf_write_frame_snapshot nf90_put_var radius_varid"  )
                     if (param%lrotation) then
                        call netcdf_io_check( nf90_put_var(nc%id, nc%Ip_varid,     pl%Ip(:,i),   start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var Ip_varid"  )
                        call netcdf_io_check( nf90_put_var(nc%id, nc%rot_varid,    pl%rot(:,i)*RAD2DEG,  start=[1, idslot, stage, eslot], count=[NDIM,1,1,1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var rotx_varid"  )
                     end if
                  end do
               end do
            end select
            end select

            if (param%lenergy) then
               call netcdf_io_check( nf90_put_var(nc%id, nc%ke_orb_varid,  collider%ke_orbit(:), start=[   1, eslot], count=[      2, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var ke_orb_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%ke_spin_varid, collider%ke_spin(:),  start=[   1, eslot], count=[      2, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var ke_spin_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%pe_varid,      collider%pe(:),       start=[   1, eslot], count=[      2, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var pe_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%be_varid,      collider%be(:),       start=[   1, eslot], count=[      2, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var pe_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%te_varid,      collider%te(:),       start=[   1, eslot], count=[      2, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var pe_varid tefore" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%L_orbit_varid,   collider%L_orbit(:,:), start=[1, 1, eslot], count=[NDIM, 2, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var L_orbit_varid before" )
               call netcdf_io_check( nf90_put_var(nc%id, nc%L_spin_varid,   collider%L_spin(:,:),  start=[1, 1, eslot], count=[NDIM, 2, 1]), "collision_io_netcdf_write_frame_snapshot nf90_put_var L_spin_varid before" )
            end if
      
            call netcdf_io_check( nf90_set_fill(nc%id, old_mode, tmp) )
         end associate
      end select
      return
   end subroutine collision_io_netcdf_write_frame_snapshot

end submodule s_collision_io