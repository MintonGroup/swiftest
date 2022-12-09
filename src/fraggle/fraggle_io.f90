!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle_classes) s_fraggle_io
   use swiftest

contains


   module subroutine fraggle_io_initialize_output(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF fragment history file system. This is a simplified version of the main simulation output NetCDF file, but with fewer variables.
      use, intrinsic :: ieee_arithmetic
      use netcdf
      implicit none
      ! Arguments
      class(fraggle_io_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),   intent(in)    :: param   
      ! Internals
      integer(I4B) :: nvar, varid, vartype
      real(DP) :: dfill
      real(SP) :: sfill
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B) :: i, ndims

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

         call check( nf90_create(nc%file_name, NF90_NETCDF4, nc%id), "fraggle_io_initialize nf90_create" )

         ! Dimensions
         call check( nf90_def_dim(nc%id, nc%time_dimname, nc%time_dimsize, nc%time_dimid), "fraggle_io_initialize nf90_def_dim time_dimid" ) ! Simulation time dimension
         call check( nf90_def_dim(nc%id, nc%space_dimname, NDIM , nc%space_dimid), "fraggle_io_initialize nf90_def_dim space_dimid" )           ! 3D space dimension
         call check( nf90_def_dim(nc%id, nc%id_dimname, param%maxid, nc%id_dimid), "fraggle_io_initialize nf90_def_dim id_dimid" )       ! dimension to store particle id numbers
         call check( nf90_def_dim(nc%id, nc%str_dimname, NAMELEN, nc%str_dimid), "fraggle_io_initialize nf90_def_dim str_dimid"  )          ! Dimension for string variables (aka character arrays)
         call check( nf90_def_dim(nc%id, nc%stage_dimname, 2, nc%stage_dimid), "fraggle_io_initialize nf90_def_dim stage_dimid"  )          ! Dimension for stage variables (aka "before" vs. "after"

         ! Dimension coordinates
         call check( nf90_def_var(nc%id, nc%time_dimname, nc%out_type, nc%time_dimid, nc%time_varid), "fraggle_io_initialize nf90_def_var time_varid"  )
         call check( nf90_def_var(nc%id, nc%space_dimname, NF90_CHAR, nc%space_dimid, nc%space_varid), "fraggle_io_initialize nf90_def_var space_varid"  )
         call check( nf90_def_var(nc%id, nc%id_dimname, NF90_INT, nc%id_dimid, nc%id_varid), "fraggle_io_initialize nf90_def_var id_varid"  )
         call check( nf90_def_var(nc%id, nc%stage_dimname, NF90_CHAR, nc%stage_dimid, nc%stage_varid), "fraggle_io_initialize nf90_def_var stage_varid"  )
      

         ! Variables

         call check( nf90_def_var(nc%id, nc%name_varname,  NF90_CHAR,  &
            [nc%str_dimid,                   nc%id_dimid                                 ], nc%name_varid), "fraggle_io_initialize nf90_def_var name_varid")

         call check( nf90_def_var(nc%id, nc%loop_varname,  NF90_INT, &
            [                                                                 nc%time_dimid], nc%loop_varid), "fraggle_io_initialize nf90_def_var loop_varid")   

         call check( nf90_def_var(nc%id, nc%rh_varname,nc%out_type,&
            [              nc%space_dimid,   nc%id_dimid,   nc%stage_dimid,   nc%time_dimid], nc%rh_varid), "fraggle_io_initialize nf90_def_var rh_varid")

         call check( nf90_def_var(nc%id, nc%vh_varname, nc%out_type,&
            [              nc%space_dimid,   nc%id_dimid,   nc%stage_dimid,   nc%time_dimid], nc%vh_varid), "fraggle_io_initialize nf90_def_var vh_varid")

         call check( nf90_def_var(nc%id, nc%Gmass_varname,  nc%out_type,&
            [                                nc%id_dimid,   nc%stage_dimid,   nc%time_dimid], nc%Gmass_varid), "fraggle_io_initialize nf90_def_var Gmass_varid")


         call check( nf90_def_var(nc%id, nc%radius_varname, nc%out_type,&
            [                                nc%id_dimid,    nc%stage_dimid,  nc%time_dimid], nc%radius_varid), "fraggle_io_initialize nf90_def_var radius_varid")

         call check( nf90_def_var(nc%id, nc%Ip_varname, nc%out_type,&
            [              nc%space_dimid,   nc%id_dimid,    nc%stage_dimid,  nc%time_dimid], nc%Ip_varid), "fraggle_io_initialize nf90_def_var Ip_varid")

         call check( nf90_def_var(nc%id, nc%rot_varname, nc%out_type,&
            [              nc%space_dimid,   nc%id_dimid,    nc%stage_dimid,  nc%time_dimid], nc%rot_varid), "fraggle_io_initialize nf90_def_var rot_varid")

         call check( nf90_def_var(nc%id, nc%ke_orb_varname, nc%out_type,&
            [                                                nc%stage_dimid,  nc%time_dimid], nc%KE_orb_varid), "netcdf_initialize_output nf90_def_var KE_orb_varid")

         call check( nf90_def_var(nc%id, nc%ke_spin_varname, nc%out_type,&
            [                                                nc%stage_dimid,  nc%time_dimid], nc%KE_spin_varid), "netcdf_initialize_output nf90_def_var KE_spin_varid"  )

         call check( nf90_def_var(nc%id, nc%pe_varname,&
                                         nc%out_type,&
            [                                                nc%stage_dimid,  nc%time_dimid], nc%PE_varid), "netcdf_initialize_output nf90_def_var PE_varid"  )

         call check( nf90_def_var(nc%id, nc%L_orb_varname, nc%out_type, &
            [              nc%space_dimid,                   nc%stage_dimid,  nc%time_dimid], nc%L_orb_varid), "netcdf_initialize_output nf90_def_var L_orb_varid"  )

         call check( nf90_def_var(nc%id, nc%L_spin_varname,  nc%out_type,&
            [              nc%space_dimid,                   nc%stage_dimid,  nc%time_dimid], nc%L_spin_varid), "netcdf_initialize_output nf90_def_var L_spin_varid"  )



         call check( nf90_inquire(nc%id, nVariables=nvar), "fraggle_io_initialize nf90_inquire nVariables"  )
         do varid = 1, nvar
            call check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), "fraggle_io_initialize nf90_inquire_variable"  )
            select case(vartype)
            case(NF90_INT)
               call check( nf90_def_var_fill(nc%id, varid, 0, NF90_FILL_INT), "fraggle_io_initialize nf90_def_var_fill NF90_INT"  )
            case(NF90_FLOAT)
               call check( nf90_def_var_fill(nc%id, varid, 0, sfill), "fraggle_io_initialize nf90_def_var_fill NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call check( nf90_def_var_fill(nc%id, varid, 0, dfill), "fraggle_io_initialize nf90_def_var_fill NF90_DOUBLE"  )
            case(NF90_CHAR)
               call check( nf90_def_var_fill(nc%id, varid, 0, 0), "fraggle_io_initialize nf90_def_var_fill NF90_CHAR"  )
            end select
         end do
         ! Take the file out of define mode
         call check( nf90_enddef(nc%id), "fraggle_io_initialize nf90_enddef"  )

         ! Add in the space and stage dimension coordinates
         call check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), "fraggle_io_initialize nf90_put_var space"  )
         call check( nf90_put_var(nc%id, nc%stage_varid, nc%stage_coords, start=[1], count=[2]), "fraggle_io_initialize nf90_put_var stage"  )

         ! Pre-fill id slots with ids
         call check( nf90_put_var(nc%id, nc%id_varid, [(-1,i=1,param%maxid)], start=[1], count=[param%maxid]), "fraggle_io_initialize nf90_put_var pl id_varid"  )
      end associate

      return

      667 continue
      write(*,*) "Error creating fragmentation output file. " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine fraggle_io_initialize_output


   module subroutine fraggle_io_write_frame(self, nc, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of a collision result
      use netcdf
      implicit none
      ! Arguments
      class(fraggle_encounter_snapshot), intent(in)    :: self   !! Swiftest encounter structure
      class(encounter_io_parameters),    intent(inout) :: nc    !! Parameters used to identify a particular encounter io NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)           :: i, tslot, idslot, old_mode, npl
      character(len=NAMELEN) :: charstring

      tslot = param%ioutput
      associate(pl => self%pl, colliders => self%colliders, fragments => self%fragments)

         call check( nf90_set_fill(nc%id, nf90_nofill, old_mode), "fraggle_io_write_frame nf90_set_fill"  )
   
         call check( nf90_put_var(nc%id, nc%time_varid, self%t, start=[tslot]), "fraggle_io_write_frame nf90_put_var time_varid"  )
         call check( nf90_put_var(nc%id, nc%loop_varid, int(self%iloop,kind=I4B), start=[tslot]), "fraggle_io_write_frame nf90_put_var pl loop_varid"  )

         ! npl = pl%nbody
         ! do i = 1, npl
         !    idslot = pl%id(i)
         !    call check( nf90_put_var(nc%id, nc%id_varid, pl%id(i),   start=[idslot]), "fraggle_io_write_frame nf90_put_var pl id_varid"  )
         !    call check( nf90_put_var(nc%id, nc%rh_varid, pl%rh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "fraggle_io_write_frame nf90_put_var pl rh_varid"  )
         !    call check( nf90_put_var(nc%id, nc%vh_varid, pl%vh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "fraggle_io_write_frame nf90_put_var pl vh_varid"  )
         !    call check( nf90_put_var(nc%id, nc%Gmass_varid, pl%Gmass(i), start=[idslot, tslot]), "fraggle_io_write_frame nf90_put_var pl Gmass_varid"  )

         !    if (param%lclose) call check( nf90_put_var(nc%id, nc%radius_varid, pl%radius(i), start=[idslot, tslot]), "fraggle_io_write_frame nf90_put_var pl radius_varid"  )

         !    if (param%lrotation) then
         !       call check( nf90_put_var(nc%id, nc%Ip_varid, pl%Ip(:,i), start=[1, idslot, tslot], count=[NDIM,1,1]), "fraggle_io_write_frame nf90_put_var pl Ip_varid"  )
         !       call check( nf90_put_var(nc%id, nc%rot_varid, pl%rot(:,i), start=[1,idslot, tslot], count=[NDIM,1,1]), "fraggle_io_write_frame nf90_put_var pl rotx_varid"  )
         !    end if

         !    charstring = trim(adjustl(pl%info(i)%name))
         !    call check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "fraggle_io_write_frame nf90_put_var pl name_varid"  )
         !    charstring = trim(adjustl(pl%info(i)%particle_type))
         !    call check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "fraggle_io_write_frame nf90_put_var pl particle_type_varid"  )
         ! end do

         call check( nf90_set_fill(nc%id, old_mode, old_mode) )
      end associate

      return
   end subroutine fraggle_io_write_frame


   module subroutine fraggle_io_log_pl(pl, param)
      !! author: David A. Minton
      !!
      !! Writes a single message to the fraggle log file
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(in) :: pl    !! Swiftest massive body object (only the new bodies generated in a collision)
      class(swiftest_parameters), intent(in) :: param !! Current swiftest run configuration parameters
      ! Internals
      integer(I4B) :: i
      character(STRMAX) :: errmsg

      return
      667 continue
      write(*,*) "Error writing Fraggle message to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_pl


   module subroutine fraggle_io_log_regime(colliders, frag)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the collisional regime determination
      implicit none
      ! Arguments
      class(fraggle_colliders),   intent(in) :: colliders !! Fraggle collider system object
      class(fraggle_fragments),   intent(in) :: frag      !! Fraggle fragment object
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Fraggle collisional regime determination results"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "True number of colliders : ",colliders%ncoll
      write(LUN, *) "Index list of true colliders  : ",colliders%idx(1:colliders%ncoll)
      select case(frag%regime) 
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
      write(LUN, *) "Energy loss                  : ", frag%Qloss
      write(LUN, *) "--------------------------------------------------------------------"
      close(LUN)

      return
      667 continue
      write(*,*) "Error writing Fraggle regime information to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_regime

end submodule s_fraggle_io