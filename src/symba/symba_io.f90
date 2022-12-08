!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (symba_classes) s_symba_io
   use swiftest
contains

   module subroutine symba_io_encounter_dump(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the time history of an encounter to file.
      implicit none
      ! Arguments
      class(symba_encounter_storage(*)),  intent(inout)        :: self   !! Encounter storage object
      class(swiftest_parameters),   intent(inout)        :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i

      ! Most of this is just temporary test code just to get something working. Eventually this should get cleaned up.
      
      do i = 1, self%nframes
         if (allocated(self%frame(i)%item)) then
            select type(snapshot => self%frame(i)%item)
            class is (symba_encounter_snapshot)
               self%nc%ienc_frame = i
               call snapshot%write_frame(self%nc,param)
            end select
         else
            exit
         end if
      end do


      return
   end subroutine symba_io_encounter_dump


   module subroutine symba_io_encounter_initialize_output(self, param)
      !! author: David A. Minton
      !!
      !! Initialize a NetCDF encounter file system. This is a simplified version of the main simulation output NetCDF file, but with fewer variables.
      use, intrinsic :: ieee_arithmetic
      use netcdf
      implicit none
      ! Arguments
      class(symba_io_encounter_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),           intent(in)    :: param   !! Current run configuration parameters
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
            self%out_type = NF90_FLOAT
         case("NETCDF_DOUBLE")
            self%out_type = NF90_DOUBLE
         end select

         ! Check if the file exists, and if it does, delete it
         inquire(file=nc%enc_file, exist=fileExists)
         if (fileExists) then
            open(unit=LUN, file=nc%enc_file, status="old", err=667, iomsg=errmsg)
            close(unit=LUN, status="delete")
         end if

         call check( nf90_create(nc%enc_file, NF90_NETCDF4, nc%id), "symba_io_encounter_initialize_output nf90_create" )

         ! Dimensions
         call check( nf90_def_dim(nc%id, nc%time_dimname, nc%time_dimsize, nc%time_dimid), "symba_io_encounter_initialize_output nf90_def_dim time_dimid" ) ! Simulation time dimension
         call check( nf90_def_dim(nc%id, nc%space_dimname, NDIM , nc%space_dimid), "symba_io_encounter_initialize_output nf90_def_dim space_dimid" )           ! 3D space dimension
         call check( nf90_def_dim(nc%id, nc%id_dimname, param%maxid, nc%id_dimid), "symba_io_encounter_initialize_output nf90_def_dim id_dimid" )       ! dimension to store particle id numbers
         call check( nf90_def_dim(nc%id, nc%str_dimname, NAMELEN, nc%str_dimid), "symba_io_encounter_initialize_output nf90_def_dim str_dimid"  )          ! Dimension for string variables (aka character arrays)

         ! Dimension coordinates
         call check( nf90_def_var(nc%id, nc%time_dimname, nc%out_type, nc%time_dimid, nc%time_varid), "symba_io_encounter_initialize_output nf90_def_var time_varid"  )
         call check( nf90_def_var(nc%id, nc%space_dimname, NF90_CHAR, nc%space_dimid, nc%space_varid), "symba_io_encounter_initialize_output nf90_def_var space_varid"  )
         call check( nf90_def_var(nc%id, nc%id_dimname, NF90_INT, nc%id_dimid, nc%id_varid), "symba_io_encounter_initialize_output nf90_def_var id_varid"  )
      
         ! Variables
         call check( nf90_def_var(nc%id, nc%name_varname, NF90_CHAR, [nc%str_dimid, nc%id_dimid], nc%name_varid), "symba_io_encounter_initialize_output nf90_def_var name_varid"  )
         call check( nf90_def_var(nc%id, nc%ptype_varname, NF90_CHAR, [nc%str_dimid, nc%id_dimid], nc%ptype_varid), "symba_io_encounter_initialize_output nf90_def_var ptype_varid"  )
         call check( nf90_def_var(nc%id, nc%rh_varname,  nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%rh_varid), "symba_io_encounter_initialize_output nf90_def_var rh_varid"  )
         call check( nf90_def_var(nc%id, nc%vh_varname,  nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%vh_varid), "symba_io_encounter_initialize_output nf90_def_var vh_varid"  )
         call check( nf90_def_var(nc%id, nc%Gmass_varname, nc%out_type, [nc%id_dimid, nc%time_dimid], nc%Gmass_varid), "symba_io_encounter_initialize_output nf90_def_var Gmass_varid"  )
         call check( nf90_def_var(nc%id, nc%level_varname, NF90_INT, [nc%id_dimid, nc%time_dimid], nc%level_varid), "symba_io_encounter_initialize_output nf90_def_var level_varid"  )
         if (param%lclose) then
            call check( nf90_def_var(nc%id, nc%radius_varname, nc%out_type, [nc%id_dimid, nc%time_dimid], nc%radius_varid), "symba_io_encounter_initialize_output nf90_def_var radius_varid"  )
         end if
         if (param%lrotation) then
            call check( nf90_def_var(nc%id, nc%Ip_varname, nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%Ip_varid), "symba_io_encounter_initialize_output nf90_def_var Ip_varid"  )
            call check( nf90_def_var(nc%id, nc%rot_varname, nc%out_type, [nc%space_dimid, nc%id_dimid, nc%time_dimid], nc%rot_varid), "symba_io_encounter_initialize_output nf90_def_var rot_varid"  )
         end if

         call check( nf90_inquire(nc%id, nVariables=nvar), "symba_io_encounter_initialize_output nf90_inquire nVariables"  )
         do varid = 1, nvar
            call check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), "symba_io_encounter_initialize_output nf90_inquire_variable"  )
            select case(vartype)
            case(NF90_INT)
               call check( nf90_def_var_fill(nc%id, varid, 0, NF90_FILL_INT), "symba_io_encounter_initialize_output nf90_def_var_fill NF90_INT"  )
            case(NF90_FLOAT)
               call check( nf90_def_var_fill(nc%id, varid, 0, sfill), "symba_io_encounter_initialize_output nf90_def_var_fill NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call check( nf90_def_var_fill(nc%id, varid, 0, dfill), "symba_io_encounter_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
            case(NF90_CHAR)
               call check( nf90_def_var_fill(nc%id, varid, 0, 0), "symba_io_encounter_initialize_output nf90_def_var_fill NF90_CHAR"  )
            end select
         end do

         ! Take the file out of define mode
         call check( nf90_enddef(nc%id), "symba_io_encounter_initialize_output nf90_enddef"  )

         ! Add in the space dimension coordinates
         call check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), "symba_io_encounter_initialize_output nf90_put_var space"  )
      end associate

      return

      667 continue
      write(*,*) "Error creating encounter output file. " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine symba_io_encounter_initialize_output


   module subroutine symba_io_encounter_write_frame(self, nc, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of an encounter trajectory.
      use netcdf
      implicit none
      ! Arguments
      class(symba_encounter_snapshot),      intent(in)    :: self  !! Swiftest encounter structure
      class(symba_io_encounter_parameters), intent(inout) :: nc    !! Parameters used to identify a particular encounter io NetCDF dataset
      class(swiftest_parameters),           intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B)           :: i,  tslot, idslot, old_mode, npl, ntp
      character(len=NAMELEN) :: charstring

      call check( nf90_set_fill(nc%id, nf90_nofill, old_mode), "symba_io_encounter_write_frame nf90_set_fill"  )

      tslot = self%tslot
      call check( nf90_put_var(nc%id, nc%time_varid, self%t, start=[tslot]), "symba_io_encounter_write_frame nf90_put_var time_varid"  )

      associate(pl => self%pl, tp => self%tp)
         npl = pl%nbody
         do i = 1, npl
            idslot = pl%id(i)
            call check( nf90_put_var(nc%id, nc%id_varid, pl%id(i), start=[idslot]), "symba_io_encounter_write_frame nf90_put_var pl id_varid"  )
            call check( nf90_put_var(nc%id, nc%rh_varid, pl%rh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "symba_io_encounter_write_frame nf90_put_var pl rh_varid"  )
            call check( nf90_put_var(nc%id, nc%vh_varid, pl%vh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "symba_io_encounter_write_frame nf90_put_var pl vh_varid"  )
            call check( nf90_put_var(nc%id, nc%Gmass_varid, pl%Gmass(i), start=[idslot, tslot]), "symba_io_encounter_write_frame nf90_put_var pl Gmass_varid"  )
            call check( nf90_put_var(nc%id, nc%level_varid, pl%levelg(i), start=[idslot, tslot]), "symba_io_encounter_write_frame nf90_put_var pl level_varid"  )

            if (param%lclose) call check( nf90_put_var(nc%id, nc%radius_varid, pl%radius(i), start=[idslot, tslot]), "symba_io_encounter_write_frame nf90_put_var pl radius_varid"  )

            if (param%lrotation) then
               call check( nf90_put_var(nc%id, nc%Ip_varid, pl%Ip(:,i), start=[1, idslot, tslot], count=[NDIM,1,1]), "symba_io_encounter_write_frame nf90_put_var pl Ip_varid"  )
               call check( nf90_put_var(nc%id, nc%rot_varid, pl%rot(:,i), start=[1,idslot, tslot], count=[NDIM,1,1]), "symba_io_encounter_write_frame nf90_put_var pl rotx_varid"  )
            end if
 
            charstring = trim(adjustl(pl%info(i)%name))
            call check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "symba_io_encounter_write_frame nf90_put_var pl name_varid"  )
            charstring = trim(adjustl(pl%info(i)%particle_type))
            call check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "symba_io_encounter_write_frame nf90_put_var pl particle_type_varid"  )
         end do

         ntp = tp%nbody
         do i = 1, ntp
            idslot = tp%id(i)
            call check( nf90_put_var(nc%id, nc%id_varid, tp%id(i), start=[idslot]), "symba_io_encounter_write_frame nf90_put_var tp id_varid"  )
            call check( nf90_put_var(nc%id, nc%rh_varid, tp%rh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "symba_io_encounter_write_frame nf90_put_var tp rh_varid"  )
            call check( nf90_put_var(nc%id, nc%vh_varid, tp%vh(:,i), start=[1,idslot,tslot], count=[NDIM,1,1]), "symba_io_encounter_write_frame nf90_put_var tp vh_varid"  )

            charstring = trim(adjustl(tp%info(i)%name))
            call check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "symba_io_encounter_write_frame nf90_put_var tp name_varid"  )
            charstring = trim(adjustl(tp%info(i)%particle_type))
            call check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "symba_io_encounter_write_frame nf90_put_var tp particle_type_varid"  )
         end do
      end associate

      call check( nf90_set_fill(nc%id, old_mode, old_mode) )

      return
   end subroutine symba_io_encounter_write_frame


   module subroutine symba_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in parameters specific to the SyMBA integrator, then calls the base io_param_reader.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
      !! Adapted from Martin Duncan's Swift routine io_init_param.f
      implicit none
      ! Arguments
      class(symba_parameters), intent(inout) :: self       !! Collection of parameters
      integer,                 intent(in)    :: unit       !! File unit number
      character(len=*),        intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                           !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
      character(len=*),        intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
      integer,                 intent(out)   :: iostat     !! IO status code
      character(len=*),        intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      ! internals
      integer(I4B)                   :: ilength, ifirst, ilast  !! Variables used to parse input file
      character(STRMAX)              :: line                    !! Line of the input file
      character (len=:), allocatable :: line_trim,param_name, param_value !! Strings used to parse the param file
      integer(I4B)                   :: nseeds, nseeds_from_file, i
      logical                        :: seed_set = .false.      !! Is the random seed set in the input file?
      character(len=*),parameter     :: linefmt = '(A)'

      associate(param => self)
         open(unit = unit, file = param%param_file_name, status = 'old', err = 667, iomsg = iomsg)
         call random_seed(size = nseeds)
         if (allocated(param%seed)) deallocate(param%seed)
         allocate(param%seed(nseeds))
         do
            read(unit = unit, fmt = linefmt, iostat = iostat, end = 1, err = 667, iomsg = iomsg) line
            line_trim = trim(adjustl(line))
            ilength = len(line_trim)
            if ((ilength /= 0)) then 
               ifirst = 1
               ! Read the pair of tokens. The first one is the parameter name, the second is the value.
               param_name = io_get_token(line_trim, ifirst, ilast, iostat)
               if (param_name == '') cycle ! No parameter name (usually because this line is commented out)
               call io_toupper(param_name)
               ifirst = ilast + 1
               param_value = io_get_token(line_trim, ifirst, ilast, iostat)
               select case (param_name)
               case ("OUT_STAT") ! We need to duplicate this from the standard io_param_reader in order to make sure that the restart flag gets set properly in SyMBA
                  call io_toupper(param_value)
                  param%out_stat = param_value 
               case ("FRAGMENTATION")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == "T") self%lfragmentation = .true.
               case ("GMTINY")
                  read(param_value, *) param%GMTINY
               case ("MIN_GMFRAG")
                  read(param_value, *) param%min_GMfrag
               case ("ENCOUNTER_SAVE")
                  call io_toupper(param_value)
                  read(param_value, *) param%encounter_save
               case ("FRAGMENTATION_SAVE")
                  call io_toupper(param_value)
                  read(param_value, *) param%fragmentation_save
               case("SEED")
                  read(param_value, *) nseeds_from_file
                  ! Because the number of seeds can vary between compilers/systems, we need to make sure we can handle cases in which the input file has a different
                  ! number of seeds than the current system. If the number of seeds in the file is smaller than required, we will use them as a source to fill in the missing elements.
                  ! If the number of seeds in the file is larger than required, we will truncate the seed array.
                  if (nseeds_from_file > nseeds) then
                     nseeds = nseeds_from_file
                     deallocate(param%seed)
                     allocate(param%seed(nseeds))
                     do i = 1, nseeds
                        ifirst = ilast + 2
                        param_value = io_get_token(line, ifirst, ilast, iostat) 
                        read(param_value, *) param%seed(i)
                     end do
                  else ! Seed array in file is too small
                     do i = 1, nseeds_from_file
                        ifirst = ilast + 2
                        param_value = io_get_token(line, ifirst, ilast, iostat) 
                        read(param_value, *) param%seed(i)
                     end do
                     param%seed(nseeds_from_file+1:nseeds) = [(param%seed(1) - param%seed(nseeds_from_file) + i, &
                                                               i=nseeds_from_file+1, nseeds)]
                  end if
                  seed_set = .true.
               end select
            end if
         end do
         1 continue
         close(unit)

         param%lrestart = (param%out_stat == "APPEND")

         if (self%GMTINY < 0.0_DP) then
            write(iomsg,*) "GMTINY invalid or not set: ", self%GMTINY
            iostat = -1
            return
         end if

         if (param%lfragmentation) then
            if (seed_set) then
               call random_seed(put = param%seed)
            else
               call random_seed(get = param%seed)
            end if
            if (param%min_GMfrag < 0.0_DP) param%min_GMfrag = param%GMTINY
         end if

         ! All reporting of collision information in SyMBA (including mergers) is now recorded in the Fraggle logfile
         call io_log_start(param, FRAGGLE_LOG_OUT, "Fraggle logfile")

         if ((param%encounter_save /= "NONE") .and. (param%encounter_save /= "TRAJECTORY") .and. (param%encounter_save /= "CLOSEST")) then
            write(iomsg,*) 'Invalid encounter_save parameter: ',trim(adjustl(param%out_type))
            write(iomsg,*) 'Valid options are NONE, TRAJECTORY, or CLOSEST'
            iostat = -1
            return
         end if

         if ((param%fragmentation_save /= "NONE") .and. (param%fragmentation_save /= "TRAJECTORY") .and. (param%fragmentation_save /= "CLOSEST")) then
            write(iomsg,*) 'Invalid fragmentation_save parameter: ',trim(adjustl(param%out_type))
            write(iomsg,*) 'Valid options are NONE, TRAJECTORY, or CLOSEST'
            iostat = -1
            return
         end if
         param%lencounter_save = (param%encounter_save == "TRAJECTORY") .or. (param%encounter_save == "CLOSEST") .or. &
                                 (param%fragmentation_save == "TRAJECTORY") .or. (param%fragmentation_save == "CLOSEST") 

         ! Call the base method (which also prints the contents to screen)
         call io_param_reader(param, unit, iotype, v_list, iostat, iomsg) 
      end associate

      iostat = 0

      return
      667 continue
      write(*,*) "Error reading SyMBA parameters in param file: ", trim(adjustl(iomsg))
   end subroutine symba_io_param_reader


   module subroutine symba_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: David A. Minton
      !!
      !! Dump integration parameters specific to SyMBA to file and then call the base io_param_writer method.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_param.f
      implicit none
      ! Arguments
      class(symba_parameters),intent(in)    :: self      !! Collection of SyMBA parameters
      integer,                intent(in)    :: unit      !! File unit number
      character(len=*),       intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                         !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
      integer,                intent(in)    :: v_list(:) !! Not used in this procedure
      integer,                intent(out)   :: iostat    !! IO status code
      character(len=*),       intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      ! Internals
      integer(I4B) :: nseeds

      associate(param => self)
         call io_param_writer(param, unit, iotype, v_list, iostat, iomsg) 

         ! Special handling is required for writing the random number seed array as its size is not known until runtime
         ! For the "SEED" parameter line, the first value will be the size of the seed array and the rest will be the seed array elements
         call io_param_writer_one("GMTINY",param%GMTINY, unit)
         call io_param_writer_one("MIN_GMFRAG",param%min_GMfrag, unit)
         call io_param_writer_one("FRAGMENTATION",param%lfragmentation, unit)
         if (param%lfragmentation) then
            nseeds = size(param%seed)
            call io_param_writer_one("SEED", [nseeds, param%seed(:)], unit)
         end if

         iostat = 0
      end associate

      return
      667 continue
      write(*,*) "Error writing parameter file for SyMBA: " // trim(adjustl(iomsg))
   end subroutine symba_io_param_writer


   module subroutine symba_io_start_encounter(self, param, t)
      !! author: David A. Minton
      !!
      !! Initializes the new encounter and/or fragmentation history
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(symba_parameters),    intent(inout) :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! Current simulation time

      if (.not. allocated(self%encounter_history)) then
         allocate(symba_encounter_storage :: self%encounter_history)
         allocate(symba_io_encounter_parameters :: self%encounter_history%nc)
      end if
      call self%encounter_history%reset()

      ! Empty out the time slot array for the next pass
      self%encounter_history%tvals(:) = huge(1.0_DP)

      ! Take the snapshot at the start of the encounter
      call self%snapshot(param, t) 

      return
   end subroutine symba_io_start_encounter


   module subroutine symba_io_stop_encounter(self, param, t)
      !! author: David A. Minton
      !!
      !! Saves the encounter and/or fragmentation data to file(s)  
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(symba_parameters),    intent(inout) :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! Current simulation time
      ! Internals
      integer(I4B) :: i

      ! Create and save the output file for this encounter
      
      ! Figure out how many time slots we need
      do i = 1, self%encounter_history%nframes
         if (self%t + param%dt <= self%encounter_history%tvals(i)) then
            self%encounter_history%nc%time_dimsize = i
            exit
         end if
      end do

      write(self%encounter_history%nc%enc_file, '("encounter_",I0.6,".nc")') param%iloop

      call self%encounter_history%nc%initialize(param)
      call self%encounter_history%dump(param)
      call self%encounter_history%nc%close()
      call self%encounter_history%reset()

      return
   end subroutine symba_io_stop_encounter


   module subroutine symba_io_write_discard(self, param)
      !! author: David A. Minton
      !!
      !! Write the metadata of the discarded body to the output file 
      implicit none
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      ! Internals

      associate(pl => self%pl, npl => self%pl%nbody, pl_adds => self%pl_adds)

         if (self%tp_discards%nbody > 0) call self%tp_discards%write_info(param%nc, param)
         select type(pl_discards => self%pl_discards)
         class is (symba_merger)
            if (pl_discards%nbody == 0) return

            call pl_discards%write_info(param%nc, param)
         end select
      end associate

      return

   end subroutine symba_io_write_discard

end submodule s_symba_io

