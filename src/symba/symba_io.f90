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
                  read(param_value, *) param%collision_save
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

         if ((param%collision_save /= "NONE") .and. (param%collision_save /= "TRAJECTORY") .and. (param%collision_save /= "CLOSEST")) then
            write(iomsg,*) 'Invalid collision_save parameter: ',trim(adjustl(param%out_type))
            write(iomsg,*) 'Valid options are NONE, TRAJECTORY, or CLOSEST'
            iostat = -1
            return
         end if
         param%lencounter_save = (param%encounter_save == "TRAJECTORY") .or. (param%encounter_save == "CLOSEST") .or. &
                                 (param%collision_save == "TRAJECTORY") .or. (param%collision_save == "CLOSEST") 

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


   module subroutine symba_io_write_discard(self, param)
      !! author: David A. Minton
      !!
      !! Write the metadata of the discarded body to the output file 
      implicit none
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      ! Internals

      associate(pl => self%pl, npl => self%pl%nbody, pl_adds => self%pl_adds)

         if (self%tp_discards%nbody > 0) call self%tp_discards%write_info(param%system_history%nc, param)
         select type(pl_discards => self%pl_discards)
         class is (symba_merger)
            if (pl_discards%nbody == 0) return

            call pl_discards%write_info(param%system_history%nc, param)
         end select
      end associate

      return

   end subroutine symba_io_write_discard

end submodule s_symba_io

