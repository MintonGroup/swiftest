!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_io
   use swiftest

contains

   module subroutine io_compact_output(self, param, timer)
      !! author: David Minton
      !!
      !! Generates the terminal output displayed when display_style is set to COMPACT. This is used by the Python driver to 
      !! make nice-looking progress reports.
      implicit none

      interface fmt
         !! author: David Minton
         !!
         !! Formats a pair of variables and corresponding values for the compact display output. Generic interface for different variable types to format.
         procedure :: fmt_I4B, fmt_I8B, fmt_DP
      end interface

      ! Arguments
      class(swiftest_nbody_system), intent(in) :: self  !! Swiftest nbody system object   
      class(swiftest_parameters),   intent(in) :: param !! Input colleciton of user-defined parameters
      class(*),                     intent(in) :: timer !! Object used for computing elapsed wall time  (must be unlimited polymorphic because the walltimer module requires swiftest_classes)
      ! Internals
      character(len=:), allocatable :: formatted_output

      select type(timer)
      class is (walltimer)
         formatted_output = fmt("ILOOP",param%iloop) // fmt("T",param%t) // fmt("NPL",self%pl%nbody) // fmt("NTP",self%tp%nbody) 
         select type(pl => self%pl)
         class is (symba_pl)
            formatted_output = formatted_output // fmt("NPLM",pl%nplm)
         end select
         if (param%lenergy) then
            formatted_output = formatted_output // fmt("LTOTERR",self%Ltot_error) // fmt("ETOTERR",self%Etot_error) // fmt("MTOTERR",self%Mtot_error) &
                             // fmt("KEOERR",self%ke_orbit_error) // fmt("PEERR",self%pe_error) // fmt("EORBERR",self%Eorbit_error) &
                             // fmt("EUNTRERR",self%Euntracked_error) // fmt("LESCERR",self%Lescape_error) // fmt("MESCERR",self%Mescape_error)
            if (param%lclose) formatted_output = formatted_output // fmt("ECOLLERR",self%Ecoll_error)
            if (param%lrotation) formatted_output = formatted_output // fmt("KESPINERR",self%ke_spin_error) // fmt("LSPINERR",self%Lspin_error) 
         end if

         if (.not. timer%main_is_started) then ! This is the start of a new run
            formatted_output =  formatted_output // fmt("WT",0.0_DP) // fmt("IWT",0.0_DP) // fmt("WTPS",0.0_DP) 
         else
            formatted_output = formatted_output // fmt("WT",timer%wall_main) // fmt("IWT",timer%wall_step) // fmt("WTPS",timer%wall_per_substep)
         end if
         write(*,*) formatted_output
      end select
      return

      contains

         function fmt_I4B(varname,val) result(pair_string)
            implicit none
            ! Arguments
            character(*), intent(in) :: varname !! The variable name of the pair
            integer(I4B), intent(in) :: val !! A 4-byte integer value
            ! Result
            character(len=:), allocatable :: pair_string
            ! Internals
            character(len=24) :: str_value
      
            write(str_value,*) val
            pair_string = trim(adjustl(varname)) // " " // trim(adjustl(str_value)) // ";"

            return 
         end function fmt_I4B

         function fmt_I8B(varname, val) result(pair_string)
            implicit none
            ! Arguments
            character(*), intent(in) :: varname !! The variable name of the pair
            integer(I8B), intent(in) :: val     !! An 8-byte integer value
            ! Result
            character(len=:), allocatable :: pair_string
            ! Internals
            character(len=24) :: str_value
      
            write(str_value,*) val
            pair_string = trim(adjustl(varname)) // " " // trim(adjustl(str_value)) // ";"

            return 
         end function fmt_I8B

         function fmt_DP(varname, val) result(pair_string)
            implicit none
            ! Arguments
            character(*), intent(in) :: varname !! The variable name of the pair
            real(DP),     intent(in) :: val     !! A double precision floating point value
            ! Result
            character(len=:), allocatable :: pair_string
            ! Internals
            character(len=24) :: str_value
      
            write(str_value,'(ES24.16)') val
            pair_string = trim(adjustl(varname)) // " " // trim(adjustl(str_value)) // ";"

            return 
         end function fmt_DP

   end subroutine io_compact_output


   module subroutine io_conservation_report(self, param, lterminal)
      !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Reports the current state of energy, mass, and angular momentum conservation in a run
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self      !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param     !! Input colleciton of user-defined parameters
      logical,                      intent(in)    :: lterminal !! Indicates whether to output information to the terminal screen
      ! Internals
      real(DP), dimension(NDIM)       :: Ltot_now,  Lorbit_now,  Lspin_now
      real(DP)                        :: ke_orbit_now,  ke_spin_now,  pe_now,  Eorbit_now
      real(DP)                        :: Eorbit_error, Etot_error, Ecoll_error
      real(DP)                        :: GMtot_now
      real(DP)                        :: Lerror, Merror
      character(len=STRMAX)           :: errmsg
      integer(I4B), parameter         :: EGYIU = 72
      character(len=*), parameter     :: EGYTERMFMT = '(" DL/L0 = ", ES12.5 &
                                                         "; DEcollisions/|E0| = ", ES12.5, &
                                                         "; D(Eorbit+Ecollisions)/|E0| = ", ES12.5, &
                                                         "; DM/M0 = ", ES12.5)'

      associate(system => self, pl => self%pl, cb => self%cb, npl => self%pl%nbody, display_unit => param%display_unit)

         call pl%vb2vh(cb)
         call pl%xh2xb(cb)

         call system%get_energy_and_momentum(param) 
         ke_orbit_now = system%ke_orbit
         ke_spin_now = system%ke_spin
         pe_now = system%pe
         Lorbit_now(:) = system%Lorbit(:)
         Lspin_now(:) = system%Lspin(:)
         Eorbit_now = ke_orbit_now + ke_spin_now + pe_now
         Ltot_now(:) = system%Ltot(:) + system%Lescape(:)
         GMtot_now = system%GMtot + system%GMescape 

         if (param%lfirstenergy) then
            system%ke_orbit_orig = ke_orbit_now
            system%ke_spin_orig = ke_spin_now
            system%pe_orig = pe_now
            system%Eorbit_orig = Eorbit_now
            system%GMtot_orig = GMtot_now
            system%Lorbit_orig(:) = Lorbit_now(:)
            system%Lspin_orig(:) = Lspin_now(:)
            system%Ltot_orig(:) = Ltot_now(:)
            param%lfirstenergy = .false.
         end if

         if (.not.param%lfirstenergy) then 
            system%ke_orbit_error = (ke_orbit_now - system%ke_orbit_orig) / abs(system%Eorbit_orig)
            system%ke_spin_error = (ke_spin_now - system%ke_spin_orig) / abs(system%Eorbit_orig)
            system%pe_error = (pe_now - system%pe_orig) / abs(system%Eorbit_orig)
            system%Eorbit_error = (Eorbit_now - system%Eorbit_orig) / abs(system%Eorbit_orig)
            system%Ecoll_error = system%Ecollisions / abs(system%Eorbit_orig)
            system%Euntracked_error = system%Euntracked / abs(system%Eorbit_orig)
            system%Etot_error = (Eorbit_now - system%Ecollisions - system%Eorbit_orig - system%Euntracked) / abs(system%Eorbit_orig)

            system%Lorbit_error = norm2(Lorbit_now(:) - system%Lorbit_orig(:)) / norm2(system%Ltot_orig(:))
            system%Lspin_error = norm2(Lspin_now(:) - system%Lspin_orig(:)) / norm2(system%Ltot_orig(:))
            system%Lescape_error = norm2(system%Lescape(:)) / norm2(system%Ltot_orig(:))
            system%Ltot_error = norm2(Ltot_now(:) - system%Ltot_orig(:)) / norm2(system%Ltot_orig(:))
            system%Mescape_error = system%GMescape / system%GMtot_orig
            system%Mtot_error = (GMtot_now - system%GMtot_orig) / system%GMtot_orig
            if (lterminal) write(display_unit, EGYTERMFMT) system%Ltot_error, system%Ecoll_error, system%Etot_error,system%Mtot_error
            if (abs(system%Mtot_error) > 100 * epsilon(system%Mtot_error)) then
               write(*,*) "Severe error! Mass not conserved! Halting!"
               ! Save the frame of data to the bin file in the slot just after the present one for diagnostics
               param%ioutput = param%ioutput + 1_I8B
               call self%write_frame(param%nciu, param)
               call param%nciu%close()
               call util_exit(FAILURE)
            end if
         end if
      end associate

      return

      667 continue
      write(*,*) "Error writing energy and momentum tracking file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_conservation_report


   module subroutine io_dump_param(self, param_file_name)
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_param.f
      implicit none
      ! Arguments
      class(swiftest_parameters),intent(in) :: self    !! Output collection of parameters
      character(len=*),          intent(in) :: param_file_name !! Parameter input file name (i.e. param.in)
      ! Internals
      character(STRMAX)        :: errmsg !! Error message in UDIO procedure
      integer(I4B)             :: ierr

      open(unit = LUN, file = param_file_name, status='replace', form = 'FORMATTED', err = 667, iomsg = errmsg)
      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    due to compiler incompatabilities
      !write(LUN,'(DT)') param
      call self%writer(LUN, iotype = "none", v_list = [0], iostat = ierr, iomsg = errmsg)
      if (ierr == 0) then
         close(LUN, err = 667, iomsg = errmsg)
         return
      end if

      667 continue
      write(*,*) "Error opening parameter dump file " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_dump_param


   module subroutine io_dump_system(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the state of the system to files in case the simulation is interrupted.
      !! As a safety mechanism, there are two dump files that are written in alternating order
      !! so that if a dump file gets corrupted during writing, the user can restart from the older one.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      ! Internals
      class(swiftest_parameters), allocatable :: dump_param !! Local parameters variable used to parameters change input file names 
                                                            !! to dump file-specific values without changing the user-defined values
      integer(I4B), save            :: idx = 1              !! Index of current dump file. Output flips between 2 files for extra security
                                                            !! in case the program halts during writing
      character(len=:), allocatable :: param_file_name
   
      allocate(dump_param, source=param)
      param_file_name    = trim(adjustl(DUMP_PARAM_FILE(idx)))
      dump_param%in_form  = "XV"
      dump_param%out_stat = 'APPEND'
      dump_param%in_type = "NETCDF_DOUBLE"
      dump_param%in_netcdf = trim(adjustl(DUMP_NC_FILE(idx)))
      dump_param%nciu%id_chunk = self%pl%nbody + self%tp%nbody
      dump_param%nciu%time_chunk = 1
      dump_param%T0 = param%t

      call dump_param%dump(param_file_name)

      dump_param%out_form = "XV"
      dump_param%outfile = trim(adjustl(DUMP_NC_FILE(idx)))
      dump_param%ioutput = 0 
      call dump_param%nciu%initialize(dump_param)
      call self%write_frame(dump_param%nciu, dump_param)
      call dump_param%nciu%close()
      ! Syncrhonize the disk and memory buffer of the NetCDF file (e.g. commit the frame files stored in memory to disk) 
      call param%nciu%flush(param)

      idx = idx + 1
      if (idx > NDUMPFILES) idx = 1

      return
   end subroutine io_dump_system


   module subroutine io_get_args(integrator, param_file_name, display_style)
      !! author: David A. Minton
      !!
      !! Reads in the name of the parameter file from command line arguments. 
      implicit none
      ! Arguments
      character(len=:), intent(inout), allocatable :: integrator      !! Symbolic code of the requested integrator  
      character(len=:), intent(inout), allocatable :: param_file_name !! Name of the input parameters file
      character(len=:), intent(inout), allocatable :: display_style   !! Style of the output display {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD"
      ! Internals
      character(len=STRMAX), dimension(:), allocatable :: arg
      integer(I4B), dimension(:), allocatable :: ierr
      integer :: i,narg
      character(len=*),parameter    :: linefmt = '(A)'

      narg = command_argument_count() 
      if (narg > 0) then
         allocate(arg(narg),ierr(narg))
         do i = 1,narg
            call get_command_argument(i, arg(i), status = ierr(i))
         end do
         if (any(ierr /= 0)) call util_exit(USAGE)
      else
         call util_exit(USAGE)
      end if
   
      if (narg == 1) then
         if (arg(1) == '-v' .or. arg(1) == '--version') then
            call util_version() 
         else if (arg(1) == '-h' .or. arg(1) == '--help') then
            call util_exit(HELP)
         else
            call util_exit(USAGE)
         end if
      else if (narg >= 2) then
         call io_toupper(arg(1))
         select case(arg(1))
         case('BS')
            integrator = BS
         case('HELIO')
            integrator = HELIO
         case('RA15')
            integrator = RA15
         case('TU4')
            integrator = TU4
         case('WHM')
            integrator = WHM
         case('RMVS')
            integrator = RMVS
         case('SYMBA')
            integrator = SYMBA
         case('RINGMOONS')
            integrator = RINGMOONS
         case default
            integrator = UNKNOWN_INTEGRATOR
            write(*,*) trim(adjustl(arg(1))) // ' is not a valid integrator.'
            call util_exit(USAGE)
         end select
         param_file_name = trim(adjustl(arg(2)))
      end if

      if (narg == 2) then
         display_style = "STANDARD"
      else if (narg == 3) then
         call io_toupper(arg(3))
         display_style = trim(adjustl(arg(3)))
      else
         call util_exit(USAGE)
      end if

      return
   end subroutine io_get_args


   module function io_get_token(buffer, ifirst, ilast, ierr) result(token)
      !! author: David A. Minton
      !!
      !! Retrieves a character token from an input string. Here a token is defined as any set of contiguous non-blank characters not 
      !! beginning with or containing "!". If "!" is present, any remaining part of the buffer including the "!" is ignored
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_get_token.f90
      implicit none
      ! Arguments
      character(len=*), intent(in)    :: buffer         !! Input string buffer
      integer(I4B),     intent(inout) :: ifirst         !! Index of the buffer at which to start the search for a token
      integer(I4B),     intent(out)   :: ilast          !! Index of the buffer at the end of the returned token
      integer(I4B),     intent(out)   :: ierr           !! Error code
      ! Result
      character(len=:), allocatable   :: token          !! Returned token string
      ! Internals
      integer(I4B) :: i,ilength
   
      ilength = len(buffer)
   
      if (ifirst > ilength) then
          ilast = ifirst
          ierr = -1 !! Bad input
          token = ''
          return
      end if
      do i = ifirst, ilength
          if (buffer(i:i) /= ' ') exit
      end do
      if ((i > ilength) .or. (buffer(i:i) == '!')) then
          ifirst = i
          ilast = i
          ierr = -2 !! No valid token
          token = ''
          return
      end if
      ifirst = i
      do i = ifirst, ilength
          if ((buffer(i:i) == ' ') .or. (buffer(i:i) == '!')) exit
      end do
      ilast = i - 1
      ierr = 0
   
      token = buffer(ifirst:ilast)

      return
   end function io_get_token


   module subroutine io_log_one_message(file, message)
      !! author: David A. Minton
      !!
      !! Writes a single message to a log file
      implicit none
      ! Arguments
      character(len=*), intent(in) :: file   !! Name of file to log
      character(len=*), intent(in) :: message
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=LUN, file=trim(adjustl(file)), status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *) trim(adjustl(message)) 
      close(LUN)

      return
      667 continue
      write(*,*) "Error writing message to log file: " // trim(adjustl(errmsg))
   end subroutine io_log_one_message


   module subroutine io_log_start(param, file, header)
      !! author: David A. Minton
      !!
      !! Checks to see if a log file needs to be created if this is a new run, or appended if this is a restarted run
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(in) :: param  !! Current Swiftest run configuration parameters
      character(len=*),           intent(in) :: file   !! Name of file to log
      character(len=*),           intent(in) :: header !! Header to print at top of log file
      ! Internals
      character(STRMAX) :: errmsg
      logical           :: fileExists

      inquire(file=trim(adjustl(file)), exist=fileExists)
      if (.not.param%lrestart .or. .not.fileExists) then
         open(unit=LUN, file=file, status="REPLACE", err = 667, iomsg = errmsg)
         write(LUN, *, err = 667, iomsg = errmsg) trim(adjustl(header))
      end if
      close(LUN)

      return

      667 continue
      write(*,*) "Error writing log file: " // trim(adjustl(errmsg))
   end subroutine io_log_start


   module subroutine io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in parameters for the integration
      !! Currently this procedure does not work in user-defined derived-type input mode 
      !!    e.g. read(unit,'(DT)') param 
      !! as the newline characters are ignored in the input file when compiled in ifort.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
      !! Adapted from Martin Duncan's Swift routine io_init_param.f
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(inout) :: self       !! Collection of parameters
      integer, intent(in)                       :: unit       !! File unit number
      character(len=*), intent(in)              :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                              !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
      character(len=*), intent(in)              :: v_list(:)  !! The first element passes the integrator code to the reader
      integer, intent(out)                      :: iostat     !! IO status code
      character(len=*), intent(inout)           :: iomsg      !! Message to pass if iostat /= 0
      ! Internals
      logical                        :: t0_set = .false.                  !! Is the initial time set in the input file?
      logical                        :: tstart_set = .false.               !! Is the final time set in the input file?
      logical                        :: tstop_set = .false.               !! Is the final time set in the input file?
      logical                        :: dt_set = .false.                  !! Is the step size set in the input file?
      integer(I4B)                   :: ilength, ifirst, ilast, i         !! Variables used to parse input file
      character(STRMAX)              :: line                              !! Line of the input file
      character (len=:), allocatable :: line_trim,param_name, param_value !! Strings used to parse the param file
      character(*),parameter         :: linefmt = '(A)'                   !! Format code for simple text string
      

      ! Parse the file line by line, extracting tokens then matching them up with known parameters if possible
      associate(param => self) 
         open(unit = unit, file = param%param_file_name, status = 'old', err = 667, iomsg = iomsg)
         do
            read(unit = unit, fmt = linefmt, end = 1, err = 667, iomsg = iomsg) line
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
               case ("T0")
                  read(param_value, *, err = 667, iomsg = iomsg) param%t0
                  t0_set = .true.
               case ("TSTART")
                  read(param_value, *, err = 667, iomsg = iomsg) param%t0
                  tstart_set = .true.                  
               case ("TSTOP")
                  read(param_value, *, err = 667, iomsg = iomsg) param%tstop
                  tstop_set = .true.
               case ("DT")
                  read(param_value, *, err = 667, iomsg = iomsg) param%dt
                  dt_set = .true.
               case ("CB_IN")
                  param%incbfile = param_value
               case ("PL_IN")
                  param%inplfile = param_value
               case ("TP_IN")
                  param%intpfile = param_value
               case ("NC_IN")
                  param%in_netcdf = param_value
               case ("IN_TYPE")
                  call io_toupper(param_value)
                  param%in_type = param_value
               case ("IN_FORM")
                  call io_toupper(param_value)
                  param%in_form = param_value
               case ("ISTEP_OUT")
                  read(param_value, *) param%istep_out
               case ("BIN_OUT")
                  param%outfile = param_value
               case ("OUT_TYPE")
                  call io_toupper(param_value)
                  param%out_type = param_value
               case ("OUT_FORM")
                  call io_toupper(param_value)
                  param%out_form = param_value
               case ("OUT_STAT")
                  call io_toupper(param_value)
                  param%out_stat = param_value
               case ("DUMP_CADENCE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%dump_cadence
               case ("CHK_CLOSE")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lclose = .true.
               case ("CHK_RMIN")
                  read(param_value, *, err = 667, iomsg = iomsg) param%rmin
               case ("CHK_RMAX")
                  read(param_value, *, err = 667, iomsg = iomsg) param%rmax
               case ("CHK_EJECT")
                  read(param_value, *, err = 667, iomsg = iomsg) param%rmaxu
               case ("CHK_QMIN")
                  read(param_value, *, err = 667, iomsg = iomsg) param%qmin
               case ("CHK_QMIN_COORD")
                  call io_toupper(param_value)
                  param%qmin_coord = param_value
               case ("CHK_QMIN_RANGE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%qmin_alo
                  ifirst = ilast + 1
                  param_value = io_get_token(line, ifirst, ilast, iostat)
                  read(param_value, *, err = 667, iomsg = iomsg) param%qmin_ahi
               case ("EXTRA_FORCE")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lextra_force = .true.
               case ("BIG_DISCARD")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T' ) param%lbig_discard = .true.
               case ("RHILL_PRESENT")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T' ) param%lrhill_present = .true.
               case ("MU2KG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%MU2KG
               case ("TU2S")
                  read(param_value, *, err = 667, iomsg = iomsg) param%TU2S
               case ("DU2M")
                  read(param_value, *, err = 667, iomsg = iomsg) param%DU2M
               case ("ENERGY")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lenergy = .true.
               case ("GR")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lgr = .true. 
               case ("ROTATION")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lrotation = .true. 
               case ("TIDES")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%ltides = .true. 
               case ("INTERACTION_LOOPS")
                  call io_toupper(param_value)
                  param%interaction_loops = param_value
               case ("ENCOUNTER_CHECK_PLPL")
                  call io_toupper(param_value)
                  param%encounter_check_plpl = param_value
               case ("ENCOUNTER_CHECK_PLTP")
                  call io_toupper(param_value)
                  param%encounter_check_pltp = param_value
               case ("ENCOUNTER_CHECK")
                  call io_toupper(param_value)
                  param%encounter_check_plpl = param_value
                  param%encounter_check_pltp = param_value
               case ("FIRSTKICK")
                  call io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstkick = .false. 
               case ("FIRSTENERGY")
                  call io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstenergy = .false. 
               case("EORBIT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Eorbit_orig 
               case("GMTOT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%GMtot_orig 
               case("LTOT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Ltot_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Ltot_orig(i)
                  end do
               case("LORBIT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Lorbit_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Lorbit_orig(i)
                  end do
               case("LSPIN_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Lspin_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Lspin_orig(i)
                  end do
               case("LESCAPE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Lescape(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Lescape(i)
                  end do
               case("GMESCAPE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%GMescape 
               case("ECOLLISIONS")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Ecollisions
               case("EUNTRACKED")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Euntracked
               case ("MAXID")
                  read(param_value, *, err = 667, iomsg = iomsg) param%maxid 
               case ("MAXID_COLLISION")
                  read(param_value, *, err = 667, iomsg = iomsg) param%maxid_collision
               case ("RESTART")
                  if (param_value == "NO" .or. param_value == 'F') then
                     param%lrestart = .false. 
                  else if (param_value == "YES" .or. param_value == 'T') then
                     param%lrestart = .true.
                  end if 
               case ("NPLMAX", "NTPMAX", "GMTINY", "MIN_GMFRAG", "FRAGMENTATION", "SEED", "YARKOVSKY", "YORP") ! Ignore SyMBA-specific, not-yet-implemented, or obsolete input parameters
               case default
                  write(*,*) "Ignoring unknown parameter -> ",param_name
               end select
            end if
         end do
         1 continue
         close(unit)
         iostat = 0

         ! Do basic sanity checks on the input values
         if ((.not. t0_set) .or. (.not. tstop_set) .or. (.not. dt_set)) then
            write(iomsg,*) 'Valid simulation time not set'
            iostat = -1
            return
         end if
         if (param%dt <= 0.0_DP) then
            write(iomsg,*) 'Invalid timestep: '
            iostat = -1
            return
         end if
         if (param%inplfile == "") then
            write(iomsg,*) 'No valid massive body file in input file'
            iostat = -1
            return
         end if
         if ((param%in_type /= "ASCII") .and. (param%in_type /= "NETCDF_FLOAT") .and. (param%in_type /= "NETCDF_DOUBLE"))  then
            write(iomsg,*) 'Invalid input file type:',trim(adjustl(param%in_type))
            iostat = -1
            return
         end if
         if (param%istep_out <= 0) then
            write(iomsg,*) 'Invalid ISTEP_OUT. Must be a positive integer'
            iostat = -1
            return
         end if
         if (param%dump_cadence < 0) then
            write(iomsg,*) 'Invalid DUMP_CADENCE. Must be a positive integer or 0.'
            iostat = -1
            return
         end if
         if ((param%istep_out > 0) .and. (param%outfile == "")) then
            write(iomsg,*) 'Invalid outfile'
            iostat = -1
            return
         end if
         param%lrestart = (param%out_stat == "APPEND")
         if (param%outfile /= "") then
            if ((param%out_type /= "NETCDF_FLOAT") .and. (param%out_type /= "NETCDF_DOUBLE")) then
               write(iomsg,*) 'Invalid out_type: ',trim(adjustl(param%out_type))
               iostat = -1
               return
            end if
            if ((param%out_form /= "EL") .and. (param%out_form /= "XV") .and. (param%out_form /= "XVEL")) then
               write(iomsg,*) 'Invalid out_form: ',trim(adjustl(param%out_form))
               iostat = -1
               return
            end if
            if ((param%out_stat /= "NEW") .and. (param%out_stat /= "REPLACE") .and. (param%out_stat /= "APPEND")  &
          .and. (param%out_stat /= "UNKNOWN")) then
               write(iomsg,*) 'Invalid out_stat: ',trim(adjustl(param%out_stat))
               iostat = -1
               return
            end if
         end if
         if (param%qmin > 0.0_DP) then
            if ((param%qmin_coord /= "HELIO") .and. (param%qmin_coord /= "BARY")) then
               write(iomsg,*) 'Invalid qmin_coord: ',trim(adjustl(param%qmin_coord))
               iostat = -1
               return
            end if
            if ((param%qmin_alo <= 0.0_DP) .or. (param%qmin_ahi <= 0.0_DP)) then
               write(iomsg,*) 'Invalid qmin vals'
               iostat = -1
               return
            end if
         end if
         if (param%ltides .and. .not. param%lrotation) then
            write(iomsg,*) 'Tides require rotation to be turned on'
            iostat = -1
            return
         end if

         if ((param%MU2KG < 0.0_DP) .or. (param%TU2S < 0.0_DP) .or. (param%DU2M < 0.0_DP)) then
            write(iomsg,*) 'Invalid unit conversion factor'
            iostat = -1
            return
         end if

         ! Calculate the G for the system units
         param%GU = GC / (param%DU2M**3 / (param%MU2KG * param%TU2S**2))

         associate(integrator => v_list(1))
            if ((integrator == RMVS) .or. (integrator == SYMBA)) then
               if (.not.param%lclose) then
                  write(iomsg,*) 'This integrator requires CHK_CLOSE to be enabled.'
                  iostat = -1
                  return
               end if
            end if
      
            ! Determine if the GR flag is set correctly for this integrator
            select case(integrator)
            case(WHM, RMVS, HELIO, SYMBA)
            case default   
               if (param%lgr) write(iomsg, *) 'GR is not yet implemented for this integrator. This parameter will be ignored.'
               param%lgr = .false.
            end select

            if (param%lgr) then
               ! Calculate the inverse speed of light in the system units
               param%inv_c2 = einsteinC * param%TU2S / param%DU2M
               param%inv_c2 = (param%inv_c2)**(-2)
            end if

         end associate

         select case(trim(adjustl(param%interaction_loops)))
         case("ADAPTIVE")
            param%ladaptive_interactions = .true.
            param%lflatten_interactions = .true.
            call io_log_start(param, INTERACTION_TIMER_LOG_OUT, "Interaction loop timer logfile")
            call io_log_one_message(INTERACTION_TIMER_LOG_OUT, "Diagnostic values: loop style, time count, nplpl, metric")
         case("TRIANGULAR")
            param%ladaptive_interactions = .false.
            param%lflatten_interactions = .false.
         case("FLAT")
            param%ladaptive_interactions = .false.
            param%lflatten_interactions = .true.
         case default
            write(*,*) "Unknown value for parameter INTERACTION_LOOPS: -> ",trim(adjustl(param%interaction_loops))
            write(*,*) "Must be one of the following: TRIANGULAR, FLAT, or ADAPTIVE"
            write(*,*) "Using default value of ADAPTIVE"
            param%interaction_loops = "ADAPTIVE"
            param%ladaptive_interactions = .true.
            param%lflatten_interactions = .true.
            call io_log_start(param, INTERACTION_TIMER_LOG_OUT, "Interaction loop timer logfile")
            call io_log_one_message(INTERACTION_TIMER_LOG_OUT, "Diagnostic values: loop style, time count, nplpl, metric")
         end select

         select case(trim(adjustl(param%encounter_check_plpl)))
         case("ADAPTIVE")
            param%ladaptive_encounters_plpl = .true.
            param%lencounter_sas_plpl = .true.
            call io_log_start(param, ENCOUNTER_PLPL_TIMER_LOG_OUT, "Encounter check loop timer logfile")
            call io_log_one_message(ENCOUNTER_PLPL_TIMER_LOG_OUT, "Diagnostic values: loop style, time count, nplpl, metric")
         case("TRIANGULAR")
            param%ladaptive_encounters_plpl = .false.
            param%lencounter_sas_plpl = .false.
         case("SORTSWEEP")
            param%ladaptive_encounters_plpl = .false.
            param%lencounter_sas_plpl = .true.
         case default
            write(*,*) "Unknown value for parameter ENCOUNTER_CHECK_PLPL: -> ",trim(adjustl(param%encounter_check_plpl))
            write(*,*) "Must be one of the following: TRIANGULAR, SORTSWEEP, or ADAPTIVE"
            write(*,*) "Using default value of ADAPTIVE"
            param%encounter_check_plpl = "ADAPTIVE"
            param%ladaptive_encounters_plpl = .true.
            param%lencounter_sas_plpl = .true.
            call io_log_start(param, ENCOUNTER_PLPL_TIMER_LOG_OUT, "Encounter check loop timer logfile")
            call io_log_one_message(ENCOUNTER_PLPL_TIMER_LOG_OUT, "Diagnostic values: loop style, time count, nplpl, metric")
         end select

         select case(trim(adjustl(param%encounter_check_pltp)))
         case("ADAPTIVE")
            param%ladaptive_encounters_pltp = .true.
            param%lencounter_sas_pltp = .true.
            call io_log_start(param, ENCOUNTER_PLTP_TIMER_LOG_OUT, "Encounter check loop timer logfile")
            call io_log_one_message(ENCOUNTER_PLTP_TIMER_LOG_OUT, "Diagnostic values: loop style, time count, npltp, metric")
         case("TRIANGULAR")
            param%ladaptive_encounters_pltp = .false.
            param%lencounter_sas_pltp = .false.
         case("SORTSWEEP")
            param%ladaptive_encounters_pltp = .false.
            param%lencounter_sas_pltp = .true.
         case default
            write(*,*) "Unknown value for parameter ENCOUNTER_CHECK_PLTP: -> ",trim(adjustl(param%encounter_check_pltp))
            write(*,*) "Must be one of the following: TRIANGULAR, SORTSWEEP, or ADAPTIVE"
            write(*,*) "Using default value of ADAPTIVE"
            param%encounter_check_pltp = "ADAPTIVE"
            param%ladaptive_encounters_pltp = .true.
            param%lencounter_sas_pltp = .true.
            call io_log_start(param, ENCOUNTER_PLTP_TIMER_LOG_OUT, "Encounter check loop timer logfile")
            call io_log_one_message(ENCOUNTER_PLTP_TIMER_LOG_OUT, "Diagnostic values: loop style, time count, npltp, metric")
         end select

         iostat = 0
         
         ! Print the contents of the parameter file to standard output
         call param%writer(unit = param%display_unit, iotype = "none", v_list = [0], iostat = iostat, iomsg = iomsg) 

      end associate

      return 
      667 continue
      write(*,*) "Error reading param file: ", trim(adjustl(iomsg))
   end subroutine io_param_reader


   module subroutine io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_param.f
      implicit none
      ! Arguments
      class(swiftest_parameters),intent(in)     :: self         !! Collection of parameters
      integer, intent(in)                       :: unit       !! File unit number
      character(len=*), intent(in)              :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                               !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
      integer, intent(in)                       :: v_list(:)  !! Not used in this procedure
      integer, intent(out)                      :: iostat     !! IO status code
      character(len=*), intent(inout)           :: iomsg      !! Message to pass if iostat /= 0
      ! Internals
      character(*),parameter :: Ifmt  = '(I0)'         !! Format label for integer values
      character(*),parameter :: Rfmt  = '(ES25.17)'    !! Format label for real values 
      character(*),parameter :: Lfmt  = '(L1)'         !! Format label for logical values 

      associate(param => self)
         call io_param_writer_one("T0", param%t0, unit)
         call io_param_writer_one("TSTOP", param%tstop, unit)
         call io_param_writer_one("DT", param%dt, unit)
         call io_param_writer_one("IN_TYPE", param%in_type, unit)
         if (param%in_type == "ASCII") then
            call io_param_writer_one("CB_IN", param%incbfile, unit)
            call io_param_writer_one("PL_IN", param%inplfile, unit)
            call io_param_writer_one("TP_IN", param%intpfile, unit)
         else 
            call io_param_writer_one("NC_IN", param%in_netcdf, unit)
         end if

         call io_param_writer_one("IN_FORM", param%in_form, unit)
         if (param%dump_cadence > 0) call io_param_writer_one("DUMP_CADENCE",param%dump_cadence, unit)
         if (param%istep_out > 0) then
            call io_param_writer_one("ISTEP_OUT", param%istep_out, unit)
            call io_param_writer_one("BIN_OUT", param%outfile, unit)
            call io_param_writer_one("OUT_TYPE", param%out_type, unit)
            call io_param_writer_one("OUT_FORM", param%out_form, unit)
            call io_param_writer_one("OUT_STAT", "APPEND", unit) 
         end if
         call io_param_writer_one("CHK_RMIN", param%rmin, unit)
         call io_param_writer_one("CHK_RMAX", param%rmax, unit)
         call io_param_writer_one("CHK_EJECT", param%rmaxu, unit)
         call io_param_writer_one("CHK_QMIN", param%qmin, unit)
         if (param%qmin >= 0.0_DP) then
            call io_param_writer_one("CHK_QMIN_COORD", param%qmin_coord, unit)
            call io_param_writer_one("CHK_QMIN_RANGE", [param%qmin_alo, param%qmin_ahi], unit)
         end if
         call io_param_writer_one("MU2KG", param%MU2KG, unit)
         call io_param_writer_one("TU2S", param%TU2S , unit)
         call io_param_writer_one("DU2M", param%DU2M, unit)
         call io_param_writer_one("RHILL_PRESENT", param%lrhill_present, unit)
         call io_param_writer_one("EXTRA_FORCE", param%lextra_force, unit)
         call io_param_writer_one("CHK_CLOSE", param%lclose, unit)
         call io_param_writer_one("ENERGY", param%lenergy, unit)
         call io_param_writer_one("GR", param%lgr, unit)
         call io_param_writer_one("ROTATION", param%lrotation, unit)
         call io_param_writer_one("TIDES", param%ltides, unit)
         call io_param_writer_one("INTERACTION_LOOPS", param%interaction_loops, unit)
         call io_param_writer_one("ENCOUNTER_CHECK_PLPL", param%encounter_check_plpl, unit)
         call io_param_writer_one("ENCOUNTER_CHECK_PLTP", param%encounter_check_pltp, unit)

         if (param%lenergy) then
            call io_param_writer_one("FIRSTENERGY", param%lfirstenergy, unit)
         end if
         call io_param_writer_one("FIRSTKICK",param%lfirstkick, unit)
         call io_param_writer_one("MAXID",param%maxid, unit)
         call io_param_writer_one("MAXID_COLLISION",param%maxid_collision, unit)
   
         iostat = 0
         iomsg = "UDIO not implemented"
      end associate

      667 continue
      return
   end subroutine io_param_writer


   module subroutine io_param_writer_one_char(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for character param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      character(len=*), intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=NAMELEN) :: param_name_fixed_width !! Parameter label converted to fixed-width string
      character(len=STRMAX)  :: iomsg             !! Message to pass if iostat /= 0

      write(param_name_fixed_width, *) param_name
      write(unit, *, err = 667, iomsg = iomsg) adjustl(param_name_fixed_width) // " " // trim(adjustl(param_value)) 

      return
      667 continue
      write(*,*) 'Error writing parameter: ',trim(adjustl(iomsg))
   end subroutine io_param_writer_one_char


   module subroutine io_param_writer_one_DP(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for real(DP) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in)    :: param_name  !! Name of parameter to print
      real(DP),         intent(in)    :: param_value !! Value of parameter to print
      integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Rfmt  = '(ES25.17)' !! Format label for real values 

      write(param_value_string,Rfmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine io_param_writer_one_DP


   module subroutine io_param_writer_one_DParr(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for real(DP) arrays () param_value type
      implicit none
      ! Arguments
      character(len=*),       intent(in) :: param_name  !! Name of parameter to print
      real(DP), dimension(:), intent(in) :: param_value !! Value of parameter to print
      integer(I4B),           intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Rfmt  = '(ES25.17)' !! Format label for real values 
      character(len=25) :: arr_val
      integer(I4B) :: i, narr

      narr = size(param_value)
      do i = 1, narr
         write(arr_val, Rfmt) param_value(i)
         if (i == 1) then
            write(param_value_string, *) arr_val
         else
            param_value_string = trim(adjustl(param_value_string)) // " " // arr_val
         end if
      end do

      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine io_param_writer_one_DParr


   module subroutine io_param_writer_one_I4B(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for integer(I4B) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      integer(I4B),     intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Ifmt  = '(I0)'      !! Format label for integer values

      write(param_value_string,Ifmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine io_param_writer_one_I4B


   module subroutine io_param_writer_one_I8B(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for integer(I8B) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in)    :: param_name  !! Name of parameter to print
      integer(I8B),     intent(in)    :: param_value !! Value of parameter to print
      integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Ifmt  = '(I0)'      !! Format label for integer values

      write(param_value_string,Ifmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine io_param_writer_one_I8B


   module subroutine io_param_writer_one_I4Barr(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for integer(I4B) arrays param_value type
      implicit none
      ! Arguments
      character(len=*),           intent(in) :: param_name  !! Name of parameter to print
      integer(I4B), dimension(:), intent(in) :: param_value !! Value of parameter to print
      integer(I4B),               intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string !! Parameter value converted to a string
      character(*),parameter :: Ifmt  = '(I0)'    !! Format label for integer values
      character(len=25) :: arr_val
      integer(I4B) :: i, narr

      narr = size(param_value)
      do i = 1, narr
         write(arr_val, Ifmt) param_value(i)
         if (i == 1) then
            write(param_value_string, *) trim(adjustl(arr_val))
         else
            param_value_string = trim(adjustl(param_value_string)) // " " // trim(adjustl(arr_val))
         end if
      end do

      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine io_param_writer_one_I4Barr


   module subroutine io_param_writer_one_logical(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for logical param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      logical,          intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Lfmt  = '(L1)'         !! Format label for logical values 

      write(param_value_string,Lfmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine io_param_writer_one_logical


   module subroutine io_param_writer_one_QP(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for real(QP) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      real(QP),         intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Rfmt  = '(ES25.17)' !! Format label for real values 

      write(param_value_string,Rfmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine io_param_writer_one_QP


   module subroutine io_read_in_base(self,param)
      !! author: Carlisle A. Wishard and David A. Minton
      !!
      !! Reads in either a central body, test particle, or massive body object. For the swiftest_body types (non-central body), it allocates array space for them
      implicit none
      class(swiftest_base),       intent(inout) :: self  !! Swiftest base object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 

      if (param%in_type /= "ASCII") return ! This method is not used in NetCDF mode, as reading is done for the whole system, not on individual particle types

      select type(self)
      class is (swiftest_body)
         call io_read_in_body(self, param)
      class is (swiftest_cb)
         call io_read_in_cb(self, param)
      end select

      return
   end subroutine io_read_in_base


   subroutine io_read_in_body(self, param) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in either test particle or massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90 and swiftest_init_tp.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f and swiftest_init_tp.f
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B)                  :: iu = LUN
      integer(I4B)                  :: i, nbody
      character(len=:), allocatable :: infile
      character(STRMAX)             :: errmsg
      ! Internals
      integer(I4B)                                :: ierr  !! Error code: returns 0 if the read is successful

      ! Select the appropriate polymorphic class (test particle or massive body)
      if (param%in_type /= "ASCII") return ! Not for NetCDF

      select type(self)
      class is (swiftest_pl)
         infile = param%inplfile
      class is (swiftest_tp)
         infile = param%intpfile
      end select

      open(unit = iu, file = infile, status = 'old', form = 'FORMATTED', err = 667, iomsg = errmsg)
      read(iu, *, err = 667, iomsg = errmsg) nbody

      call self%setup(nbody, param)
      ierr = 0
      if (nbody > 0) then
         ierr = self%read_frame(iu, param)
         self%status(:) = ACTIVE
         self%lmask(:) = .true.
         do i = 1, nbody
            call self%info(i)%set_value(status="ACTIVE")
         end do
      end if
      close(iu, err = 667, iomsg = errmsg)

      if (ierr == 0) return

      667 continue
      write(*,*) 'Error reading in initial conditions file: ',trim(adjustl(errmsg))
      return
   end subroutine io_read_in_body


   subroutine io_read_in_cb(self, param) 
      !! author: David A. Minton
      !!
      !! Reads in central body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      implicit none
      ! Arguments
      class(swiftest_cb),         intent(inout) :: self
      class(swiftest_parameters), intent(inout) :: param
      ! Internals
      integer(I4B)            :: iu = LUN
      character(len=STRMAX)   :: errmsg
      integer(I4B)            :: ierr
      character(len=NAMELEN)  :: name

      if (param%in_type /= "ASCII") return ! Not for NetCDF

      self%id = 0
      param%maxid = 0
      open(unit = iu, file = param%incbfile, status = 'old', form = 'FORMATTED', err = 667, iomsg = errmsg)
      read(iu, *, err = 667, iomsg = errmsg) name
      call self%info%set_value(name=name)
      read(iu, *, err = 667, iomsg = errmsg) self%Gmass
      self%mass = real(self%Gmass / param%GU, kind=DP)
      read(iu, *, err = 667, iomsg = errmsg) self%radius
      read(iu, *, err = 667, iomsg = errmsg) self%j2rp2
      read(iu, *, err = 667, iomsg = errmsg) self%j4rp4
      if (param%lrotation) then
         read(iu, *, err = 667, iomsg = errmsg) self%Ip(1), self%Ip(2), self%Ip(3)
         read(iu, *, err = 667, iomsg = errmsg) self%rot(1), self%rot(2), self%rot(3)
      end if
      ierr = 0
      close(iu, err = 667, iomsg = errmsg)

      if (ierr == 0) then
   
         if (param%rmin < 0.0) param%rmin = self%radius
         
         select type(cb => self)
         class is (symba_cb)
            cb%GM0 = cb%Gmass
            cb%dGM = 0.0_DP
            cb%R0 = cb%radius
            if (param%lrotation) then
               cb%L0(:) = cb%Ip(3) * cb%mass * cb%radius**2 * cb%rot(:)
               cb%dL(:) = 0.0_DP
            end if
         end select
      end if
      return

      667 continue
      write(*,*) "Error reading central body file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_read_in_cb


   module subroutine io_read_in_system(self, param)
      !! author: David A. Minton and Carlisle A. Wishard
      !!
      !! Reads in the system from input files
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self
      class(swiftest_parameters),   intent(inout) :: param
      ! Internals
      integer(I4B) :: ierr
      class(swiftest_parameters), allocatable :: tmp_param

      if (param%in_type == "ASCII") then
         call self%cb%read_in(param)
         call self%pl%read_in(param)
         call self%tp%read_in(param)
         ! Copy over param file variable inputs
         self%Eorbit_orig = param%Eorbit_orig
         self%GMtot_orig = param%GMtot_orig
         self%Ltot_orig(:) = param%Ltot_orig(:)
         self%Lorbit_orig(:) = param%Lorbit_orig(:)
         self%Lspin_orig(:) = param%Lspin_orig(:)
         self%Lescape(:) = param%Lescape(:)
         self%Ecollisions = param%Ecollisions
         self%Euntracked = param%Euntracked
      else
         allocate(tmp_param, source=param)
         tmp_param%outfile = param%in_netcdf
         tmp_param%out_form = param%in_form
         if (.not. param%lrestart) then
            ! Turn off energy computation so we don't have to feed it into the initial conditions
            tmp_param%lenergy = .false.
         end if
         ierr = self%read_frame(tmp_param%nciu, tmp_param)
         deallocate(tmp_param)
         if (ierr /=0) call util_exit(FAILURE)
      end if

      param%loblatecb = ((self%cb%j2rp2 /= 0.0_DP) .or. (self%cb%j4rp4 /= 0.0_DP))
      if (.not.param%loblatecb) then
         if (allocated(self%pl%aobl)) deallocate(self%pl%aobl)
         if (allocated(self%tp%aobl)) deallocate(self%tp%aobl)
      end if

      return
   end subroutine io_read_in_system


   module function io_read_frame_body(self, iu, param) result(ierr)
      !! author: David A. Minton
      !!
      !! Reads a frame of output of either test particle or massive body data from a binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_read_frame.f90
      !! Adapted from Hal Levison's Swift routine io_read_frame.F
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self   !! Swiftest particle object
      integer(I4B),               intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      ! Result
      integer(I4B)                              :: ierr  !! Error code: returns 0 if the read is successful
      ! Internals
      character(len=STRMAX)   :: errmsg
      character(len=NAMELEN), dimension(self%nbody)  :: name
      integer(I4B) :: i
      real(QP)                      :: val

      if (self%nbody == 0) return

      if ((param%in_form /= "EL") .and. (param%in_form /= "XV")) then
         write(errmsg, *) trim(adjustl(param%in_form)) // " is not a recognized format code for input files."
         goto 667
      end if

      associate(n => self%nbody)

         if (param%in_form == "EL") then
            if (.not.allocated(self%a))     allocate(self%a(n))
            if (.not.allocated(self%e))     allocate(self%e(n))
            if (.not.allocated(self%inc))   allocate(self%inc(n))
            if (.not.allocated(self%capom)) allocate(self%capom(n))
            if (.not.allocated(self%omega)) allocate(self%omega(n))
            if (.not.allocated(self%capm))  allocate(self%capm(n))
         end if

         select case(param%in_type)
         case ("ASCII")
            do i = 1, n
               select type(self)
               class is (swiftest_pl)
                  if (param%lrhill_present) then
                     read(iu, *, err = 667, iomsg = errmsg) name(i), val, self%rhill(i)
                  else
                     read(iu, *, err = 667, iomsg = errmsg) name(i), val
                  end if
                  self%Gmass(i) = real(val, kind=DP)
                  self%mass(i) = real(val / param%GU, kind=DP)
                  if (param%lclose) read(iu, *, err = 667, iomsg = errmsg) self%radius(i)
               class is (swiftest_tp)
                  read(iu, *, err = 667, iomsg = errmsg) name(i)
               end select
               call self%info(i)%set_value(name=name(i))

               select case(param%in_form)
               case ("XV")
                  read(iu, *, err = 667, iomsg = errmsg) self%xh(1, i), self%xh(2, i), self%xh(3, i)
                  read(iu, *, err = 667, iomsg = errmsg) self%vh(1, i), self%vh(2, i), self%vh(3, i)
               case ("EL")
                  read(iu, *, err = 667, iomsg = errmsg) self%a(i), self%e(i), self%inc(i)
                  read(iu, *, err = 667, iomsg = errmsg) self%capom(i), self%omega(i), self%capm(i)
               end select

               select type (self)
               class is (swiftest_pl)
                  if (param%lrotation) then
                     read(iu, *, err = 667, iomsg = errmsg) self%Ip(1, i), self%Ip(2, i), self%Ip(3, i)
                     read(iu, *, err = 667, iomsg = errmsg) self%rot(1, i), self%rot(2, i), self%rot(3, i)
                  end if
                  ! if (param%ltides) then
                  !    read(iu, *, err = 667, iomsg = errmsg) self%k2(i)
                  !    read(iu, *, err = 667, iomsg = errmsg) self%Q(i)
                  ! end if
               end select
               param%maxid = param%maxid + 1
               self%id(i) = param%maxid
            end do
         end select

         if (param%in_form == "EL") then
            self%inc(1:n)   = self%inc(1:n) * DEG2RAD
            self%capom(1:n) = self%capom(1:n) * DEG2RAD
            self%omega(1:n) = self%omega(1:n) * DEG2RAD
            self%capm(1:n)  = self%capm(1:n) * DEG2RAD
         end if
      end associate

      ierr = 0
      return

      667 continue
      select type (self)
      class is (swiftest_pl)
         write(*,*) "Error reading massive body file: " // trim(adjustl(errmsg))
      class is (swiftest_tp)
         write(*,*) "Error reading test particle file: " // trim(adjustl(errmsg))
      class default
         write(*,*) "Error reading body file: " // trim(adjustl(errmsg))
      end select
      call util_exit(FAILURE)
   end function io_read_frame_body


   module subroutine io_read_in_param(self, param_file_name) 
      !! author: David A. Minton
      !!
      !! Read in parameters for the integration
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
      !! Adapted from Martin Duncan's Swift routine io_init_param.f
      implicit none
      ! Arguments
      class(swiftest_parameters),intent(inout) :: self             !! Current run configuration parameters
      character(len=*), intent(in)             :: param_file_name !! Parameter input file name (i.e. param.in)
      ! Internals
      integer(I4B)      :: ierr = 0 !! Input error code
      character(STRMAX) :: errmsg   !! Error message in UDIO procedure

      ! Read in name of parameter file
      write(self%display_unit, *) 'Parameter input file is ', trim(adjustl(param_file_name))
      self%param_file_name = param_file_name

      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    as the newline characters are ignored in the input file when compiled in ifort.

      !read(LUN,'(DT)', iostat= ierr, iomsg = errmsg) self
      call self%reader(LUN, iotype= "none", v_list = [self%integrator], iostat = ierr, iomsg = errmsg)
      if (ierr == 0) return

      667 continue
      write(self%display_unit,*) "Error reading parameter file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_read_in_param


   module subroutine io_set_display_param(self, display_style)
      !! author: David A. Minton
      !!
      !! Sets the display style parameters. If display is "STANDARD" then output goes to stdout. If display is "COMPACT" 
      !! then it is redirected to a log file and a progress-bar is used for stdout
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(inout) :: self            !! Current run configuration parameters
      character(*),               intent(in)    :: display_style   !! Style of the output display 
      ! Internals
      character(STRMAX)             :: errmsg

      select case(display_style)
      case ('STANDARD')
         self%display_unit = OUTPUT_UNIT !! stdout from iso_fortran_env
         self%log_output = .false.
      case ('COMPACT', 'PROGRESS')
         open(unit=SWIFTEST_LOG_OUT, file=SWIFTEST_LOG_FILE, status='replace', err = 667, iomsg = errmsg)
         self%display_unit = SWIFTEST_LOG_OUT 
         self%log_output = .true.
      case default
         write(*,*) display_style, " is an unknown display style"
         call util_exit(USAGE)
      end select

      self%display_style = display_style

      return

      667 continue
      write(*,*) "Error opening swiftest log file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_set_display_param


   module subroutine io_toupper(string)
      !! author: David A. Minton
      !!
      !! Convert string to uppercase
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_toupper.f90
      implicit none
      ! Arguments
      character(*), intent(inout) :: string !! String to make upper case
      ! Internals
      integer(I4B) :: i, length, idx
   
      length = len(string)
      do i = 1, length
         idx = iachar(string(i:i))
         if ((idx >= lowercase_begin) .and. (idx <= lowercase_end)) then
            idx = idx + uppercase_offset
            string(i:i) = achar(idx)
         end if
      end do
   
      return
   end subroutine io_toupper


   module subroutine io_write_frame_system(self, param)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to output binary file
      !! There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      ! Internals
      logical, save                    :: lfirst = .true. !! Flag to determine if this is the first call of this method
      character(len=STRMAX)            :: errmsg
      logical                          :: fileExists

      param%nciu%id_chunk = self%pl%nbody + self%tp%nbody
      param%nciu%time_chunk = max(param%dump_cadence / param%istep_out, 1)
      if (lfirst) then
         inquire(file=param%outfile, exist=fileExists)
         
         select case(param%out_stat)
         case('APPEND')
            if (.not.fileExists) then
               errmsg = param%outfile // " not found! You must specify OUT_STAT = NEW, REPLACE, or UNKNOWN"
               goto 667
            end if
         case('NEW')
            if (fileExists) then
               errmsg = param%outfile // " Alread Exists! You must specify OUT_STAT = APPEND, REPLACE, or UNKNOWN"
               goto 667
            end if
            call param%nciu%initialize(param)
         case('REPLACE', 'UNKNOWN')
            call param%nciu%initialize(param)
         end select

         lfirst = .false.
      end if

      call self%write_frame(param%nciu, param)

      return

      667 continue
      write(*,*) "Error writing system frame: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_write_frame_system

end submodule s_io
