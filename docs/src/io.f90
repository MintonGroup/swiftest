submodule (swiftest_classes) s_io
   use swiftest
contains

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
      real(DP), dimension(NDIM), save :: Ltot_last, Lorbit_last, Lspin_last
      real(DP), save                  :: ke_orbit_last, ke_spin_last, pe_last, Eorbit_last
      real(DP)                        :: ke_orbit_now,  ke_spin_now,  pe_now,  Eorbit_now
      real(DP)                        :: Eorbit_error, Etotal_error, Ecoll_error
      real(DP)                        :: Mtot_now, Merror
      real(DP)                        :: Lmag_now, Lerror
      character(len=*), parameter     :: EGYFMT = '(ES23.16,10(",",ES23.16,:))' ! Format code for all simulation output
      character(len=*), parameter     :: EGYHEADER = '("t,Eorbit,Ecollisions,Lx,Ly,Lz,Mtot")'
      integer(I4B), parameter         :: EGYIU = 72
      character(len=*), parameter     :: EGYTERMFMT = '("  DL/L0 = ", ES12.5 &
                                                         "; DEcollisions/|E0| = ", ES12.5, &
                                                         "; D(Eorbit+Ecollisions)/|E0| = ", ES12.5, &
                                                         "; DM/M0 = ", ES12.5)'

      associate(system => self, pl => self%pl, cb => self%cb, npl => self%pl%nbody, Ecollisions => self%Ecollisions, Lescape => self%Lescape, Mescape => self%Mescape, &
               Euntracked => self%Euntracked, Eorbit_orig => param%Eorbit_orig, Mtot_orig => param%Mtot_orig, &
               Ltot_orig => param%Ltot_orig(:), Lmag_orig => param%Lmag_orig, Lorbit_orig => param%Lorbit_orig(:), Lspin_orig => param%Lspin_orig(:), &
               lfirst => param%lfirstenergy)
         if (lfirst) then
            if (param%out_stat == "OLD") then
               open(unit = EGYIU, file = ENERGY_FILE, form = "formatted", status = "old", action = "write", position = "append")
            else 
               open(unit = EGYIU, file = ENERGY_FILE, form = "formatted", status = "replace", action = "write")
               write(EGYIU,EGYHEADER)
            end if
         end if
         call system%get_energy_and_momentum(param) 
         ke_orbit_now = system%ke_orbit
         ke_spin_now = system%ke_spin
         pe_now = system%pe
         Lorbit_now = system%Lorbit
         Lspin_now = system%Lspin
         Eorbit_now = ke_orbit_now + ke_spin_now + pe_now
         Ltot_now(:) = Lorbit_now(:) + Lspin_now(:) + Lescape(:)
         Mtot_now = cb%mass + sum(pl%mass(1:npl)) + system%Mescape
         if (lfirst) then
            Eorbit_orig = Eorbit_now
            Mtot_orig = Mtot_now
            Lorbit_orig(:) = Lorbit_now(:)
            Lspin_orig(:) = Lspin_now(:)
            Ltot_orig(:) = Ltot_now(:)
            Lmag_orig = norm2(Ltot_orig(:))
            lfirst = .false.
         end if

         write(EGYIU,EGYFMT) param%t, Eorbit_now, Ecollisions, Ltot_now, Mtot_now
         flush(EGYIU)
         if (.not.lfirst .and. lterminal) then 
            Lmag_now = norm2(Ltot_now)
            Lerror = norm2(Ltot_now - Ltot_orig) / Lmag_orig
            Eorbit_error = (Eorbit_now - Eorbit_orig) / abs(Eorbit_orig)
            Ecoll_error = Ecollisions / abs(Eorbit_orig)
            Etotal_error = (Eorbit_now - Ecollisions - Eorbit_orig - Euntracked) / abs(Eorbit_orig)
            Merror = (Mtot_now - Mtot_orig) / Mtot_orig
            write(*, egytermfmt) Lerror, Ecoll_error, Etotal_error, Merror
            if (Ecoll_error > 0.0_DP) then
               write(*,*) 'Something has gone wrong! Collisional energy should not be positive!'
               write(*,*) 'dke_orbit: ',(ke_orbit_now - ke_orbit_last) / abs(Eorbit_orig)
               write(*,*) 'dke_spin : ',(ke_spin_now - ke_spin_last) / abs(Eorbit_orig)
               write(*,*) 'dpe      : ',(pe_now - pe_last) / abs(Eorbit_orig)
               write(*,*)
            end if
         end if
         ke_orbit_last = ke_orbit_now
         ke_spin_last = ke_spin_now
         pe_last = pe_now
         Eorbit_last = Eorbit_now
         Lorbit_last(:) = Lorbit_now(:)
         Lspin_last(:) = Lspin_now(:)
         Ltot_last(:) = Ltot_now(:)
      end associate
      return

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
      integer(I4B), parameter  :: LUN = 7       !! Unit number of output file
      integer(I4B)             :: ierr          !! Error code
      character(STRMAX)        :: error_message !! Error message in UDIO procedure

      open(unit = LUN, file = param_file_name, status='replace', form = 'FORMATTED', iostat =ierr)
      if (ierr /=0) then
         write(*,*) 'Swiftest error.'
         write(*,*) '   Could not open dump file: ',trim(adjustl(param_file_name))
         call util_exit(FAILURE)
      end if
      
      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    due to compiler incompatabilities
      !write(LUN,'(DT)') param
      call self%writer(LUN, iotype = "none", v_list = [0], iostat = ierr, iomsg = error_message)
      if (ierr /= 0) then
         write(*,*) trim(adjustl(error_message))
         call util_exit(FAILURE)
      end if
      close(LUN)

      return
   end subroutine io_dump_param


   module subroutine io_dump_swiftest(self, param, msg) 
      !! author: David A. Minton
      !!
      !! Dump massive body data to files
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_dump_pl.f90 and io_dump_tp.f90
      !! Adapted from Hal Levison's Swift routine io_dump_pl.f and io_dump_tp.f
      implicit none
      ! Arguments
      class(swiftest_base),       intent(inout) :: self   !! Swiftest base object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      character(*), optional,     intent(in)    :: msg  !! Message to display with dump operation
      ! Internals
      integer(I4B)                   :: ierr    !! Error code
      integer(I4B),parameter         :: LUN = 7 !! Unit number for dump file
      integer(I4B)                   :: iu = LUN
      character(len=:), allocatable  :: dump_file_name

      select type(self)
      class is(swiftest_cb)
         dump_file_name = trim(adjustl(param%incbfile)) 
      class is (swiftest_pl)
         dump_file_name = trim(adjustl(param%inplfile)) 
      class is (swiftest_tp)
         dump_file_name = trim(adjustl(param%intpfile)) 
      end select
      open(unit = iu, file = dump_file_name, form = "UNFORMATTED", status = 'replace', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Unable to open binary dump file " // dump_file_name
         call util_exit(FAILURE)
      end if
      call self%write_frame(iu, param)
      close(LUN)

      return
   end subroutine io_dump_swiftest


   module subroutine io_dump_system(self, param, msg)
      !! author: David A. Minton
      !!
      !! Dumps the state of the system to files in case the simulation is interrupted.
      !! As a safety mechanism, there are two dump files that are written in alternating order
      !! so that if a dump file gets corrupted during writing, the user can restart from the older one.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      character(*), optional,       intent(in)    :: msg   !! Message to display with dump operation
      ! Internals
      class(swiftest_parameters), allocatable :: dump_param !! Local parameters variable used to parameters change input file names 
                                                            !! to dump file-specific values without changing the user-defined values
      integer(I4B), save            :: idx = 1              !! Index of current dump file. Output flips between 2 files for extra security
                                                            !! in case the program halts during writing
      character(len=:), allocatable :: param_file_name
      real(DP) :: tfrac
     
      allocate(dump_param, source=param)
      param_file_name    = trim(adjustl(DUMP_PARAM_FILE(idx)))
      dump_param%incbfile = trim(adjustl(DUMP_CB_FILE(idx))) 
      dump_param%inplfile = trim(adjustl(DUMP_PL_FILE(idx))) 
      dump_param%intpfile = trim(adjustl(DUMP_TP_FILE(idx)))
      dump_param%out_form = XV
      dump_param%out_stat = 'APPEND'
      dump_param%T0 = param%t
      call dump_param%dump(param_file_name)

      call self%cb%dump(dump_param)
      if (self%pl%nbody > 0) call self%pl%dump(dump_param)
      if (self%tp%nbody > 0) call self%tp%dump(dump_param)

      idx = idx + 1
      if (idx > NDUMPFILES) idx = 1

      ! Print the status message (format code passed in from main driver)
      tfrac = (param%t - param%t0) / (param%tstop - param%t0)
      write(*,msg) param%t, tfrac, self%pl%nbody, self%tp%nbody
      if (param%lenergy) call self%conservation_report(param, lterminal=.true.)

      return
   end subroutine io_dump_system


   module function io_get_args(integrator, param_file_name) result(ierr)
      !! author: David A. Minton
      !!
      !! Reads in the name of the parameter file from command line arguments. 
      implicit none
      ! Arguments
      integer(I4B)                  :: integrator      !! Symbolic code of the requested integrator  
      character(len=:), allocatable :: param_file_name !! Name of the input parameters file
      ! Result
      integer(I4B)                  :: ierr             !! I/O error code
      ! Internals
      character(len=STRMAX) :: arg1, arg2
      integer :: narg,ierr_arg1, ierr_arg2
      character(len=*),parameter    :: linefmt = '(A)'

      ierr = -1 ! Default is to fail
      narg = command_argument_count() !
      if (narg == 2) then
         call get_command_argument(1, arg1, status = ierr_arg1)
         call get_command_argument(2, arg2, status = ierr_arg2)
         if ((ierr_arg1 == 0) .and. (ierr_arg2 == 0)) then
            ierr = 0
            call io_toupper(arg1)
            select case(arg1)
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
               write(*,*) trim(adjustl(arg1)) // ' is not a valid integrator.'
               ierr = -1
            end select
            param_file_name = trim(adjustl(arg2))
         end if
      else 
         call get_command_argument(1, arg1, status = ierr_arg1)
         if (ierr_arg1 == 0) then
            if (arg1 == '-v' .or. arg1 == '--version') then
               call util_version() 
            else if (arg1 == '-h' .or. arg1 == '--help') then
               call util_exit(HELP)
            end if
         end if
      end if
      if (ierr /= 0) call util_exit(USAGE) 

      return
   end function io_get_args


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
      integer, intent(in)                       :: v_list(:)  !! The first element passes the integrator code to the reader
      integer, intent(out)                      :: iostat     !! IO status code
      character(len=*), intent(inout)           :: iomsg      !! Message to pass if iostat /= 0
      ! Internals
      logical                        :: t0_set = .false.                  !! Is the initial time set in the input file?
      logical                        :: tstop_set = .false.               !! Is the final time set in the input file?
      logical                        :: dt_set = .false.                  !! Is the step size set in the input file?
      integer(I4B)                   :: ilength, ifirst, ilast, i         !! Variables used to parse input file
      character(STRMAX)              :: line                              !! Line of the input file
      character (len=:), allocatable :: line_trim,param_name, param_value !! Strings used to parse the param file
      character(*),parameter         :: linefmt = '(A)'                   !! Format code for simple text string

      ! Parse the file line by line, extracting tokens then matching them up with known parameters if possible
      associate(param => self) 
         do
            read(unit = unit, fmt = linefmt, iostat = iostat, end = 1) line
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
                  read(param_value, *) param%t0
                  t0_set = .true.
               case ("TSTOP")
                  read(param_value, *) param%tstop
                  tstop_set = .true.
               case ("DT")
                  read(param_value, *) param%dt
                  dt_set = .true.
               case ("CB_IN")
                  param%incbfile = param_value
               case ("PL_IN")
                  param%inplfile = param_value
               case ("TP_IN")
                  param%intpfile = param_value
               case ("IN_TYPE")
                  call io_toupper(param_value)
                  param%in_type = param_value
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
               case ("ISTEP_DUMP")
                  read(param_value, *) param%istep_dump
               case ("CHK_CLOSE")
                  call io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lclose = .true.
               case ("CHK_RMIN")
                  read(param_value, *) param%rmin
               case ("CHK_RMAX")
                  read(param_value, *) param%rmax
               case ("CHK_EJECT")
                  read(param_value, *) param%rmaxu
               case ("CHK_QMIN")
                  read(param_value, *) param%qmin
               case ("CHK_QMIN_COORD")
                  call io_toupper(param_value)
                  param%qmin_coord = param_value
               case ("CHK_QMIN_RANGE")
                  read(param_value, *) param%qmin_alo
                  ifirst = ilast + 1
                  param_value = io_get_token(line, ifirst, ilast, iostat)
                  read(param_value, *) param%qmin_ahi
               case ("ENC_OUT")
                  param%enc_out = param_value
               case ("DISCARD_OUT")
                  param%discard_out = param_value
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
                  read(param_value, *) param%MU2KG
               case ("TU2S")
                  read(param_value, *) param%TU2S
               case ("DU2M")
                  read(param_value, *) param%DU2M
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
               case ("FIRSTKICK")
                  call io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstkick = .false. 
               case ("FIRSTENERGY")
                  call io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstenergy = .false. 
               case("EORBIT_ORIG")
                  read(param_value, *) param%Eorbit_orig 
               case("MTOT_ORIG")
                  read(param_value, *) param%Mtot_orig 
               case("LTOT_ORIG")
                  read(param_value, *) param%Ltot_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *) param%Ltot_orig(i)
                  end do
                     param%Lmag_orig = norm2(param%Ltot_orig(:))
               case("LORBIT_ORIG")
                  read(param_value, *) param%Lorbit_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *) param%Lorbit_orig(i)
                  end do
               case("LSPIN_ORIG")
                  read(param_value, *) param%Lspin_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *) param%Lspin_orig(i)
                  end do
               case("LESCAPE")
                  read(param_value, *) param%Lescape(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *) param%Lescape(i)
                  end do
               case("MESCAPE")
                  read(param_value, *) param%Mescape 
               case("ECOLLISIONS")
                  read(param_value, *) param%Ecollisions
               case("EUNTRACKED")
                  read(param_value, *) param%Euntracked
               case ("NPLMAX", "NTPMAX", "GMTINY", "PARTICLE_FILE", "FRAGMENTATION", "SEED", "YARKOVSKY", "YORP") ! Ignore SyMBA-specific, not-yet-implemented, or obsolete input parameters
               case default
                  write(iomsg,*) "Unknown parameter -> ",param_name
                  iostat = -1
                  return
               end select
            end if
         end do
         1 continue
         iostat = 0

         !! Do basic sanity checks on the input values
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
         if ((param%in_type /= REAL8_TYPE) .and. (param%in_type /= "ASCII")) then
            write(iomsg,*) 'Invalid input file type:',trim(adjustl(param%in_type))
            iostat = -1
            return
         end if
         if ((param%istep_out <= 0) .and. (param%istep_dump <= 0)) then
            write(iomsg,*) 'Invalid istep'
            iostat = -1
            return
         end if
         if ((param%istep_out > 0) .and. (param%outfile == "")) then
            write(iomsg,*) 'Invalid outfile'
            iostat = -1
            return
         end if
         if (param%outfile /= "") then
            if ((param%out_type /= REAL4_TYPE) .and. (param%out_type /= REAL8_TYPE) .and. &
                  (param%out_type /= SWIFTER_REAL4_TYPE)  .and. (param%out_type /= SWIFTER_REAL8_TYPE)) then
               write(iomsg,*) 'Invalid out_type: ',trim(adjustl(param%out_type))
               iostat = -1
               return
            end if
            if ((param%out_form /= "EL") .and. (param%out_form /= "XV")) then
               write(iomsg,*) 'Invalid out_form: ',trim(adjustl(param%out_form))
               iostat = -1
               return
            end if
            if ((param%out_stat /= "NEW") .and. (param%out_stat /= "REPLACE") .and. (param%out_stat /= "APPEND")  .and. (param%out_stat /= "UNKNOWN")) then
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

         write(*,*) "T0             = ",param%t0
         write(*,*) "TSTOP          = ",param%tstop
         write(*,*) "DT             = ",param%dt
         write(*,*) "CB_IN          = ",trim(adjustl(param%incbfile))
         write(*,*) "PL_IN          = ",trim(adjustl(param%inplfile))
         write(*,*) "TP_IN          = ",trim(adjustl(param%intpfile))
         write(*,*) "IN_TYPE        = ",trim(adjustl(param%in_type))
         write(*,*) "ISTEP_OUT      = ",param%istep_out
         write(*,*) "BIN_OUT        = ",trim(adjustl(param%outfile))
         write(*,*) "OUT_TYPE       = ",trim(adjustl(param%out_type))
         write(*,*) "OUT_FORM       = ",trim(adjustl(param%out_form))
         write(*,*) "OUT_STAT       = ",trim(adjustl(param%out_stat))
         write(*,*) "ISTEP_DUMP     = ",param%istep_dump
         write(*,*) "CHK_CLOSE      = ",param%lclose
         write(*,*) "CHK_RMIN       = ",param%rmin
         write(*,*) "CHK_RMAX       = ",param%rmax
         write(*,*) "CHK_EJECT      = ",param%rmaxu
         write(*,*) "CHK_QMIN       = ",param%qmin
         write(*,*) "CHK_QMIN_COORD = ",trim(adjustl(param%qmin_coord))
         write(*,*) "CHK_QMIN_RANGE = ",param%qmin_alo, param%qmin_ahi
         write(*,*) "EXTRA_FORCE    = ",param%lextra_force
         write(*,*) "RHILL_PRESENT  = ",param%lrhill_present
         write(*,*) "ROTATION      = ", param%lrotation
         write(*,*) "TIDES         = ", param%ltides
         write(*,*) "ENERGY         = ",param%lenergy
         write(*,*) "MU2KG          = ",param%MU2KG       
         write(*,*) "TU2S           = ",param%TU2S        
         write(*,*) "DU2M           = ",param%DU2M        
         if (trim(adjustl(param%enc_out)) /= "") then
            write(*,*) "ENC_OUT        = ",trim(adjustl(param%enc_out))
         else
            write(*,*) "! ENC_OUT not set: Encounters will not be recorded to file"
         end if
         if (trim(adjustl(param%discard_out)) /= "") then
            write(*,*) "DISCARD_OUT    = ",trim(adjustl(param%discard_out))
            write(*,*) "BIG_DISCARD    = ",param%lbig_discard
         else
            write(*,*) "! DISCARD_OUT not set: Discards will not be recorded to file"
            write(*,*) "! BIG_DISCARD    = ",param%lbig_discard
         end if

         if ((param%MU2KG < 0.0_DP) .or. (param%TU2S < 0.0_DP) .or. (param%DU2M < 0.0_DP)) then
            write(iomsg,*) 'Invalid unit conversion factor'
            iostat = -1
            return
         end if

         ! Calculate the G for the system units
         param%GU = GC / (param%DU2M**3 / (param%MU2KG * param%TU2S**2))

         ! Calculate the inverse speed of light in the system units
         param%inv_c2 = einsteinC * param%TU2S / param%DU2M
         param%inv_c2 = (param%inv_c2)**(-2)

         associate(integrator => v_list(1))
            if (integrator == RMVS) then
               if (.not.param%lclose) then
                  write(iomsg,*) 'This integrator requires CHK_CLOSE to be enabled.'
                  iostat = -1
                  return
               end if
            end if
      
            ! Determine if the GR flag is set correctly for this integrator
            select case(integrator)
            case(WHM, RMVS, HELIO, SYMBA)
               write(*,*) "GR             = ", param%lgr
            case default   
               if (param%lgr) write(iomsg, *) 'GR is not yet implemented for this integrator. This parameter will be ignored.'
               param%lgr = .false.
            end select
         end associate

         iostat = 0
      end associate

      return 
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
      character(*),parameter :: Rarrfmt  = '(3(ES25.17,1X))'    !! Format label for real values 
      character(*),parameter :: Lfmt  = '(L1)'         !! Format label for logical values 
      character(len=*), parameter :: Afmt = '(A25,1X,64(:,A25,1X))'
      character(256)          :: param_name, param_value
      type character_array
         character(25) :: value
      end type character_array
      type(character_array), dimension(:), allocatable :: param_array
      integer(I4B) :: i

      associate(param => self)
         write(param_name, Afmt) "T0"; write(param_value,Rfmt) param%t0; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TSTOP"; write(param_value, Rfmt) param%tstop; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "DT"; write(param_value, Rfmt) param%dt; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "PL_IN"; write(param_value, Afmt) trim(adjustl(param%inplfile)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TP_in"; write(param_value, Afmt) trim(adjustl(param%intpfile)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "IN_TYPE"; write(param_value, Afmt) trim(adjustl(param%in_type)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         if (param%istep_out > 0) then
            write(param_name, Afmt) "ISTEP_OUT"; write(param_value, Ifmt) param%istep_out; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "BIN_OUT"; write(param_value, Afmt) trim(adjustl(param%outfile)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "OUT_TYPE"; write(param_value, Afmt) trim(adjustl(param%out_type)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "OUT_FORM"; write(param_value, Afmt) trim(adjustl(param%out_form)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "OUT_STAT"; write(param_value, Afmt) "APPEND"; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         end if
         write(param_name, Afmt) "ENC_OUT"; write(param_value, Afmt) trim(adjustl(param%enc_out)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         if (param%istep_dump > 0) then
            write(param_name, Afmt) "ISTEP_DUMP"; write(param_value, Ifmt) param%istep_dump; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         end if
         write(param_name, Afmt) "CHK_RMIN"; write(param_value, Rfmt) param%rmin; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_RMAX"; write(param_value, Rfmt) param%rmax; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_EJECT"; write(param_value, Rfmt) param%rmaxu; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_QMIN"; write(param_value, Rfmt) param%qmin; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         if (param%qmin >= 0.0_DP) then
            write(param_name, Afmt) "CHK_QMIN_COORD"; write(param_value, Afmt) trim(adjustl(param%qmin_coord)); write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            allocate(param_array(2))
            write(param_array(1)%value, Rfmt) param%qmin_alo
            write(param_array(2)%value, Rfmt) param%qmin_ahi
            write(param_name, Afmt) "CHK_QMIN_RANGE"; write(unit, Afmt) adjustl(param_name), adjustl(param_array(1)%value), adjustl(param_array(2)%value)
         end if
         write(param_name, Afmt) "MU2KG"; write(param_value, Rfmt) param%MU2KG; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TU2S"; write(param_value, Rfmt) param%TU2S ; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "DU2M"; write(param_value, Rfmt) param%DU2M; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "RHILL_PRESENT"; write(param_value, Lfmt) param%lrhill_present; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "EXTRA_FORCE"; write(param_value, Lfmt) param%lextra_force; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "BIG_DISCARD"; write(param_value, Lfmt) param%lbig_discard; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_CLOSE"; write(param_value, Lfmt) param%lclose; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "ENERGY"; write(param_value, Lfmt)  param%lenergy; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "GR"; write(param_value, Lfmt)  param%lgr; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "ROTATION"; write(param_value, Lfmt)  param%lrotation; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TIDES"; write(param_value, Lfmt)  param%ltides; write(unit, Afmt) adjustl(param_name), adjustl(param_value)

         if (param%lenergy) then
            write(param_name, Afmt) "FIRSTENERGY"; write(param_value, Lfmt) param%lfirstenergy; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "EORBIT_ORIG"; write(param_value, Rfmt) param%Eorbit_orig; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "MTOT_ORIG"; write(param_value, Rfmt) param%Mtot_orig; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(unit, '("LTOT_ORIG  ",3(1X,ES25.17))') param%Ltot_orig(:)
            write(unit, '("LORBIT_ORIG",3(1X,ES25.17))') param%Lorbit_orig(:)
            write(unit, '("LSPIN_ORIG ",3(1X,ES25.17))') param%Lspin_orig(:)
            write(unit, '("LESCAPE ",3(1X,ES25.17))') param%Lescape(:)
   
            write(param_name, Afmt) "MESCAPE"; write(param_value, Rfmt) param%Mescape; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "ECOLLISIONS"; write(param_value, Rfmt) param%Ecollisions; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "EUNTRACKED"; write(param_value, Rfmt) param%Euntracked; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
         end if
         write(param_name, Afmt) "FIRSTKICK"; write(param_value, Lfmt) param%lfirstkick; write(unit, Afmt) adjustl(param_name), adjustl(param_value)
   


         iostat = 0
         iomsg = "UDIO not implemented"
      end associate

      return
   end subroutine io_param_writer


   module subroutine io_read_body_in(self, param) 
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
      integer(I4B), parameter       :: LUN = 7              !! Unit number of input file
      integer(I4B)                  :: iu = LUN
      integer(I4B)                  :: i, ierr, nbody
      logical                       :: is_ascii, is_pl
      character(len=:), allocatable :: infile
      real(DP)                      :: t
      real(QP)                      :: val

      ! Select the appropriate polymorphic class (test particle or massive body)
      select type(self)
      class is (swiftest_pl)
         infile = param%inplfile
         is_pl = .true.
      class is (swiftest_tp)
         infile = param%intpfile
         is_pl = .false.
      end select

      ierr = 0
      is_ascii = (param%in_type == 'ASCII') 
      select case(param%in_type)
      case(ASCII_TYPE)
         open(unit = iu, file = infile, status = 'old', form = 'FORMATTED', iostat = ierr)
         read(iu, *, iostat = ierr) nbody
         call self%setup(nbody, param)
         if (nbody > 0) then
            do i = 1, nbody
               select type(self)
               class is (swiftest_pl)
                  if (param%lrhill_present) then
                     read(iu, *, iostat=ierr, err=100) self%id(i), val, self%rhill(i)
                  else
                     read(iu, *, iostat=ierr, err=100) self%id(i), val
                  end if
                  self%Gmass(i) = real(val, kind=DP)
                  self%mass(i) = real(val / param%GU, kind=DP)
                  if (param%lclose) read(iu, *, iostat=ierr, err=100) self%radius(i)
                  if (param%lrotation) then
                     read(iu, iostat=ierr, err=100) self%Ip(:, i)
                     read(iu, iostat=ierr, err=100) self%rot(:, i)
                  end if
                  if (param%ltides) then
                     read(iu, iostat=ierr, err=100) self%k2(i)
                     read(iu, iostat=ierr, err=100) self%Q(i)
                  end if
               class is (swiftest_tp)
                  read(iu, *, iostat=ierr, err=100) self%id(i)
               end select
               read(iu, *, iostat=ierr, err=100) self%xh(1, i), self%xh(2, i), self%xh(3, i)
               read(iu, *, iostat=ierr, err=100) self%vh(1, i), self%vh(2, i), self%vh(3, i)
               self%status(i) = ACTIVE
               self%lmask(i) = .true.
            end do
         end if
      case (REAL4_TYPE, REAL8_TYPE)  !, SWIFTER_REAL4_TYPE, SWIFTER_REAL8_TYPE)
         open(unit=iu, file=infile, status='old', form='UNFORMATTED', iostat=ierr)
         read(iu, iostat=ierr, err=100) nbody
         call self%setup(nbody, param)
         if (nbody > 0) then
            call self%read_frame(iu, param, XV, ierr)
            self%status(:) = ACTIVE
            self%lmask(:) = .true.
         end if
      case default
         write(*,*) trim(adjustl(param%in_type)) // ' is an unrecognized file type'
         ierr = -1
      end select
      close(iu)

      100 if (ierr /= 0 ) then
         write(*,*) 'Error reading in initial conditions from ',trim(adjustl(infile))
         call util_exit(FAILURE)
      end if

      return
   end subroutine io_read_body_in


   module subroutine io_read_cb_in(self, param) 
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
      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: iu = LUN
      integer(I4B)            :: ierr
      logical                 :: is_ascii 
      real(DP)                :: t
      real(QP)                :: val

      ierr = 0
      is_ascii = (param%in_type == 'ASCII') 
      if (is_ascii) then
         open(unit = iu, file = param%incbfile, status = 'old', form = 'FORMATTED', iostat = ierr)
         !read(iu, *, iostat = ierr) self%id
         read(iu, *, iostat = ierr) val 
         self%Gmass = real(val, kind=DP)
         self%mass = real(val / param%GU, kind=DP)
         read(iu, *, iostat = ierr) self%radius
         read(iu, *, iostat = ierr) self%j2rp2
         read(iu, *, iostat = ierr) self%j4rp4
      else
         open(unit = iu, file = param%incbfile, status = 'old', form = 'UNFORMATTED', iostat = ierr)
         call self%read_frame(iu, param, XV, ierr)
      end if
      close(iu)
      if (ierr /=  0) then
         write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(param%incbfile))
         call util_exit(FAILURE)
      end if
      if (self%j2rp2 /= 0.0_DP) param%loblatecb = .true.

      return
   end subroutine io_read_cb_in



   function io_read_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
      xh1, xh2, vh1, vh2, enc_out, out_type) result(ierr)
      !! author: David A. Minton
      !!
      !! Read close encounter data from input binary files
      !!     Other than time t, there is no direct file input from this function
      !!     Function returns read error status (0 = OK, nonzero = ERROR)
      !! Adapted from David E. Kaufmann's Swifter routine: io_read_encounter.f90
      implicit none
      ! Arguments
      integer(I4B),           intent(out) :: name1, name2
      real(DP),               intent(out) :: t, mass1, mass2, radius1, radius2
      real(DP), dimension(:), intent(out) :: xh1, xh2, vh1, vh2
      character(*),           intent(in)  :: enc_out, out_type
      ! Result
      integer(I4B)         :: ierr
      ! Internals
      logical , save    :: lfirst = .true.
      integer(I4B), parameter :: lun = 30
      integer(I4B), save    :: iu = lun

      if (lfirst) then
         open(unit = iu, file = enc_out, status = 'OLD', form = 'UNFORMATTED', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest Error:"
            write(*, *) "   unable to open binary encounter file"
            call util_exit(FAILURE)
         end if
         lfirst = .false.
      end if
      read(iu, iostat = ierr) t
      if (ierr /= 0) then
         close(unit = iu, iostat = ierr)
         return
      end if

      read(iu, iostat = ierr) name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), mass1, radius1
      if (ierr /= 0) then
         close(unit = iu, iostat = ierr)
         return
      end if
      read(iu, iostat = ierr) name2, xh2(2), xh2(2), xh2(3), vh2(2), vh2(2), vh2(3), mass2, radius2
      if (ierr /= 0) then
         close(unit = iu, iostat = ierr)
         return
      end if

      return
   end function io_read_encounter


   module subroutine io_read_frame_body(self, iu, param, form, ierr)
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
      character(*),               intent(in)    :: form   !! Input format code ("XV" or "EL")
      integer(I4B),               intent(out)   :: ierr   !! Error code

      associate(n => self%nbody)
         read(iu, iostat=ierr, err=100) self%id(:)
         !read(iu, iostat=ierr, err=100) self%name(1:n)
         select case (form)
         case (EL) 
            if (.not.allocated(self%a))     allocate(self%a(n))
            if (.not.allocated(self%e))     allocate(self%e(n))
            if (.not.allocated(self%inc))   allocate(self%inc(n))
            if (.not.allocated(self%capom)) allocate(self%capom(n))
            if (.not.allocated(self%omega)) allocate(self%omega(n))
            if (.not.allocated(self%capm))  allocate(self%capm(n))
            read(iu, iostat=ierr, err=100) self%a(:)
            read(iu, iostat=ierr, err=100) self%e(:)
            read(iu, iostat=ierr, err=100) self%inc(:)
            read(iu, iostat=ierr, err=100) self%capom(:)
            read(iu, iostat=ierr, err=100) self%omega(:)
            read(iu, iostat=ierr, err=100) self%capm(:)
         case (XV)
            read(iu, iostat=ierr, err=100) self%xh(1, :)
            read(iu, iostat=ierr, err=100) self%xh(2, :)
            read(iu, iostat=ierr, err=100) self%xh(3, :)
            read(iu, iostat=ierr, err=100) self%vh(1, :)
            read(iu, iostat=ierr, err=100) self%vh(2, :)
            read(iu, iostat=ierr, err=100) self%vh(3, :)
         end select
         select type(pl => self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            read(iu, iostat=ierr, err=100) pl%Gmass(:)
            pl%mass(:) = pl%Gmass(:) / param%GU 
            if (param%lrhill_present) read(iu, iostat=ierr, err=100) pl%rhill(:)
            read(iu, iostat=ierr, err=100) pl%radius(:)
            if (param%lrotation) then
               read(iu, iostat=ierr, err=100) pl%rot(1, :)
               read(iu, iostat=ierr, err=100) pl%rot(2, :)
               read(iu, iostat=ierr, err=100) pl%rot(3, :)
               read(iu, iostat=ierr, err=100) pl%Ip(1, :)
               read(iu, iostat=ierr, err=100) pl%Ip(2, :)
               read(iu, iostat=ierr, err=100) pl%Ip(3, :)
            end if
            if (param%ltides) then
               read(iu, iostat=ierr, err=100) pl%k2(1:n)
               read(iu, iostat=ierr, err=100) pl%Q(1:n)
            end if
         end select
      end associate

      100 if (ierr /=0) then
         write(*,*) 'Error reading Swiftest body data'
         call util_exit(FAILURE)
      end if

      return
   end subroutine io_read_frame_body


   module subroutine io_read_frame_cb(self, iu, param, form, ierr)
      !! author: David A. Minton
      !!
      !! Reads a frame of output of central body data to the binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_read_frame.f90
      !! Adapted from Hal Levison's Swift routine io_read_frame.F
      implicit none
      ! Arguments
      class(swiftest_cb),         intent(inout) :: self     !! Swiftest central body object
      integer(I4B),               intent(inout) :: iu       !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(inout) :: param   !! Current run configuration parameters 
      character(*),               intent(in)    :: form     !! Input format code ("XV" or "EL")
      integer(I4B),               intent(out)   :: ierr     !! Error cod

      read(iu, iostat=ierr, err=100) self%id
      !read(iu, iostat=ierr, err=100) self%name
      read(iu, iostat=ierr, err=100) self%Gmass
      self%mass = self%Gmass / param%GU
      read(iu, iostat=ierr, err=100) self%radius
      read(iu, iostat=ierr, err=100) self%j2rp2 
      read(iu, iostat=ierr, err=100) self%j4rp4 
      if (param%lrotation) then
         read(iu, iostat=ierr, err=100) self%Ip(:)
         read(iu, iostat=ierr, err=100) self%rot(:)
      end if
      if (param%ltides) then
         read(iu, iostat=ierr, err=100) self%k2
         read(iu, iostat=ierr, err=100) self%Q
      end if
      100 if (ierr /=0) then
         write(*,*) 'Error reading central body data'
         call util_exit(FAILURE)
      end if

      return
   end subroutine io_read_frame_cb


   module subroutine io_read_frame_system(self, iu, param, form, ierr)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read a frame (header plus records for each massive body and active test particle) from a output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      integer(I4B),                 intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      character(*),                 intent(in)    :: form   !! Input format code ("XV" or "EL")
      integer(I4B),                 intent(out)   :: ierr   !! Error code
      ! Internals
      logical, save             :: lfirst = .true.

      iu = BINUNIT
      if (lfirst) then
         open(unit = iu, file = param%outfile, status = 'OLD', form = 'UNFORMATTED', iostat = ierr)
         lfirst = .false.
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Binary output file already exists or cannot be accessed"
            return
         end if
      end if
      ierr = io_read_hdr(iu, param%t, self%pl%nbody, self%tp%nbody, param%out_form, param%out_type)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Binary output file already exists or cannot be accessed"
         return
      end if
      call self%cb%read_frame(iu, param, form, ierr)
      if (ierr /= 0) return
      call self%pl%read_frame(iu, param, form, ierr)
      if (ierr /= 0) return
      call self%tp%read_frame(iu, param, form, ierr)

      return
   end subroutine io_read_frame_system


   function io_read_hdr(iu, t, npl, ntp, out_form, out_type) result(ierr)
      !! author: David A. Minton
      !!
      !! Read frame header from input binary files
      !!     Function returns read error status (0 = OK, nonzero = ERROR)
      !! Adapted from David E. Kaufmann's Swifter routine: io_read_hdr.f90
      !! Adapted from Hal Levison's Swift routine io_read_hdr.f
      implicit none
      ! Arguments
      integer(I4B), intent(in)   :: iu
      integer(I4B), intent(out)  :: npl, ntp
      character(*), intent(out)  ::  out_form
      real(DP),     intent(out)  :: t
      character(*), intent(in)   :: out_type
      ! Result
      integer(I4B)      :: ierr
      ! Internals
      real(SP)             :: ttmp

      select case (out_type)
      case (REAL4_TYPE, SWIFTER_REAL4_TYPE)
         read(iu, iostat = ierr) ttmp, npl, ntp, out_form
         if (ierr /= 0) return
         t = ttmp
      case (REAL8_TYPE, SWIFTER_REAL8_TYPE)
         read(iu, iostat = ierr) t
         read(iu, iostat = ierr) npl
         read(iu, iostat = ierr) ntp
         read(iu, iostat = ierr) out_form
      case default
         write(*,*) trim(adjustl(out_type)) // ' is an unrecognized file type'
         ierr = -1
      end select

      return
   end function io_read_hdr

   module subroutine io_read_param_in(self, param_file_name) 
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
      integer(I4B), parameter :: LUN = 7                 !! Unit number of input file
      integer(I4B)            :: ierr = 0                !! Input error code
      character(STRMAX)       :: error_message           !! Error message in UDIO procedure

      ! Read in name of parameter file
      write(*, *) 'Parameter input file is ', trim(adjustl(param_file_name))
      write(*, *) ' '
      100 format(A)
      open(unit = LUN, file = param_file_name, status = 'old', iostat = ierr)
      if (ierr /= 0) then
         write(*,*) 'Swiftest error: ', ierr
         write(*,*) '   Unable to open file ',trim(adjustl(param_file_name))
         call util_exit(FAILURE)
      end if

      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    as the newline characters are ignored in the input file when compiled in ifort.

      !read(LUN,'(DT)', iostat= ierr, iomsg = error_message) param
      call self%reader(LUN, iotype= "none", v_list = [self%integrator], iostat = ierr, iomsg = error_message)
      if (ierr /= 0) then
         write(*,*) 'Swiftest error reading ', trim(adjustl(param_file_name))
         write(*,*) ierr,trim(adjustl(error_message))
         call util_exit(FAILURE)
      end if

      return 
   end subroutine io_read_param_in


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


   module subroutine io_write_discard(self, param)
      !! author: David A. Minton
      !!
      !! Write out information about discarded test particle
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_discard_write.f90
      !! Adapted from Hal Levison's Swift routine io_discard_write.f
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self     !! Swiftest system object
      class(swiftest_parameters),   intent(in)    :: param   !! Current run configuration parameters 
      ! Internals
      integer(I4B), parameter   :: LUN = 40
      integer(I4B)          :: i, ierr
      logical, save :: lfirst = .true. 
      real(DP), dimension(:,:), allocatable :: vh
      character(*), parameter :: HDRFMT    = '(E23.16, 1X, I8, 1X, L1)'
      character(*), parameter :: NAMEFMT   = '(A, 2(1X, I8))'
      character(*), parameter :: VECFMT    = '(3(E23.16, 1X))'
      character(*), parameter :: NPLFMT    = '(I8)'
      character(*), parameter :: PLNAMEFMT = '(I8, 2(1X, E23.16))'
      class(swiftest_body), allocatable :: pltemp

      associate(tp_discards => self%tp_discards, nsp => self%tp_discards%nbody, pl => self%pl, npl => self%pl%nbody)
         if (nsp == 0) return
         if (lfirst) then
            select case(param%out_stat)
            case('APPEND')
               open(unit = LUN, file = param%discard_out, status = 'OLD', position = 'APPEND', form = 'FORMATTED', iostat = ierr)
            case('NEW', 'REPLACE', 'UNKNOWN')
               open(unit = LUN, file = param%discard_out, status = param%out_stat, form = 'FORMATTED', iostat = ierr)
            case default
               write(*,*) 'Invalid status code for OUT_STAT: ',trim(adjustl(param%out_stat))
               call util_exit(FAILURE)
            end select
            lfirst = .false.
         else
            open(unit = LUN, file = param%discard_out, status = 'OLD', position = 'APPEND', form = 'FORMATTED', iostat = ierr)
         end if
         if (param%lgr) call tp_discards%pv2v(param) 

         write(LUN, HDRFMT) param%t, nsp, param%lbig_discard
         do i = 1, nsp
            write(LUN, NAMEFMT) SUB, tp_discards%id(i), tp_discards%status(i)
            write(LUN, VECFMT) tp_discards%xh(1, i), tp_discards%xh(2, i), tp_discards%xh(3, i)
            write(LUN, VECFMT) tp_discards%vh(1, i), tp_discards%vh(2, i), tp_discards%vh(3, i)
         end do
         if (param%lbig_discard) then
               if (param%lgr) then
                  allocate(pltemp, source = pl)
                  call pltemp%pv2v(param)
                  allocate(vh, source = pltemp%vh)
                  deallocate(pltemp)
               else
                  allocate(vh, source = pl%vh)
               end if

               write(LUN, NPLFMT) npl
               do i = 1, npl
                  write(LUN, PLNAMEFMT) pl%id(i), pl%Gmass(i), pl%radius(i)
                  write(LUN, VECFMT) pl%xh(1, i), pl%xh(2, i), pl%xh(3, i)
                  write(LUN, VECFMT) vh(1, i), vh(2, i), vh(3, i)
               end do
               deallocate(vh)
         end if
         close(LUN)
      end associate

      return
   end subroutine io_write_discard


   module subroutine io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
                                 xh1, xh2, vh1, vh2, enc_out, out_type)
      !! author: David A. Minton
      !!
      !! Write close encounter data to output binary files
      !!  There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_write_encounter.f90
      !! Adapted from Hal Levison's Swift routine io_write_encounter.f
      implicit none
      ! Arguments
      integer(I4B),           intent(in) :: name1, name2
      real(DP),               intent(in) :: t, mass1, mass2, radius1, radius2
      real(DP), dimension(:), intent(in) :: xh1, xh2, vh1, vh2
      character(*),           intent(in) :: enc_out, out_type
      ! Internals
      logical , save    :: lfirst = .true.
      integer(I4B), parameter :: lun = 30
      integer(I4B)        :: ierr
      integer(I4B), save    :: iu = lun

      open(unit = iu, file = enc_out, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
      if ((ierr /= 0) .and. lfirst) then
         open(unit = iu, file = enc_out, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
      end if
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to open binary encounter file"
         call util_exit(FAILURE)
      end if
      lfirst = .false.
      write(iu, iostat = ierr) t
      if (ierr < 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to write binary file record"
         call util_exit(FAILURE)
      end if
      write(iu) name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), mass1, radius1
      write(iu) name2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), mass2, radius2
      close(unit = iu, iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to close binary encounter file"
         call util_exit(FAILURE)
      end if

      return
   end subroutine io_write_encounter


   module subroutine io_write_frame_body(self, iu, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      implicit none
      ! Arguments
      class(swiftest_body),       intent(in)    :: self   !! Swiftest particle object
      integer(I4B),               intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 

      associate(n => self%nbody)
         if (n == 0) return
         write(iu) self%id(1:n)
         !write(iu) self%name(1:n)
         select case (param%out_form)
         case (EL) 
            write(iu) self%a(1:n)
            write(iu) self%e(1:n)
            write(iu) self%inc(1:n)
            write(iu) self%capom(1:n)
            write(iu) self%omega(1:n)
            write(iu) self%capm(1:n)
         case (XV)
            write(iu) self%xh(1, 1:n)
            write(iu) self%xh(2, 1:n)
            write(iu) self%xh(3, 1:n)
            write(iu) self%vh(1, 1:n)
            write(iu) self%vh(2, 1:n)
            write(iu) self%vh(3, 1:n)
         end select
         select type(pl => self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            write(iu) pl%Gmass(1:n)
            write(iu) pl%rhill(1:n)
            write(iu) pl%radius(1:n)
            if (param%lrotation) then
               write(iu) pl%rot(1, 1:n)
               write(iu) pl%rot(2, 1:n)
               write(iu) pl%rot(3, 1:n)
               write(iu) pl%Ip(1, 1:n)
               write(iu) pl%Ip(2, 1:n)
               write(iu) pl%Ip(3, 1:n)
            end if
            if (param%ltides) then
               write(iu) pl%k2(1:n)
               write(iu) pl%Q(1:n)
            end if
         end select
      end associate

      return
   end subroutine io_write_frame_body


   module subroutine io_write_frame_cb(self, iu, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of central body data to the binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      implicit none
      ! Arguments
      class(swiftest_cb),         intent(in)    :: self   !! Swiftest central body object 
      integer(I4B),               intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 

      associate(cb => self)
         !write(iu) cb%name
         write(iu) cb%id
         write(iu) cb%Gmass
         write(iu) cb%radius
         write(iu) cb%j2rp2 
         write(iu) cb%j4rp4 
         if (param%lrotation) then
            write(iu) cb%rot(1)
            write(iu) cb%rot(2)
            write(iu) cb%rot(3)
            write(iu) cb%Ip(1)
            write(iu) cb%Ip(2)
            write(iu) cb%Ip(3)
         end if
         if (param%ltides) then
            write(iu) cb%k2
            write(iu) cb%Q
         end if
      end associate

      return
   end subroutine io_write_frame_cb


   module subroutine io_write_frame_system(self, iu, param)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to output binary file
      !! There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(in)    :: self   !! Swiftest system object
      integer(I4B),                 intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
      ! Internals
      logical, save                    :: lfirst = .true. !! Flag to determine if this is the first call of this method
      integer(I4B)                     :: ierr            !! I/O error code

      class(swiftest_cb), allocatable  :: cb         !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_pl), allocatable  :: pl         !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_tp), allocatable  :: tp          !! Temporary local version of pl structure used for non-destructive conversions

      allocate(cb, source = self%cb)
      allocate(pl, source = self%pl)
      allocate(tp, source = self%tp)
      iu = BINUNIT

      if (lfirst) then
         select case(param%out_stat)
         case('APPEND')
            open(unit = iu, file = param%outfile, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
         case('NEW', 'REPLACE', 'UNKNOWN')
            open(unit = iu, file = param%outfile, status = param%out_stat, form = 'UNFORMATTED', iostat = ierr)
         case default
            write(*,*) 'Invalid status code for OUT_STAT: ',trim(adjustl(param%out_stat))
            call util_exit(FAILURE)
         end select
         if (ierr /= 0) then
            write(*, *) "Swiftest error: io_write_frame_system - first", ierr
            write(*, *) "   Binary output file " // trim(adjustl(param%outfile)) // " already exists or cannot be accessed"
            write(*, *) "   OUT_STAT: " // trim(adjustl(param%out_stat))
            call util_exit(FAILURE)
         end if
         lfirst = .false.
      else
         open(unit = iu, file = param%outfile, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest error: io_write_frame_system"
            write(*, *) "   Unable to open binary output file for APPEND"
            call util_exit(FAILURE)
         end if
      end if
      call io_write_hdr(iu, param%t, pl%nbody, tp%nbody, param%out_form, param%out_type)

      if (param%lgr) then
         call pl%pv2v(param)
         call tp%pv2v(param)
      end if

      if (param%out_form == EL) then ! Do an orbital element conversion prior to writing out the frame, as we have access to the central body here
         call pl%xv2el(cb)
         call tp%xv2el(cb)
      end if
      
      ! Write out each data type frame
      call cb%write_frame(iu, param)
      call pl%write_frame(iu, param)
      call tp%write_frame(iu, param)

      deallocate(cb, pl, tp)

      close(iu)

      return
   end subroutine io_write_frame_system


   subroutine io_write_hdr(iu, t, npl, ntp, out_form, out_type)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write frame header to output binary file
      !!
      !! Adapted from David Adapted from David E. Kaufmann's Swifter routine io_write_hdr.f90
      !! Adapted from Hal Levison's Swift routine io_write_hdr.F
      implicit none
      ! Arguments
      integer(I4B), intent(in) :: iu       !! Output file unit number
      real(DP),     intent(in) :: t        !! Current time of simulation
      integer(I4B), intent(in) :: npl      !! Number of massive bodies
      integer(I4B), intent(in) :: ntp      !! Number of test particles
      character(*), intent(in) :: out_form !! Output format type ("EL" or  "XV")
      character(*), intent(in) :: out_type !! Output file format type (REAL4, REAL8 - see swiftest module for symbolic name definitions)
      ! Internals
      integer(I4B)               :: ierr !! Error code
   
      select case (out_type)
      case (REAL4_TYPE,SWIFTER_REAL4_TYPE)
         write(iu, iostat = ierr) real(t, kind=SP)
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Unable to write binary file header"
            call util_exit(FAILURE)
         end if
      case (REAL8_TYPE,SWIFTER_REAL8_TYPE)
         write(iu, iostat = ierr) t
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Unable to write binary file header"
            call util_exit(FAILURE)
         end if
      end select
      write(iu, iostat = ierr) npl
      write(iu, iostat = ierr) ntp
      write(iu, iostat = ierr) out_form
   
      return
   end subroutine io_write_hdr

end submodule s_io
