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
      character(len=STRMAX)           :: errmsg
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
         if (param%energy_out /= "") then
            if (lfirst .and. (param%out_stat /= "OLD")) then
               open(unit = EGYIU, file = param%energy_out, form = "formatted", status = "replace", action = "write", err = 667, iomsg = errmsg)
               write(EGYIU,EGYHEADER, err = 667, iomsg = errmsg)
            else
               open(unit = EGYIU, file = param%energy_out, form = "formatted", status = "old", action = "write", position = "append", err = 667, iomsg = errmsg)
            end if
         end if
         call pl%vb2vh(cb)
         call pl%xh2xb(cb)
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

         if (param%energy_out /= "") then
            write(EGYIU,EGYFMT, err = 667, iomsg = errmsg) param%t, Eorbit_now, Ecollisions, Ltot_now, Mtot_now
            close(EGYIU, err = 667, iomsg = errmsg)
         end if
         if (.not.lfirst .and. lterminal) then 
            Lmag_now = norm2(Ltot_now)
            Lerror = norm2(Ltot_now - Ltot_orig) / Lmag_orig
            Eorbit_error = (Eorbit_now - Eorbit_orig) / abs(Eorbit_orig)
            Ecoll_error = Ecollisions / abs(Eorbit_orig)
            Etotal_error = (Eorbit_now - Ecollisions - Eorbit_orig - Euntracked) / abs(Eorbit_orig)
            Merror = (Mtot_now - Mtot_orig) / Mtot_orig
            write(*, egytermfmt) Lerror, Ecoll_error, Etotal_error, Merror
            ! if (Ecoll_error > 0.0_DP) then
            !    write(*,*) 'Something has gone wrong! Collisional energy should not be positive!'
            !    write(*,*) 'dke_orbit: ',(ke_orbit_now - ke_orbit_last) / abs(Eorbit_orig)
            !    write(*,*) 'dke_spin : ',(ke_spin_now - ke_spin_last) / abs(Eorbit_orig)
            !    write(*,*) 'dpe      : ',(pe_now - pe_last) / abs(Eorbit_orig)
            !    write(*,*)
            ! end if
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
      integer(I4B), parameter  :: LUN = 7       !! Unit number of output file
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


   module subroutine io_dump_swiftest(self, param)
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
      ! Internals
      integer(I4B)                   :: ierr    !! Error code
      integer(I4B),parameter         :: LUN = 7 !! Unit number for dump file
      integer(I4B)                   :: iu = LUN
      character(len=:), allocatable  :: dump_file_name
      character(STRMAX)              :: errmsg 

      select type(self)
      class is(swiftest_cb)
         dump_file_name = trim(adjustl(param%incbfile)) 
      class is (swiftest_pl)
         dump_file_name = trim(adjustl(param%inplfile)) 
      class is (swiftest_tp)
         dump_file_name = trim(adjustl(param%intpfile)) 
      end select
      open(unit = iu, file = dump_file_name, form = "UNFORMATTED", status = 'replace', err = 667, iomsg = errmsg)
      call self%write_frame(iu, param)
      close(LUN, err = 667, iomsg = errmsg)

      return

      667 continue
      write(*,*) "Error dumping body data to file " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_dump_swiftest


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
      real(DP)                      :: deltawall, wallperstep, tfrac
      integer(I8B)                  :: clock_count, count_rate, count_max
      character(*),     parameter   :: statusfmt   = '("Time = ", ES12.5, "; fraction done = ", F6.3, "; Number of active pl, tp = ", I5, ", ", I5)'
      character(len=*), parameter   :: walltimefmt = '("      Wall time (s): ", es12.5, "; Wall time/step in this interval (s):  ", es12.5)'
      logical, save                 :: lfirst = .true.
      real(DP), save                :: start, finish
    
      if (lfirst) then
         call system_clock(clock_count, count_rate, count_max)
         start = clock_count / (count_rate * 1.0_DP)
         finish = start
         lfirst = .false.
         if (param%lenergy) call self%conservation_report(param, lterminal=.false.)
      else
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

         tfrac = (param%t - param%t0) / (param%tstop - param%t0)
         
         call system_clock(clock_count, count_rate, count_max)
         deltawall = clock_count / (count_rate * 1.0_DP) - finish
         wallperstep = deltawall / param%istep_dump
         finish = clock_count / (count_rate * 1.0_DP)
         write(*, statusfmt) param%t, tfrac, self%pl%nbody, self%tp%nbody
         write(*, walltimefmt) finish - start, wallperstep
         if (param%lenergy) call self%conservation_report(param, lterminal=.true.)
      end if

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
               case ("ENC_OUT")
                  param%enc_out = param_value
               case ("DISCARD_OUT")
                  param%discard_out = param_value
               case ("ENERGY_OUT")
                  param%energy_out = param_value
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
               case ("FIRSTKICK")
                  call io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstkick = .false. 
               case ("FIRSTENERGY")
                  call io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstenergy = .false. 
               case("EORBIT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Eorbit_orig 
               case("MTOT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Mtot_orig 
               case("LTOT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Ltot_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Ltot_orig(i)
                  end do
                     param%Lmag_orig = norm2(param%Ltot_orig(:))
               case("LORBIT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Lorbit_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Lorbit_orig(i)
                  end do
               case("LSPIN_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Lspin_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Lspin_orig(i)
                  end do
               case("LESCAPE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Lescape(1)
                  do i = 2, NDIM
                     ifirst = ilast + 1
                     param_value = io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%Lescape(i)
                  end do
               case("MESCAPE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Mescape 
               case("ECOLLISIONS")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Ecollisions
               case("EUNTRACKED")
                  read(param_value, *, err = 667, iomsg = iomsg) param%Euntracked
               case ("NPLMAX", "NTPMAX", "GMTINY", "PARTICLE_OUT", "FRAGMENTATION", "SEED", "YARKOVSKY", "YORP") ! Ignore SyMBA-specific, not-yet-implemented, or obsolete input parameters
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
         param%lrestart = (param%out_stat == "APPEND")
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
         if (param%rmin > 0.0) then
            write(*,*) "CHK_RMIN       = ",param%rmin
         else
            write(*,*) "! CHK_RMIN value will be the central body radius"
         end if
         if (param%rmax > 0.0_DP) write(*,*) "CHK_RMAX       = ",param%rmax
         if (param%rmaxu > 0.0_DP) write(*,*) "CHK_EJECT      = ",param%rmaxu
         if ((param%qmin > 0.0_DP) .or. (param%qmin_alo > 0.0_DP) .or. (param%qmin_ahi > 0.0_DP)) write(*,*) "CHK_QMIN_COORD = ",trim(adjustl(param%qmin_coord))
         if (param%qmin > 0.0_DP) write(*,*) "CHK_QMIN       = ",param%qmin
         if ((param%qmin_alo > 0.0_DP) .or. (param%qmin_ahi > 0.0_DP))  write(*,*) "CHK_QMIN_RANGE = ",param%qmin_alo, param%qmin_ahi
         write(*,*) "EXTRA_FORCE    = ",param%lextra_force
         write(*,*) "RHILL_PRESENT  = ",param%lrhill_present
         write(*,*) "ROTATION      = ", param%lrotation
         write(*,*) "TIDES         = ", param%ltides
         write(*,*) "ENERGY         = ",param%lenergy
         if (param%lenergy) write(*,*) "ENERGY_OUT     = ",trim(adjustl(param%energy_out))
         write(*,*) "MU2KG          = ",param%MU2KG       
         write(*,*) "TU2S           = ",param%TU2S        
         write(*,*) "DU2M           = ",param%DU2M        
         if (trim(adjustl(param%outfile)) == "") then
            param%outfile = BIN_OUTFILE
         end if
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
            param%lbig_discard = .false.
            write(*,*) "! BIG_DISCARD    = ",param%lbig_discard
         end if

         if ((param%MU2KG < 0.0_DP) .or. (param%TU2S < 0.0_DP) .or. (param%DU2M < 0.0_DP)) then
            write(iomsg,*) 'Invalid unit conversion factor'
            iostat = -1
            return
         end if

         ! Calculate the G for the system units
         param%GU = GC / (param%DU2M**3 / (param%MU2KG * param%TU2S**2))
         if (param%lgr) then
            ! Calculate the inverse speed of light in the system units
            param%inv_c2 = einsteinC * param%TU2S / param%DU2M
            param%inv_c2 = (param%inv_c2)**(-2)
         end if

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
               write(*,*) "GR             = ", param%lgr
            case default   
               if (param%lgr) write(iomsg, *) 'GR is not yet implemented for this integrator. This parameter will be ignored.'
               param%lgr = .false.
            end select
         end associate

         iostat = 0
      end associate

      667 continue
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
         write(param_name, Afmt) "T0"; write(param_value,Rfmt) param%t0; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TSTOP"; write(param_value, Rfmt) param%tstop; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "DT"; write(param_value, Rfmt) param%dt; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "PL_IN"; write(param_value, Afmt) trim(adjustl(param%inplfile)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TP_in"; write(param_value, Afmt) trim(adjustl(param%intpfile)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "IN_TYPE"; write(param_value, Afmt) trim(adjustl(param%in_type)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         if (param%istep_out > 0) then
            write(param_name, Afmt) "ISTEP_OUT"; write(param_value, Ifmt) param%istep_out; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "BIN_OUT"; write(param_value, Afmt) trim(adjustl(param%outfile)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "OUT_TYPE"; write(param_value, Afmt) trim(adjustl(param%out_type)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "OUT_FORM"; write(param_value, Afmt) trim(adjustl(param%out_form)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "OUT_STAT"; write(param_value, Afmt) "APPEND"; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         end if
         write(param_name, Afmt) "ENC_OUT"; write(param_value, Afmt) trim(adjustl(param%enc_out)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         if (param%istep_dump > 0) then
            write(param_name, Afmt) "ISTEP_DUMP"; write(param_value, Ifmt) param%istep_dump; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         end if
         write(param_name, Afmt) "CHK_RMIN"; write(param_value, Rfmt) param%rmin; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_RMAX"; write(param_value, Rfmt) param%rmax; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_EJECT"; write(param_value, Rfmt) param%rmaxu; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_QMIN"; write(param_value, Rfmt) param%qmin; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         if (param%qmin >= 0.0_DP) then
            write(param_name, Afmt) "CHK_QMIN_COORD"; write(param_value, Afmt) trim(adjustl(param%qmin_coord)); write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            allocate(param_array(2))
            write(param_array(1)%value, Rfmt) param%qmin_alo
            write(param_array(2)%value, Rfmt) param%qmin_ahi
            write(param_name, Afmt) "CHK_QMIN_RANGE"; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_array(1)%value), adjustl(param_array(2)%value)
         end if
         write(param_name, Afmt) "MU2KG"; write(param_value, Rfmt) param%MU2KG; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TU2S"; write(param_value, Rfmt) param%TU2S ; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "DU2M"; write(param_value, Rfmt) param%DU2M; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "RHILL_PRESENT"; write(param_value, Lfmt) param%lrhill_present; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "EXTRA_FORCE"; write(param_value, Lfmt) param%lextra_force; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "BIG_DISCARD"; write(param_value, Lfmt) param%lbig_discard; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "CHK_CLOSE"; write(param_value, Lfmt) param%lclose; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "ENERGY"; write(param_value, Lfmt)  param%lenergy; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "GR"; write(param_value, Lfmt)  param%lgr; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "ROTATION"; write(param_value, Lfmt)  param%lrotation; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         write(param_name, Afmt) "TIDES"; write(param_value, Lfmt)  param%ltides; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)

         if (param%lenergy) then
            write(param_name, Afmt) "FIRSTENERGY"; write(param_value, Lfmt) param%lfirstenergy; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "EORBIT_ORIG"; write(param_value, Rfmt) param%Eorbit_orig; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "MTOT_ORIG"; write(param_value, Rfmt) param%Mtot_orig; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(unit, '("LTOT_ORIG  ",3(1X,ES25.17))') param%Ltot_orig(:)
            write(unit, '("LORBIT_ORIG",3(1X,ES25.17))') param%Lorbit_orig(:)
            write(unit, '("LSPIN_ORIG ",3(1X,ES25.17))') param%Lspin_orig(:)
            write(unit, '("LESCAPE ",3(1X,ES25.17))') param%Lescape(:)
   
            write(param_name, Afmt) "MESCAPE"; write(param_value, Rfmt) param%Mescape; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "ECOLLISIONS"; write(param_value, Rfmt) param%Ecollisions; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
            write(param_name, Afmt) "EUNTRACKED"; write(param_value, Rfmt) param%Euntracked; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
         end if
         write(param_name, Afmt) "FIRSTKICK"; write(param_value, Lfmt) param%lfirstkick; write(unit, Afmt, err = 667, iomsg = iomsg) adjustl(param_name), adjustl(param_value)
   
         iostat = 0
         iomsg = "UDIO not implemented"
      end associate

      667 continue

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
      integer(I4B)                  :: i, nbody
      logical                       :: is_ascii, is_pl
      character(len=:), allocatable :: infile
      real(DP)                      :: t
      real(QP)                      :: val
      character(STRMAX)             :: errmsg

      ! Select the appropriate polymorphic class (test particle or massive body)
      select type(self)
      class is (swiftest_pl)
         infile = param%inplfile
         is_pl = .true.
      class is (swiftest_tp)
         infile = param%intpfile
         is_pl = .false.
      end select

      is_ascii = (param%in_type == 'ASCII') 
      select case(param%in_type)
      case(ASCII_TYPE)
         open(unit = iu, file = infile, status = 'old', form = 'FORMATTED', err = 667, iomsg = errmsg)
         read(iu, *, err = 667, iomsg = errmsg) nbody
         call self%setup(nbody, param)
         if (nbody > 0) then
            do i = 1, nbody
               select type(self)
               class is (swiftest_pl)
                  if (param%lrhill_present) then
                     read(iu, *, err = 667, iomsg = errmsg) self%id(i), val, self%rhill(i)
                  else
                     read(iu, *, err = 667, iomsg = errmsg) self%id(i), val
                  end if
                  self%Gmass(i) = real(val, kind=DP)
                  self%mass(i) = real(val / param%GU, kind=DP)
                  if (param%lclose) read(iu, *, err = 667, iomsg = errmsg) self%radius(i)
               class is (swiftest_tp)
                  read(iu, *, err = 667, iomsg = errmsg) self%id(i)
               end select
               read(iu, *, err = 667, iomsg = errmsg) self%xh(1, i), self%xh(2, i), self%xh(3, i)
               read(iu, *, err = 667, iomsg = errmsg) self%vh(1, i), self%vh(2, i), self%vh(3, i)
               select type (self)
               class is (swiftest_pl)
                  if (param%lrotation) then
                     read(iu, *, err = 667, iomsg = errmsg) self%Ip(1, i), self%Ip(2, i), self%Ip(3, i)
                     read(iu, *, err = 667, iomsg = errmsg) self%rot(1, i), self%rot(2, i), self%rot(3, i)
                  end if
                  if (param%ltides) then
                     read(iu, *, err = 667, iomsg = errmsg) self%k2(i)
                     read(iu, *, err = 667, iomsg = errmsg) self%Q(i)
                  end if
               end select
               self%status(i) = ACTIVE
               self%lmask(i) = .true.
            end do
         end if
      case (REAL4_TYPE, REAL8_TYPE)  !, SWIFTER_REAL4_TYPE, SWIFTER_REAL8_TYPE)
         open(unit=iu, file=infile, status='old', form='UNFORMATTED', err = 667, iomsg = errmsg)
         read(iu, err = 667, iomsg = errmsg) nbody
         call self%setup(nbody, param)
         if (nbody > 0) then
            call self%read_frame(iu, param, XV)
            self%status(:) = ACTIVE
            self%lmask(:) = .true.
         end if
      case default
         write(errmsg,*) trim(adjustl(param%in_type)) // ' is an unrecognized file type'
         goto 667
      end select
      close(iu, err = 667, iomsg = errmsg)

      return

      667 continue
      write(*,*) 'Error reading in initial conditions file: ',trim(adjustl(errmsg))
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
      character(len=STRMAX)   :: errmsg

      if (param%in_type == 'ASCII') then
         open(unit = iu, file = param%incbfile, status = 'old', form = 'FORMATTED', err = 667, iomsg = errmsg)
         read(iu, *, err = 667, iomsg = errmsg) self%id
         read(iu, *, err = 667, iomsg = errmsg) self%Gmass
         self%mass = real(self%Gmass / param%GU, kind=DP)
         read(iu, *, err = 667, iomsg = errmsg) self%radius
         read(iu, *, err = 667, iomsg = errmsg) self%j2rp2
         read(iu, *, err = 667, iomsg = errmsg) self%j4rp4
         if (param%lrotation) then
            read(iu, *, err = 667, iomsg = errmsg) self%Ip(1), self%Ip(2), self%Ip(3)
            read(iu, *, err = 667, iomsg = errmsg) self%rot(1), self%rot(2), self%rot(3)
         end if
      else
         open(unit = iu, file = param%incbfile, status = 'old', form = 'UNFORMATTED', err = 667, iomsg = errmsg)
         call self%read_frame(iu, param, XV)
      end if
      close(iu, err = 667, iomsg = errmsg)

      if (self%j2rp2 /= 0.0_DP) param%loblatecb = .true.
      if (param%rmin < 0.0) param%rmin = self%radius

      return

      667 continue
      write(*,*) "Error reading central body file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_read_cb_in


   function io_read_encounter(t, id1, id2, Gmass1, Gmass2, radius1, radius2, &
      xh1, xh2, vh1, vh2, enc_out, out_type) result(ierr)
      !! author: David A. Minton
      !!
      !! Read close encounter data from input binary files
      !!     Other than time t, there is no direct file input from this function
      !!     Function returns read error status (0 = OK, nonzero = ERROR)
      !! Adapted from David E. Kaufmann's Swifter routine: io_read_encounter.f90
      implicit none
      ! Arguments
      integer(I4B),           intent(out) :: id1, id2
      real(DP),               intent(out) :: t, Gmass1, Gmass2, radius1, radius2
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

      read(iu, iostat = ierr) id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), vh1(3), Gmass1, radius1
      if (ierr /= 0) then
         close(unit = iu, iostat = ierr)
         return
      end if
      read(iu, iostat = ierr) id2, xh2(2), xh2(2), xh2(3), vh2(2), vh2(2), vh2(3), Gmass2, radius2
      if (ierr /= 0) then
         close(unit = iu, iostat = ierr)
         return
      end if

      return
   end function io_read_encounter


   module subroutine io_read_frame_body(self, iu, param, form)
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
      ! Internals
      character(len=STRMAX)   :: errmsg

      associate(n => self%nbody)
         read(iu, err = 667, iomsg = errmsg) self%id(:)
         !read(iu, err = 667, iomsg = errmsg) self%name(1:n)
         select case (form)
         case (EL) 
            if (.not.allocated(self%a))     allocate(self%a(n))
            if (.not.allocated(self%e))     allocate(self%e(n))
            if (.not.allocated(self%inc))   allocate(self%inc(n))
            if (.not.allocated(self%capom)) allocate(self%capom(n))
            if (.not.allocated(self%omega)) allocate(self%omega(n))
            if (.not.allocated(self%capm))  allocate(self%capm(n))
            read(iu, err = 667, iomsg = errmsg)  self%a(:)
            read(iu, err = 667, iomsg = errmsg)  self%e(:)
            read(iu, err = 667, iomsg = errmsg)  self%inc(:)
            read(iu, err = 667, iomsg = errmsg) self%capom(:)
            read(iu, err = 667, iomsg = errmsg) self%omega(:)
            read(iu, err = 667, iomsg = errmsg) self%capm(:)
         case (XV)
            read(iu, err = 667, iomsg = errmsg) self%xh(1, :)
            read(iu, err = 667, iomsg = errmsg) self%xh(2, :)
            read(iu, err = 667, iomsg = errmsg) self%xh(3, :)
            read(iu, err = 667, iomsg = errmsg) self%vh(1, :)
            read(iu, err = 667, iomsg = errmsg) self%vh(2, :)
            read(iu, err = 667, iomsg = errmsg) self%vh(3, :)
         end select
         select type(pl => self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            read(iu, err = 667, iomsg = errmsg) pl%Gmass(:)
            pl%mass(:) = pl%Gmass(:) / param%GU 
            if (param%lrhill_present) read(iu, err = 667, iomsg = errmsg) pl%rhill(:)
            if (param%lclose) read(iu, err = 667, iomsg = errmsg) pl%radius(:)
            if (param%lrotation) then
               read(iu, err = 667, iomsg = errmsg) pl%rot(1, :)
               read(iu, err = 667, iomsg = errmsg) pl%rot(2, :)
               read(iu, err = 667, iomsg = errmsg) pl%rot(3, :)
               read(iu, err = 667, iomsg = errmsg) pl%Ip(1, :)
               read(iu, err = 667, iomsg = errmsg) pl%Ip(2, :)
               read(iu, err = 667, iomsg = errmsg) pl%Ip(3, :)
            end if
            if (param%ltides) then
               read(iu, err = 667, iomsg = errmsg) pl%k2(1:n)
               read(iu, err = 667, iomsg = errmsg) pl%Q(1:n)
            end if
         end select
      end associate

      return

      667 continue
      write(*,*) "Error reading central body file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_read_frame_body


   module subroutine io_read_frame_cb(self, iu, param, form)
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
      ! Internals
      character(len=STRMAX)   :: errmsg

      read(iu, err = 667, iomsg = errmsg) self%id
      !read(iu, err = 667, iomsg = errmsg) self%name
      read(iu, err = 667, iomsg = errmsg) self%Gmass
      self%mass = self%Gmass / param%GU
      read(iu, err = 667, iomsg = errmsg) self%radius
      read(iu, err = 667, iomsg = errmsg) self%j2rp2 
      read(iu, err = 667, iomsg = errmsg) self%j4rp4 
      if (param%lrotation) then
         read(iu, err = 667, iomsg = errmsg) self%Ip(:)
         read(iu, err = 667, iomsg = errmsg) self%rot(:)
      end if
      if (param%ltides) then
         read(iu, err = 667, iomsg = errmsg) self%k2
         read(iu, err = 667, iomsg = errmsg) self%Q
      end if
      return

      667 continue
      write(*,*) "Error reading central body file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_read_frame_cb


   module subroutine io_read_frame_system(self, iu, param, form)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read a frame (header plus records for each massive body and active test particle) from a output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      integer(I4B),                 intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      character(*),                 intent(in)    :: form   !! Input format code ("XV" or "EL")
      ! Internals
      logical, save         :: lfirst = .true.
      character(len=STRMAX) :: errmsg
      integer(I4B)          :: ierr

      iu = BINUNIT
      if (lfirst) then
         open(unit = iu, file = param%outfile, status = 'OLD', form = 'UNFORMATTED', err = 667, iomsg = errmsg) 
         lfirst = .false.
      end if
      ierr = io_read_hdr(iu, param%t, self%pl%nbody, self%tp%nbody, param%out_form, param%out_type)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Binary output file already exists or cannot be accessed"
         return
      end if
      call self%cb%read_frame(iu, param, form)
      call self%pl%read_frame(iu, param, form)
      call self%tp%read_frame(iu, param, form)

      return

      667 continue
      write(*,*) "Error reading system frame: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
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
      real(SP)              :: ttmp
      character(len=STRMAX) :: errmsg

      select case (out_type)
      case (REAL4_TYPE, SWIFTER_REAL4_TYPE)
         read(iu, iostat = ierr, err = 667, iomsg = errmsg) ttmp, npl, ntp, out_form
         if (ierr /= 0) return
         t = ttmp
      case (REAL8_TYPE, SWIFTER_REAL8_TYPE)
         read(iu, iostat = ierr, err = 667, iomsg = errmsg) t
         read(iu, iostat = ierr, err = 667, iomsg = errmsg) npl
         read(iu, iostat = ierr, err = 667, iomsg = errmsg) ntp
         read(iu, iostat = ierr, err = 667, iomsg = errmsg) out_form
      case default
         write(errmsg,*) trim(adjustl(out_type)) // ' is an unrecognized file type'
         ierr = -1
      end select

      return

      667 continue
      write(*,*) "Error reading header: " // trim(adjustl(errmsg))

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
      character(STRMAX)       :: errmsg           !! Error message in UDIO procedure

      ! Read in name of parameter file
      write(*, *) 'Parameter input file is ', trim(adjustl(param_file_name))
      write(*, *) ' '
      100 format(A)
      open(unit = LUN, file = param_file_name, status = 'old', iostat = ierr, err = 667, iomsg = errmsg)

      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    as the newline characters are ignored in the input file when compiled in ifort.

      !read(LUN,'(DT)', iostat= ierr, iomsg = errmsg) param
      call self%reader(LUN, iotype= "none", v_list = [self%integrator], iostat = ierr, iomsg = errmsg)
      if (ierr == 0) return

      667 continue
      write(*,*) "Error reading parameter file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
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
      integer(I4B)          :: i
      logical, save :: lfirst = .true. 
      real(DP), dimension(:,:), allocatable :: vh
      character(*), parameter :: HDRFMT    = '(E23.16, 1X, I8, 1X, L1)'
      character(*), parameter :: NAMEFMT   = '(A, 2(1X, I8))'
      character(*), parameter :: VECFMT    = '(3(E23.16, 1X))'
      character(*), parameter :: NPLFMT    = '(I8)'
      character(*), parameter :: PLNAMEFMT = '(I8, 2(1X, E23.16))'
      class(swiftest_body), allocatable :: pltemp
      character(len=STRMAX)   :: errmsg, out_stat

      if (param%discard_out == "") return

      associate(tp_discards => self%tp_discards, nsp => self%tp_discards%nbody, pl => self%pl, npl => self%pl%nbody)
         if (nsp == 0) return
         if (lfirst) then
            out_stat = param%out_stat
         else
            out_stat = 'APPEND'
         end if
         select case(out_stat)
         case('APPEND')
            open(unit = LUN, file = param%discard_out, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
         case('NEW', 'REPLACE', 'UNKNOWN')
            open(unit = LUN, file = param%discard_out, status = param%out_stat, form = 'FORMATTED', err = 667, iomsg = errmsg)
         case default
            write(*,*) 'Invalid status code for OUT_STAT: ',trim(adjustl(param%out_stat))
            call util_exit(FAILURE)
         end select
         lfirst = .false.
         if (param%lgr) call tp_discards%pv2v(param) 

         write(LUN, HDRFMT) param%t, nsp, param%lbig_discard
         do i = 1, nsp
            write(LUN, NAMEFMT, err = 667, iomsg = errmsg) SUB, tp_discards%id(i), tp_discards%status(i)
            write(LUN, VECFMT, err = 667, iomsg = errmsg) tp_discards%xh(1, i), tp_discards%xh(2, i), tp_discards%xh(3, i)
            write(LUN, VECFMT, err = 667, iomsg = errmsg) tp_discards%vh(1, i), tp_discards%vh(2, i), tp_discards%vh(3, i)
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
                  write(LUN, PLNAMEFMT, err = 667, iomsg = errmsg) pl%id(i), pl%Gmass(i), pl%radius(i)
                  write(LUN, VECFMT, err = 667, iomsg = errmsg) pl%xh(1, i), pl%xh(2, i), pl%xh(3, i)
                  write(LUN, VECFMT, err = 667, iomsg = errmsg) vh(1, i), vh(2, i), vh(3, i)
               end do
               deallocate(vh)
         end if
         close(LUN)
      end associate

      return

      667 continue
      write(*,*) "Error writing discard file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_write_discard


   module subroutine io_write_encounter(self, pl, encbody, param)
      implicit none
      ! Arguments
      class(swiftest_encounter),  intent(in) :: self    !! Swiftest encounter list object
      class(swiftest_pl),         intent(in) :: pl      !! Swiftest massive body object
      class(swiftest_body),       intent(in) :: encbody !! Encountering body - Swiftest generic body object (pl or tp) 
      class(swiftest_parameters), intent(in) :: param   !! Current run configuration parameters 
      ! Internals
      logical , save          :: lfirst = .true.
      integer(I4B), parameter :: LUN = 30
      integer(I4B)            :: k, ierr
      character(len=STRMAX)   :: errmsg

      if (param%enc_out == "" .or. self%nenc == 0) return

      open(unit = LUN, file = param%enc_out, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr, iomsg = errmsg)
      if (ierr /= 0) then
         if (lfirst) then
            open(unit = LUN, file = param%enc_out, status = 'NEW', form = 'UNFORMATTED', err = 667, iomsg = errmsg)
         else
            goto 667
         end if
      end if
      lfirst = .false.

      associate(ind1 => self%index1, ind2 => self%index2)
         select type(encbody)
         class is (swiftest_pl)
            do k = 1, self%nenc
               call io_write_frame_encounter(LUN, self%t(k), &
                                             pl%id(ind1(k)),     encbody%id(ind2(k)), &
                                             pl%Gmass(ind1(k)),  encbody%Gmass(ind2(k)), &
                                             pl%radius(ind1(k)), encbody%radius(ind2(k)), &
                                             self%x1(:,k),       self%x2(:,k), &
                                             self%v1(:,k),       self%v2(:,k))
            end do
         class is (swiftest_tp)
            do k = 1, self%nenc
               call io_write_frame_encounter(LUN, self%t(k), &
                                             pl%id(ind1(k)),     encbody%id(ind2(k)), &
                                             pl%Gmass(ind1(k)),  0.0_DP, &
                                             pl%radius(ind1(k)), 0.0_DP, &
                                             self%x1(:,k),       self%x2(:,k), &
                                             self%v1(:,k),       self%v2(:,k))
            end do 
         end select
      end associate

      close(unit = LUN, err = 667, iomsg = errmsg)

      return
      667 continue
      write(*,*) "Error writing encounter file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
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
      ! Internals
      character(len=STRMAX)   :: errmsg

      associate(n => self%nbody)
         if (n == 0) return
         write(iu, err = 667, iomsg = errmsg) self%id(1:n)
         !write(iu, err = 667, iomsg = errmsg) self%name(1:n)
         select case (param%out_form)
         case (EL) 
            write(iu, err = 667, iomsg = errmsg) self%a(1:n)
            write(iu, err = 667, iomsg = errmsg) self%e(1:n)
            write(iu, err = 667, iomsg = errmsg) self%inc(1:n)
            write(iu, err = 667, iomsg = errmsg) self%capom(1:n)
            write(iu, err = 667, iomsg = errmsg) self%omega(1:n)
            write(iu, err = 667, iomsg = errmsg) self%capm(1:n)
         case (XV)
            write(iu, err = 667, iomsg = errmsg) self%xh(1, 1:n)
            write(iu, err = 667, iomsg = errmsg) self%xh(2, 1:n)
            write(iu, err = 667, iomsg = errmsg) self%xh(3, 1:n)
            write(iu, err = 667, iomsg = errmsg) self%vh(1, 1:n)
            write(iu, err = 667, iomsg = errmsg) self%vh(2, 1:n)
            write(iu, err = 667, iomsg = errmsg) self%vh(3, 1:n)
         end select
         select type(pl => self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            write(iu, err = 667, iomsg = errmsg) pl%Gmass(1:n)
            if (param%lrhill_present) write(iu, err = 667, iomsg = errmsg) pl%rhill(1:n)
            if (param%lclose) write(iu, err = 667, iomsg = errmsg) pl%radius(1:n)
            if (param%lrotation) then
               write(iu, err = 667, iomsg = errmsg) pl%Ip(1, 1:n)
               write(iu, err = 667, iomsg = errmsg) pl%Ip(2, 1:n)
               write(iu, err = 667, iomsg = errmsg) pl%Ip(3, 1:n)
               write(iu, err = 667, iomsg = errmsg) pl%rot(1, 1:n)
               write(iu, err = 667, iomsg = errmsg) pl%rot(2, 1:n)
               write(iu, err = 667, iomsg = errmsg) pl%rot(3, 1:n)
            end if
            if (param%ltides) then
               write(iu, err = 667, iomsg = errmsg) pl%k2(1:n)
               write(iu, err = 667, iomsg = errmsg) pl%Q(1:n)
            end if
         end select
      end associate

      return
      667 continue
      write(*,*) "Error writing body frame: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_write_frame_body

   module subroutine io_netcdf_write_frame_body(self, iu, param)
      !! author: David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      use netcdf
      implicit none
      ! Arguments
      class(swiftest_body),       intent(in)    :: self   !! Swiftest particle object
      integer(I4B),               intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters

      integer(I4B)                              :: i, j
      integer(I4B)                              :: ncid         !! NetCDF ID for the output file
      integer(I4B)                              :: dimids(2)    !! Dimensions of the NetCDF file
      integer(I4B)                              :: time_dimid   !! NetCDF ID for the time dimension 
      integer(I4B)                              :: name_dimid   !! NetCDF ID for the particle name dimension
      integer(I4B)                              :: ioutput      !! The current output number
      integer(I4B)                              :: a_varid      !! NetCDF ID for the semimajor axis variable 
      integer(I4B)                              :: e_varid      !! NetCDF ID for the eccentricity variable 
      integer(I4B)                              :: inc_varid    !! NetCDF ID for the inclination variable 
      integer(I4B)                              :: capom_varid  !! NetCDF ID for the long. asc. node variable 
      integer(I4B)                              :: omega_varid  !! NetCDF ID for the arg. periapsis variable 
      integer(I4B)                              :: capm_varid   !! NetCDF ID for the mean anomaly variable 
      integer(I4B)                              :: xhx_varid    !! NetCDF ID for the heliocentric position x variable 
      integer(I4B)                              :: xhy_varid    !! NetCDF ID for the heliocentric position y variable 
      integer(I4B)                              :: xhz_varid    !! NetCDF ID for the heliocentric position z variable 
      integer(I4B)                              :: vhx_varid    !! NetCDF ID for the heliocentric velocity x variable 
      integer(I4B)                              :: vhy_varid    !! NetCDF ID for the heliocentric velocity y variable 
      integer(I4B)                              :: vhz_varid    !! NetCDF ID for the heliocentric velocity z variable 
      integer(I4B)                              :: Gmass_varid  !! NetCDF ID for the mass variable
      integer(I4B)                              :: rhill_varid  !! NetCDF ID for the hill radius variable
      integer(I4B)                              :: radius_varid !! NetCDF ID for the radius variable
      integer(I4B)                              :: Ip1_varid    !! NetCDF ID for the axis 1 principal moment of inertial variable
      integer(I4B)                              :: Ip2_varid    !! NetCDF ID for the axis 2 principal moment of inertial variable
      integer(I4B)                              :: Ip3_varid    !! NetCDF ID for the axis 3 principal moment of inertial variable
      integer(I4B)                              :: rotx_varid   !! NetCDF ID for the rotation x variable
      integer(I4B)                              :: roty_varid   !! NetCDF ID for the rotation y variable
      integer(I4B)                              :: rotz_varid   !! NetCDF ID for the rotation z variable
      integer(I4B)                              :: k2_varid     !! NetCDF ID for the Love number variable
      integer(I4B)                              :: Q_varid      !! NetCDF ID for the energy dissipation variable

      !! Open the netCDF file
      call check( nf90_open(param%outfile, nf90_write, ncid) )

      associate(n => self%nbody)
         if (n == 0) return

      !! Calculate the output number that we are currently on
      ioutput = (param%t / param%dt) / param%istep_out
      !call check( nf90_inq_varid(ncid, "Time", ioutput))
      !call check( nf90_inquire_dimension(ncid, time_dimid, len=ioutput))

         select case (param%out_form)
         case (EL) 
            do j = 1, n
               do i = 1, n
                  if (self%id(i) == j) then

                     !! Reassign all variable IDs
                     call check( nf90_inq_varid(ncid, "a", a_varid))
                     call check( nf90_inq_varid(ncid, "e", e_varid))
                     call check( nf90_inq_varid(ncid, "inc", inc_varid))
                     call check( nf90_inq_varid(ncid, "capom", capom_varid))
                     call check( nf90_inq_varid(ncid, "omega", omega_varid))
                     call check( nf90_inq_varid(ncid, "capm", capm_varid))

                     call check( nf90_put_var(ncid, a_varid, self%a(j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, e_varid, self%e(j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, inc_varid, self%inc(j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, capom_varid, self%capom(j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, omega_varid, self%omega(j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, capm_varid, self%capm(j), start=(/ioutput + 1, j/)) )
                  end if
               end do 
            end do 
         case (XV)
            do j = 1, n
               do i = 1, n 
                  if (self%id(i) == j) then

                     !! Reassign all variable IDs
                     call check( nf90_inq_varid(ncid, "xhx", xhx_varid))
                     call check( nf90_inq_varid(ncid, "xhy", xhy_varid))
                     call check( nf90_inq_varid(ncid, "xhz", xhz_varid))
                     call check( nf90_inq_varid(ncid, "vhx", vhx_varid))
                     call check( nf90_inq_varid(ncid, "vhy", vhy_varid))
                     call check( nf90_inq_varid(ncid, "vhz", vhz_varid))

                     call check( nf90_put_var(ncid, xhx_varid, self%xh(1, j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, xhy_varid, self%xh(2, j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, xhz_varid, self%xh(3, j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, vhx_varid, self%vh(1, j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, vhy_varid, self%vh(2, j), start=(/ioutput + 1, j/)) )
                     call check( nf90_put_var(ncid, vhz_varid, self%vh(3, j), start=(/ioutput + 1, j/)) )
                  end if
               end do 
            end do 
         end select
         select type(pl => self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            do j = 1, n
               do i = 1, n
                  if (self%id(i) == j) then

                     !! Reassign all variable IDs
                     call check( nf90_inq_varid(ncid, "Gmass", Gmass_varid))
                     call check( nf90_put_var(ncid, Gmass_varid, pl%Gmass(j), start=(/ioutput + 1, j/)) )
                     if (param%lrhill_present) then 
                        !! Reassign all variable IDs
                        call check( nf90_inq_varid(ncid, "rhill", rhill_varid))
                        call check( nf90_put_var(ncid, rhill_varid, pl%rhill(j), start=(/ioutput + 1, j/)) )
                     end if
                     if (param%lclose) then
                        !! Reassign all variable IDs
                        call check( nf90_inq_varid(ncid, "radius", radius_varid))
                        call check( nf90_put_var(ncid, radius_varid, pl%radius(j), start=(/ioutput + 1, j/)) )
                     end if
                     if (param%lrotation) then

                        !! Reassign all variable IDs
                        call check( nf90_inq_varid(ncid, "Ip1", Ip1_varid))
                        call check( nf90_inq_varid(ncid, "Ip2", Ip2_varid))
                        call check( nf90_inq_varid(ncid, "Ip3", Ip3_varid))
                        call check( nf90_inq_varid(ncid, "rotx", rotx_varid))
                        call check( nf90_inq_varid(ncid, "roty", roty_varid))
                        call check( nf90_inq_varid(ncid, "rotz", rotz_varid))

                        call check( nf90_put_var(ncid, Ip1_varid, pl%Ip(1, j), start=(/ioutput + 1, j/)) )
                        call check( nf90_put_var(ncid, Ip2_varid, pl%Ip(2, j), start=(/ioutput + 1, j/)) )
                        call check( nf90_put_var(ncid, Ip3_varid, pl%Ip(3, j), start=(/ioutput + 1, j/)) )
                        call check( nf90_put_var(ncid, rotx_varid, pl%rot(1, j), start=(/ioutput + 1, j/)) )
                        call check( nf90_put_var(ncid, roty_varid, pl%rot(2, j), start=(/ioutput + 1, j/)) )
                        call check( nf90_put_var(ncid, rotz_varid, pl%rot(3, j), start=(/ioutput + 1, j/)) )
                     end if
                     if (param%ltides) then

                        !! Reassign all variable IDs
                        call check( nf90_inq_varid(ncid, "k2", k2_varid))
                        call check( nf90_inq_varid(ncid, "Q", Q_varid))

                        call check( nf90_put_var(ncid, k2_varid, pl%k2(j), start=(/ioutput + 1, j/)) )
                        call check( nf90_put_var(ncid, Q_varid, pl%Q(j), start=(/ioutput + 1, j/)) )
                     end if
                  end if
               end do
            end do
         end select
      end associate

      !! Close the netCDF file
      call check( nf90_close(ncid) )

      contains

      !! Checks the status of all NetCDF operations to catch errors
      subroutine check(status)
         integer, intent ( in) :: status

         if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "NetCDF Error: Stopped"
         end if
      end subroutine check

      !return
   end subroutine io_netcdf_write_frame_body


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
      ! Internals
      character(len=STRMAX)   :: errmsg

      associate(cb => self)
         !write(iu, err = 667, iomsg = errmsg) cb%name
         write(iu, err = 667, iomsg = errmsg) cb%id
         write(iu, err = 667, iomsg = errmsg) cb%Gmass
         write(iu, err = 667, iomsg = errmsg) cb%radius
         write(iu, err = 667, iomsg = errmsg) cb%j2rp2 
         write(iu, err = 667, iomsg = errmsg) cb%j4rp4 
         if (param%lrotation) then
            write(iu, err = 667, iomsg = errmsg) cb%Ip(1)
            write(iu, err = 667, iomsg = errmsg) cb%Ip(2)
            write(iu, err = 667, iomsg = errmsg) cb%Ip(3)
            write(iu, err = 667, iomsg = errmsg) cb%rot(1)
            write(iu, err = 667, iomsg = errmsg) cb%rot(2)
            write(iu, err = 667, iomsg = errmsg) cb%rot(3)
         end if
         if (param%ltides) then
            write(iu, err = 667, iomsg = errmsg) cb%k2
            write(iu, err = 667, iomsg = errmsg) cb%Q
         end if
      end associate

      return
      667 continue
      write(*,*) "Error writing central body frame: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_write_frame_cb


  module subroutine io_write_frame_encounter(iu, t, id1, id2, Gmass1, Gmass2, radius1, radius2, xh1, xh2, vh1, vh2)
      !! author: David A. Minton
      !!
      !! Write a single frame of close encounter data to output binary files
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_write_encounter.f90
      !! Adapted from Hal Levison's Swift routine io_write_encounter.f
      implicit none
      ! Arguments
      integer(I4B),           intent(in) :: iu               !! Open file unit number
      real(DP),               intent(in) :: t                !! Time of encounter
      integer(I4B),           intent(in) :: id1, id2         !! ids of the two encountering bodies
      real(DP),               intent(in) :: Gmass1, Gmass2   !! G*mass of the two encountering bodies
      real(DP),               intent(in) :: radius1, radius2 !! Radii of the two encountering bodies
      real(DP), dimension(:), intent(in) :: xh1, xh2         !! Heliocentric position vectors of the two encountering bodies 
      real(DP), dimension(:), intent(in) :: vh1, vh2         !! Heliocentric velocity vectors of the two encountering bodies  
      ! Internals
      character(len=STRMAX)   :: errmsg

      write(iu, err = 667, iomsg = errmsg) t
      write(iu, err = 667, iomsg = errmsg) id1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), Gmass1, radius1
      write(iu, err = 667, iomsg = errmsg) id2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), Gmass2, radius2

      return
      667 continue
      write(*,*) "Error writing encounter file: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine

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
      class(swiftest_cb), allocatable  :: cb         !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_pl), allocatable  :: pl         !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_tp), allocatable  :: tp          !! Temporary local version of pl structure used for non-destructive conversions
      character(len=STRMAX)            :: errmsg

      allocate(cb, source = self%cb)
      allocate(pl, source = self%pl)
      allocate(tp, source = self%tp)
      iu = BINUNIT

      if (lfirst) then
         select case(param%out_stat)
         case('APPEND')
            open(unit = iu, file = param%outfile, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', err = 667, iomsg = errmsg)
         case('NEW', 'REPLACE', 'UNKNOWN')
            open(unit = iu, file = param%outfile, status = param%out_stat, form = 'UNFORMATTED', err = 667, iomsg = errmsg)
         case default
            write(*,*) 'Invalid status code for OUT_STAT: ',trim(adjustl(param%out_stat))
            call util_exit(FAILURE)
         end select
         lfirst = .false.
      else
         open(unit = iu, file = param%outfile, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', err = 667, iomsg = errmsg)
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

      close(iu, err = 667, iomsg = errmsg)

      return
      667 continue
      write(*,*) "Error writing system frame: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_write_frame_system

   module subroutine io_netcdf_write_frame_system(self, iu, param)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to output binary file
      !! There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      use netcdf
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(in)    :: self            !! Swiftest system object
      integer(I4B),                 intent(inout) :: iu              !! Unit number for the output file to write frame to
      class(swiftest_parameters),   intent(in)    :: param           !! Current run configuration parameters 
      ! Internals
      logical, save                               :: lfirst = .true. !! Flag to determine if this is the first call of this method
      integer(I4B)                                :: ierr            !! I/O error code

      class(swiftest_cb), allocatable             :: cb              !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_pl), allocatable             :: pl              !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_tp), allocatable             :: tp              !! Temporary local version of pl structure used for non-destructive conversions
      
      integer(I4B)                                :: ncid            !! NetCDF ID for the output file
      integer(I4B)                                :: dimids(2)       !! Dimensions of the NetCDF file
      integer(I4B)                                :: time_dimid      !! NetCDF ID for the time dimension 
      integer(I4B)                                :: name_dimid      !! NetCDF ID for the particle name dimension
      integer(I4B)                                :: noutput         !! Number of output events covering the total simulation time
      integer(I4B)                                :: a_varid         !! NetCDF ID for the semimajor axis variable 
      integer(I4B)                                :: e_varid         !! NetCDF ID for the eccentricity variable 
      integer(I4B)                                :: inc_varid       !! NetCDF ID for the inclination variable 
      integer(I4B)                                :: capom_varid     !! NetCDF ID for the long. asc. node variable 
      integer(I4B)                                :: omega_varid     !! NetCDF ID for the arg. periapsis variable 
      integer(I4B)                                :: capm_varid      !! NetCDF ID for the mean anomaly variable 
      integer(I4B)                                :: xhx_varid       !! NetCDF ID for the heliocentric position x variable 
      integer(I4B)                                :: xhy_varid       !! NetCDF ID for the heliocentric position y variable 
      integer(I4B)                                :: xhz_varid       !! NetCDF ID for the heliocentric position z variable 
      integer(I4B)                                :: vhx_varid       !! NetCDF ID for the heliocentric velocity x variable 
      integer(I4B)                                :: vhy_varid       !! NetCDF ID for the heliocentric velocity y variable 
      integer(I4B)                                :: vhz_varid       !! NetCDF ID for the heliocentric velocity z variable 
      integer(I4B)                                :: Gmass_varid     !! NetCDF ID for the mass variable
      integer(I4B)                                :: rhill_varid     !! NetCDF ID for the hill radius variable
      integer(I4B)                                :: radius_varid    !! NetCDF ID for the radius variable
      integer(I4B)                                :: Ip1_varid       !! NetCDF ID for the axis 1 principal moment of inertial variable
      integer(I4B)                                :: Ip2_varid       !! NetCDF ID for the axis 2 principal moment of inertial variable
      integer(I4B)                                :: Ip3_varid       !! NetCDF ID for the axis 3 principal moment of inertial variable
      integer(I4B)                                :: rotx_varid      !! NetCDF ID for the rotation x variable
      integer(I4B)                                :: roty_varid      !! NetCDF ID for the rotation y variable
      integer(I4B)                                :: rotz_varid      !! NetCDF ID for the rotation z variable
      integer(I4B)                                :: k2_varid        !! NetCDF ID for the Love number variable
      integer(I4B)                                :: Q_varid         !! NetCDF ID for the energy dissipation variable

      allocate(cb, source = self%cb)
      allocate(pl, source = self%pl)
      allocate(tp, source = self%tp)
      iu = BINUNIT

      if (param%lgr) then
         call pl%pv2v(param)
         call tp%pv2v(param)
      end if

      if (param%out_form == EL) then ! Do an orbital element conversion prior to writing out the frame, as we have access to the central body here
         call pl%xv2el(cb)
         call tp%xv2el(cb)
      end if

      if (lfirst) then
         select case(param%out_stat)
         case('APPEND')

            call cb%write_frame(iu, param)
            call pl%write_frame(iu, param)
            call tp%write_frame(iu, param)

         case('NEW', 'REPLACE', 'UNKNOWN')
         
            !! Create the new output file, deleting any previously existing output file of the same name
            call check( nf90_create(param%outfile, NF90_HDF5, ncid) )

            !! Calculate the number of outputs needed to cover the entire simulation time
            noutput = ((param%tstop / param%dt) / param%istep_out) + 2 !! +2 because t=0 gets put in spot 1 and need a stop for the final output

            !! Define the NetCDF dimensions with particle name as the record dimension
            call check( nf90_def_dim(ncid, "Name", NF90_UNLIMITED, name_dimid) )     !! 'x' dimension
            call check( nf90_def_dim(ncid, "Time", NF90_UNLIMITED, time_dimid) )     !! 'y' dimension
            dimids = (/ time_dimid, name_dimid /)

            !! Define the variables
            select case (param%out_form)
            case (EL)
               call check( nf90_def_var(ncid, "a", NF90_FLOAT, dimids, a_varid) )
               call check( nf90_def_var(ncid, "e", NF90_FLOAT, dimids, e_varid) )
               call check( nf90_def_var(ncid, "inc", NF90_FLOAT, dimids, inc_varid) )
               call check( nf90_def_var(ncid, "capom", NF90_FLOAT, dimids, capom_varid) )
               call check( nf90_def_var(ncid, "omega", NF90_FLOAT, dimids, omega_varid) )
               call check( nf90_def_var(ncid, "capm", NF90_FLOAT, dimids, capm_varid) )
            case (XV)
               call check( nf90_def_var(ncid, "xhx", NF90_FLOAT, dimids, xhx_varid) )
               call check( nf90_def_var(ncid, "xhy", NF90_FLOAT, dimids, xhy_varid) )
               call check( nf90_def_var(ncid, "xhz", NF90_FLOAT, dimids, xhz_varid) )
               call check( nf90_def_var(ncid, "vhx", NF90_FLOAT, dimids, vhx_varid) )
               call check( nf90_def_var(ncid, "vhy", NF90_FLOAT, dimids, vhy_varid) )
               call check( nf90_def_var(ncid, "vhz", NF90_FLOAT, dimids, vhz_varid) )
            end select

            call check( nf90_def_var(ncid, "Gmass", NF90_FLOAT, dimids, Gmass_varid) )
            if (param%lrhill_present) call check( nf90_def_var(ncid, "rhill", NF90_FLOAT, dimids, rhill_varid) )
            if (param%lclose) call check( nf90_def_var(ncid, "radius", NF90_FLOAT, dimids, radius_varid) )
            if (param%lrotation) then
               call check( nf90_def_var(ncid, "Ip1", NF90_FLOAT, dimids, Ip1_varid) )
               call check( nf90_def_var(ncid, "Ip2", NF90_FLOAT, dimids, Ip3_varid) )
               call check( nf90_def_var(ncid, "Ip3", NF90_FLOAT, dimids, Ip3_varid) )
               call check( nf90_def_var(ncid, "rotx", NF90_FLOAT, dimids, rotx_varid) )
               call check( nf90_def_var(ncid, "roty", NF90_FLOAT, dimids, roty_varid) )
               call check( nf90_def_var(ncid, "rotz", NF90_FLOAT, dimids, rotz_varid) )
            end if
            if (param%ltides) then
               call check( nf90_def_var(ncid, "k2", NF90_FLOAT, dimids, k2_varid) )
               call check( nf90_def_var(ncid, "Q", NF90_FLOAT, dimids, Q_varid) )
            end if

            !! Exit define mode 
            call check( nf90_enddef(ncid) )

            !! Close the netCDF file
            call check( nf90_close(ncid) )

            !! Write the first frame of the output 
            call cb%write_frame(iu, param)
            call pl%write_frame(iu, param)
            call tp%write_frame(iu, param)

         case default
            write(*,*) 'Invalid status code for OUT_STAT: ',trim(adjustl(param%out_stat))
            call util_exit(FAILURE)
         end select

         lfirst = .false.
         
      else

         call cb%write_frame(iu, param)
         call pl%write_frame(iu, param)
         call tp%write_frame(iu, param)

      end if
      
      deallocate(cb, pl, tp)

   contains

   !! Checks the status of all NetCDF operations to catch errors
   subroutine check(status)
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "NetCDF Error: Stopped"
      end if
   end subroutine check

      !return
   end subroutine io_netcdf_write_frame_system


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
      character(len=STRMAX) :: errmsg
   
      select case (out_type)
      case (REAL4_TYPE,SWIFTER_REAL4_TYPE)
         write(iu, err = 667, iomsg = errmsg) real(t, kind=SP)
      case (REAL8_TYPE,SWIFTER_REAL8_TYPE)
         write(iu, err = 667, iomsg = errmsg) t
      end select
      write(iu, err = 667, iomsg = errmsg) npl
      write(iu, err = 667, iomsg = errmsg) ntp
      write(iu, err = 667, iomsg = errmsg) out_form
   
      return

      667 continue
      write(*,*) "Error writing header: " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine io_write_hdr

end submodule s_io
