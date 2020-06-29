submodule (swiftest_classes) io_read
   !! author: David A. Minton
   !! 
   !! This submodule contains implementations of the following procedures:
   !!    io_read_config_in
   !!    io_config_reader
   !!    io_get_token
   !!    io_read_initialize_system
   !!    io_read_cb_in
   !!    io_read_body_in
   !!    io_read_hdr
   !!    io_read_frame_system
   !!    io_read_frame_cb
   !!    io_read_frame_body
   !!    io_read_encounter
contains
   module procedure io_read_config_in
      !! author: David A. Minton
      !!
      !! Read in parameters for the integration
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_config.f90
      !! Adapted from Martin Duncan's Swift routine io_init_config.f
      use swiftest
      implicit none

      integer(I4B), parameter :: LUN = 7                 !! Unit number of input file
      integer(I4B)            :: ierr = 0                !! Input error code
      character(STRMAX)       :: error_message           !! Error message in UDIO procedure

      ! Read in name of parameter file
      write(*, *) 'Configuration data file is ', trim(adjustl(config_file_name))
      write(*, *) ' '
      100 format(A)
      open(unit = LUN, file = config_file_name, status = 'old', iostat = ierr)
      if (ierr /= 0) then
         write(*,*) 'Swiftest error: ', ierr
         write(*,*) '   Unable to open file ',trim(adjustl(config_file_name))
         call util_exit(FAILURE)
      end if

      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    as the newline characters are ignored in the input file when compiled in ifort.

      !read(LUN,'(DT)', iostat= ierr, iomsg = error_message) config
      call self%config_reader(LUN, iotype= "none", v_list = [integrator], iostat = ierr, iomsg = error_message)
      if (ierr /= 0) then
         write(*,*) 'Swiftest error reading ', trim(adjustl(config_file_name))
         write(*,*) ierr,trim(adjustl(error_message))
         call util_exit(FAILURE)
      end if

      return 
   end procedure io_read_config_in

   module procedure io_config_reader
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in parameters for the integration
      !! Currently this procedure does not work in user-defined derived-type input mode 
      !!    e.g. read(unit,'(DT)') param 
      !! as the newline characters are ignored in the input file when compiled in ifort.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_config.f90
      !! Adapted from Martin Duncan's Swift routine io_init_config.f
      use swiftest
      !$ use omp_lib
      !use util, only: util_exit ! IMPLEMENTATION TBD
      implicit none

      logical                        :: t0_set = .false.        !! Is the initial time set in the input file?
      logical                        :: tstop_set = .false.     !! Is the final time set in the input file?
      logical                        :: dt_set = .false.        !! Is the step size set in the input file?
      logical                        :: mtiny_set = .false.     !! Is the mtiny value set?
      integer(I4B)                   :: ilength, ifirst, ilast  !! Variables used to parse input file
      character(STRMAX)              :: line                    !! Line of the input file
      character (len=:), allocatable :: line_trim,config_name, config_value !! Strings used to parse the config file
      character(*),parameter         :: linefmt = '(A)'         !! Format code for simple text string
      integer(I4B)                   :: integrator              !! Symbolic name of integrator being used

      integrator = v_list(1)
      ! Parse the file line by line, extracting tokens then matching them up with known parameters if possible
      do
         read(unit = unit, fmt = linefmt, iostat = iostat, end = 1) line
         line_trim = trim(adjustl(line))
         ilength = len(line_trim)
         if ((ilength /= 0)) then 
            ifirst = 1
            ! Read the pair of tokens. The first one is the parameter name, the second is the value.
            config_name = io_get_token(line_trim, ifirst, ilast, iostat)
            if (config_name == '') cycle ! No parameter name (usually because this line is commented out)
            call util_toupper(config_name)
            ifirst = ilast + 1
            config_value = io_get_token(line_trim, ifirst, ilast, iostat)
            select case (config_name)
            case ("NPLMAX")
               read(config_value, *) self%nplmax
            case ("NTPMAX")
               read(config_value, *) self%ntpmax
            case ("T0")
               read(config_value, *) self%t0
               t0_set = .true.
            case ("TSTOP")
               read(config_value, *) self%tstop
               tstop_set = .true.
            case ("DT")
               read(config_value, *) self%dt
            case ("CB_IN")
               self%incbfile = config_value
            case ("PL_IN")
               self%inplfile = config_value
            case ("TP_IN")
               self%intpfile = config_value
            case ("IN_TYPE")
               call util_toupper(config_value)
               self%in_type = config_value
            case ("ISTEP_OUT")
               read(config_value, *) self%istep_out
            case ("BIN_OUT")
               self%outfile = config_value
            case ("OUT_TYPE")
               call util_toupper(config_value)
               self%out_type = config_value
            case ("OUT_FORM")
               call util_toupper(config_value)
               self%out_form = config_value
            case ("OUT_STAT")
               call util_toupper(config_value)
               self%out_stat = config_value
            case ("ISTEP_DUMP")
               read(config_value, *) self%istep_dump
            case ("CHK_CLOSE")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%lclose = .true.
            case ("CHK_RMIN")
               read(config_value, *) self%rmin
            case ("CHK_RMAX")
               read(config_value, *) self%rmax
            case ("CHK_EJECT")
               read(config_value, *) self%rmaxu
            case ("CHK_QMIN")
               read(config_value, *) self%qmin
            case ("CHK_QMIN_COORD")
               call util_toupper(config_value)
               self%qmin_coord = config_value
            case ("CHK_QMIN_RANGE")
               read(config_value, *) self%qmin_alo
               ifirst = ilast + 1
               config_value = io_get_token(line, ifirst, ilast, iostat)
               read(config_value, *) self%qmin_ahi
            case ("ENC_OUT")
               self%encounter_file = config_value
            case ("EXTRA_FORCE")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%lextra_force = .true.
            case ("BIG_DISCARD")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T' ) self%lbig_discard = .true.
            case ("FRAGMENTATION")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == "T") self%lfragmentation = .true.
            case ("MU2KG")
               read(config_value, *) self%MU2KG
            case ("TU2S")
               read(config_value, *) self%TU2S
            case ("DU2M")
               read(config_value, *) self%DU2M
            case ("MTINY")
               read(config_value, *) self%mtiny
               mtiny_set = .true.
            case ("ENERGY")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%lenergy = .true.
            case ("ROTATION")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%lrotation = .true. 
            case ("TIDES")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%ltides = .true. 
            case ("GR")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%lgr = .true. 
            case ("YARKOVSKY")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%lyarkovsky = .true. 
            case ("YORP")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') self%lyorp = .true. 
            case default
               write(iomsg,*) "Unknown parameter -> ",config_name
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
      if (self%dt <= 0.0_DP) then
         write(iomsg,*) 'Invalid timestep: '
         iostat = -1
         return
      end if
      if (self%inplfile == "") then
         write(iomsg,*) 'No valid massive body file in input file'
         iostat = -1
         return
      end if
      if ((self%in_type /= REAL8_TYPE) .and. (self%in_type /= "ASCII")) then
         write(iomsg,*) 'Invalid input file type:',trim(adjustl(self%in_type))
         iostat = -1
         return
      end if
      if ((self%istep_out <= 0) .and. (self%istep_dump <= 0)) then
         write(iomsg,*) 'Invalid istep'
         iostat = -1
         return
      end if
      if ((self%istep_out > 0) .and. (self%outfile == "")) then
         write(iomsg,*) 'Invalid outfile'
         iostat = -1
         return
      end if
      if (self%outfile /= "") then
         if ((self%out_type /= REAL4_TYPE) .and. (self%out_type /= REAL8_TYPE) .and. &
               (self%out_type /= SWIFTER_REAL4_TYPE)  .and. (self%out_type /= SWIFTER_REAL8_TYPE)) then
            write(iomsg,*) 'Invalid out_type: ',trim(adjustl(self%out_type))
            iostat = -1
            return
         end if
         if ((self%out_form /= "EL") .and. (self%out_form /= "XV")) then
            write(iomsg,*) 'Invalid out_form: ',trim(adjustl(self%out_form))
            iostat = -1
            return
         end if
         if ((self%out_stat /= "NEW") .and. (self%out_stat /= "REPLACE") .and. (self%out_stat /= "APPEND")) then
            write(iomsg,*) 'Invalid out_stat: ',trim(adjustl(self%out_stat))
            iostat = -1
            return
         end if
      end if
      if (self%qmin > 0.0_DP) then
         if ((self%qmin_coord /= "HELIO") .and. (self%qmin_coord /= "BARY")) then
            write(iomsg,*) 'Invalid qmin_coord: ',trim(adjustl(self%qmin_coord))
            iostat = -1
            return
         end if
         if ((self%qmin_alo <= 0.0_DP) .or. (self%qmin_ahi <= 0.0_DP)) then
            write(iomsg,*) 'Invalid qmin vals'
            iostat = -1
            return
         end if
      end if

      write(*,*) "NPLMAX         = ",self%nplmax
      write(*,*) "NTPMAX         = ",self%ntpmax
      write(*,*) "T0             = ",self%t0
      write(*,*) "TSTOP          = ",self%tstop
      write(*,*) "DT             = ",self%dt
      write(*,*) "PL_IN          = ",trim(adjustl(self%inplfile))
      write(*,*) "TP_IN          = ",trim(adjustl(self%intpfile))
      write(*,*) "IN_TYPE        = ",trim(adjustl(self%in_type))
      write(*,*) "ISTEP_OUT      = ",self%istep_out
      write(*,*) "BIN_OUT        = ",trim(adjustl(self%outfile))
      write(*,*) "OUT_TYPE       = ",trim(adjustl(self%out_type))
      write(*,*) "OUT_FORM       = ",trim(adjustl(self%out_form))
      write(*,*) "OUT_STAT       = ",trim(adjustl(self%out_stat))
      write(*,*) "ISTEP_DUMP     = ",self%istep_dump
      write(*,*) "CHK_CLOSE      = ",self%lclose
      write(*,*) "CHK_RMIN       = ",self%rmin
      write(*,*) "CHK_RMAX       = ",self%rmax
      write(*,*) "CHK_EJECT      = ",self%rmaxu
      write(*,*) "CHK_QMIN       = ",self%qmin
      write(*,*) "CHK_QMIN_COORD = ",trim(adjustl(self%qmin_coord))
      write(*,*) "CHK_QMIN_RANGE = ",self%qmin_alo, self%qmin_ahi
      write(*,*) "ENC_OUT        = ",trim(adjustl(self%encounter_file))
      write(*,*) "EXTRA_FORCE    = ",self%lextra_force
      write(*,*) "BIG_DISCARD    = ",self%lbig_discard
      if (self%lenergy) write(*,*) "ENERGY         = ",self%lenergy

      if ((self%MU2KG < 0.0_DP) .or. (self%TU2S < 0.0_DP) .or. (self%DU2M < 0.0_DP)) then
         write(iomsg,*) 'Invalid unit conversion factor'
         iostat = -1
         return
      end if

      ! Calculate the G for the system units
      self%GU = GC / (self%DU2M**3 / (self%MU2KG * self%TU2S**2))

      ! Calculate the inverse speed of light in the system units
      self%inv_c2 = einstinC * self%TU2S / self%DU2M
      self%inv_c2 = (self%inv_c2)**(-2)

      ! The fragmentation model requires the user to set the unit system explicitly.
      if ((integrator == SYMBA) .or. (integrator == RINGMOONS)) then 
         write(*,*) "FRAGMENTATION  = ",self%lfragmentation
         if (.not.mtiny_set) then
            write(iomsg,*) 'SyMBA requres an MTINY value'
            iostat = -1
         end if
      else
         if (self%lfragmentation) then
            write(iomsg,*) 'This integrator does not support fragmentation.'
            iostat = -1
            return
         end if
         if (mtiny_set) then
            write(iomsg,*) 'This integrator does not support MTINY'
            iostat = -1
            return
         end if
      end if

      if ((integrator == SYMBA) .or. (integrator == RINGMOONS) .or. (integrator == RMVS)) then
         if (.not.self%lclose) then
            write(iomsg,*) 'This integrator requires CHK_CLOSE to be enabled'
            iostat = -1
            return
         end if
      end if

      if (mtiny_set) then
         if (self%mtiny < 0.0_DP) then
            write(iomsg,*) "Invalid MTINY: ", self%mtiny
            iostat = -1
            return
         else
            write(*,*) "MTINY          = ", self%mtiny   
         end if
      end if

      ! Determine if the GR flag is set correctly for this integrator
      select case(integrator)
      case(WHM)
         write(*,*) "GR             = ", self%lgr
      case default   
         write(iomsg, *) 'GR is implemented compatible with this integrator'
         iostat = -1
      end select

      iostat = 0

      return 
   end procedure io_config_reader

   module procedure io_get_token
      !! author: David A. Minton
      !!
      !! Retrieves a character token from an input string. Here a token is defined as any set of contiguous non-blank characters not 
      !! beginning with or containing "!". If "!" is present, any remaining part of the buffer including the "!" is ignored
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_get_token.f90
      use swiftest
      implicit none
   
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
   
   end procedure io_get_token

   module procedure io_read_initialize_system
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      use swiftest
      implicit none
   
      call self%cb%initialize(config)
      call self%pl%initialize(config)
      call self%tp%initialize(config)
      call self%set_msys()
      call self%pl%set_vec(config%dt)
      call self%pl%set_vec(self%cb) 
      call self%tp%set_vec(config%dt)
      call self%tp%set_vec(self%cb) 
   
   end procedure io_read_initialize_system

   module procedure io_read_cb_in
      !! author: David A. Minton
      !!
      !! Readsin central body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      use swiftest
      implicit none

      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: iu = LUN
      integer(I4B)            :: i, ierr
      logical                 :: is_ascii 
      real(DP)                :: t

      ierr = 0
      is_ascii = (config%in_type == 'ASCII') 
      if (is_ascii) then
         open(unit = iu, file = config%incbfile, status = 'old', form = 'FORMATTED', iostat = ierr)
         read(iu, *, iostat = ierr) self%mass
         read(iu, *, iostat = ierr) self%radius
         read(iu, *, iostat = ierr) self%j2rp2
         read(iu, *, iostat = ierr) self%j4rp4
         if (config%lrotation) then
            read(iu, *, iostat = ierr) self%Ip(:)
            read(iu, *, iostat = ierr) self%rot(:)
         end if
         if (config%ltides) then
            read(iu, *, iostat = ierr) self%k2
            read(iu, *, iostat = ierr) self%Q
         end if
            
      else
         open(unit = iu, file = config%incbfile, status = 'old', form = 'UNFORMATTED', iostat = ierr)
         call self%read_frame(iu, config, XV, t, ierr)
      end if
      close(iu)
      if (ierr /=  0) then
         write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(config%incbfile))
         call util_exit(FAILURE)
      end if

      self%Gmass = config%GU * self%mass
      return
   end procedure io_read_cb_in

   module procedure io_read_body_in
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in either test particle or massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90 and swiftest_init_tp.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f and swiftest_init_tp.f
      use swiftest
      implicit none


      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: iu = LUN
      integer(I4B)            :: i, ierr, nbody
      logical                 :: is_ascii, is_pl
      character(len=:), allocatable :: infile
      real(DP)               :: t


      ! Select the appropriate polymorphic class (test particle or massive body)
      select type(self)
      class is (swiftest_pl)
         infile = config%inplfile
         is_pl = .true.
      class is (swiftest_tp)
         infile = config%intpfile
         is_pl = .false.
      end select

      ierr = 0
      is_ascii = (config%in_type == 'ASCII') 
      select case(config%in_type)
      case(ASCII_TYPE)
         open(unit = iu, file = infile, status = 'old', form = 'FORMATTED', iostat = ierr)
         read(iu, *, iostat = ierr) nbody
         call self%setup(nbody)
         if (nbody > 0) then
            do i = 1, nbody
               select type(self)
               class is (swiftest_pl)
                  read(iu, *, iostat = ierr) self%name(i), self%mass(i)
                  if (config%lclose) then
                     read(iu, *, iostat = ierr) self%radius(i)
                     if (ierr /= 0 ) exit
                  else
                     self%radius(i) = 0.0_DP
                  end if
                  if (config%lrotation) then
                     read(iu, iostat = ierr) self%Ip(:,i)
                     read(iu, iostat = ierr) self%rot(:,i)
                  end if
                  if (config%ltides) then
                     read(iu, iostat = ierr) self%k2(i)
                     read(iu, iostat = ierr) self%Q(i)
                  end if
               class is (swiftest_tp)
                  read(iu, *, iostat = ierr) self%name(i)
               end select
               if (ierr /= 0 ) exit
               read(iu, *, iostat = ierr) self%xh(:,i)
               read(iu, *, iostat = ierr) self%vh(:,i)
               if (ierr /= 0 ) exit
               self%status(i) = ACTIVE
            end do
         end if
      case (REAL4_TYPE, REAL8_TYPE)  !, SWIFTER_REAL4_TYPE, SWIFTER_REAL8_TYPE)
         open(unit = iu, file = infile, status = 'old', form = 'UNFORMATTED', iostat = ierr)
         read(iu, iostat = ierr) nbody
         call self%setup(nbody)
         if (nbody > 0) then
            call self%read_frame(iu, config, XV, t, ierr)
            self%status(:) = ACTIVE
         end if
      case default
         write(*,*) trim(adjustl(config%in_type)) // ' is an unrecognized file type'
         ierr = -1
      end select
      close(iu)
      if (ierr /= 0 ) then
         write(*,*) 'Error reading in massive body initial conditions from ',trim(adjustl(infile))
         call util_exit(FAILURE)
      end if

      select type(self)
      class is (swiftest_pl)
         self%Gmass(:) = config%GU * self%mass(:)
      end select

      if (config%lgr) then
         select type(self)
         class is (whm_pl)
            call self%gr_vh2pv(config)
         class is (whm_tp)
            call self%gr_vh2pv(config)
         end select
      end if

      return
   end procedure io_read_body_in

   module procedure io_read_hdr
      !! author: David A. Minton
      !!
      !! Read frame header from input binary files
      !!     Function returns read error status (0 = OK, nonzero = ERROR)
      !! Adapted from David E. Kaufmann's Swifter routine: io_read_hdr.f90
      !! Adapted from Hal Levison's Swift routine io_read_hdr.f
      use swiftest
      implicit none
      integer(I4B)         :: ierr
      real(SP)             :: ttmp

      select case (out_type)
      case (REAL4_TYPE, SWIFTER_REAL4_TYPE)
         read(iu, iostat = ierr) ttmp, npl, ntp, out_form
         io_read_hdr = ierr
         if (ierr /= 0) return
         t = ttmp
      case (REAL8_TYPE, SWIFTER_REAL8_TYPE)
         read(iu, iostat = ierr) t
         read(iu, iostat = ierr) npl
         read(iu, iostat = ierr) ntp
         read(iu, iostat = ierr) out_form
         io_read_hdr = ierr
      case default
         write(*,*) trim(adjustl(out_type)) // ' is an unrecognized file type'
         io_read_hdr = -1
      end select

      return
   end procedure io_read_hdr

   module procedure io_read_frame_system
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read a frame (header plus records for each massive body and active test particle) from a output binary file
      use swiftest
      implicit none

      logical, save             :: lfirst = .true.
      integer(I4B)              :: i, j
      real(DP),dimension(:),allocatable :: a, e, inc, capom, omega, capm
      real(DP), dimension(NDIM) :: xtmp, vtmp

      iu = BINUNIT
      if (lfirst) then
         open(unit = iu, file = config%outfile, status = 'OLD', form = 'UNFORMATTED', iostat = ierr)
         lfirst = .false.
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Binary output file already exists or cannot be accessed"
            return
         end if
      end if
      ierr =  io_read_hdr(iu, t, self%pl%nbody, self%tp%nbody, config%out_form, config%out_type)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Binary output file already exists or cannot be accessed"
         return
      end if
      call self%cb%read_frame(iu, config, form, t, ierr)
      if (ierr /= 0) return
      call self%pl%read_frame(iu, config, form, t, ierr)
      if (ierr /= 0) return
      call self%tp%read_frame(iu, config, form, t, ierr)
      return
   end procedure io_read_frame_system

   module procedure io_read_frame_cb
      !! author: David A. Minton
      !!
      !! Reads a frame of output of central body data to the binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_read_frame.f90
      !! Adapted from Hal Levison's Swift routine io_read_frame.F
      use swiftest
      implicit none

      read(iu, iostat = ierr) self%mass
      read(iu, iostat = ierr) self%radius
      read(iu, iostat = ierr) self%j2rp2 
      read(iu, iostat = ierr) self%j4rp4 
      if (config%lrotation) then
         read(iu, iostat = ierr) self%Ip(:)
         read(iu, iostat = ierr) self%rot(:)
      end if
      if (config%ltides) then
         read(iu, iostat = ierr) self%k2
         read(iu, iostat = ierr) self%Q
      end if
      if (ierr /=0) then
         write(*,*) 'Error reading central body data'
         call util_exit(FAILURE)
      end if

      return
   end procedure io_read_frame_cb

   module procedure io_read_frame_body
      !! author: David A. Minton
      !!
      !! Reads a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_read_frame.f90
      !! Adapted from Hal Levison's Swift routine io_read_frame.F
      use swiftest
      implicit none

      associate(n => self%nbody)
         select case (form)
         case (EL) 
            read(iu, iostat = ierr) self%a(1:n)
            read(iu, iostat = ierr) self%e(1:n)
            read(iu, iostat = ierr) self%inc(1:n)
            read(iu, iostat = ierr) self%capom(:)
            read(iu, iostat = ierr) self%omega(:)
            read(iu, iostat = ierr) self%capm(:)
         case (XV)
            read(iu, iostat = ierr) self%xh(1:n, 1)
            read(iu, iostat = ierr) self%xh(1:n, 2)
            read(iu, iostat = ierr) self%xh(1:n, 3)
            read(iu, iostat = ierr) self%vh(1:n, 1)
            read(iu, iostat = ierr) self%vh(1:n, 2)
            read(iu, iostat = ierr) self%vh(1:n, 3)
         end select
         select type(self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            read(iu, iostat = ierr) self%mass(1:n)
            read(iu, iostat = ierr) self%radius(1:n)
            if (config%lrotation) then
               read(iu, iostat = ierr) self%Ip(1:n, 1)
               read(iu, iostat = ierr) self%Ip(1:n, 2)
               read(iu, iostat = ierr) self%Ip(1:n, 3)
               read(iu, iostat = ierr) self%rot(1:n, 1)
               read(iu, iostat = ierr) self%rot(1:n, 2)
               read(iu, iostat = ierr) self%rot(1:n, 3)
            end if
            if (config%ltides) then
               read(iu, iostat = ierr) self%k2(1:n)
               read(iu, iostat = ierr) self%Q(1:n)
            end if
         end select
      end associate

      if (ierr /=0) then
         write(*,*) 'Error reading Swiftest body data'
         call util_exit(FAILURE)
      end if

      return
   end procedure io_read_frame_body


   module procedure io_read_encounter
      !! author: David A. Minton
      !!
      !! Read close encounter data from input binary files
      !!     Other than time t, there is no direct file input from this function
      !!     Function returns read error status (0 = OK, nonzero = ERROR)
      !! Adapted from David E. Kaufmann's Swifter routine: io_read_encounter.f90
      use swiftest
      implicit none
      logical         :: lxdr
      logical , save    :: lfirst = .true.
      integer(I4B), parameter :: lun = 30
      integer(I4B)        :: ierr
      integer(I4B), save    :: iu = lun

      if (lfirst) then
         open(unit = iu, file = encounter_file, status = 'OLD', form = 'UNFORMATTED', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest Error:"
            write(*, *) "   unable to open binary encounter file"
            call util_exit(FAILURE)
         end if
         lfirst = .false.
      end if
      read(iu, iostat = ierr) t
      io_read_encounter = ierr
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
   end procedure io_read_encounter


end submodule io_read
