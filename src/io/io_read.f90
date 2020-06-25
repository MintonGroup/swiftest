submodule (swiftest_classes) io_read
   !! author: David A. Minton
   !! 
   !! This submodule contains implementations of the following procedures:
   !!    io_read_config_in
   !!    io_config_reader
   !!    io_read_cb_in
   !!    io_read_pl_in
   !!    io_read_tp_in
   !!    io_read_line_swifter
contains
   module procedure io_read_config_in
      !! author: David A. Minton
      !!
      !! Read in parameters for the integration
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_config.f90
      !! Adapted from Martin Duncan's Swift routine io_init_config.f
      
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
      call config_reader(LUN, iotype= "none", v_list = [integrator], iostat = ierr, iomsg = error_message)
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

      integrator = v_list(0)
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
               read(config_value, *) config%nplmax
            case ("NTPMAX")
               read(config_value, *) config%ntpmax
            case ("T0")
               read(config_value, *) config%t0
               t0_set = .true.
            case ("TSTOP")
               read(config_value, *) config%tstop
               tstop_set = .true.
            case ("DT")
               read(config_value, *) config%dt
            case ("CB_IN")
               config%incbfile = config_value
            case ("PL_IN")
               config%inplfile = config_value
            case ("TP_IN")
               config%intpfile = config_value
            case ("IN_TYPE")
               call util_toupper(config_value)
               config%in_type = config_value
            case ("ISTEP_OUT")
               read(config_value, *) config%istep_out
            case ("BIN_OUT")
               config%outfile = config_value
            case ("OUT_TYPE")
               call util_toupper(config_value)
               config%out_type = config_value
            case ("OUT_FORM")
               call util_toupper(config_value)
               config%out_form = config_value
            case ("OUT_STAT")
               call util_toupper(config_value)
               config%out_stat = config_value
            case ("ISTEP_DUMP")
               read(config_value, *) config%istep_dump
            case ("CHK_CLOSE")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lclose = .true.
            case ("CHK_RMIN")
               read(config_value, *) config%rmin
            case ("CHK_RMAX")
               read(config_value, *) config%rmax
            case ("CHK_EJECT")
               read(config_value, *) config%rmaxu
            case ("CHK_QMIN")
               read(config_value, *) config%qmin
            case ("CHK_QMIN_COORD")
               call util_toupper(config_value)
               config%qmin_coord = config_value
            case ("CHK_QMIN_RANGE")
               read(config_value, *) config%qmin_alo
               ifirst = ilast + 1
               config_value = io_get_token(line, ifirst, ilast, iostat)
               read(config_value, *) config%qmin_ahi
            case ("ENC_OUT")
               config%encounter_file = config_value
            case ("EXTRA_FORCE")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lextra_force = .true.
            case ("BIG_DISCARD")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T' ) config%lbig_discard = .true.
            case ("FRAGMENTATION")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == "T") config%lfragmentation = .true.
            case ("MU2KG")
               read(config_value, *) config%MU2KG
            case ("TU2S")
               read(config_value, *) config%TU2S
            case ("DU2M")
               read(config_value, *) config%DU2M
            case ("MTINY")
               read(config_value, *) config%mtiny
               mtiny_set = .true.
            case ("ENERGY")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lenergy = .true.
            case ("RING_OUTFILE")
               config%ring_outfile = config_value
            case ("ROTATION")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lrotation = .true. 
            case ("TIDES")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%ltides = .true. 
            case ("GR")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lgr = .true. 
            case ("YARKOVSKY")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lyarkovsky = .true. 
            case ("YORP")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lyorp = .true. 
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
      if (config%dt <= 0.0_DP) then
         write(iomsg,*) 'Invalid timestep: '
         iostat = -1
         return
      end if
      if (config%inplfile == "") then
         write(iomsg,*) 'No valid massive body file in input file'
         iostat = -1
         return
      end if
      if ((config%in_type /= REAL8_TYPE) .and. (config%in_type /= "ASCII")) then
         write(iomsg,*) 'Invalid input file type:',trim(adjustl(config%in_type))
         iostat = -1
         return
      end if
      if ((config%istep_out <= 0) .and. (config%istep_dump <= 0)) then
         write(iomsg,*) 'Invalid istep'
         iostat = -1
         return
      end if
      if ((config%istep_out > 0) .and. (config%outfile == "")) then
         write(iomsg,*) 'Invalid outfile'
         iostat = -1
         return
      end if
      if (config%outfile /= "") then
         if ((config%out_type /= REAL4_TYPE) .and. (config%out_type /= REAL8_TYPE) .and. &
               (config%out_type /= SWIFTER_REAL4_TYPE)  .and. (config%out_type /= SWIFTER_REAL8_TYPE)) then
            write(iomsg,*) 'Invalid out_type: ',trim(adjustl(config%out_type))
            iostat = -1
            return
         end if
         if ((config%out_form /= "EL") .and. (config%out_form /= "XV")) then
            write(iomsg,*) 'Invalid out_form: ',trim(adjustl(config%out_form))
            iostat = -1
            return
         end if
         if ((config%out_stat /= "NEW") .and. (config%out_stat /= "REPLACE") .and. (config%out_stat /= "APPEND")) then
            write(iomsg,*) 'Invalid out_stat: ',trim(adjustl(config%out_stat))
            iostat = -1
            return
         end if
      end if
      if (config%qmin > 0.0_DP) then
         if ((config%qmin_coord /= "HELIO") .and. (config%qmin_coord /= "BARY")) then
            write(iomsg,*) 'Invalid qmin_coord: ',trim(adjustl(config%qmin_coord))
            iostat = -1
            return
         end if
         if ((config%qmin_alo <= 0.0_DP) .or. (config%qmin_ahi <= 0.0_DP)) then
            write(iomsg,*) 'Invalid qmin vals'
            iostat = -1
            return
         end if
      end if

      write(*,*) "NPLMAX         = ",config%nplmax
      write(*,*) "NTPMAX         = ",config%ntpmax
      write(*,*) "T0             = ",config%t0
      write(*,*) "TSTOP          = ",config%tstop
      write(*,*) "DT             = ",config%dt
      write(*,*) "PL_IN          = ",trim(adjustl(config%inplfile))
      write(*,*) "TP_IN          = ",trim(adjustl(config%intpfile))
      write(*,*) "IN_TYPE        = ",trim(adjustl(config%in_type))
      write(*,*) "ISTEP_OUT      = ",config%istep_out
      write(*,*) "BIN_OUT        = ",trim(adjustl(config%outfile))
      write(*,*) "OUT_TYPE       = ",trim(adjustl(config%out_type))
      write(*,*) "OUT_FORM       = ",trim(adjustl(config%out_form))
      write(*,*) "OUT_STAT       = ",trim(adjustl(config%out_stat))
      write(*,*) "ISTEP_DUMP     = ",config%istep_dump
      write(*,*) "CHK_CLOSE      = ",config%lclose
      write(*,*) "CHK_RMIN       = ",config%rmin
      write(*,*) "CHK_RMAX       = ",config%rmax
      write(*,*) "CHK_EJECT      = ",config%rmaxu
      write(*,*) "CHK_QMIN       = ",config%qmin
      write(*,*) "CHK_QMIN_COORD = ",trim(adjustl(config%qmin_coord))
      write(*,*) "CHK_QMIN_RANGE = ",config%qmin_alo, config%qmin_ahi
      write(*,*) "ENC_OUT        = ",trim(adjustl(config%encounter_file))
      write(*,*) "EXTRA_FORCE    = ",config%lextra_force
      write(*,*) "BIG_DISCARD    = ",config%lbig_discard
      if (config%lenergy) write(*,*) "ENERGY         = ",config%lenergy

      if ((MU2KG < 0.0_DP) .or. (TU2S < 0.0_DP) .or. (DU2M < 0.0_DP)) then
         write(iomsg,*) 'Invalid unit conversion factor'
         iostat = -1
         return
      end if

      ! The fragmentation model requires the user to set the unit system explicitly.
      if ((integrator == SYMBA) .or. (integrator == RINGMOONS)) then 
         write(*,*) "FRAGMENTATION  = ",config%lfragmentation
         if (.not.mtiny_set) then
            write(iomsg,*) 'SyMBA requres an MTINY value'
            iostat = -1
         end if
      else
         if (config%lfragmentation) then
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
         if (.not.config%lclose) then
            write(iomsg,*) 'This integrator requires CHK_CLOSE to be enabled'
            iostat = -1
            return
         end if
      end if

      if (mtiny_set) then
         if (config%mtiny < 0.0_DP) then
            write(iomsg,*) "Invalid MTINY: ",config%mtiny
            iostat = -1
            return
         else
            write(*,*) "MTINY          = ",config%mtiny   
         end if
      end if

      iostat = 0

      return 
   end procedure io_config_reader

   module procedure io_read_cb_in
      !! author: David A. Minton
      !!
      !! Read in central body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      implicit none

      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: i, iu, ierr, npl
      logical                 :: is_ascii 

      ierr = 0
      is_ascii = (config%in_type == 'ASCII') 
      if (is_ascii) then
         open(unit = LUN, file = config%incbfile, status = 'old', form = 'formatted', iostat = ierr)
         read(LUN, *, iostat = ierr) self%mass
         read(LUN, *, iostat = ierr) self%radius
         read(LUN, *, iostat = ierr) self%j2rp2
         read(LUN, *, iostat = ierr) self%j4rp4
         if (config%lrotation) then
            read(LUN, *, iostat = ierr) self%Ip(:)
            read(LUN, *, iostat = ierr) self%rot(:)
         end if
         if (config%tides) then
            read(LUN, *, iostat = ierr) self%k2
            read(LUN, *, iostat = ierr) self%Q
         end if
            
      else
         open(unit = LUN, file = config%incbfile, status = 'old', form = 'unformatted', iostat = ierr)
      end if
      close(LUN)
      if (ierr /=  0) then
         write(*,*) 'Error opening massive body initial conditions file ',trim(adjustl(config%inplfile))
         call util_exit(FAILURE)
      end if
      return
   end procedure io_read_cb_in

   module procedure io_read_pl_in
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      implicit none

      integer(I4B), parameter :: LUN = 7              !! Unit number of input file
      integer(I4B)            :: i, iu, ierr, npl
      logical                 :: is_ascii 

      ierr = 0
      is_ascii = (config%in_type == 'ASCII') 
      select case(config%in_type)
      case(ASCII_TYPE)
         open(unit = LUN, file = config%inplfile, status = 'old', form = 'formatted', iostat = ierr)
         read(LUN, *, iostat = ierr) npl
         call self%alloc(npl)
         if (npl > 0) then
            do i = 1, npl
               read(LUN, *, iostat = ierr) self%name(i), self%mass(i)
               if (ierr /= 0 ) exit
               if (config%lclose) then
                  read(LUN, *, iostat = ierr) self%radius(i)
                  if (ierr /= 0 ) exit
               else
                  self%radius(i) = 0.0_DP
               end if
               read(LUN, *, iostat = ierr) self%xh(:,i)
               read(LUN, *, iostat = ierr) self%vh(:,i)
               if (ierr /= 0 ) exit
               self%status(i) = ACTIVE
            end do
         end if
      case (REAL4_TYPE, REAL8_TYPE, SWIFTER_REAL4_TYPE, SWIFTER_REAL8_TYPE)
         open(unit = LUN, file = config%inplfile, status = 'old', form = 'unformatted', iostat = ierr)
         read(LUN, iostat = ierr) npl
         call self%alloc(npl)
         if (npl > 0) then
            read(LUN, iostat = ierr) self%name(1:npl)
            read(LUN, iostat = ierr) self%mass(1:npl)
            if (config%lclose) then
               read(LUN, iostat = ierr) self%radius(1:npl)
            else
               self%radius(:) = 0.0_DP
            end if
            read(LUN, iostat = ierr) self%xh(1,1:npl)
            read(LUN, iostat = ierr) self%xh(2,1:npl)
            read(LUN, iostat = ierr) self%xh(3,1:npl)
            read(LUN, iostat = ierr) self%vh(1,1:npl)
            read(LUN, iostat = ierr) self%vh(2,1:npl)
            read(LUN, iostat = ierr) self%vh(3,1:npl)
            if (ierr /= 0 ) exit
            self%status(:) = ACTIVE
         end if
      case default
         write(*,*) trim(adjustl(config%in_type)) // ' is an unrecognized file type'
         ierr = -1
      end select
      close(LUN)
      if (ierr /= 0 ) then
         write(*,*) 'Error reading in massive body initial conditions from ',trim(adjustl(config%inplfile))
         call util_exit(FAILURE)
      end if

      return
   end procedure io_read_pl_in

   module procedure io_read_tp_in
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in test particle data
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine io_init_tp.f
      implicit none
   
      integer(I4B), parameter  :: LUN = 7              !! Unit number of input file
      integer(I4B)             :: i, iu, ierr, ntp
   
      ierr = 0
      select case(config%in_type)
      case(ASCII_TYPE)
         open(unit = LUN, file = config%intpfile, status = 'old', form = 'formatted', iostat = ierr)
         read(LUN, *, iostat = ioerr) ntp
         call self%alloc(ntp)
         if (ntp > 0) then 
            do i = 1, self%nbody
               read(LUN, *) self%name(i)
               read(LUN, *) self%xh(:,i)
               read(LUN, *) self%vh(:,i)
               self%status(i) = ACTIVE
            end do
         end if
      case (REAL4_TYPE, REAL8_TYPE, SWIFTER_REAL4_TYPE, SWIFTER_REAL8_TYPE)
         open(unit = LUN, file = config%intpfile, status = 'old', form = 'unformatted', iostat = ierr)
         read(LUN, iostat = ioerr) ntp
         call self%alloc(ntp)
         if (ntp > 0) then 
            read(LUN) self%name(1:ntp)
            read(LUN) self%xh(1,1:ntp)
            read(LUN) self%xh(2,1:ntp)
            read(LUN) self%xh(3,1:ntp)
            read(LUN) self%vh(1,1:ntp)
            read(LUN) self%vh(2,1:ntp)
            read(LUN) self%vh(3,1:ntp)
            self%status(:) = ACTIVE
         end if
      case default
         write(*,*) trim(adjustl(config%in_type)) // ' is an unrecognized file type'
         ierr = -1
      end select
      close(LUN)
      if (ierr /=  0) then
         write(*,*) 'Error opening test particle initial conditions file ',trim(adjustl(config%intpfile))
         call util_exit(FAILURE)
      end if
   
      return
   end procedure io_read_tp_in

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

      select case (config%out_type)
      case (REAL4_TYPE, SWIFTER_REAL4_TYPE)
         read(iu, iostat = ierr) ttmp, npl, ntp, iout_form
         io_read_hdr = ierr
         if (ierr /= 0) return
         t = ttmp
      case (REAL8_TYPE, SWIFTER_REAL8_TYPE)
         read(iu, iostat = ierr) t
         read(iu, iostat = ierr) npl
         read(iu, iostat = ierr) ntp
         read(iu, iostat = ierr) iout_form
         io_read_hdr = ierr
      case default
         write(*,*) trim(adjustl(config%in_type)) // ' is an unrecognized file type'
         io_read_hder = -1
      end select

      return
   end procedure io_read_hdr

   module procedure io_read_line
      !! author: David A. Minton
      !!
      !! Read a line (record) from input binary files
      !!     Function returns read error status (0 = OK, nonzero = ERROR)
      !! Adapted from David E. Kaufmann's Swifter routine: io_read_line.f90
      !! Adapted from Hal Levison's Swift routine io_read_line.f
      use swiftest
      implicit none
      logical                    :: lmass, lradius
      integer( I4B)              :: ierr
      real(SP)                   :: smass, sradius
      real(SP), dimension(NDIM2) :: svec
      real(DP), dimension(NDIM2) :: dvec

      lmass = present(mass)
      if (lmass) then
         lradius = present(radius)
         if (.not. lradius) then
            write(*, *) "Swiftest Error:"
            write(*, *) "   function io_read_line called with optional mass but without optional radius"
            call util_exit(FAILURE)
         end if
      end if
      select case (config%out_type)
         case (REAL4_TYPE, SWIFTER_REAL4_TYPE)
            if (lmass) then
               read(iu, iostat = ierr) name, smass, sradius, svec
            else
               read(iu, iostat = ierr) name, svec
            end if
            io_read_line = ierr
            if (ierr /= 0) return
            if (lmass) mass = smass
            d1 = svec(1); d2 = svec(2); d3 = svec(3); d4 = svec(4); d5 = svec(5); d6 = svec(6)
         case (REAL8_TYPE, SWIFTER_REAL8_TYPE)
            if (lmass) then
               read(iu, iostat = ierr) name, mass, radius, dvec
            else
               read(iu, iostat = ierr) name, dvec
            end if
            io_read_line = ierr
            if (ierr /= 0) return
            d1 = dvec(1); d2 = dvec(2); d3 = dvec(3); d4 = dvec(4); d5 = dvec(5); d6 = dvec(6)
      end select

      return
   end procedure io_read_line

   module procedure io_read_line_swifter
      !! author: David A. Minton
      !!
      !! Read a line (record) from input binary files using the old Swifter style
      !!     Function returns read error status (0 = OK, nonzero = ERROR)
      !! Adapted from David E. Kaufmann's Swifter routine: io_read_line.f90
      !! Adapted from Hal Levison's Swift routine io_read_line.f
      use swiftest
      implicit none
      logical                    :: lmass, lradius
      integer( I4B)              :: ierr
      real(SP)                   :: smass, sradius
      real(SP), dimension(NDIM2) :: svec
      real(DP), dimension(NDIM2) :: dvec

      lmass = present(MASS)
      if (lmass) then
         lradius = present(RADIUS)
         if (.not. lradius) then
            write(*, *) "Swiftest Error:"
            write(*, *) "   Function io_read_line_swifter called with optional mass but without optional radius"
            call util_exit(FAILURE)
         end if
      end if
      select case (config%out_type)
      case (SWIFTER_REAL4_TYPE)
         if (lmass) then
            read(iu, iostat = ierr) name, smass, sradius, svec
         else
            read(iu, iostat = ierr) name, svec
         end if
         io_read_line = ierr
         if (ierr /= 0) return
         if (lmass) mass = smass
         d1 = svec(1); d2 = svec(2); d3 = svec(3); d4 = svec(4); d5 = svec(5); d6 = svec(6)
      case (SWIFTER_REAL8_TYPE)
         if (lmass) then
            read(iu, iostat = ierr) name, mass, radius, dvec
         else
            read(iu, iostat = ierr) name, dvec
         end if
         io_read_line = ierr
         if (ierr /= 0) return
         d1 = dvec(1); d2 = dvec(2); d3 = dvec(3); d4 = dvec(4); d5 = dvec(5); d6 = dvec(6)
      end select

      return
   end procedure io_read_line_swifter

end submodule io_read
