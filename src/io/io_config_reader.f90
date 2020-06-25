submodule (swiftest_classes) s_io_config_reader
contains
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

      logical                 :: t0_set = .false.        !! Is the initial time set in the input file?
      logical                 :: tstop_set = .false.     !! Is the final time set in the input file?
      logical                 :: dt_set = .false.        !! Is the step size set in the input file?
      integer(I4B)            :: ilength, ifirst, ilast  !! Variables used to parse input file
      character(STRMAX)       :: line                    !! Line of the input file
      character (len=:), allocatable :: line_trim,config_name, config_value
      character(*),parameter :: linefmt = '(A)'

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
            case ("J2")
               read(config_value, *) config%j2rp2
            case ("J4")
               read(config_value, *) config%j4rp4
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
            case ("RHILL_PRESENT")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == "T") config%lrhill_present = .true.
            ! Added by the Purdue Swiftest development group (Minton, Wishard, Populin, and Elliott)
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
            case ("ENERGY")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lenergy = .true.
            ! The following are not yet implemented
            case ("RINGMOONS")
               call util_toupper(config_value)
               if (config_value == "YES" .or. config_value == 'T') config%lringmoons = .true.
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
      if ((config%j2rp2 == 0.0_DP) .and. (config%j4rp4 /= 0.0_DP)) then
         write(iomsg,*) 'Cannot have j4 without j2'
         return
         iostat = -1
      end if
      if (config%qmin > 0.0_DP) then
         if ((config%qmin_coord /= "HELIO") .and. (config%qmin_coord /= "BARY")) then
            write(iomsg,*) 'Invalid qmin_coord: ',trim(adjustl(config%qmin_coord))
            return
            iostat = -1
         end if
         if ((config%qmin_alo <= 0.0_DP) .or. (config%qmin_ahi <= 0.0_DP)) then
            write(iomsg,*) 'Invalid qmin vals'
            return
            iostat = -1
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
      write(*,*) "J2             = ",config%j2rp2
      write(*,*) "J4             = ",config%j4rp4
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
      write(*,*) "RHILL_PRESENT  = ",config%lrhill_present

      ! Added by D. Minton
      MU2KG = config%MU2KG
      TU2S  = config%TU2S 
      DU2M = config%DU2M
      ! The fragmentation model requires the user to set the unit system explicitly.
      write(*,*) "FRAGMENTATION  = ",config%lfragmentation
      if (config%lfragmentation) then
         write(*,*) "MU2KG          = ",MU2KG
         write(*,*) "TU2S           = ",TU2S 
         write(*,*) "DU2M          = ",DU2M
         if ((MU2KG < 0.0_DP) .or. (TU2S < 0.0_DP) .or. (DU2M < 0.0_DP)) then
            write(*,*) 'Invalid unit conversion factor'
            write(*,*) 'MU2KG: ',MU2KG
            write(*,*) 'TU2S: ',TU2S
            write(*,*) 'DU2M: ',DU2M
            iostat = -1
         end if
      end if 
      !Added mtiny to the argument list rather than from the terminal
      if (config%mtiny < 0.0_DP) then
         write(*,*) "Invalid MTINY: ",config%mtiny
         iostat = -1
      else
         write(*,*) "MTINY          = ",config%mtiny   
      end if
      if (config%lenergy) write(*,*) "ENERGY         = ",config%lenergy
      if (config%lringmoons) write(*,*) "RINGMOONS      = ",config%lringmoons

      return 

   end procedure io_config_reader

end submodule s_io_config_reader
