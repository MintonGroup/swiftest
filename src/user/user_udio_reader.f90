submodule (user) s_user_udio_reader
contains
   module procedure user_udio_reader
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Read in parameters for the integration
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
   !! Adapted from Martin Duncan's Swift routine io_init_param.f
   !$ use omp_lib
   !use util, only: util_exit ! IMPLEMENTATION TBD
   implicit none

   logical                 :: t0_set = .false.        !! Is the initial time set in the input file?
   logical                 :: tstop_set = .false.     !! Is the final time set in the input file?
   logical                 :: dt_set = .false.        !! Is the step size set in the input file?
   integer(I4B)            :: ilength, ifirst, ilast  !! Variables used to parse input file
   character(STRMAX)       :: line                    !! Line of the input file
   character (len=:), allocatable :: line_trim,param_name, param_value

   100 format(A)
   ! Parse the file line by line, extracting tokens then matching them up with known parameters if possible
   do
      read(unit = unit, fmt = 100, iostat = iostat, end = 1) line
      line_trim = trim(adjustl(line))
      ilength = len(line_trim)
      if ((ilength /= 0)) then 
         ifirst = 1
         ! Read the pair of tokens. The first one is the parameter name, the second is the value.
         param_name = user_get_token(line_trim, ifirst, ilast, iostat)
         if (param_name == '') cycle ! No parameter name (usually because this line is commented out)
         call util_toupper(param_name)
         ifirst = ilast + 1
         param_value = user_get_token(line_trim, ifirst, ilast, iostat)
         select case (param_name)
         case ("NPLMAX")
            read(param_value, *) param%nplmax
         case ("NTPMAX")
            read(param_value, *) param%ntpmax
         case ("T0")
            read(param_value, *) param%t0
            t0_set = .true.
         case ("TSTOP")
            read(param_value, *) param%tstop
            tstop_set = .true.
         case ("DT")
            read(param_value, *) param%dt
            dt_set = .true.
         case ("PL_IN")
            param%inplfile = param_value
         case ("TP_IN")
            param%intpfile = param_value
         case ("IN_TYPE")
            call util_toupper(param_value)
            param%in_type = param_value
         case ("ISTEP_OUT")
            read(param_value, *) param%istep_out
         case ("BIN_OUT")
            param%outfile = param_value
         case ("OUT_TYPE")
            call util_toupper(param_value)
            param%out_type = param_value
         case ("OUT_FORM")
            call util_toupper(param_value)
            param%out_form = param_value
         case ("OUT_STAT")
            call util_toupper(param_value)
            param%out_stat = param_value
         case ("ISTEP_DUMP")
            read(param_value, *) param%istep_dump
         case ("J2")
            read(param_value, *) param%j2rp2
         case ("J4")
            read(param_value, *) param%j4rp4
         case ("CHK_CLOSE")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lclose = .true.
         case ("CHK_RMIN")
            read(param_value, *) param%rmin
         case ("CHK_RMAX")
            read(param_value, *) param%rmax
         case ("CHK_EJECT")
            read(param_value, *) param%rmaxu
         case ("CHK_QMIN")
            read(param_value, *) param%qmin
         case ("CHK_QMIN_COORD")
            call util_toupper(param_value)
            param%qmin_coord = param_value
         case ("CHK_QMIN_RANGE")
            read(param_value, *) param%qmin_alo
            ifirst = ilast + 1
            param_value = user_get_token(line, ifirst, ilast, iostat)
            read(param_value, *) param%qmin_ahi
         case ("ENC_OUT")
            param%encounter_file = param_value
         case ("EXTRA_FORCE")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lextra_force = .true.
         case ("BIG_DISCARD")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T' ) param%feature%lbig_discard = .true.
         case ("RHILL_PRESENT")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == "T") param%feature%lrhill_present = .true.
         ! Added by the Purdue Swiftest development group (Minton, Wishard, Populin, and Elliott)
         case ("FRAGMENTATION")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == "T") param%feature%lfragmentation = .true.
         case ("MU2KG")
            read(param_value, *) param%MU2KG
         case ("TU2S")
            read(param_value, *) param%TU2S
         case ("DU2M")
            read(param_value, *) param%DU2M
         case ("MTINY")
            read(param_value, *) param%mtiny
         case ("ENERGY")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lenergy = .true.

         ! The following are not yet implemented
         case ("RINGMOONS")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lringmoons = .true.
         case ("RING_OUTFILE")
            param%ring_outfile = param_value
         case ("ROTATION")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lrotation = .true. 
         case ("TIDES")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%ltides = .true. 
         case ("GR")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lgr = .true. 
         case ("YARKOVSKY")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lyarkovsky = .true. 
         case ("YORP")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lyorp = .true. 
         case default
            write(iomsg,*) "Unknown parameter -> ",param_name
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
   if (param%dt <= 0.0_DP) then
      write(iomsg,*) 'Invalid timestep: '
      iostat = -1
      return
   end if
   if (param%inplfile == "") then
      write(iomsg,*) 'No valid planet file in input file'
      iostat = -1
      return
   end if
   if ((param%in_type /= SWIFTER_REAL8_TYPE) .and. (param%in_type /= "ASCII")) then
      write(iomsg,*) 'Invalid input file type:',param%in_type
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
         write(iomsg,*) 'Invalid out_type: ',param%out_type
         iostat = -1
         return
      end if
      if ((param%out_form /= "EL") .and. (param%out_form /= "XV")) then
         write(iomsg,*) 'Invalid out_form: ',param%out_form
         iostat = -1
         return
      end if
      if ((param%out_stat /= "NEW") .and. (param%out_stat /= "REPLACE") .and. (param%out_stat /= "APPEND")) then
         write(iomsg,*) 'Invalid out_stat: ',param%out_stat
         iostat = -1
         return
      end if
   end if
   if ((param%j2rp2 == 0.0_DP) .and. (param%j4rp4 /= 0.0_DP)) then
      write(iomsg,*) 'Cannot have j4 without j2'
      return
      iostat = -1
   end if
   if (param%qmin > 0.0_DP) then
      if ((param%qmin_coord /= "HELIO") .and. (param%qmin_coord /= "BARY")) then
         write(iomsg,*) 'Invalid qmin_coord: ',param%qmin_coord
         return
         iostat = -1
      end if
      if ((param%qmin_alo <= 0.0_DP) .or. (param%qmin_ahi <= 0.0_DP)) then
         write(iomsg,*) 'Invalid qmin vals'
         return
         iostat = -1
      end if
   end if

   return 

   end procedure user_udio_reader

end submodule s_user_udio_reader
