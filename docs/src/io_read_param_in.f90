submodule (io) s_io_read_param_in
contains
   module procedure io_read_param_in
   !! Reads in the param.in file that sets up the run parameters
   !$ use omp_lib
   implicit none

   ! Internals
   integer(I4B), parameter :: LUN = 7                 !! Unit number of input file
   logical                 :: t0_set = .false.        !! Is the initial time set in the input file?
   logical                 :: tstop_set = .false.     !! Is the final time set in the input file?
   logical                 :: dt_set = .false.        !! Is the step size set in the input file?
   integer(I4B)            :: ierr = 0                !! Input error code
   integer(I4B)            :: ilength, ifirst, ilast  !! Variables used to parse input file
   character(STRMAX)       :: line                    !! Line of the input file
   character(STRMAX)       :: token                   !! Input file token

   ! Executable code

   ! Read in name of parameter file
   write(*, 100, advance = "NO") "Parameter data file is "
   write(*, 100) inparfile
   write(*, *) " "
   100 format(A)
   call io_open(LUN, inparfile, "OLD", "FORMATTED", ierr)
   if (ierr /= 0) then
      write(*, 100, advance = "NO") "Unable to open file "
      write(*, 100) inparfile
      call util_exit(FAILURE)
   end if

   ! Parse the file line by line, extracting tokens then matching them up with known parameters if possible
   do
      read(LUN, 100, IOSTAT = ierr, end = 1) line
      line = adjustl(line)
      ilength = len_trim(line)
      if ((ilength /= 0) .and. (line(1:1) /= "!")) then
         ifirst = 1
         call io_get_token(line, ilength, ifirst, ilast, ierr)
         token = line(ifirst:ilast)
         call util_toupper(token)
         select case (token)
           case ("NPLMAX")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%nplmax
           case ("NTPMAX")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%ntpmax
           case ("T0")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%t0
             t0_set = .true.
           case ("TSTOP")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%tstop
             tstop_set = .true.
           case ("DT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%dt
             dt_set = .true.
           case ("PL_IN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             param%inplfile = token
           case ("TP_IN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             param%intpfile = token
           case ("IN_TYPE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             param%in_type = token
           case ("ISTEP_OUT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%istep_out
           case ("BIN_OUT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             param%outfile = token
           case ("OUT_TYPE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             param%out_type = token
           case ("OUT_FORM")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             param%out_form = token
           case ("OUT_STAT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             param%out_stat = token
           case ("ISTEP_DUMP")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%istep_dump
           case ("J2")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%j2rp2
           case ("J4")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%j4rp4
           case ("CHK_CLOSE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == 'T') param%feature%lclose = .true.
           case ("CHK_RMIN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%rmin
           case ("CHK_RMAX")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%rmax
           case ("CHK_EJECT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%rmaxu
           case ("CHK_QMIN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%qmin
           case ("CHK_QMIN_COORD")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             param%qmin_coord = token
           case ("CHK_QMIN_RANGE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%qmin_alo
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%qmin_ahi
           case ("ENC_OUT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             param%encounter_file = token
           case ("EXTRA_FORCE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == 'T') param%feature%lextra_force = .true.
           case ("BIG_DISCARD")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == 'T' ) param%feature%lbig_discard = .true.
           case ("RHILL_PRESENT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == "T") param%feature%lrhill_present = .true.

           ! Added by the Purdue Swiftest development group (Minton, Wishard, Populin, and Elliott)
           case ("FRAGMENTATION")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == "T") param%feature%lfragmentation = .true.
           case ("MU2GM")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%MU2GM
           case ("TU2S")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%TU2S
           case ("DU2CM")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%DU2CM
           case ("MTINY")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) param%mtiny
           case ("PYTHON")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == 'T') param%feature%lpython = .true.
           case ("ENERGY")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == 'T') param%feature%lenergy = .true.
           case ("RINGMOONS")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == 'T') param%feature%lringmoons = .true.
           case ("RING_OUTFILE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             param%ring_outfile = token
           case ("ROTATION")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES" .or. token == 'T') param%feature%lrotation = .true. 
           case default
             write(*, 100, advance = "NO") "Unknown parameter -> "
             write(*, *) token
             call util_exit(FAILURE)
         end select
      end if
   end do
   1 close(LUN)
   write(*, 100, advance = "NO") "NPLMAX      = "
   write(*, *) param%nplmax
   write(*, 100, advance = "NO") "NTPMAX      = "
   write(*, *) param%ntpmax
   write(*, 100, advance = "NO") "T0       = "
   write(*, *) param%t0
   write(*, 100, advance = "NO") "TSTOP     = "
   write(*, *) param%tstop
   write(*, 100, advance = "NO") "DT       = "
   write(*, *) param%dt
   write(*, 100, advance = "NO") "PL_IN     = "
   write(*, *) trim(adjustl(param%inplfile))
   write(*, 100, advance = "NO") "TP_IN     = "
   write(*, *) trim(adjustl(param%intpfile))
   write(*, 100, advance = "NO") "IN_TYPE     = "
   write(*, *) trim(adjustl(param%in_type))
   write(*, 100, advance = "NO") "ISTEP_OUT   = "
   write(*, *) param%istep_out
   write(*, 100, advance = "NO") "BIN_OUT     = "
   write(*, *) trim(adjustl(param%outfile))
   write(*, 100, advance = "NO") "OUT_TYPE    = "
   write(*, *) trim(adjustl(param%out_type))
   write(*, 100, advance = "NO") "OUT_FORM    = "
   write(*, *) trim(adjustl(param%out_form))
   write(*, 100, advance = "NO") "OUT_STAT    = "
   write(*, *) trim(adjustl(param%out_stat))
   write(*, 100, advance = "NO") "ISTEP_DUMP   = "
   write(*, *) param%istep_dump
   write(*, 100, advance = "NO") "J2       = "
   write(*, *) param%j2rp2
   write(*, 100, advance = "NO") "J4       = "
   write(*, *) param%j4rp4
   write(*, 100, advance = "NO") "CHK_CLOSE   = "
   write(*, *) param%feature%lclose
   write(*, 100, advance = "NO") "CHK_RMIN    = "
   write(*, *) param%rmin
   write(*, 100, advance = "NO") "CHK_RMAX    = "
   write(*, *) param%rmax
   write(*, 100, advance = "NO") "CHK_EJECT   = "
   write(*, *) param%rmaxu
   write(*, 100, advance = "NO") "CHK_QMIN    = "
   write(*, *) param%qmin
   write(*, 100, advance = "NO") "CHK_QMIN_COORD = "
   write(*, *) trim(adjustl(param%qmin_coord))
   write(*, 100, advance = "NO") "CHK_QMIN_RANGE = "
   write(*, *) param%qmin_alo, param%qmin_ahi
   write(*, 100, advance = "NO") "ENC_OUT     = "
   write(*, *) trim(adjustl(param%encounter_file))
   write(*, 100, advance = "NO") "EXTRA_FORCE   = "
   write(*, *) param%feature%lextra_force
   write(*, 100, advance = "NO") "BIG_DISCARD   = "
   write(*, *) param%feature%lbig_discard
   write(*, 100, advance = "NO") "RHILL_PRESENT  = "
   write(*, *) param%feature%lrhill_present
   write(*, *) " "
   ierr = 0
   if ((.not. t0_set) .or. (.not. tstop_set) .or. (.not. dt_set)) then
      write(*,*) 'Valid simulation time not set'
      write(*,*) 't0_set:   ',t0_set
      write(*,*) 'tstop_set:',tstop_set
      write(*,*) 'dt_set:   ',dt_set
      ierr = -1
   end if
   if (param%dt <= 0.0_DP) then
      write(*,*) 'Invalid timestep: '
      write(*,*) 'dt: ',param%dt
      ierr = -1
   end if
   if (param%inplfile == "") then
      write(*,*) 'No valid planet file in input file'
      ierr = -1
   end if
   if ((param%in_type /= XDR8_TYPE) .and. (param%in_type /= "ASCII")) then
      write(*,*) 'Invalid input file type:',param%in_type
      ierr = -1
   end if
   if ((param%istep_out <= 0) .and. (param%istep_dump <= 0)) then
      write(*,*) 'Invalid istep'
      write(*,*) 'istep_out:  ',param%istep_out
      write(*,*) 'istep_dump: ',param%istep_dump
      ierr = -1
   end if
   if ((param%istep_out > 0) .and. (param%outfile == "")) then
      write(*,*) 'invalid outfile'
      ierr = -1
   end if
   if (param%outfile /= "") then
      if ((param%out_type /= REAL4_TYPE) .and. (param%out_type /= REAL8_TYPE) .and. &
            (param%out_type /= XDR4_TYPE)  .and. (param%out_type /= XDR8_TYPE)) then
         write(*,*) 'Invalid out_type: ',param%out_type
         ierr = -1
      end if
      if ((param%out_form /= "EL") .and. (param%out_form /= "XV") .and. (param%out_form /= "FILT")) then
         write(*,*) 'Invalid out_form: ',param%out_form
         ierr = -1
      end if
      if ((param%out_stat /= "NEW") .and. (param%out_stat /= "UNKNOWN") .and. (param%out_stat /= "APPEND")) then
         write(*,*) 'Invalid out_stat: ',param%out_stat
         ierr = -1
      end if
   end if
   if ((param%j2rp2 == 0.0_DP) .and. (param%j4rp4 /= 0.0_DP)) then
      write(*,*) 'Cannot have j4 without j2'
      ierr = -1
   end if
   if (param%qmin > 0.0_DP) then
      if ((param%qmin_coord /= "HELIO") .and. (param%qmin_coord /= "BARY")) then
         write(*,*) 'Invalid qmin_coord: ',param%qmin_coord
         ierr = -1
      end if
      if ((param%qmin_alo <= 0.0_DP) .or. (param%qmin_ahi <= 0.0_DP)) then
         write(*,*) 'Invalid qmin vals'
         write(*,*) 'qmin_alo: ',param%qmin_alo
         write(*,*) 'qmin_ahi: ',param%qmin_ahi
         ierr = -1
      end if
   end if

   ! Added by D. Minton
   MU2GM = param%MU2GM
   TU2S  = param%TU2S 
   DU2CM = param%DU2CM
   ! The fragmentation model requires the user to set the unit system explicitly.
   write(*, 100, advance = "NO") "FRAGMENTATION  = "
   write(*, *) param%feature%lfragmentation
   if (param%feature%lfragmentation) then
      write(*, 100, advance = "NO") "MU2GM     = "
      write(*, *) MU2GM
      write(*, 100, advance = "NO") "TU2S      = "
      write(*, *) TU2S 
      write(*, 100, advance = "NO") "DU2CM     = "
      write(*, *) DU2CM
      if ((MU2GM < 0.0_DP) .or. (TU2S < 0.0_DP) .or. (DU2CM < 0.0_DP)) then
         write(*,*) 'Invalid unit conversion factor'
         write(*,*) 'MU2GM: ',MU2GM
         write(*,*) 'TU2S: ',TU2S
         write(*,*) 'DU2CM: ',DU2CM
         ierr = -1
      end if
   end if 
   !Added mtiny to the argument list rather than from the terminal
   if (param%mtiny < 0.0_DP) then
      write(*,*) "Invalid MTINY: ",param%mtiny
      ierr = -1
   else
      write(*, 100, advance = "NO") "MTINY     = "
      write(*, *) param%mtiny   
   end if
   if (param%feature%lpython) then
      write(*, 100, advance = "NO") "PYTHON   = "
      write(*, *) param%feature%lpython
   end if
   if (param%feature%lenergy) then
      write(*, 100, advance = "NO") "ENERGY   = "
      write(*, *) param%feature%lenergy
   end if
   if (param%feature%lringmoons) then
      write(*, 100, advance = "NO") "RINGMOONS   = "
      write(*, *) param%feature%lringmoons
   end if


   if (ierr < 0) then
      write(*, 100) "Input parameter(s) failed check"
      call util_exit(FAILURE)
   end if

   !> Define the maximum number of threads
   nthreads = 1            ! In the *serial* case
   !$ nthreads = omp_get_max_threads() ! In the *parallel* case
   !$ write(*,'(a)')   ' OpenMP parameters:'
   !$ write(*,'(a)')   ' ------------------'
   !$ write(*,'(a,i3,/)') ' Number of threads  = ', nthreads   

   return 

   end procedure io_read_param_in

end submodule s_io_read_param_in
