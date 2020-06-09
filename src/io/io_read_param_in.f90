submodule (io) s_io_read_param_in
contains
   module procedure io_read_param_in
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Read in parameters for the integration
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
   !! Adapted from Martin Duncan's Swift routine io_init_param.f
   !$ use omp_lib
   use io
   implicit none

   integer(I4B), parameter :: LUN = 7                 !! Unit number of input file
   logical                 :: t0_set = .false.        !! Is the initial time set in the input file?
   logical                 :: tstop_set = .false.     !! Is the final time set in the input file?
   logical                 :: dt_set = .false.        !! Is the step size set in the input file?
   integer(I4B)            :: ierr = 0                !! Input error code
   integer(I4B)            :: ilength, ifirst, ilast  !! Variables used to parse input file
   character(STRMAX)       :: line                    !! Line of the input file
   character (len=:), allocatable :: line_trim,param_name, param_value

   ! Read in name of parameter file
   write(*, *) "Parameter data file is ", trim(adjustl(inparfile))
   write(*, *) " "
   100 format(A)
   open(unit = LUN, file = inparfile, status = 'old', iostat = ierr)
   if (ierr /= 0) then
      write(*,*) "Unable to open file ",ierr
      write(*,*) inparfile
      call util_exit(FAILURE)
   end if

   ! Parse the file line by line, extracting tokens then matching them up with known parameters if possible
   do
      read(LUN, 100, IOSTAT = ierr, end = 1) line
      line_trim = trim(adjustl(line))
      ilength = len(line_trim)
      if ((ilength /= 0)) then 
         ifirst = 1
         ! Read the pair of tokens. The first one is the parameter name, the second is the value.
         param_name = io_get_token(line_trim, ifirst, ilast, ierr)
         if (param_name == '') cycle ! No parameter name (usually because this line is commented out)
         call util_toupper(param_name)
         ifirst = ilast + 1
         param_value = io_get_token(line_trim, ifirst, ilast, ierr)
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
            param_value = io_get_token(line, ifirst, ilast, ierr)
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
         case ("PYTHON")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lpython = .true.
         case ("ENERGY")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lenergy = .true.
         case ("RINGMOONS")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lringmoons = .true.
         case ("RING_OUTFILE")
            param%ring_outfile = param_value
         case ("ROTATION")
            call util_toupper(param_value)
            if (param_value == "YES" .or. param_value == 'T') param%feature%lrotation = .true. 
         case default
            write(*,*) "Unknown parameter -> ",param_name
            call util_exit(FAILURE)
         end select
      end if
   end do
   1 close(LUN)
   write(*,*) "NPLMAX         = ",param%nplmax
   write(*,*) "NTPMAX         = ",param%ntpmax
   write(*,*) "T0             = ",param%t0
   write(*,*) "TSTOP          = ",param%tstop
   write(*,*) "DT             = ",param%dt
   write(*,*) "PL_IN          = ",trim(adjustl(param%inplfile))
   write(*,*) "TP_IN          = ",trim(adjustl(param%intpfile))
   write(*,*) "IN_TYPE        = ",trim(adjustl(param%in_type))
   write(*,*) "ISTEP_OUT      = ",param%istep_out
   write(*,*) "BIN_OUT        = ",trim(adjustl(param%outfile))
   write(*,*) "OUT_TYPE       = ",trim(adjustl(param%out_type))
   write(*,*) "OUT_FORM       = ",trim(adjustl(param%out_form))
   write(*,*) "OUT_STAT       = ",trim(adjustl(param%out_stat))
   write(*,*) "ISTEP_DUMP     = ",param%istep_dump
   write(*,*) "J2             = ",param%j2rp2
   write(*,*) "J4             = ",param%j4rp4
   write(*,*) "CHK_CLOSE      = ",param%feature%lclose
   write(*,*) "CHK_RMIN       = ",param%rmin
   write(*,*) "CHK_RMAX       = ",param%rmax
   write(*,*) "CHK_EJECT      = ",param%rmaxu
   write(*,*) "CHK_QMIN       = ",param%qmin
   write(*,*) "CHK_QMIN_COORD = ",trim(adjustl(param%qmin_coord))
   write(*,*) "CHK_QMIN_RANGE = ",param%qmin_alo, param%qmin_ahi
   write(*,*) "ENC_OUT        = ",trim(adjustl(param%encounter_file))
   write(*,*) "EXTRA_FORCE    = ",param%feature%lextra_force
   write(*,*) "BIG_DISCARD    = ",param%feature%lbig_discard
   write(*,*) "RHILL_PRESENT  = ",param%feature%lrhill_present
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
   if ((param%in_type /= SWIFTER_REAL8_TYPE) .and. (param%in_type /= "ASCII")) then
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
      write(*,*) 'Invalid outfile'
      ierr = -1
   end if
   if (param%outfile /= "") then
      if ((param%out_type /= REAL4_TYPE) .and. (param%out_type /= REAL8_TYPE) .and. &
            (param%out_type /= SWIFTER_REAL4_TYPE)  .and. (param%out_type /= SWIFTER_REAL8_TYPE)) then
         write(*,*) 'Invalid out_type: ',param%out_type
         ierr = -1
      end if
      if ((param%out_form /= "EL") .and. (param%out_form /= "XV")) then
         write(*,*) 'Invalid out_form: ',param%out_form
         ierr = -1
      end if
      if ((param%out_stat /= "NEW") .and. (param%out_stat /= "REPLACE") .and. (param%out_stat /= "APPEND")) then
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
   MU2KG = param%MU2KG
   TU2S  = param%TU2S 
   DU2M = param%DU2M
   ! The fragmentation model requires the user to set the unit system explicitly.
   write(*,*) "FRAGMENTATION  = ",param%feature%lfragmentation
   if (param%feature%lfragmentation) then
      write(*,*) "MU2KG          = ",MU2KG
      write(*,*) "TU2S           = ",TU2S 
      write(*,*) "DU2M          = ",DU2M
      if ((MU2KG < 0.0_DP) .or. (TU2S < 0.0_DP) .or. (DU2M < 0.0_DP)) then
         write(*,*) 'Invalid unit conversion factor'
         write(*,*) 'MU2KG: ',MU2KG
         write(*,*) 'TU2S: ',TU2S
         write(*,*) 'DU2M: ',DU2M
         ierr = -1
      end if
   end if 
   !Added mtiny to the argument list rather than from the terminal
   if (param%mtiny < 0.0_DP) then
      write(*,*) "Invalid MTINY: ",param%mtiny
      ierr = -1
   else
      write(*,*) "MTINY          = ",param%mtiny   
   end if
   if (param%feature%lpython) write(*,*) "PYTHON         = ",param%feature%lpython
   if (param%feature%lenergy) write(*,*) "ENERGY         = ",param%feature%lenergy
   if (param%feature%lringmoons) write(*,*) "RINGMOONS      = ",param%feature%lringmoons

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
