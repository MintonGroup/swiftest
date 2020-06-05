submodule (io) s_io_read_param_in
contains
   module procedure io_read_param_in
   !! Reads in the param.in file that sets up the run parameters
   !$ use omp_lib
   implicit none

   ! Internals
   logical(LGT)      :: t0_set, tstop_set, dt_set
   integer(I4B), parameter :: LUN = 7
   integer(I4B)      :: ierr = 0, ilength, ifirst, ilast
   character(STRMAX)    :: line, token

   ! Executable code
   nplmax = -1
   ntpmax = -1
   t0_set = .false.
   tstop_set = .false.
   dt_set = .false.
   t0 = 0.0_DP
   tstop = 0.0_DP
   dt = 0.0_DP
   inplfile = ""
   intpfile = ""
   in_type = "ASCII"
   istep_out = -1
   outfile = ""
   out_type = XDR4_TYPE
   out_form = "XV"
   out_stat = "NEW"
   istep_dump = -1
   j2rp2 = 0.0_DP
   j4rp4 = 0.0_DP

   rmin = -1.0_DP
   rmax = -1.0_DP
   rmaxu = -1.0_DP
   qmin = -1.0_DP
   qmin_coord = "HELIO"
   qmin_alo = -1.0_DP
   qmin_ahi = -1.0_DP
   encounter_file = ""

   mtiny = -1.0_DP

   write(*, 100, advance = "NO") "Parameter data file is "
   write(*, 100) inparfile
   write(*, *) " "
   100 format(A)
   call io_open(LUN, inparfile, "OLD", "formatTED", ierr)
   if (ierr /= 0) then
      write(*, 100, advance = "NO") "Unable to open file "
      write(*, 100) inparfile
      call util_exit(FAILURE)
   end if
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
             read(token, *) nplmax
           case ("NTPMAX")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) ntpmax
           case ("T0")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) t0
             t0_set = .true.
           case ("TSTOP")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) tstop
             tstop_set = .true.
           case ("DT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) dt
             dt_set = .true.
           case ("PL_IN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             inplfile = token
           case ("TP_IN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             intpfile = token
           case ("IN_TYPE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             in_type = token
           case ("ISTEP_OUT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) istep_out
           case ("BIN_OUT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             outfile = token
           case ("OUT_TYPE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             out_type = token
           case ("OUT_FORM")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             out_form = token
           case ("OUT_STAT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             out_stat = token
           case ("ISTEP_DUMP")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) istep_dump
           case ("J2")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) j2rp2
           case ("J4")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) j4rp4
           case ("CHK_CLOSE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lclose = .true.
           case ("CHK_RMIN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) rmin
           case ("CHK_RMAX")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) rmax
           case ("CHK_EJECT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) rmaxu
           case ("CHK_QMIN")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) qmin
           case ("CHK_QMIN_COORD")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             qmin_coord = token
           case ("CHK_QMIN_RANGE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) qmin_alo
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) qmin_ahi
           case ("ENC_OUT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             encounter_file = token
           case ("EXTRA_FORCE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lextra_force = .true.
           case ("BIG_DISCARD")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lbig_discard = .true.
           case ("RHILL_PRESENT")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lrhill_present = .true.

           ! Added by the Purdue Swiftest development group (Minton, Wishard, Populin, and Elliott)
           case ("FRAGMENTATION")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lfragmentation = .true.
           case ("MU2GM")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) MU2GM
           case ("TU2S")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) TU2S
           case ("DU2CM")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             read(token, *) DU2CM
           case ("MTINY")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             if (present(mtiny)) read(token, *) mtiny
           case ("PYTHON")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lpython = .true.
           case ("ENERGY")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lenergy = .true.
           case ("RINGMOONS")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lringmoons = .true.
           case ("RING_OUTFILE")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             if (present(ring_outfile)) ring_outfile = token
           case ("ROTATION")
             ifirst = ilast + 1
             call io_get_token(line, ilength, ifirst, ilast, ierr)
             token = line(ifirst:ilast)
             call util_toupper(token)
             if (token == "YES") feature%lrotation = .true. 
           case default
             write(*, 100, advance = "NO") "Unknown parameter -> "
             write(*, *) token
             call util_exit(FAILURE)
         end select
      end if
   end do
   1  close(unit = LUN)
   write(*, 100, advance = "NO") "NPLMAX      = "
   write(*, *) nplmax
   write(*, 100, advance = "NO") "NTPMAX      = "
   write(*, *) ntpmax
   write(*, 100, advance = "NO") "T0       = "
   write(*, *) t0
   write(*, 100, advance = "NO") "TSTOP     = "
   write(*, *) tstop
   write(*, 100, advance = "NO") "DT       = "
   write(*, *) dt
   write(*, 100, advance = "NO") "PL_IN     = "
   ilength = len_trim(inplfile)
   write(*, *) inplfile(1:ilength)
   write(*, 100, advance = "NO") "TP_IN     = "
   ilength = len_trim(intpfile)
   write(*, *) intpfile(1:ilength)
   write(*, 100, advance = "NO") "IN_TYPE     = "
   ilength = len_trim(in_type)
   write(*, *) in_type(1:ilength)
   write(*, 100, advance = "NO") "ISTEP_OUT   = "
   write(*, *) istep_out
   write(*, 100, advance = "NO") "BIN_OUT     = "
   ilength = len_trim(outfile)
   write(*, *) outfile(1:ilength)
   write(*, 100, advance = "NO") "OUT_TYPE    = "
   ilength = len_trim(out_type)
   write(*, *) out_type(1:ilength)
   write(*, 100, advance = "NO") "OUT_FORM    = "
   ilength = len_trim(out_form)
   write(*, *) out_form(1:ilength)
   write(*, 100, advance = "NO") "OUT_STAT    = "
   ilength = len_trim(out_stat)
   write(*, *) out_stat(1:ilength)
   write(*, 100, advance = "NO") "ISTEP_DUMP   = "
   write(*, *) istep_dump
   write(*, 100, advance = "NO") "J2       = "
   write(*, *) j2rp2
   write(*, 100, advance = "NO") "J4       = "
   write(*, *) j4rp4
   write(*, 100, advance = "NO") "CHK_CLOSE   = "
   write(*, *) feature%lclose
   write(*, 100, advance = "NO") "CHK_RMIN    = "
   write(*, *) rmin
   write(*, 100, advance = "NO") "CHK_RMAX    = "
   write(*, *) rmax
   write(*, 100, advance = "NO") "CHK_EJECT   = "
   write(*, *) rmaxu
   write(*, 100, advance = "NO") "CHK_QMIN    = "
   write(*, *) qmin
   write(*, 100, advance = "NO") "CHK_QMIN_COORD = "
   ilength = len_trim(qmin_coord)
   write(*, *) qmin_coord(1:ilength)
   write(*, 100, advance = "NO") "CHK_QMIN_RANGE = "
   write(*, *) qmin_alo, qmin_ahi
   write(*, 100, advance = "NO") "ENC_OUT     = "
   ilength = len_trim(encounter_file)
   write(*, *) encounter_file(1:ilength)
   write(*, 100, advance = "NO") "EXTRA_FORCE   = "
   write(*, *) feature%lextra_force
   write(*, 100, advance = "NO") "BIG_DISCARD   = "
   write(*, *) feature%lbig_discard
   write(*, 100, advance = "NO") "RHILL_PRESENT  = "
   write(*, *) feature%lrhill_present
   write(*, *) " "
   ierr = 0
   if ((.not. t0_set) .or. (.not. tstop_set) .or. (.not. dt_set)) then
      write(*,*) 'Valid simulation time not set'
      write(*,*) 't0_set:   ',t0_set
      write(*,*) 'tstop_set:',tstop_set
      write(*,*) 'dt_set:   ',dt_set
      ierr = -1
   end if
   if (dt <= 0.0_DP) then
      write(*,*) 'Invalid timestep: '
      write(*,*) 'dt: ',dt
      ierr = -1
   end if
   if (inplfile == "") then
      write(*,*) 'No valid planet file in input file'
      ierr = -1
   end if
   if ((in_type /= XDR8_TYPE) .and. (in_type /= "ASCII")) then
      write(*,*) 'Invalid input file type:',in_type
      ierr = -1
   end if
   if ((istep_out <= 0) .and. (istep_dump <= 0)) then
      write(*,*) 'Invalid istep'
      write(*,*) 'istep_out:  ',istep_out
      write(*,*) 'istep_dump: ',istep_dump
      ierr = -1
   end if
   if ((istep_out > 0) .and. (outfile == "")) then
      write(*,*) 'invalid outfile'
      ierr = -1
   end if
   if (outfile /= "") then
      if ((out_type /= REAL4_TYPE) .and. (out_type /= REAL8_TYPE) .and. &
            (out_type /= XDR4_TYPE)  .and. (out_type /= XDR8_TYPE)) then
         write(*,*) 'Invalid out_type: ',out_type
         ierr = -1
      end if
      if ((out_form /= "EL") .and. (out_form /= "XV") .and. (out_form /= "FILT")) then
         write(*,*) 'Invalid out_form: ',out_form
         ierr = -1
      end if
      if ((out_stat /= "NEW") .and. (out_stat /= "UNKNOWN") .and. (out_stat /= "APPEND")) then
         write(*,*) 'Invalid out_stat: ',out_stat
         ierr = -1
      end if
   end if
   if ((j2rp2 == 0.0_DP) .and. (j4rp4 /= 0.0_DP)) then
      write(*,*) 'Cannot have j4 without j2'
      ierr = -1
   end if
   if (qmin > 0.0_DP) then
      if ((qmin_coord /= "HELIO") .and. (qmin_coord /= "BARY")) then
         write(*,*) 'Invalid qmin_coord: ',qmin_coord
         ierr = -1
      end if
      if ((qmin_alo <= 0.0_DP) .or. (qmin_ahi <= 0.0_DP)) then
         write(*,*) 'Invalid qmin vals'
         write(*,*) 'qmin_alo: ',qmin_alo
         write(*,*) 'qmin_ahi: ',qmin_ahi
         ierr = -1
      end if
   end if

   ! Added by D. Minton
   ! The fragmentation model requires the user to set the unit system explicitly.
   write(*, 100, advance = "NO") "FRAGMENTATION  = "
   write(*, *) feature%lfragmentation
   if (feature%lfragmentation) then
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
   if (present(mtiny)) then
      if (mtiny < 0.0_DP) then
         write(*,*) "MTINY not set or invalid value"
         ierr = -1
      else
         write(*, 100, advance = "NO") "MTINY     = "
         write(*, *) mtiny   
      end if
   end if
   if (feature%lpython) then
      write(*, 100, advance = "NO") "PYTHON   = "
      write(*, *) feature%lpython
   end if
   if (feature%lenergy) then
      write(*, 100, advance = "NO") "ENERGY   = "
      write(*, *) feature%lenergy
   end if
   if (feature%lringmoons) then
      write(*, 100, advance = "NO") "RINGMOONS   = "
      write(*, *) feature%lringmoons
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
