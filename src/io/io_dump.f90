submodule (swiftest_classes) io_dump
   !! author: David A. Minton
   !! 
   !! This submodule contains implementations of the following procedures:
   !!    io_dump_system
   !!    io_dump_config
   !!    io_dump_config_writer
   !!    io_dump_swiftest
contains
   module procedure io_dump_system
      !! author: David A. Minton
      !!
      !! Dumps the state of the system to files in case the simulation is interrupted.
      !! As a safety mechanism, there are two dump files that are written in alternating order
      !! so that if a dump file gets corrupted during writing, the user can restart from the older one.
      use swiftest
      implicit none
      type(swiftest_configuration) :: dump_config   !! Local configuration variable used to configuration change input file names 
                                                    !!    to dump file-specific values without changing the user-defined values
      integer(I4B), save           :: idx = 1       !! Index of current dump file. Output flips between 2 files for extra security
                                                    !!    in case the program halts during writing
      character(len=:), allocatable :: config_file_name

      dump_config = config
      config_file_name = trim(adjustl(DUMP_CONFIG_FILE(idx)))
      dump_config%incbfile = trim(adjustl(DUMP_CB_FILE(idx))) 
      dump_config%inplfile = trim(adjustl(DUMP_PL_FILE(idx))) 
      dump_config%intpfile = trim(adjustl(DUMP_TP_FILE(idx)))
      dump_config%out_form = XV
      dump_config%out_stat = 'APPEND'
      call dump_config%dump(config_file_name,t,dt)

      call self%cb%dump(dump_config, t, dt, tfrac)
      call self%pl%dump(dump_config, t, dt, tfrac)
      call self%tp%dump(dump_config, t, dt, tfrac)

      idx = idx + 1
      if (idx > NDUMPFILES) idx = 1

      return
   end procedure io_dump_system

   module procedure io_dump_config
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_config.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_config.f
      use swiftest
      implicit none

      integer(I4B), parameter      :: LUN = 7       !! Unit number of output file
      integer(I4B)                 :: ierr          !! Error code
      character(STRMAX)            :: error_message !! Error message in UDIO procedure

 
      open(unit = LUN, file = config_file_name, status='replace', form = 'FORMATTED', iostat =ierr)
      if (ierr /=0) then
         write(*,*) 'Swiftest error.'
         write(*,*) '   Could not open dump file: ',trim(adjustl(config_file_name))
         call util_exit(FAILURE)
      end if
      
      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    due to compiler incompatabilities
      !write(LUN,'(DT)') config
      call self%config_writer(LUN, iotype = "none", v_list = [0], iostat = ierr, iomsg = error_message)
      if (ierr /= 0) then
         write(*,*) trim(adjustl(error_message))
         call util_exit(FAILURE)
      end if
      close(LUN)

      return
   end procedure io_dump_config

   module procedure io_config_writer
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_config.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_config.f
      use swiftest
      implicit none
                                                               !! In user-defined derived-type output, we need newline characters at the end of each format statement
      !character(*),parameter :: Ifmt  = '(A20,1X,I0/)'         !! Format label for integer values
      !character(*),parameter :: Rfmt  = '(A20,1X,ES25.17/)'    !! Format label for real values 
      !character(*),parameter :: R2fmt = '(A20,2(1X,ES25.17)/)'  !! Format label for 2x real values 
      !character(*),parameter :: Sfmt  = '(A20,1X,A/)'          !! Format label for string values 
      !character(*),parameter :: Lfmt  = '(A20,1X,L1/)'         !! Format label for logical values 
      !character(*),parameter :: Pfmt  = '(A20/)'               !! Format label for single parameter string
      character(*),parameter :: Ifmt  = '(A20,1X,I0)'         !! Format label for integer values
      character(*),parameter :: Rfmt  = '(A20,1X,ES25.17)'    !! Format label for real values 
      character(*),parameter :: R2fmt = '(A20,2(1X,ES25.17))'  !! Format label for 2x real values 
      character(*),parameter :: Sfmt  = '(A20,1X,A)'          !! Format label for string values 
      character(*),parameter :: Lfmt  = '(A20,1X,L1)'         !! Format label for logical values 
      character(*),parameter :: Pfmt  = '(A20)'               !! Format label for single parameter string

      write(unit, Ifmt) "NPLMAX",                   self%nplmax
      write(unit, Ifmt) "NTPMAX",                   self%ntpmax
      write(unit, Rfmt) "T0",                       self%t0
      write(unit, Rfmt) "TSTOP",                    self%tstop
      write(unit, Rfmt) "DT",                       self%dt
      write(unit, Sfmt) "CB_IN",                    trim(adjustl(self%incbfile))
      write(unit, Sfmt) "PL_IN",                    trim(adjustl(self%inplfile))
      write(unit, Sfmt) "TP_IN",                    trim(adjustl(self%intpfile))
      write(unit, Sfmt) "IN_TYPE",                  trim(adjustl(self%out_type))
      if (self%istep_out > 0) then
         write(unit, Ifmt) "ISTEP_OUT",             self%istep_out
         write(unit, Sfmt) "BIN_OUT",               trim(adjustl(self%outfile))
         write(unit, Sfmt) "OUT_TYPE",              trim(adjustl(self%out_type))
         write(unit, Sfmt) "OUT_FORM",              trim(adjustl(self%out_form))
         write(unit, Sfmt) "OUT_STAT",              "APPEND"
      else
         write(unit, Pfmt) "!ISTEP_OUT "
         write(unit, Pfmt) "!BIN_OUT"
         write(unit, Pfmt) "!OUT_TYPE"
         write(unit, Pfmt) "!OUT_FORM"
         write(unit, Pfmt) "!OUT_STAT"
      end if
      write(unit, Sfmt) "ENC_OUT",                  trim(adjustl(self%encounter_file))
      if (self%istep_dump > 0) then
         write(unit, Ifmt) "ISTEP_DUMP",            self%istep_dump
      else
         write(unit, Pfmt) "!ISTEP_DUMP" 
      end if
      write(unit, Rfmt) "CHK_RMIN",                 self%rmin
      write(unit, Rfmt) "CHK_RMAX",                 self%rmax
      write(unit, Rfmt) "CHK_EJECT",                self%rmaxu
      write(unit, Rfmt) "CHK_QMIN",                 self%qmin
      if (self%qmin >= 0.0_DP) then
         write(unit, Sfmt) "CHK_QMIN_COORD",        trim(adjustl(self%qmin_coord))
         write(unit, R2fmt) "CHK_QMIN_RANGE",       self%qmin_alo, self%qmin_ahi
      else
         write(unit, Pfmt) "!CHK_QMIN_COORD"
         write(unit, Pfmt) "!CHK_QMIN_RANGE"
      end if
      if (self%lmtiny) write(unit, Rfmt) "MTINY", self%mtiny
      write(unit, Rfmt) "MU2KG",                    MU2KG
      write(unit, Rfmt) "TU2S",                     TU2S 
      write(unit, Rfmt) "DU2M",                     DU2M
      
      write(unit, Lfmt) "EXTRA_FORCE",              self%lextra_force
      write(unit, Lfmt) "BIG_DISCARD",              self%lbig_discard
      write(unit, Lfmt) "CHK_CLOSE",                self%lclose
      write(unit, Lfmt) "FRAGMENTATION",            self%lfragmentation
      write(unit, Lfmt) "ROTATION",                 self%lrotation
      write(unit, Lfmt) "TIDES",                    self%ltides
      write(unit, Lfmt) "GR",                       self%lgr
      write(unit, Lfmt) "ENERGY",                   self%lenergy
      !write(unit, Lfmt) "YARKOVSKY", self%lyarkovsky
      !write(unit, Lfmt) "YORP", self%lyorp

      return
   end procedure io_config_writer

   module procedure io_dump_swiftest
      !! author: David A. Minton
      !!
      !! Dump massive body data to files
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_dump_pl.f90 and io_dump_tp.f90
      !! Adapted from Hal Levison's Swift routine io_dump_pl.f and io_dump_tp.f
      use swiftest
      implicit none
      integer(I4B)                   :: ierr    !! Error code
      integer(I4B),parameter         :: LUN = 7 !! Unit number for dump file
      integer(I4B)                   :: iu = LUN
      character(len=:), allocatable  :: dump_file_name

      select type(self)
      class is(swiftest_central_body)
         dump_file_name = trim(adjustl(config%incbfile)) 
      class is (swiftest_pl)
         dump_file_name = trim(adjustl(config%inplfile)) 
      class is (swiftest_tp)
         dump_file_name = trim(adjustl(config%intpfile)) 
      end select
      open(unit = iu, file = dump_file_name, form = "UNFORMATTED", status = 'replace', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Unable to open binary dump file " // dump_file_name
         call util_exit(FAILURE)
      end if
      call self%write_frame(iu, config, t, dt)
      close(LUN)

      return
   end procedure io_dump_swiftest
end submodule io_dump
