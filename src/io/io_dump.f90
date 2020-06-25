submodule (swiftest_classes) io_dump
   !! author: David A. Minton
   !! 
   !! This submodule contains implementations of the following procedures:
   !!    io_dump_system
   !!    io_dump_config
   !!    io_dump_config_writer
   !!    io_dump_cb
   !!    io_dump_pl
   !!    io_dump_tp
contains
   module procedure io_dump_system
      !! author: David A. Minton
      !!
      !! Dumps the state of the system to files in case the simulation is interrupted.
      !! As a safety mechanism, there are two dump files that are written in alternating order
      !! so that if a dump file gets corrupted during writing, the user can restart from the older one.
      use swiftest
      implicit none

      config%dump(config,t,dt)
      self%cb%dump(config,t,dt)
      self%pl%dump(config,t,dt)
      self%tp%dump(config,t,dt)

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

      type(swiftest_configuration) :: dump_config   !! Local configuration variable used to configuration change input file names 
                                                    !!    to dump file-specific values without changing the user-defined values
      integer(I4B), parameter      :: LUN = 7       !! Unit number of output file
      integer(I4B)                 :: ierr          !! Error code
      integer(I4B), save           :: idx = 1       !! Index of current dump file. Output flips between 2 files for extra security
                                                    !!    in case the program halts during writing
      character(STRMAX)            :: error_message !! Error message in UDIO procedure

      dump_config = config
      dump_config%t0 = t
      dump_config%incbfile = trim(adjustl(DUMP_CB_FILE(idx))) 
      dump_config%inplfile = trim(adjustl(DUMP_PL_FILE(idx))) 
      dump_config%intpfile = trim(adjustl(DUMP_TP_FILE(idx))) 
      open(unit = LUN, file = DUMP_CONFIG_FILE(idx), status='replace', form = 'formatted', iostat =ierr)
      if (ierr /=0) then
         write(*,*) 'Swiftest error.'
         write(*,*) '   Could not open dump file: ',trim(adjustl(DUMP_CONFIG_FILE(idx)))
         call util_exit(FAILURE)
      end if
      
      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    due to compiler incompatabilities
      !write(LUN,'(DT)') config
      call dump_config_writer(LUN, iotype = "none", v_list = [0], iostat = ierr, iomsg = error_message)
      if (ierr /= 0) then
         write(*,*) trim(adjustl(error_message))
         call util_exit(FAILURE)
      end if

      idx = idx + 1
      if (idx > NDUMPFILES) idx = 1

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

      write(unit, Ifmt) "NPLMAX",                   config%nplmax
      write(unit, Ifmt) "NTPMAX",                   config%ntpmax
      write(unit, Rfmt) "T0",                       config%t0
      write(unit, Rfmt) "TSTOP",                    config%tstop
      write(unit, Rfmt) "DT",                       config%dt
      write(unit, Sfmt) "CB_IN",                    trim(adjustl(config%incbfile))
      write(unit, Sfmt) "PL_IN",                    trim(adjustl(config%inplfile))
      write(unit, Sfmt) "TP_IN",                    trim(adjustl(config%intpfile))
      write(unit, Sfmt) "IN_TYPE",                  trim(adjustl(config%out_type))
      if (config%istep_out > 0) then
         write(unit, Ifmt) "ISTEP_OUT",             config%istep_out
         write(unit, Sfmt) "BIN_OUT",               trim(adjustl(config%outfile))
         write(unit, Sfmt) "OUT_TYPE",              trim(adjustl(config%out_type))
         write(unit, Sfmt) "OUT_FORM",              trim(adjustl(config%out_form))
         write(unit, Sfmt) "OUT_STAT",              "APPEND"
      else
         write(unit, Pfmt) "!ISTEP_OUT "
         write(unit, Pfmt) "!BIN_OUT"
         write(unit, Pfmt) "!OUT_TYPE"
         write(unit, Pfmt) "!OUT_FORM"
         write(unit, Pfmt) "!OUT_STAT"
      end if
      write(unit, Sfmt) "ENC_OUT",                  trim(adjustl(config%encounter_file))
      if (config%istep_dump > 0) then
         write(unit, Ifmt) "ISTEP_DUMP",            config%istep_dump
      else
         write(unit, Pfmt) "!ISTEP_DUMP" 
      end if
      write(unit, Rfmt) "CHK_RMIN",                 config%rmin
      write(unit, Rfmt) "CHK_RMAX",                 config%rmax
      write(unit, Rfmt) "CHK_EJECT",                config%rmaxu
      write(unit, Rfmt) "CHK_QMIN",                 config%qmin
      if (config%qmin >= 0.0_DP) then
         write(unit, Sfmt) "CHK_QMIN_COORD",        trim(adjustl(config%qmin_coord))
         write(unit, R2fmt) "CHK_QMIN_RANGE",       config%qmin_alo, config%qmin_ahi
      else
         write(unit, Pfmt) "!CHK_QMIN_COORD"
         write(unit, Pfmt) "!CHK_QMIN_RANGE"
      end if
      if (config%lmtiny) write(unit, Rfmt) "MTINY", config%mtiny
      write(unit, Rfmt) "MU2KG",                    MU2KG
      write(unit, Rfmt) "TU2S",                     TU2S 
      write(unit, Rfmt) "DU2M",                     DU2M
      
      write(unit, Lfmt) "EXTRA_FORCE",              config%lextra_force
      write(unit, Lfmt) "BIG_DISCARD",              config%lbig_discard
      write(unit, Lfmt) "RHILL_PRESENT",            config%lrhill_present
      write(unit, Lfmt) "CHK_CLOSE",                config%lclose
      write(unit, Lfmt) "FRAGMENTATION",            config%lfragmentation
      write(unit, Lfmt) "ROTATION",                 config%lrotation
      write(unit, Lfmt) "TIDES",                    config%ltides
      write(unit, Lfmt) "GR",                       config%lgr
      write(unit, Lfmt) "ENERGY",                   config%lenergy
      !write(unit, Lfmt) "YARKOVSKY", config%lyarkovsky
      !write(unit, Lfmt) "YORP", config%lyorp

      return
   end procedure io_config_writer

   module procedure io_dump_cb
      !! author: David A. Minton
      !!
      !! Dump massive body data to files
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_dump_pl.f90
      !! Adapted from Hal Levison's Swift routine io_dump_pl.f
      use swiftest
      implicit none
      integer(I4B)                :: i, iu, ierr
      integer(I4B), save          :: idx = 1
      integer(I4B), parameter     :: LUN = 7

      open(unit = LUN, file = DUMP_CB_FILE(idx), form = "unformatted", status = 'replace', iostat = ierr)

      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Unable to open binary dump file ", trim(DUMP_CP_FILE(idx))
         call util_exit(FAILURE)
      end if
      write(LUN) self%mass
      write(LUN) self%radius
      write(LUN) self%j2rp2 
      write(LUN) self%j4rp4 
      if (config%lrotation) then
         write(LUN) self%Ip(:)
         write(LUN) self%rot(:)
      end if
      if (config%tides) then
         write(LUN) self%k2
         write(LUN) self%Q
      end if
      
      close(LUN)
      idx = idx + 1
      if (idx > 2) idx = 1

      return

   end procedure io_dump_cb

   module procedure io_dump_pl
      !! author: David A. Minton
      !!
      !! Dump planet data to files
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_dump_pl.f90
      !! Adapted from Hal Levison's Swift routine io_dump_pl.f
      use swiftest
      implicit none
      integer(I4B)             :: i, iu, ierr
      integer(I4B), save         :: idx = 1
      integer(I4B),parameter         :: LUN = 7

      open(unit = LUN, file = DUMP_PL_FILE(idx), form = "unformatted", status = 'replace', iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Unable to open binary dump file ", trim(DUMP_PL_FILE(idx))
         call util_exit(FAILURE)
      end if
      associate(npl => self%nbody)
         write(LUN) npl
         if (npl > 0) then
            write(LUN) self%name(1:npl)
            write(LUN) self%mass(1:npl)
            if (config%lrhill_present) write(LUN) self%rhill(1:npl) 
            if (config%lclose) write(LUN) self%radius(1:npl) 
            write(LUN) self%xh(:,1:npl)
            write(LUN) self%vh(:,1:npl)
            if (config%lrotation) then
               write(LUN) self%Ip(:,1:npl)
               write(LUN) self%rot(:,1:npl)
            end if
            if (config%ltides) then
               write(LUN) self%k2(1:npl)
               write(LUN) self%Q(1:npl)
            end if
         end if
      end associate
      close(LUN)
      idx = idx + 1
      if (idx > 2) idx = 1

      return

   end procedure io_dump_pl

   module procedure io_dump_tp
      !! author: David A. Minton
      !!
      !! Dump test particle data to files
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_dump_tp.f90
      !! Adapted from Hal Levison's Swift routine io_dump_tp.f
      use swiftest
      implicit none
      integer(I4B)            :: i, iu, ierr
      integer(I4B), save      :: idx = 1
      integer(I4B), parameter :: LUN = 7
   
      open(unit = LUN, file = DUMP_TP_FILE(idx), form = "unformatted", status = 'replace', iostat = ierr)
   
      if (ierr /= 0) then
         write(*, *) "Swiftest error:"
         write(*, *) "   Unable to open binary dump file ", trim(DUMP_TP_FILE(idx))
         call util_exit(FAILURE)
      end if
      associate(ntp => self%nbody)
         write(LUN) ntp
         if (ntp > 0) then
            write(LUN) self%name(1:ntp)
            write(LUN) self%xh(:,1:ntp)
            write(LUN) self%vh(:,1:ntp)
         end if
      end associate
      close(LUN)
      idx = idx + 1
      if (idx > 2) idx = 1
   
      return
   
    end procedure io_dump_tp
end submodule io_dump
