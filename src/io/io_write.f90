submodule (swiftest_classes) s_io_write
   !! author: David A. Minton
   !! 
   !! This submodule contains implementations of the following procedures:
   !!    io_write_frame_system
   !!    io_write_hdr
   !!    io_write_frame_cb
   !!    io_write_frame_body
   !!    io_write_write_encounter
   !!    io_dump_system
   !!    io_dump_config
   !!    io_config_writer
   !!    io_dump_swiftest
   !!    io_write_discards
contains
   module procedure io_write_frame_system
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to output binary file
      !! There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      use swiftest
      implicit none

      logical, save                         :: lfirst = .true. !! Flag to determine if this is the first call of this method
      integer(I4B)                          :: ierr            !! I/O error code

      class(swiftest_cb), allocatable       :: cb         !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_pl), allocatable       :: pl              !! Temporary local version of pl structure used for non-destructive conversions
      class(swiftest_tp), allocatable       :: tp              !! Temporary local version of pl structure used for non-destructive conversions

      allocate(cb, source = self%cb)
      allocate(pl, source = self%pl)
      allocate(tp, source = self%tp)
      iu = BINUNIT

      if (lfirst) then
         select case(config%out_stat)
         case('APPEND')
            open(unit = iu, file = config%outfile, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
         case('NEW')
            open(unit = iu, file = config%outfile, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
         case ('REPLACE')
            open(unit = iu, file = config%outfile, status = 'REPLACE', form = 'UNFORMATTED', iostat = ierr)
         case default
            write(*,*) 'Invalid status code',trim(adjustl(config%out_stat))
            call util_exit(FAILURE)
         end select
         if (ierr /= 0) then
            write(*, *) "Swiftest error: io_write_frame_system - first", ierr
            write(*, *) "   Binary output file " // trim(adjustl(config%outfile)) // " already exists or cannot be accessed"
            write(*, *) "   out_stat: " // trim(adjustl(config%out_stat))
            call util_exit(FAILURE)
         end if
         lfirst = .false.
      else
         open(unit = iu, file = config%outfile, status = 'OLD', position =  'APPEND', form = 'UNFORMATTED', iostat = ierr)
         if (ierr /= 0) then
            write(*, *) "Swiftest error: io_write_frame_system"
            write(*, *) "   Unable to open binary output file for APPEND"
            call util_exit(FAILURE)
         end if
      end if
      call io_write_hdr(iu, t, pl%nbody, tp%nbody, config%out_form, config%out_type)

      if (config%lgr) then
         associate(vh => pl%vh, vht => tp%vh)
            select type(pl)
            class is (whm_pl)
               call pl%gr_pv2vh(config)
            end select
            select type(tp) 
            class is (whm_tp)
               call tp%gr_pv2vh(config)
            end select
         end associate
      end if

      if (config%out_form == EL) then ! Do an orbital element conversion prior to writing out the frame, as we have access to the central body here
         call pl%xv2el(cb)
         call tp%xv2el(cb)
      end if
      
      ! Write out each data type frame
      call cb%write_frame(iu, config, t, dt)
      call pl%write_frame(iu, config, t, dt)
      call tp%write_frame(iu, config, t, dt)

      deallocate(cb, pl, tp)

      close(iu)

      return
   end procedure io_write_frame_system

   module procedure io_write_hdr
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write frame header to output binary file
      !!
      !! Adapted from David Adapted from David E. Kaufmann's Swifter routine io_write_hdr.f90
      !! Adapted from Hal Levison's Swift routine io_write_hdr.F
      use swiftest
      implicit none
   
      integer(I4B)               :: ierr !! Error code
   
      select case (out_type)
      case (REAL4_TYPE,SWIFTER_REAL4_TYPE)
         write(iu, iostat = ierr) real(t, kind=SP)
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Unable to write binary file header"
            call util_exit(FAILURE)
         end if
      case (REAL8_TYPE,SWIFTER_REAL8_TYPE)
         write(iu, iostat = ierr) t
         if (ierr /= 0) then
            write(*, *) "Swiftest error:"
            write(*, *) "   Unable to write binary file header"
            call util_exit(FAILURE)
         end if
      end select
      write(iu, iostat = ierr) npl
      write(iu, iostat = ierr) ntp
      write(iu, iostat = ierr) out_form
   
      return
   
   end procedure io_write_hdr

   module procedure io_write_frame_cb
      !! author: David A. Minton
      !!
      !! Write a frame of output of central body data to the binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      use swiftest
      implicit none

      write(iu) self%Gmass
      write(iu) self%radius
      write(iu) self%j2rp2 
      write(iu) self%j4rp4 
      if (config%lrotation) then
         write(iu) self%Ip(1)
         write(iu) self%Ip(2)
         write(iu) self%Ip(3)
         write(iu) self%rot(1)
         write(iu) self%rot(2)
         write(iu) self%rot(3)
      end if
      if (config%ltides) then
         write(iu) self%k2
         write(iu) self%Q
      end if

      return
   end procedure io_write_frame_cb

   module procedure io_write_frame_body
      !! author: David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      use swiftest
      implicit none

      associate(n => self%nbody)
         if (n == 0) return
         write(iu) self%name(1:n)
         select case (config%out_form)
         case (EL) 
            write(iu) self%a(1:n)
            write(iu) self%e(1:n)
            write(iu) self%inc(1:n)
            write(iu) self%capom(1:n)
            write(iu) self%omega(1:n)
            write(iu) self%capm(1:n)
         case (XV)
            write(iu) self%xh(1, 1:n)
            write(iu) self%xh(2, 1:n)
            write(iu) self%xh(3, 1:n)
            write(iu) self%vh(1, 1:n)
            write(iu) self%vh(2, 1:n)
            write(iu) self%vh(3, 1:n)
         end select
         select type(self)  
         class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
            write(iu) self%Gmass(1:n)
            write(iu) self%radius(1:n)
            if (config%lrotation) then
               write(iu) self%Ip(1, 1:n)
               write(iu) self%Ip(2, 1:n)
               write(iu) self%Ip(3, 1:n)
               write(iu) self%rot(1, 1:n)
               write(iu) self%rot(2, 1:n)
               write(iu) self%rot(3, 1:n)
            end if
            if (config%ltides) then
               write(iu) self%k2(1:n)
               write(iu) self%Q(1:n)
            end if
         end select
      end associate

      return
   end procedure io_write_frame_body

   module procedure io_write_encounter
      !! author: David A. Minton
      !!
      !! Write close encounter data to output binary files
      !!  There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: io_write_encounter.f90
      !! Adapted from Hal Levison's Swift routine io_write_encounter.f
      use swiftest
      implicit none
      logical         :: lxdr
      logical , save    :: lfirst = .true.
      integer(I4B), parameter :: lun = 30
      integer(I4B)        :: ierr
      integer(I4B), save    :: iu = lun

      open(unit = iu, file = encounter_file, status = 'OLD', position = 'APPEND', form = 'UNFORMATTED', iostat = ierr)
      if ((ierr /= 0) .and. lfirst) then
         open(unit = iu, file = encounter_file, status = 'NEW', form = 'UNFORMATTED', iostat = ierr)
      end if
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to open binary encounter file"
         call util_exit(FAILURE)
      end if
      lfirst = .false.
      write(iu, iostat = ierr) t
      if (ierr < 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to write binary file record"
         call util_exit(FAILURE)
      end if
      write(iu) name1, xh1(1), xh1(2), xh1(3), vh1(1), vh1(2), mass1, radius1
      write(iu) name2, xh2(1), xh2(2), xh2(3), vh2(1), vh2(2), mass2, radius2
      close(unit = iu, iostat = ierr)
      if (ierr /= 0) then
         write(*, *) "Swiftest Error:"
         write(*, *) "   Unable to close binary encounter file"
         call util_exit(FAILURE)
      end if

      return

   end procedure io_write_encounter

   module procedure io_dump_system
      !! author: David A. Minton
      !!
      !! Dumps the state of the system to files in case the simulation is interrupted.
      !! As a safety mechanism, there are two dump files that are written in alternating order
      !! so that if a dump file gets corrupted during writing, the user can restart from the older one.
      use swiftest
      implicit none
      class(swiftest_configuration), allocatable :: dump_config   !! Local configuration variable used to configuration change input file names 
                                                    !!    to dump file-specific values without changing the user-defined values
      integer(I4B), save           :: idx = 1       !! Index of current dump file. Output flips between 2 files for extra security
                                                    !!    in case the program halts during writing
      character(len=:), allocatable :: config_file_name
      real(DP) :: tfrac
     

      allocate(dump_config, source=config)
      config_file_name = trim(adjustl(DUMP_CONFIG_FILE(idx)))
      dump_config%incbfile = trim(adjustl(DUMP_CB_FILE(idx))) 
      dump_config%inplfile = trim(adjustl(DUMP_PL_FILE(idx))) 
      dump_config%intpfile = trim(adjustl(DUMP_TP_FILE(idx)))
      dump_config%out_form = XV
      dump_config%out_stat = 'APPEND'
      call dump_config%dump(config_file_name,t,dt)

      call self%cb%dump(dump_config, t, dt)
      if (self%pl%nbody > 0) call self%pl%dump(dump_config, t, dt)
      if (self%tp%nbody > 0) call self%tp%dump(dump_config, t, dt)

      idx = idx + 1
      if (idx > NDUMPFILES) idx = 1

      ! Print the status message (format code passed in from main driver)
      tfrac = (t - config%t0) / (config%tstop - config%t0)
      write(*,msg) t, tfrac, self%pl%nbody, self%tp%nbody

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
      call self%writer(LUN, iotype = "none", v_list = [0], iostat = ierr, iomsg = error_message)
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
      if (self%lmtiny) write(unit, Rfmt) "MTINY",   self%mtiny
      write(unit, Rfmt) "MU2KG",                    self%MU2KG
      write(unit, Rfmt) "TU2S",                     self%TU2S 
      write(unit, Rfmt) "DU2M",                     self%DU2M
      
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
      class is(swiftest_cb)
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

   module procedure io_write_discard 
      !!
      !! Write out information about discarded test particle
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_discard_write.f90
      !! Adapted from Hal Levison's Swift routine io_discard_write.f
      use swiftest
      implicit none
      integer(I4B), parameter   :: LUN = 40
      integer(I4B)          :: i, ierr
      logical, save :: lfirst = .true. 
      real(DP), dimension(:,:), allocatable :: vh
      character(*), parameter :: HDRFMT    = '(E23.16, 1X, I8, 1X, L1)'
      character(*), parameter :: NAMEFMT   = '(A, 2(1X, I8))'
      character(*), parameter :: VECFMT    = '(3(E23.16, 1X))'
      character(*), parameter :: NPLFMT    = '(I8)'
      character(*), parameter :: PLNAMEFMT = '(I8, 2(1X, E23.16))'
      class(swiftest_body), allocatable :: pltemp

      associate(t => config%t, config => config, nsp => discards%nbody, dxh => discards%xh, dvh => discards%vh, &
                dname => discards%name, dstatus => discards%status) 
         
         if (config%out_stat == 'APPEND' .or. (.not.lfirst)) then
            open(unit = LUN, file = DISCARD_FILE, status = 'OLD', position = 'APPEND', form = 'FORMATTED', iostat = ierr)
         else if (config%out_stat == 'NEW') then
            open(unit = LUN, file = DISCARD_FILE, status = 'NEW', form = 'FORMATTED', iostat = ierr)
         else if (config%out_stat == 'REPLACE') then
            open(unit = LUN, file = DISCARD_FILE, status = 'REPLACE', form = 'FORMATTED', iostat = ierr)
         else
            write(*,*) 'Invalid status code',trim(adjustl(config%out_stat))
            call util_exit(FAILURE)
         end if
         lfirst = .false.
         if (config%lgr) then
            select type(discards)
            class is (whm_tp)
               call discards%gr_pv2vh(config)
            end select
         end if
         write(LUN, HDRFMT) t, nsp, config%lbig_discard
         do i = 1, nsp
            write(LUN, NAMEFMT) sub, dname(i), dstatus(i)
            write(LUN, VECFMT) dxh(1, i), dxh(2, i), dxh(3, i)
            write(LUN, VECFMT) dvh(1, i), dvh(2, i), dvh(3, i)
         end do
         if (config%lbig_discard) then
            associate(npl => self%pl%nbody, pl => self%pl, GMpl => self%pl%Gmass, &
                     Rpl => self%pl%radius, name => self%pl%name, xh => self%pl%xh)

               if (config%lgr) then
                  allocate(pltemp, source = pl)
                  select type(pltemp)
                  class is (whm_pl)
                     call pltemp%gr_pv2vh(config)
                     allocate(vh, source = pltemp%vh)
                  end select
                  deallocate(pltemp)
               else
                  allocate(vh, source = pl%vh)
               end if

               write(LUN, NPLFMT) npl
               do i = 1, npl
                  write(LUN, PLNAMEFMT) name(i), GMpl(i), Rpl(i)
                  write(LUN, VECFMT) xh(1, i), xh(2, i), xh(3, i)
  
                  write(LUN, VECFMT) vh(1, i), vh(2, i), vh(3, i)
               end do
               deallocate(vh)
            end associate
         end if
         close(LUN)
      end associate
      return
   
   end procedure io_write_discard


end submodule s_io_write
