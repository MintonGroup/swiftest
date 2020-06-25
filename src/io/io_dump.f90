submodule (swiftest_classes) io_dump
contains

   module procedure io_dump_config
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_config.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_config.f
      implicit none

      type(swiftest_configuration)  :: dump_config !! Data type of dumped parameter file
      integer(I4B), parameter :: LUN = 7 !! Unit number of output file
      integer(I4B)            :: ierr     !! Error code
      integer(I4B), save      :: idx = 1  !! Index of current dump file. Output flips between 2 files for extra security
                                          !!    in case the program halts during writing
      character(*), dimension(2), parameter :: DUMP_CONFIG_FILE = (/ "dump_config1.dat", "dump_config2.dat" /) !! Dump file names
      character(STRMAX)       :: error_message           !! Error message in UDIO procedure

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
      call dump_config_writer(LUN, iotype="none",v_list=(/0/),iostat=ierr,iomsg=error_message)

      idx = idx + 1
      if (idx > 2) idx = 1

      close(LUN)

      return

   end procedure io_dump_config


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
