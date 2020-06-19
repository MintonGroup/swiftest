submodule(nbody_data_structures) s_io_dump_config
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
end submodule s_io_dump_config
