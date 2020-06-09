submodule(user) s_user_dump_param
contains
   module procedure user_dump_param
   !! author: David A. Minton
   !!
   !! Dump integration parameters to file
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
   !! Adapted from Martin Duncan's Swift routine io_dump_param.f
   implicit none

   integer(I4B), parameter :: LUN = 7 !! Unit number of output file
   integer(I4B)            :: ierr     !! Error code
   integer(I4B), save      :: idx = 1  !! Index of current dump file. Output flips between 2 files for extra security
                                       !!    in case the program halts during writing
   character(*), dimension(2), parameter :: DUMP_PARAM_FILE = (/ "dump_param1.dat", "dump_param2.dat" /) !! Dump file names

   open(unit = LUN, file = DUMP_PARAM_FILE(idx), status='replace', form = 'formatted', iostat =ierr)
   if (ierr /=0) then
      write(*,*) 'Swiftest error.'
      write(*,*) '   Could not open dump file: ',trim(adjustl(DUMP_PARAM_FILE(idx)))
      call util_exit(FAILURE)
   end if

   100 format(A20) 
   200 format(A20,1X,I0)      ! Format label for integer values
   210 format(A20,2(1X,F0.16,:))   ! Format label for real values
   220 format(A20,1X,A)       ! Format label for string values
   230 format(A20,1X,L1)      ! Format label for logical values
   write(LUN, 200) "NPLMAX ",param%nplmax
   write(LUN, 200) "NTPMAX ",param%ntpmax
   write(LUN, 210) "T0", t
   write(LUN, 210) "TSTOP ",param%tstop
   write(LUN, 210) "DT",param%dt
   write(LUN, 220) "PL_IN ",trim(adjustl(DUMP_PL_FILE(idx)))
   write(LUN, 220) "IN_TYPE ",SWIFTER_REAL8_TYPE
   if (param%istep_out > 0) then
      write(LUN, 200) "ISTEP_OUT",param%istep_out
      write(LUN, 220) "BIN_OUT",trim(adjustl(param%outfile))
      write(LUN, 220) "OUT_TYPE",trim(adjustl(param%out_type))
      write(LUN, 220) "OUT_FORM",trim(adjustl(param%out_form))
      write(LUN, 220) "OUT_STAT","APPEND"
   else
      write(LUN, 100) "!ISTEP_OUT "
      write(LUN, 100) "!BIN_OUT "
      write(LUN, 100) "!OUT_TYPE "
      write(LUN, 100) "!OUT_FORM "
      write(LUN, 100) "!OUT_STAT "
   end if
   if (param%istep_dump > 0) then
      write(LUN, 200) "ISTEP_DUMP",param%istep_dump
   else
      write(LUN, 100) "!ISTEP_DUMP" 
   end if
   if (param%j2rp2 > TINY) then
      write(LUN, 210) "J2 ",param%j2rp2
      if (param%j4rp4 > TINY) then
         write(LUN, 210) "J4 ",param%j4rp4
      else
         write(LUN, 100) "!J4 "
      end if
   else
      write(LUN, 100) "!J2 "
      write(LUN, 100) "!J4 "
   end if
   write(LUN, 230) "CHK_CLOSE ",param%feature%lclose
   write(LUN, 210) "CHK_RMIN ",param%rmin
   write(LUN, 210) "CHK_RMAX ",param%rmax
   write(LUN, 210) "CHK_EJECT ",param%rmaxu
   write(LUN, 210) "CHK_QMIN ",param%qmin
   if (param%qmin >= 0.0_DP) then
      write(LUN, 220) "CHK_QMIN_COORD ",trim(adjustl(param%qmin_coord))
      write(LUN, 210) "CHK_QMIN_RANGE ",param%qmin_alo, param%qmin_ahi
   else
      write(LUN, 100) "!CHK_QMIN_COORD "
      write(LUN, 100) "!CHK_QMIN_RANGE "
   end if
   write(LUN, 220) "ENC_OUT ",trim(adjustl(param%encounter_file))
   write(LUN, 230) "EXTRA_FORCE ",param%feature%lextra_force
   write(LUN, 230) "BIG_DISCARD ",param%feature%lbig_discard
   write(LUN, 230) "RHILL_PRESENT ",param%feature%lrhill_present
   if (param%feature%lmtiny) write(LUN, 210) "MTINY ",param%mtiny
   idx = idx + 1
   if (idx > 2) idx = 1

   ! the fragmentation model requires the user to set the unit system explicitly.
   write(LUN, 230) "FRAGMENTATION", param%feature%lfragmentation
   write(LUN, 210) "MU2KG",MU2KG
   write(LUN, 210) "TU2S ",TU2S 
   write(LUN, 210) "DU2M",DU2M

   close(LUN)

   return

   end procedure user_dump_param
end submodule s_user_dump_param
