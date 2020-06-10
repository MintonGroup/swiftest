submodule(user) s_user_udio_writer
contains
   module procedure user_udio_writer
   !! author: David A. Minton
   !!
   !! Dump integration parameters to file
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
   !! Adapted from Martin Duncan's Swift routine io_dump_param.f
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

   write(unit, Ifmt) "NPLMAX",param%nplmax
   write(unit, Ifmt) "NTPMAX",param%ntpmax
   write(unit, Rfmt) "T0", param%t0
   write(unit, Rfmt) "TSTOP",param%tstop
   write(unit, Rfmt) "DT",param%dt
   write(unit, Sfmt) "PL_IN",trim(adjustl(param%inplfile))
   write(unit, Sfmt) "TP_IN",trim(adjustl(param%intpfile))
   write(unit, Sfmt) "IN_TYPE",trim(adjustl(param%out_type))
   if (param%istep_out > 0) then
      write(unit, Ifmt) "ISTEP_OUT",param%istep_out
      write(unit, Sfmt) "BIN_OUT",trim(adjustl(param%outfile))
      write(unit, Sfmt) "OUT_TYPE",trim(adjustl(param%out_type))
      write(unit, Sfmt) "OUT_FORM",trim(adjustl(param%out_form))
      write(unit, Sfmt) "OUT_STAT","APPEND"
   else
      write(unit, Pfmt) "!ISTEP_OUT "
      write(unit, Pfmt) "!BIN_OUT"
      write(unit, Pfmt) "!OUT_TYPE"
      write(unit, Pfmt) "!OUT_FORM"
      write(unit, Pfmt) "!OUT_STAT"
   end if
   write(unit, Sfmt) "ENC_OUT",trim(adjustl(param%encounter_file))
   if (param%istep_dump > 0) then
      write(unit, Ifmt) "ISTEP_DUMP",param%istep_dump
   else
      write(unit, Pfmt) "!ISTEP_DUMP" 
   end if
   if (param%j2rp2 > TINY) then
      write(unit, Rfmt) "J2 ",param%j2rp2
      if (param%j4rp4 > TINY) then
         write(unit, Rfmt) "J4 ",param%j4rp4
      else
         write(unit, Pfmt) "!J4 "
      end if
   else
      write(unit, Pfmt) "!J2 "
      write(unit, Pfmt) "!J4 "
   end if
   write(unit, Rfmt) "CHK_RMIN",param%rmin
   write(unit, Rfmt) "CHK_RMAX",param%rmax
   write(unit, Rfmt) "CHK_EJECT",param%rmaxu
   write(unit, Rfmt) "CHK_QMIN",param%qmin
   if (param%qmin >= 0.0_DP) then
      write(unit, Sfmt) "CHK_QMIN_COORD",trim(adjustl(param%qmin_coord))
      write(unit, R2fmt) "CHK_QMIN_RANGE",param%qmin_alo, param%qmin_ahi
   else
      write(unit, Pfmt) "!CHK_QMIN_COORD"
      write(unit, Pfmt) "!CHK_QMIN_RANGE"
   end if
   if (param%lmtiny) write(unit, Rfmt) "MTINY",param%mtiny
   write(unit, Rfmt) "MU2KG",MU2KG
   write(unit, Rfmt) "TU2S",TU2S 
   write(unit, Rfmt) "DU2M",DU2M
   
   write(unit, Lfmt) "EXTRA_FORCE",param%lextra_force
   write(unit, Lfmt) "BIG_DISCARD",param%lbig_discard
   write(unit, Lfmt) "RHILL_PRESENT",param%lrhill_present
   write(unit, Lfmt) "CHK_CLOSE",param%lclose
   write(unit, Lfmt) "FRAGMENTATION", param%lfragmentation
   !write(unit, Lfmt) "ROTATION", param%lrotation
   !write(unit, Lfmt) "TIDES", param%ltides
   !write(unit, Lfmt) "GR", param%lgr
   !write(unit, Lfmt) "YARKOVSKY", param%lyarkovsky
   !write(unit, Lfmt) "YORP", param%lyorp
   !write(unit, Lfmt) "ENERGY", param%lenergy

   !write(unit, Lfmt) "RINGMOONS", param%lringmoons
   !if (param%lringmoons) write(unit, Sfmt) "RING_OUTFILE",trim(adjustl(param%ring_outfile))



   return

   end procedure user_udio_writer
end submodule s_user_udio_writer
