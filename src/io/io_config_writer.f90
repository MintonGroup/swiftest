submodule(nbody_data_structures) s_io_config_writer
contains
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

   write(unit, Ifmt) "NPLMAX",config%nplmax
   write(unit, Ifmt) "NTPMAX",config%ntpmax
   write(unit, Rfmt) "T0", config%t0
   write(unit, Rfmt) "TSTOP",config%tstop
   write(unit, Rfmt) "DT",config%dt
   write(unit, Sfmt) "PL_IN",trim(adjustl(config%inplfile))
   write(unit, Sfmt) "TP_IN",trim(adjustl(config%intpfile))
   write(unit, Sfmt) "IN_TYPE",trim(adjustl(config%out_type))
   if (config%istep_out > 0) then
      write(unit, Ifmt) "ISTEP_OUT",config%istep_out
      write(unit, Sfmt) "BIN_OUT",trim(adjustl(config%outfile))
      write(unit, Sfmt) "OUT_TYPE",trim(adjustl(config%out_type))
      write(unit, Sfmt) "OUT_FORM",trim(adjustl(config%out_form))
      write(unit, Sfmt) "OUT_STAT","APPEND"
   else
      write(unit, Pfmt) "!ISTEP_OUT "
      write(unit, Pfmt) "!BIN_OUT"
      write(unit, Pfmt) "!OUT_TYPE"
      write(unit, Pfmt) "!OUT_FORM"
      write(unit, Pfmt) "!OUT_STAT"
   end if
   write(unit, Sfmt) "ENC_OUT",trim(adjustl(config%encounter_file))
   if (config%istep_dump > 0) then
      write(unit, Ifmt) "ISTEP_DUMP",config%istep_dump
   else
      write(unit, Pfmt) "!ISTEP_DUMP" 
   end if
   if (config%j2rp2 > VSMALL) then
      write(unit, Rfmt) "J2 ",config%j2rp2
      if (config%j4rp4 > VSMALL) then
         write(unit, Rfmt) "J4 ",config%j4rp4
      else
         write(unit, Pfmt) "!J4 "
      end if
   else
      write(unit, Pfmt) "!J2 "
      write(unit, Pfmt) "!J4 "
   end if
   write(unit, Rfmt) "CHK_RMIN",config%rmin
   write(unit, Rfmt) "CHK_RMAX",config%rmax
   write(unit, Rfmt) "CHK_EJECT",config%rmaxu
   write(unit, Rfmt) "CHK_QMIN",config%qmin
   if (config%qmin >= 0.0_DP) then
      write(unit, Sfmt) "CHK_QMIN_COORD",trim(adjustl(config%qmin_coord))
      write(unit, R2fmt) "CHK_QMIN_RANGE",config%qmin_alo, config%qmin_ahi
   else
      write(unit, Pfmt) "!CHK_QMIN_COORD"
      write(unit, Pfmt) "!CHK_QMIN_RANGE"
   end if
   if (config%lmtiny) write(unit, Rfmt) "MTINY",config%mtiny
   write(unit, Rfmt) "MU2KG",MU2KG
   write(unit, Rfmt) "TU2S",TU2S 
   write(unit, Rfmt) "DU2M",DU2M
   
   write(unit, Lfmt) "EXTRA_FORCE",config%lextra_force
   write(unit, Lfmt) "BIG_DISCARD",config%lbig_discard
   write(unit, Lfmt) "RHILL_PRESENT",config%lrhill_present
   write(unit, Lfmt) "CHK_CLOSE",config%lclose
   write(unit, Lfmt) "FRAGMENTATION", config%lfragmentation
   !write(unit, Lfmt) "ROTATION", config%lrotation
   !write(unit, Lfmt) "TIDES", config%ltides
   !write(unit, Lfmt) "GR", config%lgr
   !write(unit, Lfmt) "YARKOVSKY", config%lyarkovsky
   !write(unit, Lfmt) "YORP", config%lyorp
   !write(unit, Lfmt) "ENERGY", config%lenergy

   !write(unit, Lfmt) "RINGMOONS", config%lringmoons
   !if (config%lringmoons) write(unit, Sfmt) "RING_OUTFILE",trim(adjustl(config%ring_outfile))



   return

   end procedure io_config_writer
end submodule s_io_config_writer
