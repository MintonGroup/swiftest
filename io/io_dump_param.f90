!**********************************************************************************************************************************
!
!  Unit Name   : io_dump_param
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Dump integration parameters to file
!
!  Input
!    Arguments : nplmax         : maximum allowed number of planets
!                ntpmax         : maximum allowed number of test particles
!                ntp            : number of active test particles
!                t              : time
!                tstop          : integration stop time
!                dt             : time step
!                in_type        : format of input data files
!                istep_out      : number of time steps between binary outputs
!                outfile        : name of output binary file
!                out_type       : binary format of output file
!                out_form       : data to write to output file
!                istep_dump     : number of time steps between dumps
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                lclose         : logical flag indicating whether to check for planet-test particle encounters
!                rmin           : minimum heliocentric radius for test particle
!                rmax           : maximum heliocentric radius for test particle
!                rmaxu          : maximum unbound heliocentric radius for test particle
!                qmin           : minimum pericenter distance for test particle
!                qmin_coord     : coordinate frame to use for qmin
!                qmin_alo       : minimum semimajor axis for qmin
!                qmin_ahi       : maximum semimajor axis for qmin
!                encounter_file : name of output file for encounters
!                lextra_force   : logical flag indicating whether to use user-supplied accelerations
!                lbig_discard   : logical flag indicating whether to dump planet data with discards
!                lrhill_present : logical flag indicating whether Hill's sphere radii are present in planet data
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : none
!    File      : (same quantities listed above as input from argument are dumped to file, with the following exceptions)
!                ntp            : number of active test particles                    (NOT WRITTEN)
!                inplfile       : name of corresponding planet data dump file        (WRITTEN)
!                intpfile       : name of corresponding test particle data dump file (WRITTEN)
!                out_stat       : open status for output binary file                 (WRITTEN)
!
!  Invocation  : CALL io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,
!                                   istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,
!                                   encounter_file, lextra_force, lbig_discard, lrhill_present)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_dump_param.f
!
!**********************************************************************************************************************************
SUBROUTINE io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form, istep_dump, j2rp2,   &
     j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file, lextra_force, lbig_discard,          &
     lrhill_present)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => io_dump_param
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lclose, lextra_force, lbig_discard, lrhill_present
     INTEGER(I4B), INTENT(IN) :: nplmax, ntpmax, ntp, istep_out, istep_dump
     REAL(DP), INTENT(IN)     :: t, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN) :: qmin_coord, encounter_file, in_type, outfile, out_type, out_form

! Internals
     INTEGER(I4B), PARAMETER :: LUN = 7
     INTEGER(I4B)            :: ierr
     INTEGER(I4B), SAVE      :: idx = 1

! Executable code

     CALL io_open(LUN, DUMP_PARAM_FILE(idx), "REPLACE", "FORMATTED", ierr)
 100 FORMAT(A20)
     WRITE(LUN, 100, ADVANCE = "NO") "NPLMAX "
     WRITE(LUN, *) nplmax
     WRITE(LUN, 100, ADVANCE = "NO") "NTPMAX "
     WRITE(LUN, *) ntpmax
     WRITE(LUN, 100, ADVANCE = "NO") "T0 "
     WRITE(LUN, *) t
     WRITE(LUN, 100, ADVANCE = "NO") "TSTOP "
     WRITE(LUN, *) tstop
     WRITE(LUN, 100, ADVANCE = "NO") "DT "
     WRITE(LUN, *) dt
     WRITE(LUN, 100, ADVANCE = "NO") "PL_IN "
     WRITE(LUN, *) TRIM(DUMP_PL_FILE(idx))
     IF (ntp > 0) THEN
          WRITE(LUN, 100, ADVANCE = "NO") "TP_IN "
          WRITE(LUN, *) TRIM(DUMP_TP_FILE(idx))
     ELSE
          WRITE(LUN, 100) "!TP_IN "
     END IF
     WRITE(LUN, 100, ADVANCE = "NO") "IN_TYPE "
     WRITE(LUN, *) XDR8_TYPE
     IF (istep_out > 0) THEN
          WRITE(LUN, 100, ADVANCE = "NO") "ISTEP_OUT "
          WRITE(LUN, *) istep_out
          WRITE(LUN, 100, ADVANCE = "NO") "BIN_OUT "
          WRITE(LUN, *) TRIM(outfile)
          WRITE(LUN, 100, ADVANCE = "NO") "OUT_TYPE "
          WRITE(LUN, *) TRIM(out_type)
          WRITE(LUN, 100, ADVANCE = "NO") "OUT_FORM "
          WRITE(LUN, *) TRIM(out_form)
          WRITE(LUN, 100, ADVANCE = "NO") "OUT_STAT "
          WRITE(LUN, *) "APPEND"
     ELSE
          WRITE(LUN, 100) "!ISTEP_OUT "
          WRITE(LUN, 100) "!BIN_OUT "
          WRITE(LUN, 100) "!OUT_TYPE "
          WRITE(LUN, 100) "!OUT_FORM "
          WRITE(LUN, 100) "!OUT_STAT "
     END IF
     IF (istep_dump > 0) THEN
          WRITE(LUN, 100, ADVANCE = "NO") "ISTEP_DUMP "
          WRITE(LUN, *) istep_dump
     ELSE
          WRITE(LUN, 100) "!ISTEP_DUMP "
     END IF
     IF (j2rp2 /= 0.0_DP) THEN
          WRITE(LUN, 100, ADVANCE = "NO") "J2 "
          WRITE(LUN, *) j2rp2
          IF (j4rp4 /= 0.0_DP) THEN
               WRITE(LUN, 100, ADVANCE= "NO") "J4 "
               WRITE(LUN, *) j4rp4
          ELSE
               WRITE(LUN, 100) "!J4 "
          END IF
     ELSE
          WRITE(LUN, 100) "!J2 "
          WRITE(LUN, 100) "!J4 "
     END IF
     WRITE(LUN, 100, ADVANCE = "NO") "CHK_CLOSE "
     IF (lclose) THEN
          WRITE(LUN, *) "YES"
     ELSE
          WRITE(LUN, *) "NO"
     END IF
     WRITE(LUN, 100, ADVANCE = "NO") "CHK_RMIN "
     WRITE(LUN, *) rmin
     WRITE(LUN, 100, ADVANCE = "NO") "CHK_RMAX "
     WRITE(LUN, *) rmax
     WRITE(LUN, 100, ADVANCE = "NO") "CHK_EJECT "
     WRITE(LUN, *) rmaxu
     WRITE(LUN, 100, ADVANCE = "NO") "CHK_QMIN "
     WRITE(LUN, *) qmin
     IF (qmin >= 0.0_DP) THEN
          WRITE(LUN, 100, ADVANCE = "NO") "CHK_QMIN_COORD "
          WRITE(LUN, *) TRIM(qmin_coord)
          WRITE(LUN, 100, ADVANCE = "NO") "CHK_QMIN_RANGE "
          WRITE(LUN, *) qmin_alo, qmin_ahi
     ELSE
          WRITE(LUN, 100) "!CHK_QMIN_COORD "
          WRITE(LUN, 100) "!CHK_QMIN_RANGE "
     END IF
     WRITE(LUN, 100, ADVANCE = "NO") "ENC_OUT "
     WRITE(LUN, *) TRIM(encounter_file)
     WRITE(LUN, 100, ADVANCE = "NO") "EXTRA_FORCE "
     IF (lextra_force) THEN
          WRITE(LUN, *) "YES"
     ELSE
          WRITE(LUN, *) "NO"
     END IF
     WRITE(LUN, 100, ADVANCE = "NO") "BIG_DISCARD "
     IF (lbig_discard) THEN
          WRITE(LUN, *) "YES"
     ELSE
          WRITE(LUN, *) "NO"
     END IF
     WRITE(LUN, 100, ADVANCE = "NO") "RHILL_PRESENT "
     IF (lrhill_present) THEN
          WRITE(LUN, *) "YES"
     ELSE
          WRITE(LUN, *) "NO"
     END IF
     CLOSE(UNIT = LUN)
     idx = idx + 1
     IF (idx > 2) idx = 1

     ! The fragmentation model requires the user to set the unit system explicitly.
     WRITE(*, 100, ADVANCE = "NO") "FRAGMENTATION  = "
     WRITE(*, *) lfragmentation
     IF (lfragmentation) THEN
          WRITE(*, 100, ADVANCE = "NO") "MU2GM          = "
          WRITE(*, *) MU2GM
          WRITE(*, 100, ADVANCE = "NO") "TU2S           = "
          WRITE(*, *) TU2S 
          WRITE(*, 100, ADVANCE = "NO") "DU2CM          = "
          WRITE(*, *) DU2CM
          IF ((MU2GM < 0.0_DP) .OR. (TU2S < 0.0_DP) .OR. (DU2CM < 0.0_DP)) ierr = -1
     END IF 


     RETURN

END SUBROUTINE io_dump_param
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
