!**********************************************************************************************************************************
!
!  Unit Name   : io_discard_write
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write out information about discarded test particles
!
!  Input
!    Arguments : t             : time
!                npl           : number of planets
!                nsp           : number of spilled test particles
!                swifter_pl1P  : pointer to head of Swifter planet structure linked-list
!                swifter_tpd1P : pointer to head of discard Swifter test particle structure linked-list
!                fname         : name of file to write
!                lbig_discard  : logical flag indicating whether to dump planet data with discards
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : discard data to discard file
!
!  Invocation  : CALL io_discard_write(t, npl, nsp, swifter_pl1P, swifter_tpd1P, fname, lbig_discard)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_discard_write.f
!
!**********************************************************************************************************************************
SUBROUTINE io_discard_write(t, npl, nsp, swifter_pl1P, swifter_tpd1P, fname, lbig_discard)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => io_discard_write
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lbig_discard
     INTEGER(I4B), INTENT(IN)  :: npl, nsp
     REAL(DP), INTENT(IN)      :: t
     CHARACTER(*), INTENT(IN)  :: fname
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tpd1P

! Internals
     INTEGER(I4B), PARAMETER   :: LUN = 40
     INTEGER(I4B)              :: i, ierr
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     CALL io_open(LUN, fname, "APPEND", "FORMATTED", ierr)
     IF (ierr /= 0) THEN
          CALL io_open(LUN, fname, "NEW", "FORMATTED", ierr)
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to open discard output file, ", fname
               CALL util_exit(FAILURE)
          END IF
     END IF
     WRITE(LUN, 100) t, nsp, lbig_discard
 100 FORMAT(E23.16, 1X, I8, 1X, L1)
     swifter_tpP => swifter_tpd1P
     DO i = 1, nsp
          WRITE(LUN, 200) SUB, swifter_tpP%id, swifter_tpP%status
 200      FORMAT(A, 2(1X, I8))
          WRITE(LUN, 300) swifter_tpP%xh(:)
 300      FORMAT(3(E23.16, 1X))
          WRITE(LUN, 300) swifter_tpP%vh(:)
          swifter_tpP => swifter_tpP%nextP
     END DO
     IF (lbig_discard) THEN
          WRITE(LUN, 400) npl
 400      FORMAT(I8)
          swifter_plP => swifter_pl1P
          DO i = 2, npl
               swifter_plP => swifter_plP%nextP
               WRITE(LUN, 500) swifter_plP%id, swifter_plP%mass, swifter_plP%radius
 500           FORMAT(I8, 2(1X, E23.16))
               WRITE(LUN, 300) swifter_plP%xh(:)
               WRITE(LUN, 300) swifter_plP%vh(:)
          END DO
     END IF
     CLOSE(LUN)

     RETURN

END SUBROUTINE io_discard_write
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
