!**********************************************************************************************************************************
!
!  Unit Name   : io_dump_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Dump test particle data to file
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : to dump file
!                ntp          : number of active test particles
!                id           : test particle identifier        (from test particle structure, for each test particle)
!                xh           : heliocentric position           (from test particle structure, for each test particle)
!                vh           : heliocentric velocity           (from test particle structure, for each test particle)
!
!  Invocation  : CALL io_dump_tp(ntp, swifter_tp1P)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_dump_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE io_dump_tp(ntp, swiftest_tpA)

! Modules
     USE swiftest
     USE module_swiftest
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_dump_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)            :: ntp
     TYPE(swiftest_tp), INTENT(INOUT)    :: swiftest_tpA

! Internals
     INTEGER(I4B)                        :: i, iu, ierr
     INTEGER(I4B), SAVE                  :: idx = 1

! Executable code
     CALL io_open_fxdr(DUMP_TP_FILE(idx), "W", .TRUE., iu, ierr)
     IF (ierr /= 0) THEN
          WRITE(*, *) "SWIFTEST Error:"
          WRITE(*, *) "   Unable to open binary dump file ", TRIM(DUMP_TP_FILE(idx))
          CALL util_exit(FAILURE)
     END IF
     ierr = ixdrint(iu, ntp)
     DO i = 1, ntp
          ierr = ixdrint(iu, swiftest_tpA%name(i))
          ierr = ixdrdmat(iu, NDIM, swiftest_tpA%xh(:,i))
          ierr = ixdrdmat(iu, NDIM, swiftest_tpA%vh(:,i))
     END DO
     ierr = ixdrclose(iu)
     idx = idx + 1
     IF (idx > 2) idx = 1

     RETURN

END SUBROUTINE io_dump_tp
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
