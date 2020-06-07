!**********************************************************************************************************************************
!
!  Unit Name   : io_init_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read in test particle data
!
!  Input
!    Arguments : intpfile     : name of input file for test particles
!                in_type      : format of input data file
!                ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : id           : test particle identifier (all test particles)
!                xh           : heliocentric position    (all test particles)
!                vh           : heliocentric velocity    (all test particles)
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL io_init_tp(intpfile, in_type, ntp, swifter_tp1P)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_init_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE io_init_tp(intpfile, in_type, ntp, symba_tpA)

! Modules
     USE swiftest
     USE module_symba
     USE module_helio
     USE module_swiftest
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_init_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)            :: ntp
     CHARACTER(*), INTENT(IN)            :: intpfile, in_type
     TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA

! Internals
     INTEGER(I4B), PARAMETER             :: LUN = 7
     INTEGER(I4B)                        :: i, iu, ierr, intp
    

! Executable code
     IF (ntp == 0) RETURN
     IF (in_type == "ASCII") THEN
          CALL io_open(LUN, intpfile, "OLD", "FORMATTED", ierr)
          READ(LUN, *) intp
          DO i = 1, ntp
               READ(LUN, *) symba_tpA%helio%swiftest%name(i)
               READ(LUN, *) symba_tpA%helio%swiftest%xh(:,i)
               READ(LUN, *) symba_tpA%helio%swiftest%vh(:,i)
               symba_tpA%helio%swiftest%status(i) = ACTIVE
          END DO
          CLOSE(UNIT = LUN)
     ELSE
          CALL io_open_fxdr(intpfile, "R", .TRUE., iu, ierr)
          ierr = ixdrint(iu, intp)
          DO i = 1, ntp
               ierr = ixdrint(iu, symba_tpA%helio%swiftest%name(i))
               ierr = ixdrdmat(iu, NDIM, symba_tpA%helio%swiftest%xh(:,i))
               ierr = ixdrdmat(iu, NDIM, symba_tpA%helio%swiftest%vh(:,i))
               symba_tpA%helio%swiftest%status(i) = ACTIVE
          END DO
          ierr = ixdrclose(iu)
     END IF

     RETURN

END SUBROUTINE io_init_tp 
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann (Checked by Jennifer Pouplin & Carlisle Wishard)
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
