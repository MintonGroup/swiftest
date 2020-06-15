!**********************************************************************************************************************************
!
!  Unit Name   : coord_vh2vb_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from heliocentric to barycentric coordinates, active test particle velocities only
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!                vs           : barycentric velocity of the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL coord_vh2vb_tp(ntp, swifter_tp1P, vs)
!
!  Notes       : Adapted from Hal Levison's Swift routine coord_vh2b_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_vh2vb_tp(ntp, swiftest_tpA, vs)

! Modules
     USE swiftest, EXCEPT_THIS_ONE => coord_vh2vb_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)              :: ntp
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: vs
     TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA

! Internals
     INTEGER(I4B)              :: i

! Executable code
     DO i = 1, ntp
          IF (swiftest_tpA%status(i) == ACTIVE) swiftest_tpA%vb(:,i) = swiftest_tpA%vh(:,i) + vs(:)
     END DO

     RETURN

END SUBROUTINE coord_vh2vb_tp
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
