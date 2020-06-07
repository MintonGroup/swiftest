!**********************************************************************************************************************************
!
!  Unit Name   : helio_lindrift_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Perform linear drift of test particles due to barycentric momentum of Sun
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!                dt           : time step
!                pt           : negative barycentric velocity of the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_lindrift_tp(ntp, swifter_tp1P, dt, pt)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_lindrift_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_lindrift_tp(ntp, swiftest_tpA, dt, pt)

! Modules
     USE swiftest
     USE module_swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => helio_lindrift_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)              :: ntp
     REAL(DP), INTENT(IN)                  :: dt
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: pt
     TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA

! Internals
     INTEGER(I4B)              :: i

! Executable code
     DO i = 1, ntp
          IF (swiftest_tpA%status(i) == ACTIVE) THEN 
               swiftest_tpA%xh(:,i) = swiftest_tpA%xh(:,i) + pt(:)*dt
          END IF
     END DO

     RETURN

END SUBROUTINE helio_lindrift_tp
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
