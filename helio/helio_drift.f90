!**********************************************************************************************************************************
!
!  Unit Name   : helio_drift
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Loop through planets and call Danby drift routine
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of Swifter planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL helio_drift(npl, swifter_pl1P, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine drift.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_drift(npl, swiftest_plA, dt)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => helio_drift
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     REAL(DP), INTENT(IN)      :: dt
     TYPE(swifter_pl), DIMENSION(:),INTENT(INOUT) :: swifter_plA

! Internals
     INTEGER(I4B)              :: i, iflag
     REAL(DP)                  :: mu

! Executable code
     mu = swiftest_plA%mass(1)
     DO i = 2, npl
          CALL drift_one(mu, swiftest_plA%xh(:,i), swiftest_plA%vb(:,i), dt, iflag)
          IF (iflag /= 0) THEN
               WRITE(*, *) " Planet ", swiftest_plA%id(i), " is lost!!!!!!!!!!"
               WRITE(*, *) mu, dt
               WRITE(*, *) swiftest_plA%xh(:,i)
               WRITE(*, *) swiftest_plA%vb(:,i)
               WRITE(*, *) " STOPPING "
               CALL util_exit(FAILURE)
          END IF
     END DO

     RETURN

END SUBROUTINE helio_drift
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
