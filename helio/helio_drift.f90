!**********************************************************************************************************************************
!
!  Unit Name   : helio_drift
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE helio_drift(npl, swifter_pl1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => helio_drift
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     REAL(DP), INTENT(IN)      :: dt
     TYPE(swifter_pl), POINTER :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i, iflag
     REAL(DP)                  :: mu
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     mu = swifter_pl1P%mass
     swifter_plP => swifter_pl1P
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          CALL drift_one(mu, swifter_plP%xh(:), swifter_plP%vb(:), dt, iflag)
          IF (iflag /= 0) THEN
               WRITE(*, *) " Planet ", swifter_plP%id, " is lost!!!!!!!!!!"
               WRITE(*, *) mu, dt
               WRITE(*, *) swifter_plP%xh(:)
               WRITE(*, *) swifter_plP%vb(:)
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
