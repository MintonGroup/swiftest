!**********************************************************************************************************************************
!
!  Unit Name   : discard_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : discard
!  Language    : Fortran 90/95
!
!  Description : Check to see if test particles should be discarded based on their positions relative to the planets
!
!  Input
!    Arguments : t            : time
!                dt           : time step
!                npl          : number of planets
!                ntp          : number of active test particles
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL discard_pl(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE discard_pl(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => discard_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     REAL(DP), INTENT(IN)      :: t, dt
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     INTEGER(I4B)              :: i, j, isp
     REAL(DP)                  :: r2min, radius
     REAL(DP), DIMENSION(NDIM) :: dx, dv
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     swifter_tpP => swifter_tp1P
     DO i = 1, ntp
          IF (swifter_tpP%status == ACTIVE) THEN
               swifter_plP => swifter_pl1P
               DO j = 2, npl
                    swifter_plP => swifter_plP%nextP
                    dx(:) = swifter_tpP%xh(:) - swifter_plP%xh(:)
                    dv(:) = swifter_tpP%vh(:) - swifter_plP%vh(:)
                    radius = swifter_plP%radius
                    CALL discard_pl_close(dx(:), dv(:), dt, radius*radius, isp, r2min)
                    IF (isp /= 0) THEN
                         swifter_tpP%status = DISCARDED_PLR
                         WRITE(*, *) "Particle ", swifter_tpP%id, " too close to Planet ", swifter_plP%id, " at t = ", t
                         EXIT
                    END IF
               END DO
          END IF
          swifter_tpP => swifter_tpP%nextP
     END DO

     RETURN

END SUBROUTINE discard_pl
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
