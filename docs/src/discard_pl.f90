!**********************************************************************************************************************************
!
!  Unit Name   : discard_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
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
SUBROUTINE discard_pl(t, dt, npl, ntp, swiftest_plA, swiftest_tpA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => discard_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)         :: npl, ntp
     REAL(DP), INTENT(IN)             :: t, dt
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
     TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA

! Internals
     INTEGER(I4B)              :: i, j, isp
     REAL(DP)                  :: r2min, radius
     REAL(DP), DIMENSION(NDIM) :: dx, dv

! Executable code
     DO i = 1, ntp
          IF (swiftest_tpA%status(i) == ACTIVE) THEN
               DO j = 2, npl
                    dx(:) = swiftest_tpA%xh(:,i) - swiftest_plA%xh(:,i)
                    dv(:) = swiftest_tpA%vh(:,i) - swiftest_plA%vh(:,i)
                    radius = swiftest_plA%radius(i)
                    CALL discard_pl_close(dx(:), dv(:), dt, radius*radius, isp, r2min)
                    IF (isp /= 0) THEN
                         swiftest_tpA%status(i) = DISCARDED_PLR
                         ldiscard = .TRUE.
                         WRITE(*, *) "Particle ", swiftest_tpA%name(i), " too close to Planet ", swiftest_plA%name(i), " at t = ", t
                         ldiscard_tp = .TRUE.
                         EXIT
                    END IF
               END DO
          END IF
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
