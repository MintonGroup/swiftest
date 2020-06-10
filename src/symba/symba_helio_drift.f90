!**********************************************************************************************************************************
!
!  Unit Name   : symba_helio_drift
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Loop through planets and call Danby drift routine
!
!  Input
!    Arguments : irec       : input recursion level
!                npl        : number of planets
!                symba_pl1P : pointer to head of SyMBA planet structure linked-list
!                dt         : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P : pointer to head of SyMBA planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_helio_drift(irec, npl, symba_pl1P, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_helio_drift.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_helio_drift(irec, npl, symba_plA, dt)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_helio_drift
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)       :: irec, npl
     REAL(DP), INTENT(IN)           :: dt
     TYPE(symba_pl), INTENT(INOUT)  :: symba_plA

! Internals
     INTEGER(I4B)              :: i, iflag
     REAL(DP)                  :: mu

! Executable code
     mu = symba_plA%helio%swiftest%mass(1)
     DO i = 2, npl
          IF ((symba_plA%levelg(i) == irec) .AND. (symba_plA%helio%swiftest%status(i) == ACTIVE)) THEN
               CALL drift_one(mu, symba_plA%helio%swiftest%xh(:,i), symba_plA%helio%swiftest%vb(:,i), dt, iflag)
               IF (iflag /= 0) THEN
                    WRITE(*, *) " Planet ", symba_plA%helio%swiftest%name(i), " is lost!!!!!!!!!!"
                    WRITE(*, *) mu, dt
                    WRITE(*, *) symba_plA%helio%swiftest%xh(:,i)
                    WRITE(*, *) symba_plA%helio%swiftest%vb(:,i)
                    WRITE(*, *) " STOPPING "
                    CALL util_exit(FAILURE)
               END IF
          END IF
     END DO

     RETURN

END SUBROUTINE symba_helio_drift
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
