!**********************************************************************************************************************************
!
!  Unit Name   : symba_chk
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check for an encounter
!
!  Input
!    Arguments : xr         : position of body 2 relative to body 1
!                vr         : velocity of body 2 relative to body 1
!                rhill1     : Hill sphere radius of body 1
!                rhill2     : Hill sphere radius of body 2
!                dt         : time step
!                irec       : recursion level
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lencounter : logical flag indicating whether there is an encounter this time step
!                lvdotr     : logical flag indicating whether the two bodies are approaching
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_chk.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_chk
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(OUT)          :: lencounter, lvdotr
     INTEGER(I4B), INTENT(IN)           :: irec
     REAL(DP), INTENT(IN)               :: rhill1, rhill2, dt
     REAL(DP), DIMENSION(:), INTENT(IN) :: xr, vr

! Internals
     ! LOGICAL(LGT) :: iflag lvdotr_flag
     REAL(DP)     :: rcrit, r2crit, vdotr, r2, v2, tmin, r2min

! Executable code
     lencounter = .FALSE.
     lvdotr = .FALSE.

     rcrit = (rhill1 + rhill2)*RHSCALE*(RSHELL**(irec)) 
     r2crit = rcrit*rcrit 

     r2 = DOT_PRODUCT(xr(:), xr(:)) 
     vdotr = DOT_PRODUCT(vr(:), xr(:))

     lvdotr = (vdotr < 0.0_DP)

     IF (r2 < r2crit) THEN
          lencounter = .TRUE.
     ELSE
          IF (vdotr < 0.0_DP) THEN
               v2 = DOT_PRODUCT(vr(:), vr(:))
               tmin = -vdotr/v2
               IF (tmin < dt) THEN
                    r2min = r2 - vdotr*vdotr/v2
               ELSE
                    r2min = r2 + 2.0_DP*vdotr*dt + v2*dt*dt
               END IF
               r2min = MIN(r2min, r2)
               IF (r2min <= r2crit) lencounter = .TRUE.
          END IF
     END IF


     RETURN

END SUBROUTINE symba_chk
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
