!**********************************************************************************************************************************
!
!  Unit Name   : discard_pl_close
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : discard
!  Language    : Fortran 90/95
!
!  Description : Check to see if a test particle and planet are having, or will have within the next time step, an encounter such
!                that the separation distance r is less than some critical radius rcrit (or r**2 < rcrit**2 = r2crit)
!
!  Input
!    Arguments : dx     : relative position of test particle with respect to planet
!                dv     : relative velocity of test particle with respect to planet
!                dt     : time step
!                r2crit : square of the boundary of the encounter region
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : iflag  : flag indicating encounter (0 = NO, 1 = YES)
!                r2min  : square of the smallest predicted separation distance
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_pl_close.f
!
!**********************************************************************************************************************************
SUBROUTINE discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => discard_pl_close
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT)             :: iflag
     REAL(DP), INTENT(IN)                  :: dt, r2crit
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: dx, dv
     REAL(DP), INTENT(OUT)                 :: r2min

! Internals
     REAL(DP) :: r2, v2, vdotr, tmin

! Executable code
     r2 = DOT_PRODUCT(dx(:), dx(:))
     IF (r2 <= r2crit) THEN
          iflag = 1
     ELSE
          vdotr = DOT_PRODUCT(dx(:), dv(:))
          IF (vdotr > 0.0_DP) THEN
               iflag = 0
          ELSE
               v2 = DOT_PRODUCT(dv(:), dv(:))
               tmin = -vdotr/v2
               IF (tmin < dt) THEN
                    r2min = r2 - vdotr*vdotr/v2
               ELSE
                    r2min = r2 + 2.0_DP*vdotr*dt + v2*dt*dt
               END IF
               r2min = MIN(r2min, r2)
               IF (r2min <= r2crit) THEN
                    iflag = 1
               ELSE
                    iflag = 0
               END IF
          END IF
     END IF

     RETURN

END SUBROUTINE discard_pl_close
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
