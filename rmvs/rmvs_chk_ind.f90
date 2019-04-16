!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_chk_ind
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Determine whether a test particle and planet are having or will have an encounter within the next time step
!
!  Input
!    Arguments : xr     : relative position
!                vr     : relative velocity
!                dt     : time step
!                r2crit : square of the radius of the encounter region
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : iflag  : flag indicating encounter ( 1 ) or no encounter ( 0 )
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs_chk_ind.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_chk_ind(xr, vr, dt, r2crit, iflag)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_chk_ind
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)                  :: dt, r2crit
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: xr, vr
     INTEGER(I4B), INTENT(OUT)             :: iflag

! Internals
     REAL(DP) :: r2, v2, vdotr, tmin, r2min

! Executable code
     iflag = 0
     r2 = DOT_PRODUCT(xr(:), xr(:))
     IF (r2 < r2crit) THEN
          iflag = 1
     ELSE
          vdotr = DOT_PRODUCT(vr(:), xr(:))
          IF (vdotr < 0.0_DP) THEN
               v2 = DOT_PRODUCT(vr(:), vr(:))
               tmin = -vdotr/v2
               IF (tmin < dt) THEN
                    r2min = r2 - vdotr*vdotr/v2
               ELSE
                    r2min = r2 + 2.0_DP*vdotr*dt + v2*dt*dt
               END IF
               r2min = MIN(r2min, r2)
               IF (r2min <= r2crit) iflag = 1
          END IF
     END IF

     RETURN

END SUBROUTINE rmvs_chk_ind
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
