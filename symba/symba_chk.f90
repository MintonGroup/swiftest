!**********************************************************************************************************************************
!
!  Unit Name   : symba_chk
!  Unit Type   : subroutine
!  Project     : Swifter
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
     USE module_swifter
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
     INTEGER(I4B) :: iflag
     REAL(DP)     :: rcrit, r2crit, vdotr

! Executable code
     lencounter = .FALSE.
     rcrit = (rhill1 + rhill2)*RHSCALE*(RSHELL**(irec))
     r2crit = rcrit*rcrit
     CALL rmvs_chk_ind(xr(:), vr(:), dt, r2crit, iflag)
     IF (iflag /= 0) lencounter = .TRUE.
     vdotr = DOT_PRODUCT(vr(:), xr(:))
     lvdotr = (vdotr < 0.0_DP)

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
