!**********************************************************************************************************************************
!
!  Unit Name   : orbel_scget
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : orbel
!  Language    : Fortran 90/95
!
!  Description : Efficiently compute the sine and cosine of an input angle
!
!  Input
!    Arguments : angle : input angle
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : sx    : sine of the angle
!                cx    : cosine of the angle
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL orbel_scget(angle, sx, cx)
!
!  Notes       : Adapted from Martin Duncan's Swift routine orbel_scget.f
!
!                Input angle must be in radians
!
!**********************************************************************************************************************************
SUBROUTINE orbel_scget(angle, sx, cx)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => orbel_scget
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)  :: angle
     REAL(DP), INTENT(OUT) :: sx, cx

! Internals
     INTEGER(I4B) :: nper
     REAL(DP)     :: x

! Executable code
     nper = angle/TWOPI
     x = angle - nper*TWOPI
     IF (x < 0.0_DP) x = x + TWOPI
     sx = SIN(x)
     cx = SQRT(1.0_DP - sx*sx)
     IF ((x > PIBY2) .AND. (x < PI3BY2)) cx = -cx

     RETURN

END SUBROUTINE orbel_scget
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
