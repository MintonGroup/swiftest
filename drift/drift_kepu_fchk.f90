!**********************************************************************************************************************************
!
!  Unit Name   : drift_kepu_fchk
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Computes the value of f, the function whose root we are trying to find in universal variables
!
!  Input
!    Arguments : dt    : time step
!                r0    : distance between two bodies
!                mu    : G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
!                alpha : twice the binding energy
!                u     : dot product of position and velocity vectors
!                s     : universal variable (approximate root of f)
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : f     : function value
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
!
!  Notes       : Adapted from Martin Duncan's Swift routine drift_kepu_fchk.f
!
!**********************************************************************************************************************************
SUBROUTINE drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_kepu_fchk
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u, s
     REAL(DP), INTENT(OUT) :: f

! Internals
     REAL(DP) :: x, c0, c1, c2, c3

! Executable code
     x = s*s*alpha
     CALL drift_kepu_stumpff(x, c0, c1, c2, c3)
     c1 = c1*s
     c2 = c2*s*s
     c3 = c3*s*s*s
     f = r0*c1 + u*c2 + mu*c3 - dt

     RETURN

END SUBROUTINE drift_kepu_fchk
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
