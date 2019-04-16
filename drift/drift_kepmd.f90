!**********************************************************************************************************************************
!
!  Unit Name   : drift_kepmd
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Solve Kepler's equation in difference form for an ellipse for small input dm and eccentricity
!
!  Input
!    Arguments : dm : increment in mean anomaly
!                es : eccentricity times the sine of eccentric anomaly
!                ec : eccentricity times the cosine of eccentric anomaly
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : x  : solution to Kepler's equation in difference form (x = dE)
!                s  : sine of x
!                c  : cosine of x
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_kepmd(dm, es, ec, x, s, c)
!
!  Notes       : Adapted from Martin Duncan's Swift routine drift_kepmd.f
!
!                Original disclaimer: built for speed, does not check how well the original equation is solved
!                Can do that in calling routine by checking how close (x - ec*s + es*(1.0 - c) - dm) is to zero
!
!**********************************************************************************************************************************
SUBROUTINE drift_kepmd(dm, es, ec, x, s, c)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_kepmd
     IMPLICIT NONE

! Arguments
     REAL(DP), INTENT(IN)  :: dm, es, ec
     REAL(DP), INTENT(OUT) :: x, s, c

! Internals
     REAL(DP), PARAMETER :: A0 = 39916800.0_DP, A1 = 6652800.0_DP, A2 = 332640.0_DP, A3 = 7920.0_DP, A4 = 110.0_DP
     REAL(DP)            :: dx, fac1, fac2, q, y, f, fp, fpp, fppp

! Executable code
     fac1 = 1.0_DP/(1.0_DP - ec)
     q = fac1*dm
     fac2 = es*es*fac1 - ec/3.0_DP
     x = q*(1.0_DP - 0.5_DP*fac1*q*(es - q*fac2))
     y = x*x
     s = x*(A0 - y*(A1 - y*(A2 - y*(A3 - y*(A4 - y)))))/A0
     c = SQRT(1.0_DP - s*s)
     f = x - ec*s + es*(1.0_DP - c) - dm
     fp = 1.0_DP - ec*c + es*s
     fpp = ec*s + es*c
     fppp = ec*c - es*s
     dx = -f/fp
     dx = -f/(fp + dx*fpp/2.0_DP)
     dx = -f/(fp + dx*fpp/2.0_DP + dx*dx*fppp/6.0_DP)
     x = x + dx
     y = x*x
     s = x*(A0 - y*(A1 - y*(A2 - y*(A3 - y*(A4 - y)))))/A0
     c = SQRT(1.0_DP - s*s)

     RETURN

END SUBROUTINE drift_kepmd
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
