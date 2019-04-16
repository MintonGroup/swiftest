!**********************************************************************************************************************************
!
!  Unit Name   : drift_dan
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Perform Kepler drift, solving Kepler's equation in appropriate variables
!
!  Input
!    Arguments : mu    : G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
!                x0    : position of body to drift
!                v0    : velocity of body to drift
!                dt0   : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : x0    : position of body to drift
!                v0    : velocity of body to drift
!                iflag : error status flag for Kepler drift (0 = OK, nonzero = NO CONVERGENCE)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_dan(mu, x0, v0, dt0, iflag)
!
!  Notes       : Adapted from Hal Levison and Martin Duncan's Swift routine drift_dan.f
!
!**********************************************************************************************************************************
SUBROUTINE drift_dan(mu, x0, v0, dt0, iflag)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => drift_dan
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT)                :: iflag
     REAL(DP), INTENT(IN)                     :: mu, dt0
     REAL(DP), DIMENSION(NDIM), INTENT(INOUT) :: x0, v0

! Internals
     REAL(DP)                  :: dt, f, g, fdot, gdot, c1, c2, c3, u, alpha, fp, r0
     REAL(DP)                  :: v0s, a, asq, en, dm, ec, es, esq, xkep, fchk, s, c
     REAL(DP), DIMENSION(NDIM) :: x, v

! Executable code
     iflag = 0
     dt = dt0
     r0 = SQRT(DOT_PRODUCT(x0(:), x0(:)))
     v0s = DOT_PRODUCT(v0(:), v0(:))
     u = DOT_PRODUCT(x0(:), v0(:))
     alpha = 2.0_DP*mu/r0 - v0s
     IF (alpha > 0.0_DP) THEN
          a = mu/alpha
          asq = a*a
          en = SQRT(mu/(a*asq))
          ec = 1.0_DP - r0/a
          es = u/(en*asq)
          esq = ec*ec + es*es
          dm = dt*en - INT(dt*en/TWOPI)*TWOPI
          dt = dm/en
          IF ((esq < E2MAX) .AND. (dm*dm < DM2MAX) .AND. (esq*dm*dm < E2DM2MAX)) THEN
               CALL drift_kepmd(dm, es, ec, xkep, s, c)
               fchk = (xkep - ec*s + es*(1.0_DP - c) - dm)
! DEK - original code compared fchk*fchk with DANBYB, but I think it should
! DEK - be compared with DANBYB*DANBYB, and I changed it accordingly - please
! DEK - check with Hal and/or Martin about this
               IF (fchk*fchk > DANBYB*DANBYB) THEN
                    iflag = 1
                    RETURN
               END IF
               fp = 1.0_DP - ec*c + es*s
               f = a/r0*(c - 1.0_DP) + 1.0_DP
               g = dt + (s - xkep)/en
               fdot = -(a/(r0*fp))*en*s
               gdot = (c - 1.0_DP)/fp + 1.0_DP
               x(:) = x0(:)*f + v0(:)*g
               v(:) = x0(:)*fdot + v0(:)*gdot
               x0(:) = x(:)
               v0(:) = v(:)
               iflag = 0
               RETURN
          END IF
     END IF
     CALL drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
     IF (iflag == 0) THEN
          f = 1.0_DP - mu/r0*c2
          g = dt - mu*c3
          fdot = -mu/(fp*r0)*c1
          gdot = 1.0_DP - mu/fp*c2
          x(:) = x0(:)*f + v0(:)*g
          v(:) = x0(:)*fdot + v0(:)*gdot
          x0(:) = x(:)
          v0(:) = v(:)
     END IF

     RETURN

END SUBROUTINE drift_dan
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
