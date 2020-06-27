!**********************************************************************************************************************************
!
!  Unit Name   : gr_vel2pseudovel
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Converts the heliocentric velocity into a pseudovelocity with
!  relativistic corrections. Uses Newton-Raphson method with direct inversion
!  of the Jacobian (yeah, it's slow, but this is only done once per run).
!
!  Input
!    Arguments : xh, vh  : heliocentric velocity and position vectors
!                c2      : inverse speed of light squared
!                mu      : Gravitational constant G*(Msun+Mp)
!    Terminal  : none
!
!  Output
!    Arguments : pvh      : pseudovelocity 
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL gr_vel2pseudovel(xh,vh,mu,c2,pvh)
!
!  Notes       : Based on Saha & Tremaine (1994) Eq. 32
!
!
!**********************************************************************************************************************************
SUBROUTINE gr_vel2pseudovel(xh, vh, mu, c2, pvh)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => gr_vel2pseudovel
     IMPLICIT NONE

! Arguments
     REAL(DP), DIMENSION(NDIM), INTENT(IN)  :: xh,vh
     REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: pvh
     REAL(DP), INTENT(IN)                   :: mu, c2

! Internals
     REAL(DP)                       :: G,pv2,rterm,det
     REAL(DP), DIMENSION(NDIM,NDIM) :: J,Jinv
     REAL(DP), DIMENSION(NDIM)      :: F
     INTEGER(I4B)                   :: n,i,k
     INTEGER(I4B)                   :: maxiter=50
     REAL(DP),PARAMETER             :: tol=1.0e-12_DP

     pvh = vh ! Initial guess
     rterm = 3*mu/SQRT(DOT_PRODUCT(xh, xh))
     DO n = 1, maxiter
        pv2 = DOT_PRODUCT(pvh, pvh)
        G = 1.0_DP - c2 *(0.5_DP*pv2 + rterm)
        F(:) = pvh(:) * G - vh(:)
        if (ABS(SUM(F)/DOT_PRODUCT(vh, vh))<tol) exit ! Root found

        ! Calculate the Jacobian
        DO k = 1, NDIM
            DO i = 1, NDIM
               IF (i==k) THEN
                  J(i,k) = G - c2*pvh(k)
               ELSE
                  J(i,k) = -c2*pvh(k)
               END IF
            END DO
        END DO

        ! Inverse of the Jacobian
        det = J(1,1)*(J(3,3)*J(2,2)-J(3,2)*J(2,3))
        det = det - J(2,1)*(J(3,3)*J(1,2)-J(3,2)*J(1,3))
        det = det + J(3,1)*(J(2,3)*J(1,2)-J(2,2)*J(1,3))

        Jinv(1,1) =   J(3,3)*J(2,2)-J(3,2)*J(2,3)
        Jinv(1,2) = -(J(3,3)*J(1,2)-J(3,2)*J(1,3))
        Jinv(1,3) =   J(2,3)*J(1,2)-J(2,2)*J(1,3)

        Jinv(2,1) = -(J(3,3)*J(2,1)-J(3,1)*J(2,3))
        Jinv(2,2) =   J(3,3)*J(1,1)-J(3,1)*J(1,3)
        Jinv(2,3) = -(J(2,3)*J(1,1)-J(2,1)*J(1,3))

        Jinv(3,1) =   J(3,2)*J(2,1)-J(3,1)*J(2,2)
        Jinv(3,2) = -(J(3,2)*J(1,1)-J(3,1)*J(1,2))
        Jinv(3,3) =   J(2,2)*J(1,1)-J(2,1)*J(1,2)

        Jinv = Jinv * det

        DO i = 1, NDIM
           pvh(i) = pvh(i) - DOT_PRODUCT(Jinv(:, i) ,F(:))
        END DO
        
     END DO 

     RETURN

END SUBROUTINE gr_vel2pseudovel
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
