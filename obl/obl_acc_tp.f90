!**********************************************************************************************************************************
!
!  Unit Name   : obl_acc_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : obl
!  Language    : Fortran 90/95
!
!  Description : Compute the barycentric accelerations of test particles due to the oblateness of the central body
!
!  Input
!    Arguments : ntp   : number of active test particles
!                xht   : heliocentric positions of test particles
!                j2rp2 : J2 * R**2 for the Sun
!                j4rp4 : J4 * R**4 for the Sun
!                irht  : inverse heliocentric radii of test particles
!                msun  : mass of the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : aoblt : barycentric accelerations of test particles due to central body oblateness
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)
!
!  Notes       : Adapted from Hal Levison's Swift routine obl_acc_tp.f
!
!                Returned values do not include monopole term or terms higher than J4
!
!**********************************************************************************************************************************
SUBROUTINE obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => obl_acc_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                    :: ntp
     REAL(DP), INTENT(IN)                        :: j2rp2, j4rp4, msun
     REAL(DP), DIMENSION(ntp), INTENT(IN)        :: irht
     REAL(DP), DIMENSION(NDIM, ntp), INTENT(IN)  :: xht
     REAL(DP), DIMENSION(NDIM, ntp), INTENT(OUT) :: aoblt

! Internals
     INTEGER(I4B) :: i
     REAL(DP)     :: rinv2, t0, t1, t2, t3, r2, fac1, fac2

! Executable code
     DO i = 1, ntp
          rinv2 = irht(i)**2
          t0 = -msun*rinv2*rinv2*irht(i)
          t1 = 1.5_DP*j2rp2
          t2 = xht(3, i)*xht(3, i)*rinv2
          t3 = 1.875_DP*j4rp4*rinv2
          fac1 = t0*(t1 - t3 - (5.0_DP*t1 - (14.0_DP - 21.0_DP*t2)*t3)*t2)
          fac2 = 2.0_DP*t0*(t1 - (2.0_DP - (14.0_DP*t2/3.0_DP))*t3)
          aoblt(:, i) = fac1*xht(:, i)
          aoblt(3, i) = aoblt(3, i) + fac2*xht(3, i)
     END DO

     RETURN

END SUBROUTINE obl_acc_tp
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
