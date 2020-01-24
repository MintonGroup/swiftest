!**********************************************************************************************************************************
!
!  Unit Name   : obl_acc
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : obl
!  Language    : Fortran 90/95
!
!  Description : Compute the barycentric accelerations of planets due to the oblateness of the central body
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                xh           : heliocentric positions of planets
!                irh          : inverse heliocentric radii of planets
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : aobl         : barycentric accelerations of planets due to central body oblateness
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
!
!  Notes       : Adapted from Martin Duncan's Swift routine obl_acc.f
!
!                Returned values do not include monopole term or terms higher than J4
!
!**********************************************************************************************************************************
SUBROUTINE obl_acc(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, aobl)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => obl_acc
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                     :: npl
     REAL(DP), INTENT(IN)                         :: j2rp2, j4rp4
     REAL(DP), DIMENSION(npl), INTENT(IN)         :: irh
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN)   :: xh
     REAL(DP), DIMENSION(NDIM, npl), INTENT(OUT)  :: aobl
     TYPE(swiftest_pl), INTENT(INOUT)             :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: rinv2, t0, t1, t2, t3, fac1, fac2, msun

! Executable code
     msun = swiftest_plA%mass(1)
     DO i = 2, npl
          rinv2 = irh(i)**2
          t0 = -msun*rinv2*rinv2*irh(i)
          t1 = 1.5_DP*j2rp2
          t2 = xh(3, i)*xh(3, i)*rinv2
          t3 = 1.875_DP*j4rp4*rinv2
          fac1 = t0*(t1 - t3 - (5.0_DP*t1 - (14.0_DP - 21.0_DP*t2)*t3)*t2)
          fac2 = 2.0_DP*t0*(t1 - (2.0_DP - (14.0_DP*t2/3.0_DP))*t3)
          aobl(:, i) = fac1*xh(:, i)
          aobl(3, i) = fac2*xh(3, i) + aobl(3, i)
     END DO
     aobl(:, 1) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     DO i = 2, npl
          aobl(:, 1) = aobl(:, 1) - swiftest_plA%mass(i)*aobl(:, i)/msun
     END DO

     RETURN

END SUBROUTINE obl_acc
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
