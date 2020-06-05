!**********************************************************************************************************************************
!
!  Unit Name   : util_hills
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Compute Hill sphere radii of planets
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_hills(npl, swifter_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine util_hills.f
!
!**********************************************************************************************************************************
SUBROUTINE util_hills(npl, swiftest_plA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => util_hills
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: msun, mp, mu, energy, ap, r, v2

! Executable code
     msun = swiftest_plA%mass(1)
     DO i = 2, npl
          mp = swiftest_plA%mass(i)
          IF (mp > 0.0_DP) THEN
               mu = msun + mp
               r = SQRT(DOT_PRODUCT(swiftest_plA%xh(:,i), swiftest_plA%xh(:,i)))
               v2 = DOT_PRODUCT(swiftest_plA%vh(:,i), swiftest_plA%vh(:,i))
               energy = 0.5_DP*v2 - mu/r
               ap = -0.5_DP*mu/energy
               swiftest_plA%rhill(i) = ap*(((mp/mu)/3.0_DP)**(1.0_DP/3.0_DP))
          ELSE
               swiftest_plA%rhill(i) = 0.0_DP
          END IF
     END DO

     RETURN

END SUBROUTINE util_hills
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
