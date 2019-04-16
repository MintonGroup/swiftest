!**********************************************************************************************************************************
!
!  Unit Name   : util_hills
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE util_hills(npl, swifter_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => util_hills
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     TYPE(swifter_pl), POINTER :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: msun, mp, mu, energy, ap, r, v2
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     msun = swifter_pl1P%mass
     swifter_plP => swifter_pl1P
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          mp = swifter_plP%mass
          IF (mp > 0.0_DP) THEN
               mu = msun + mp
               r = SQRT(DOT_PRODUCT(swifter_plP%xh(:), swifter_plP%xh(:)))
               v2 = DOT_PRODUCT(swifter_plP%vh(:), swifter_plP%vh(:))
               energy = 0.5_DP*v2 - mu/r
               ap = -0.5_DP*mu/energy
               swifter_plP%rhill = ap*(((mp/mu)/3.0_DP)**(1.0_DP/3.0_DP))
          ELSE
               swifter_plP%rhill = 0.0_DP
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
