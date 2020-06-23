!**********************************************************************************************************************************
!
!  Unit Name   : whm_getacch_ah2
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Compute second term heliocentric accelerations of planets
!
!  Input
!    Arguments : npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!                ir3j     : inverse cubed Jacobi radii of planets
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_getacch_ah2(npl, whm_pl1P, ir3j)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_ah2.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_getacch_ah2(npl, whm_pl1P, ir3j)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_getacch_ah2
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)           :: npl
     REAL(DP), DIMENSION(:), INTENT(IN) :: ir3j
     TYPE(whm_pl), POINTER              :: whm_pl1P

! Internals
     INTEGER(I4B)          :: i
     REAL(DP)              :: etaj, fac, msun
     TYPE(whm_pl), POINTER :: whm_plP, whm_ploP

! Executable code
     msun = whm_pl1P%swifter%mass
     whm_pl1P%ah2(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     IF (npl > 1) THEN
          whm_plP => whm_pl1P%nextP
          whm_plP%ah2(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          etaj = msun
     ENDIF
     DO i = 3, npl
          whm_ploP => whm_plP
          whm_plP => whm_plP%nextP
          etaj = etaj + whm_ploP%swifter%mass
          fac = whm_plP%swifter%mass*msun*ir3j(i)/etaj
          whm_plP%ah2(:) = whm_ploP%ah2(:) + fac*whm_plP%xj(:)
     END DO

     RETURN

END SUBROUTINE whm_getacch_ah2
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
