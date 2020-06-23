!**********************************************************************************************************************************
!
!  Unit Name   : whm_getacch_ah3
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Compute direct cross (third) term heliocentric accelerations of planets
!
!  Input
!    Arguments : npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_getacch_ah3(npl, whm_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_ah3.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_getacch_ah3(npl, whm_pl1P)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_getacch_ah3
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     TYPE(whm_pl), POINTER    :: whm_pl1P

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: rji2, irij3, faci, facj
     REAL(DP), DIMENSION(NDIM) :: dx
     TYPE(whm_pl), POINTER     :: whm_pliP, whm_pljP

! Executable code
     whm_pliP => whm_pl1P
     DO i = 1, npl
          whm_pliP%ah3(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          whm_pliP => whm_pliP%nextP
     END DO
     whm_pliP => whm_pl1P
     DO i = 2, npl - 1
          whm_pliP => whm_pliP%nextP
          whm_pljP => whm_pliP
          DO j = i + 1, npl
               whm_pljP => whm_pljP%nextP
               dx(:) = whm_pljP%swifter%xh(:) - whm_pliP%swifter%xh(:)
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP/(rji2*SQRT(rji2))
               faci = whm_pliP%swifter%mass*irij3
               facj = whm_pljP%swifter%mass*irij3
               whm_pliP%ah3(:) = whm_pliP%ah3(:) + facj*dx(:)
               whm_pljP%ah3(:) = whm_pljP%ah3(:) - faci*dx(:)
          END DO
     END DO

     RETURN

END SUBROUTINE whm_getacch_ah3
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
