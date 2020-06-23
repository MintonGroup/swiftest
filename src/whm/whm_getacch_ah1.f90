!**********************************************************************************************************************************
!
!  Unit Name   : whm_getacch_ah1
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Compute first term heliocentric accelerations of planets
!
!  Input
!    Arguments : npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!                ir3h     : inverse cubed heliocentric radii of planets
!                ir3j     : inverse cubed Jacobi radii of planets
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_getacch_ah1(npl, whm_pl1P, ir3h, ir3j)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_ah1.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_getacch_ah1(npl, whm_pl1P, ir3h, ir3j)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_getacch_ah1
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)           :: npl
     REAL(DP), DIMENSION(:), INTENT(IN) :: ir3h, ir3j
     TYPE(whm_pl), POINTER              :: whm_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: msun
     REAL(DP), DIMENSION(NDIM) :: ah1h, ah1j
     TYPE(whm_pl), POINTER     :: whm_plP

! Executable code
     msun = whm_pl1P%swifter%mass
     whm_pl1P%ah1(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     IF (npl > 1) THEN
          whm_plP => whm_pl1P%nextP
          whm_plP%ah1(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END IF
     DO i = 3, npl
          whm_plP => whm_plP%nextP
          ah1j(:) = whm_plP%xj(:)*ir3j(i)
          ah1h(:) = whm_plP%swifter%xh(:)*ir3h(i)
          whm_plP%ah1(:) = msun*(ah1j(:) - ah1h(:))
     END DO

     RETURN

END SUBROUTINE whm_getacch_ah1
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
