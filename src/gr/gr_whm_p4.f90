!**********************************************************************************************************************************
!
!  Unit Name   : gr_whm_p4
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Position kick due to p**4 term in the post-Newtonian correction
!
!  Input
!    Arguments : npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!                dt       : time step
!                c2       : inverse speed of light squared
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL gr_whm_p4(npl, whm_pl1P, dt, c2)
!
!  Notes       : Based on Saha & Tremaine (1994) Eq. 28
!
!**********************************************************************************************************************************
SUBROUTINE gr_whm_p4(npl, whm_pl1P, dt, c2)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => gr_whm_p4
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     TYPE(whm_pl), POINTER    :: whm_pl1P
     REAL(DP), INTENT(IN)     :: dt, c2

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: dr
     TYPE(whm_pl), POINTER     :: whm_plP
     REAL(DP)                  :: vjmag2

! Executable code
     whm_plP => whm_pl1P
     DO i = 2, npl
          whm_plP => whm_plP%nextP
          vjmag2 = DOT_PRODUCT(whm_plP%vj(:), whm_plP%vj(:))
          dr(:) = - c2 * whm_plP%vj(:) * vjmag2
          whm_plP%xj(:) = whm_plP%xj(:) + dr(:)*dt
     END DO

     RETURN

END SUBROUTINE gr_whm_p4
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton
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
