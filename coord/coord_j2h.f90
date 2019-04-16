!**********************************************************************************************************************************
!
!  Unit Name   : coord_j2h
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from Jacobi to heliocentric coordinates, planets only
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
!  Invocation  : CALL coord_j2h(npl, whm_pl1P)
!
!  Notes       : Adapted from Martin Duncan's Swift routine coord_j2h.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_j2h(npl, whm_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => coord_j2h
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     TYPE(whm_pl), POINTER    :: whm_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: eta
     REAL(DP), DIMENSION(NDIM) :: sum, sumv
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(whm_pl), POINTER     :: whm_plP

! Executable code
     eta = whm_pl1P%swifter%mass
     whm_pl1P%eta = eta
     whm_plP => whm_pl1P
     DO i = 2, npl
          whm_plP => whm_plP%nextP
          eta = eta + whm_plP%swifter%mass
          whm_plP%eta = eta
     END DO
     whm_plP => whm_pl1P
     swifter_plP => whm_plP%swifter
     swifter_plP%xh(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     swifter_plP%vh(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     whm_plP => whm_plP%nextP
     swifter_plP => whm_plP%swifter
     swifter_plP%xh(:) = whm_plP%xj(:)
     swifter_plP%vh(:) = whm_plP%vj(:)
     sum(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     sumv(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     DO i = 3, npl
          sum(:) = sum(:) + swifter_plP%mass*whm_plP%xj(:)/whm_plP%eta
          sumv(:) = sumv(:) + swifter_plP%mass*whm_plP%vj(:)/whm_plP%eta
          whm_plP => whm_plP%nextP
          swifter_plP => whm_plP%swifter
          swifter_plP%xh(:) = whm_plP%xj(:) + sum(:)
          swifter_plP%vh(:) = whm_plP%vj(:) + sumv(:)
     END DO

     RETURN

END SUBROUTINE coord_j2h
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
