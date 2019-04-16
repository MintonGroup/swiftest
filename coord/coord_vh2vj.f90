!**********************************************************************************************************************************
!
!  Unit Name   : coord_vh2vj
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from heliocentric to Jacobi coordinates, planet velocities only
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
!  Invocation  : CALL coord_vh2vj(npl, whm_pl1P)
!
!  Notes       : Adapted from Martin Duncan's Swift routine coord_vh2vj.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_vh2vj(npl, whm_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => coord_vh2vj
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     TYPE(whm_pl), POINTER    :: whm_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: eta
     REAL(DP), DIMENSION(NDIM) :: sumv, capv
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
     whm_plP%vj(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     whm_plP => whm_plP%nextP
     swifter_plP => whm_plP%swifter
     whm_plP%vj(:) = swifter_plP%vh(:)
     sumv(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     DO i = 3, npl
          sumv(:) = sumv(:) + swifter_plP%mass*swifter_plP%vh(:)
          capv(:) = sumv(:)/whm_plP%eta
          whm_plP => whm_plP%nextP
          swifter_plP => whm_plP%swifter
          whm_plP%vj(:) = swifter_plP%vh(:) - capv(:)
     END DO

     RETURN

END SUBROUTINE coord_vh2vj
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
