!**********************************************************************************************************************************
!
!  Unit Name   : coord_b2h
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from barycentric to heliocentric coordinates, planets only
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of Swifter planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL coord_b2h(npl, swifter_pl1P)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine coord_b2h.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_b2h(npl, swifter_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => coord_b2h
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     TYPE(swifter_pl), POINTER :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     xtmp(:) = swifter_pl1P%xb(:)
     vtmp(:) = swifter_pl1P%vb(:)
     swifter_plP => swifter_pl1P
     DO i = 1, npl
          swifter_plP%xh(:) = swifter_plP%xb(:) - xtmp(:)
          swifter_plP%vh(:) = swifter_plP%vb(:) - vtmp(:)
          swifter_plP => swifter_plP%nextP
     END DO

     RETURN

END SUBROUTINE coord_b2h
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
