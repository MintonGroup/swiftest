!**********************************************************************************************************************************
!
!  Unit Name   : coord_h2b_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from heliocentric to barycentric coordinates, active test particles only
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL coord_h2b_tp(ntp, swifter_tp1P, swifter_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine coord_h2b_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_h2b_tp(ntp, swifter_tp1P, swifter_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => coord_h2b_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: ntp
     TYPE(swifter_tp), POINTER :: swifter_tp1P
     TYPE(swifter_pl), POINTER :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     xtmp(:) = swifter_pl1P%xb(:)
     vtmp(:) = swifter_pl1P%vb(:)
     swifter_tpP => swifter_tp1P
     DO i = 1, ntp
          swifter_tpP%xb(:) = swifter_tpP%xh(:) + xtmp(:)
          swifter_tpP%vb(:) = swifter_tpP%vh(:) + vtmp(:)
          swifter_tpP => swifter_tpP%nextP
     END DO

     RETURN

END SUBROUTINE coord_h2b_tp
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
