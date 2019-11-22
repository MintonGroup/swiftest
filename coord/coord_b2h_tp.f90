!**********************************************************************************************************************************
!
!  Unit Name   : coord_b2h_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from barycentric to heliocentric coordinates, active test particles only
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
!  Invocation  : CALL coord_b2h_tp(ntp, swifter_tp1P, swifter_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine coord_b2h_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_b2h_tp(ntp, swiftest_tpA, swiftest_plA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => coord_b2h_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: ntp
     TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp

! Executable code
     xtmp(:) = swiftest_plA%xb(:,1)
     vtmp(:) = swiftest_plA%vb(:,1)
     DO i = 1, ntp
          swiftest_tpA%xh(:,1) = swiftest_tpA%xb(:,1) - xtmp(:)
          swiftest_tpA%vh(:,1) = swiftest_tpA%vb(:,1) - vtmp(:)
     END DO

     RETURN

END SUBROUTINE coord_b2h_tp
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
