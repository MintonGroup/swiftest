!**********************************************************************************************************************************
!
!  Unit Name   : coord_b2h
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from barycentric to heliocentric coordinates, planets only
!
!  Input
!    Arguments : npl          : number of planets
!                swiftest_pl1P : pointer to head of Swiftest planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swiftest_pl1P : pointer to head of Swiftest planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL coord_b2h(npl, swiftest_pl1P)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine coord_b2h.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_b2h(npl, swiftest_plA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => coord_b2h
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp

! Executable code
     xtmp(:) = swiftest_plA%xb(:,1)
     vtmp(:) = swiftest_plA%vb(:,1)
     DO i = 1, npl
          swiftest_plA%xh(:,i) = swiftest_plA%xb(:,i) - xtmp(:)
          swiftest_plA%vh(:,i) = swiftest_plA%vb(:,i) - vtmp(:)
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
