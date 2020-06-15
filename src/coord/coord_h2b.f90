!**********************************************************************************************************************************
!
!  Unit Name   : coord_h2b
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from heliocentric to barycentric coordinates, planets only
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                msys         : total system mass
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL coord_h2b(npl, swifter_pl1P, msys)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine coord_h2b.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_h2b(npl, swiftest_plA, msys)

! Modules
     USE swiftest, EXCEPT_THIS_ONE => coord_h2b
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     REAL(DP), INTENT(OUT)     :: msys
     TYPE(swiftest_pl),INTENT(INOUT) :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp

! Executable code
     msys = swiftest_plA%mass(1)
     xtmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     vtmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     DO i = 2, npl
          msys = msys + swiftest_plA%mass(i)
          xtmp(:) = xtmp(:) + swiftest_plA%mass(i)*swiftest_plA%xh(:,i)
          vtmp(:) = vtmp(:) + swiftest_plA%mass(i)*swiftest_plA%vh(:,i)
     END DO
     swiftest_plA%xb(:,1) = -xtmp(:)/msys                                  
     swiftest_plA%vb(:,1) = -vtmp(:)/msys                                  
     xtmp(:) = swiftest_plA%xb(:,1)
     vtmp(:) = swiftest_plA%vb(:,1)
     DO i = 2, npl
          swiftest_plA%xb(:,i) = swiftest_plA%xh(:,i) + xtmp(:)
          swiftest_plA%vb(:,i) = swiftest_plA%vh(:,i) + vtmp(:)
     END DO

     RETURN

END SUBROUTINE coord_h2b
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
