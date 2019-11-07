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
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => coord_h2b
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
     !Removed by D. minton
     !swifter_plP => swifter_pl1P
     !^^^^^^^^^^^^^^^^^^^^^
     xtmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     vtmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     !^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_plP) &
     !$OMP SHARED(npl,swifter_pl1P) &
     !$OMP REDUCTION(+:msys,xtmp,vtmp)     
     DO i = 2, npl
          ! Removed by D. Minton
          !swifter_plP => swifter_plP%nextP
          !^^^^^^^^^^^^^^^^^^^^^
          ! Added by D. Minton
          !swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^^
          msys = msys + swiftest_plA%mass(i)
          xtmp(:) = xtmp(:) + swiftest_plA%mass(i)*swiftest_plA%xh(:,i)
          vtmp(:) = vtmp(:) + swiftest_plA%mass(i)*swiftest_plA%vh(:,i)
     END DO
     !$OMP END PARALLEL DO
     swiftest_plA%xb(:,1) = -xtmp(:)/msys
     swiftest_plA%vb(:,1) = -vtmp(:)/msys
     xtmp(:) = swiftest_plA%xb(:,1)
     vtmp(:) = swiftest_plA%vb(:,1)
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_plP) &
     !$OMP SHARED(npl,swifter_pl1P,xtmp,vtmp) 
     DO i = 2, npl
          !Removed by D. Minton
          !swifter_plP => swifter_plP%nextP
          !^^^^^^^^^^^^^^^^^^
          !Added by D. Minton
          !swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^
          swiftest_plA%xb(:,i) = swiftest_plA%xh(:,i) + xtmp(:)
          swiftest_plA%vb(:,i) = swiftest_plA%vh(:,i) + vtmp(:)
     END DO
     !$OMP END PARALLEL DO

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
