!**********************************************************************************************************************************
!
!  Unit Name   : coord_vb2vh
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from barycentric to heliocentric coordinates, planet velocities only
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
!  Invocation  : CALL coord_vb2vh(npl, swifter_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine coord_vb2h.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_vb2vh(npl, swiftest_plA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => coord_vb2vh
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: vtmp

! Executable code
     vtmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     !Removed by D. Minton
     !swifter_plP => swifter_pl1P
     !^^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_plP) &
     !$OMP SHARED(npl,swifter_pl1P) &
     !$OMP REDUCTION(+:vtmp)
     DO i = 2, npl
          !Removed by D. Minton
          !swifter_plP => swifter_plP%nextP
          !^^^^^^^^^^^^^^^^^^^^
          !Added by D. Minton
          !swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^
          vtmp(:) = vtmp(:) - swiftest_plA%mass(i)*swiftest_plA%vb(:,i)
     END DO
     !$OMP END PARALLEL DO
     vtmp(:) = vtmp(:)/swiftest_plA%mass(1)
     swiftest_plA%vb(:,1) = vtmp(:)
     !Removed by D. Minton
     !swifter_plP => swifter_pl1P
     !^^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_plP) &
     !$OMP SHARED(npl,swifter_pl1P,vtmp) 
     DO i = 2, npl
          !Removed by D. Minton
          !swifter_plP => swifter_plP%nextP
          !^^^^^^^^^^^^^^^^^^^^
          !Added by D. Minton
          !swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^
          swiftest_plA%vh(:,i) = swiftest_plA%vb(:,i) - vtmp(:)
     END DO
     !$OMP END PARALLEL DO

     RETURN

END SUBROUTINE coord_vb2vh
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
