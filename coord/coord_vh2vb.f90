!**********************************************************************************************************************************
!
!  Unit Name   : coord_vh2vb
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from heliocentric to barycentric coordinates, planet velocities only
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
!  Invocation  : CALL coord_vh2vb(npl, swifter_pl1P, msys)
!
!  Notes       : Adapted from Hal Levison's Swift routine coord_vh2b.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_vh2vb(npl, swiftest_plA, msys)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => coord_vh2vb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)         :: npl
     REAL(DP), INTENT(OUT)            :: msys
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: vtmp

! Executable code
     ! Removed by D. Minton
     !swifter_plP => swifter_pl1P
     !^^^^^^^^^^^^^^^^^^^^^
     vtmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     msys = swiftest_plA%mass(1)
     !^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton

! EDIT FOR PARALLELIZATION

     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_plP) &
     !$OMP SHARED(npl,swifter_pl1P) &
     !$OMP REDUCTION(+:vtmp,msys)
     DO i = 2, npl
          ! Removed by D. Minton
          !swifter_plP => swifter_plP%nextP
          !^^^^^^^^^^^^^^^^^^^^^
          ! Added by D. Minton
          !swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^^
          msys = msys + swiftest_plA%mass(i)
          vtmp(:) = vtmp(:) + swiftest_plA%mass(i)*swiftest_plA%vh(:,i)
     END DO
     !$OMP END PARALLEL DO
     swiftest_plA%vb(:,1) = -vtmp(:)/msys
     vtmp(:) = swiftest_plA%vb(:,1)
     !^^^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_plP) &
     !$OMP SHARED(npl,swifter_pl1P,vtmp) 
     DO i = 2, npl
          ! Removed by D. Minton
          !swifter_plP => swifter_plP%nextP
          !^^^^^^^^^^^^^^^^^^^^^
          ! Added by D. Minton
          !swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^^^^
          swiftest_plA%vb(:,i) = swiftest_plA%vh(:,i) + vtmp(:)
     END DO
     !$OMP END PARALLEL DO

     RETURN

END SUBROUTINE coord_vh2vb
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann (Checked by Jennifer Pouplin & Carlisle Wishard)
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
