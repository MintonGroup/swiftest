!**********************************************************************************************************************************
!
!  Unit Name   : coord_h2b
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE coord_h2b(npl, swifter_pl1P, msys)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => coord_h2b
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     REAL(DP), INTENT(OUT)     :: msys
     TYPE(swifter_pl), POINTER :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     msys = swifter_pl1P%mass
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
          swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^^
          msys = msys + swifter_plP%mass
          xtmp(:) = xtmp(:) + swifter_plP%mass*swifter_plP%xh(:)
          vtmp(:) = vtmp(:) + swifter_plP%mass*swifter_plP%vh(:)
     END DO
     !$OMP END PARALLEL DO
     swifter_plP => swifter_pl1P
     swifter_plP%xb(:) = -xtmp(:)/msys
     swifter_plP%vb(:) = -vtmp(:)/msys
     xtmp(:) = swifter_plP%xb(:)
     vtmp(:) = swifter_plP%vb(:)
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) & 
     !$OMP PRIVATE(i,swifter_plP) &
     !$OMP SHARED(npl,swifter_pl1P,xtmp,vtmp) 
     DO i = 2, npl
          !Removed by D. Minton
          !swifter_plP => swifter_plP%nextP
          !^^^^^^^^^^^^^^^^^^
          !Added by D. Minton
          swifter_plP => swifter_pl1P%swifter_plPA(i)%thisP
          !^^^^^^^^^^^^^^^^^^
          swifter_plP%xb(:) = swifter_plP%xh(:) + xtmp(:)
          swifter_plP%vb(:) = swifter_plP%vh(:) + vtmp(:)
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
