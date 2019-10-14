!**********************************************************************************************************************************
!
!  Unit Name   : coord_vb2vh
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE coord_vb2vh(npl, swifter_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_random_access, EXCEPT_THIS_ONE => coord_vb2vh
     !USE module_interfaces, EXCEPT_THIS_ONE => coord_vb2vh
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     TYPE(swifter_pl), POINTER :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP), DIMENSION(NDIM) :: vtmp = 0.0_DP
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(AUTO) & 
     !$OMP SHARED(npl) &
     !$OMP REDUCTION(-:vtmp)
     DO i = 2, npl
          CALL get_point(i,swifter_plP)
          vtmp(:) = vtmp(:) -  swifter_plP%vb(:) * swifter_plP%mass 
     END DO
     !$OMP END PARALLEL DO
     vtmp(:) = vtmp(:) / swifter_pl1P%mass
     swifter_pl1P%vb(:) = vtmp(:)
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(AUTO) & 
     !$OMP SHARED(npl,vtmp) 
     DO i = 2, npl
          CALL get_point(i,swifter_plP)
          swifter_plP%vh(:) = swifter_plP%vb(:) - vtmp(:)
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
