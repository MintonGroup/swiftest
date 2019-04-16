!**********************************************************************************************************************************
!
!  Unit Name   : symba_reorder_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Rearrange SyMBA planet structure linked-list in order of decreasing mass
!
!  Input
!    Arguments : npl        : number of planets
!                symba_pl1P : pointer to head of SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P : pointer to head of SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_reorder_pl(npl, symba_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_reorder_pl(npl, symba_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_reorder_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     TYPE(symba_pl), POINTER  :: symba_pl1P

! Internals
     INTEGER(I4B)                              :: i
     INTEGER(I4B), DIMENSION(:), ALLOCATABLE   :: index
     REAL(DP), DIMENSION(:), ALLOCATABLE       :: mass
     TYPE(symba_pl), DIMENSION(:), ALLOCATABLE :: symba_plwkspA
     TYPE(symba_pl), POINTER                   :: symba_plP

! Executable code
     ALLOCATE(index(npl), mass(npl), symba_plwkspA(npl))
     symba_plP => symba_pl1P
     DO i = 1, npl
          symba_plwkspA(i) = symba_plP
          mass(i) = symba_plP%helio%swifter%mass
          symba_plP => symba_plP%nextP
     END DO
     CALL util_index(mass, index)
     symba_plP => symba_pl1P
     DO i = 1, npl
          symba_plP%helio%swifter%id = symba_plwkspA(index(npl-i+1))%helio%swifter%id
          symba_plP%helio%swifter%status = symba_plwkspA(index(npl-i+1))%helio%swifter%status
          symba_plP%helio%swifter%mass = symba_plwkspA(index(npl-i+1))%helio%swifter%mass
          symba_plP%helio%swifter%radius = symba_plwkspA(index(npl-i+1))%helio%swifter%radius
          symba_plP%helio%swifter%rhill = symba_plwkspA(index(npl-i+1))%helio%swifter%rhill
          symba_plP%helio%swifter%xh(:) = symba_plwkspA(index(npl-i+1))%helio%swifter%xh(:)
          symba_plP%helio%swifter%vh(:) = symba_plwkspA(index(npl-i+1))%helio%swifter%vh(:)
          symba_plP => symba_plP%nextP
     END DO
     IF (ALLOCATED(symba_plwkspA)) DEALLOCATE(symba_plwkspA)
     IF (ALLOCATED(mass)) DEALLOCATE(mass)
     IF (ALLOCATED(index)) DEALLOCATE(index)

     RETURN

END SUBROUTINE symba_reorder_pl
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
