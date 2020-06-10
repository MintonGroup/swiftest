!**********************************************************************************************************************************
!
!  Unit Name   : symba_reorder_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Rearrange SyMBA planet arrays in order of decreasing mass
!
!  Input
!    Arguments : npl        : number of planets
!                symba_plA  : structure of arrays of SyMBA planet 
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_plA : structure of arrays of SyMBA planet 
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_reorder_pl(npl, symba_plA)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_reorder_pl(npl, symba_plA)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_reorder_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     TYPE(symba_pl), INTENT(INOUT)  :: symba_plA

! Internals
     INTEGER(I4B)                              :: i
     INTEGER(I4B), DIMENSION(:), ALLOCATABLE   :: index
     REAL(DP), DIMENSION(:), ALLOCATABLE       :: mass
     REAL(DP), DIMENSION(:,:), allocatable     :: symba_plwkspA
     INTEGER(I4B), DIMENSION(:,:), allocatable :: symba_plwkspA_id_status

! Executable code
     ALLOCATE(index(npl), mass(npl))
     ALLOCATE(symba_plwkspA(12,npl))
     ALLOCATE(symba_plwkspA_id_status(2,npl))

     DO i = 1, npl
          mass(i) = symba_plA%helio%swiftest%mass(i)
          symba_plwkspA_id_status(1,i) = symba_plA%helio%swiftest%name(i)
          symba_plwkspA_id_status(2,i) = symba_plA%helio%swiftest%status(i)
          symba_plwkspA(1,i) = symba_plA%helio%swiftest%mass(i)
          symba_plwkspA(2,i) = symba_plA%helio%swiftest%radius(i)
          symba_plwkspA(3,i) = symba_plA%helio%swiftest%xh(1,i)
          symba_plwkspA(4,i) = symba_plA%helio%swiftest%xh(2,i)
          symba_plwkspA(5,i) = symba_plA%helio%swiftest%xh(3,i)
          symba_plwkspA(6,i) = symba_plA%helio%swiftest%vh(1,i)
          symba_plwkspA(7,i) = symba_plA%helio%swiftest%vh(2,i)
          symba_plwkspA(8,i) = symba_plA%helio%swiftest%vh(3,i)
          symba_plwkspA(9,i) = symba_plA%helio%swiftest%rhill(i)
          symba_plwkspA(10,i) = symba_plA%helio%ah(1,i)
          symba_plwkspA(11,i) = symba_plA%helio%ah(2,i)
          symba_plwkspA(12,i) = symba_plA%helio%ah(3,i)
     END DO
     CALL util_index(mass, index)
     WRITE(*,*) "************ REORDER ***************"
     DO i = 1, npl
          symba_plA%helio%swiftest%name(i) = symba_plwkspA_id_status(1,index(npl-i+1))
          symba_plA%helio%swiftest%status(i) = symba_plwkspA_id_status(2,index(npl-i+1))
          symba_plA%helio%swiftest%mass(i) = symba_plwkspA(1,index(npl-i+1))
          symba_plA%helio%swiftest%radius(i) = symba_plwkspA(2,index(npl-i+1))
          symba_plA%helio%swiftest%xh(1,i) = symba_plwkspA(3,index(npl-i+1))
          symba_plA%helio%swiftest%xh(2,i) = symba_plwkspA(4,index(npl-i+1))
          symba_plA%helio%swiftest%xh(3,i) = symba_plwkspA(5,index(npl-i+1))
          symba_plA%helio%swiftest%vh(1,i) = symba_plwkspA(6,index(npl-i+1))
          symba_plA%helio%swiftest%vh(2,i) = symba_plwkspA(7,index(npl-i+1))
          symba_plA%helio%swiftest%vh(3,i) = symba_plwkspA(8,index(npl-i+1))
          symba_plA%helio%swiftest%rhill(i) = symba_plwkspA(9,index(npl-i+1))
          symba_plA%helio%ah(1,i) = symba_plwkspA(10,index(npl-i+1))
          symba_plA%helio%ah(2,i) = symba_plwkspA(11,index(npl-i+1))
          symba_plA%helio%ah(3,i) = symba_plwkspA(12,index(npl-i+1))


     END DO
     IF (ALLOCATED(symba_plwkspA)) DEALLOCATE(symba_plwkspA)
     IF (ALLOCATED(symba_plwkspA_id_status)) DEALLOCATE(symba_plwkspA_id_status)
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
