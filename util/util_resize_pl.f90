!**********************************************************************************************************************************
!
!  Unit Name   : util_resize_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Index input real array into ascending numerical order using Quicksort algorithm
!
!  Input
!    Arguments : arr   : array to index
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : index : index table for sorted array
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL util_resize_pl(symba_plA, npl)
!
!  Notes       : Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
!                Vetterling, and Flannery, 2nd ed., pp. 1173-4
!
!**********************************************************************************************************************************
SUBROUTINE util_resize_pl(symba_plA, npl_new, npl_old)

! Modules
     USE module_parameters
     USE module_symba
     USE module_swiftest
     USE module_helio
     USE module_nrutil
     USE module_swiftestalloc
     USE module_interfaces, EXCEPT_THIS_ONE => util_resize_pl
     IMPLICIT NONE

! Arguments
     TYPE(symba_pl), INTENT(INOUT) :: symba_plA
     INTEGER(I4B), INTENT(IN)      :: npl_old, npl_new

! Internals
     TYPE(symba_pl)                :: new_symba_plA

! Executable code
     IF (npl_new >= npl_old) THEN 
          CALL symba_pl_allocate(new_symba_plA, npl_new)
          new_symba_plA%helio%swiftest%name(1:npl_old) = symba_plA%helio%swiftest%name(1:npl_old)
          new_symba_plA%helio%swiftest%status(1:npl_old) = symba_plA%helio%swiftest%status(1:npl_old)
          new_symba_plA%helio%swiftest%mass(1:npl_old) = symba_plA%helio%swiftest%mass(1:npl_old)
          new_symba_plA%helio%swiftest%radius(1:npl_old) = symba_plA%helio%swiftest%radius(1:npl_old)
          new_symba_plA%helio%swiftest%xh(1,1:npl_old) = symba_plA%helio%swiftest%xh(1,1:npl_old)
          new_symba_plA%helio%swiftest%xh(2,1:npl_old) = symba_plA%helio%swiftest%xh(2,1:npl_old)
          new_symba_plA%helio%swiftest%xh(3,1:npl_old) = symba_plA%helio%swiftest%xh(3,1:npl_old)
          new_symba_plA%helio%swiftest%vh(1,1:npl_old) = symba_plA%helio%swiftest%vh(1,1:npl_old)
          new_symba_plA%helio%swiftest%vh(2,1:npl_old) = symba_plA%helio%swiftest%vh(2,1:npl_old)
          new_symba_plA%helio%swiftest%vh(3,1:npl_old) = symba_plA%helio%swiftest%vh(3,1:npl_old)
          new_symba_plA%helio%swiftest%rhill(1:npl_old) = symba_plA%helio%swiftest%rhill(1:npl_old)
          new_symba_plA%helio%swiftest%xb(1,1:npl_old) = symba_plA%helio%swiftest%xb(1,1:npl_old)
          new_symba_plA%helio%swiftest%xb(2,1:npl_old) = symba_plA%helio%swiftest%xb(2,1:npl_old)
          new_symba_plA%helio%swiftest%xb(3,1:npl_old) = symba_plA%helio%swiftest%xb(3,1:npl_old)
          new_symba_plA%helio%swiftest%vb(1,1:npl_old) = symba_plA%helio%swiftest%vb(1,1:npl_old)
          new_symba_plA%helio%swiftest%vb(2,1:npl_old) = symba_plA%helio%swiftest%vb(2,1:npl_old)
          new_symba_plA%helio%swiftest%vb(3,1:npl_old) = symba_plA%helio%swiftest%vb(3,1:npl_old)
          new_symba_plA%helio%ah(1,1:npl_old) = symba_plA%helio%ah(1,1:npl_old)
          new_symba_plA%helio%ah(2,1:npl_old) = symba_plA%helio%ah(2,1:npl_old)
          new_symba_plA%helio%ah(3,1:npl_old) = symba_plA%helio%ah(3,1:npl_old)

     END IF



     IF (npl_new < npl_old) THEN 
          CALL symba_pl_allocate(new_symba_plA, npl_new)
          new_symba_plA%helio%swiftest%name(1:npl_new) = symba_plA%helio%swiftest%name(1:npl_new)
          new_symba_plA%helio%swiftest%status(1:npl_new) = symba_plA%helio%swiftest%status(1:npl_new)
          new_symba_plA%helio%swiftest%mass(1:npl_new) = symba_plA%helio%swiftest%mass(1:npl_new)
          new_symba_plA%helio%swiftest%radius(1:npl_new) = symba_plA%helio%swiftest%radius(1:npl_new)
          new_symba_plA%helio%swiftest%xh(1,1:npl_new) = symba_plA%helio%swiftest%xh(1,1:npl_new)
          new_symba_plA%helio%swiftest%xh(2,1:npl_new) = symba_plA%helio%swiftest%xh(2,1:npl_new)
          new_symba_plA%helio%swiftest%xh(3,1:npl_new) = symba_plA%helio%swiftest%xh(3,1:npl_new)
          new_symba_plA%helio%swiftest%vh(1,1:npl_new) = symba_plA%helio%swiftest%vh(1,1:npl_new)
          new_symba_plA%helio%swiftest%vh(2,1:npl_new) = symba_plA%helio%swiftest%vh(2,1:npl_new)
          new_symba_plA%helio%swiftest%vh(3,1:npl_new) = symba_plA%helio%swiftest%vh(3,1:npl_new)
          new_symba_plA%helio%swiftest%rhill(1:npl_new) = symba_plA%helio%swiftest%rhill(1:npl_new)
          new_symba_plA%helio%swiftest%xb(1,1:npl_new) = symba_plA%helio%swiftest%xb(1,1:npl_new)
          new_symba_plA%helio%swiftest%xb(2,1:npl_new) = symba_plA%helio%swiftest%xb(2,1:npl_new)
          new_symba_plA%helio%swiftest%xb(3,1:npl_new) = symba_plA%helio%swiftest%xb(3,1:npl_new)
          new_symba_plA%helio%swiftest%vb(1,1:npl_new) = symba_plA%helio%swiftest%vb(1,1:npl_new)
          new_symba_plA%helio%swiftest%vb(2,1:npl_new) = symba_plA%helio%swiftest%vb(2,1:npl_new)
          new_symba_plA%helio%swiftest%vb(3,1:npl_new) = symba_plA%helio%swiftest%vb(3,1:npl_new)
          new_symba_plA%helio%ah(1,1:npl_new) = symba_plA%helio%ah(1,1:npl_old)
          new_symba_plA%helio%ah(2,1:npl_new) = symba_plA%helio%ah(2,1:npl_old)
          new_symba_plA%helio%ah(3,1:npl_new) = symba_plA%helio%ah(3,1:npl_old)

     END IF
     CALL symba_pl_deallocate(symba_plA)
     CALL symba_pl_allocate(symba_plA, npl_new)
     symba_plA%helio%swiftest%name = new_symba_plA%helio%swiftest%name
     symba_plA%helio%swiftest%status = new_symba_plA%helio%swiftest%status
     symba_plA%helio%swiftest%mass = new_symba_plA%helio%swiftest%mass
     symba_plA%helio%swiftest%radius = new_symba_plA%helio%swiftest%radius
     symba_plA%helio%swiftest%xh = new_symba_plA%helio%swiftest%xh
     symba_plA%helio%swiftest%vh = new_symba_plA%helio%swiftest%vh
     symba_plA%helio%swiftest%rhill = new_symba_plA%helio%swiftest%rhill
     symba_plA%helio%swiftest%xb = new_symba_plA%helio%swiftest%xb 
     symba_plA%helio%swiftest%vb = new_symba_plA%helio%swiftest%vb
     symba_plA%helio%ah = new_symba_plA%helio%ah
     CALL symba_pl_deallocate(new_symba_plA)


     RETURN

END SUBROUTINE util_resize_pl
!**********************************************************************************************************************************
!
!  Author(s)   : C.Wishard and J.Pouplin
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
