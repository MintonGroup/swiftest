!**********************************************************************************************************************************
!
!  Unit Name   : util_dist_eucl_plpl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Calculates the Euclidean distance matrix (but in array form)
!
!  Input
!    Arguments : npl          : number of planets
!              : swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
!              : ik
!              : jk
!              : l
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : l            : length of the distance array
!              : ik           : 
!              : jk
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_dist_index(npl, swifter_pl1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_dist_eucl_plpl(npl, invar, l, ik, jk, outvar)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_eucl_plpl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     INTEGER(I4B), DIMENSION(:),INTENT(IN) :: ik, jk
     INTEGER(I4B), INTENT(IN) :: l
     REAL(DP),DIMENSION(NDIM,npl),INTENT(IN) :: invar
     REAL(DP), DIMENSION(NDIM,l),INTENT(INOUT) :: outvar

! Internals
     INTEGER(I4B)              :: i
     
! Executable code

!$omp parallel do schedule(static) default(none) &
!$omp shared (outvar, invar, jk, ik, l) &
!$omp private(i)
     do i = 1,l
          outvar(:,i) = invar(:,jk(i)) - invar(:,ik(i))
     enddo
!$omp end parallel do

     RETURN

END SUBROUTINE util_dist_eucl_plpl
!**********************************************************************************************************************************
!
!  Author(s)   : Jacob R. Elliott 
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
