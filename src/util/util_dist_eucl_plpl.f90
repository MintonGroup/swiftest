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
!    Arguments : npl : number of massive bodies
!              : invar : variable that we want to make comparisons between
!              : num_comparisons : number of comparisons to make
!              : k_plpl : matrix to convert linear index k into i,j indices
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : outvar : results of comparisons
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_dist_index(invar, num_comparisons, k_plpl_outvar)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_dist_eucl_plpl(npl, invar, num_comparisons, k_plpl, outvar)

! Modules
     use swiftest, EXCEPT_THIS_ONE => util_dist_eucl_plpl
     !$ USE omp_lib
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     INTEGER(I4B), DIMENSION(:,:),INTENT(IN) :: k_plpl
     INTEGER(I4B), INTENT(IN) :: num_comparisons
     REAL(DP),DIMENSION(:,:),INTENT(IN) :: invar
     REAL(DP), DIMENSION(:,:),INTENT(INOUT) :: outvar

! Internals
     INTEGER(I4B) :: k
     
! Executable code

!$omp parallel do schedule(static) default(none) &
!$omp num_threads(min(omp_get_max_threads(),ceiling(num_comparisons/10000.))) &
!$omp shared (outvar, invar, num_comparisons, k_plpl) &
!$omp private(k)
      do k = 1,num_comparisons
           outvar(:,k) = invar(:,k_plpl(2,k)) - invar(:,k_plpl(1,k))
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
