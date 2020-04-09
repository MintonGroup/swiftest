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
!    Arguments : npl : number of planets
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
     USE module_parameters
     USE module_swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_eucl_plpl
     USE omp_lib
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     INTEGER(I4B), DIMENSION(num_comparisons,2),INTENT(IN) :: k_plpl
     INTEGER(I4B), INTENT(IN) :: num_comparisons
     REAL(DP),DIMENSION(NDIM,npl),INTENT(IN) :: invar
     REAL(DP), DIMENSION(NDIM,num_comparisons),INTENT(INOUT) :: outvar

! Internals
     INTEGER(I4B) :: k
     
! Executable code

!$omp parallel do schedule(static) default(none) &
!$omp num_threads(min(omp_get_max_threads(),ceiling(num_comparisons/10000.))) &
!$omp shared (outvar, invar, num_comparisons, k_plpl) &
!$omp private(k)
      do k = 1,num_comparisons
           outvar(:,k) = invar(:,k_plpl(k,2)) - invar(:,k_plpl(k,1))
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
