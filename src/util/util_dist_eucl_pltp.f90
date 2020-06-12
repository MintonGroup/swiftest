! !**********************************************************************************************************************************
! !
! !  Unit Name   : util_dist_eucl_pltp
! !  Unit Type   : subroutine
! !  Project     : Swiftest
! !  Package     : util
! !  Language    : Fortran 90/95
! !
! !  Description : Calculates the Euclidean distance matrix (but in array form)
! !
! !  Input
! !    Arguments : npl          : number of planets
! !              : swifter_pl1P : pointer to head of SWIFTER planet structure linked-list
! !              : ik
! !              : jk
! !              : l
! !    Terminal  : none
! !    File      : none
! !
! !  Output
! !    Arguments : l            : length of the distance array
! !              : ik           : 
! !              : jk
! !    Terminal  : none
! !    File      : none
! !
! !  Invocation  : CALL util_dist_index(npl, swifter_pl1P)
! !
! !  Notes       : 
! !
! !**********************************************************************************************************************************
SUBROUTINE util_dist_eucl_pltp(npl, ntp, planets, test_particles, num_pltp_comparisons, k_pltp, outvar)

! Modules
     USE swiftest_globals
     USE swiftest_data_structures
     USE symba
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_eucl_pltp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl, ntp
     INTEGER(I4B), DIMENSION(2,num_pltp_comparisons),INTENT(IN) :: k_pltp
     INTEGER(I4B), INTENT(IN) :: num_pltp_comparisons
     REAL(DP),DIMENSION(NDIM,npl),INTENT(IN) :: planets
     REAL(DP),DIMENSION(NDIM,ntp),INTENT(IN) :: test_particles
     REAL(DP), DIMENSION(NDIM,num_pltp_comparisons),INTENT(INOUT) :: outvar

! Internals
     INTEGER(I4B)              :: k
     
! Executable code

!$omp parallel do default(none) schedule(static) &
!$omp shared (num_pltp_comparisons, test_particles, planets, outvar, k_pltp) &
!$omp private (k)
     do k = 1,num_pltp_comparisons
          outvar(:,k) = test_particles(:,k_pltp(2,k)) - planets(:,k_pltp(1,k))
     enddo
!$omp end parallel do
     RETURN

END SUBROUTINE util_dist_eucl_pltp
! !**********************************************************************************************************************************
! !
! !  Author(s)   : Jacob R. Elliott 
! !
! !  Revision Control System (RCS) Information
! !
! !  Source File : $RCSfile$
! !  Full Path   : $Source$
! !  Revision    : $Revision$
! !  Date        : $Date$
! !  Programmer  : $Author$
! !  Locked By   : $Locker$
! !  State       : $State$
! !
! !  Modification History:
! !
! !  $Log$
! !**********************************************************************************************************************************
