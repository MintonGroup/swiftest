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
SUBROUTINE util_dist_eucl_pltp(npl, ntp, planets, test_particles, ik_pltp, jk_pltp, outvar)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_eucl_pltp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     INTEGER(I4B), DIMENSION((npl-1)*ntp),INTENT(IN) :: ik_pltp, jk_pltp
     REAL(DP),DIMENSION(NDIM,npl),INTENT(IN) :: planets
     REAL(DP),DIMENSION(NDIM,ntp),INTENT(IN) :: test_particles
     REAL(DP), DIMENSION(NDIM,(npl-1)*ntp),INTENT(INOUT) :: outvar

! Internals
     INTEGER(I4B)              :: i
     
! Executable code

     do i = 1,(npl-1)*ntp
          outvar(:,i) = test_particles(:,jk_pltp(i)) - planets(:,ik_pltp(i))
     enddo

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
