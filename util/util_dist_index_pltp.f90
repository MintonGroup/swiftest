!**********************************************************************************************************************************
!
!  Unit Name   : util_dist_index_plpl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Turns i,j indices into k index for use in the Euclidean distance matrix
!
!  Input
!    Arguments : npl          : number of planets
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
SUBROUTINE util_dist_index_pltp(npl, ntp, ik_pltp, jk_pltp)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_index_pltp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     INTEGER(I4B), DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: ik_pltp, jk_pltp

! Internals
     INTEGER(I4B)              :: i,j,l,k_count

! Executable code
     l = (npl - 1) * ntp ! number of entries in our distance array

     allocate(ik_pltp(l))
     allocate(jk_pltp(l)) 

     k_count = 1

     do i = 2, npl
          do j = 1, ntp
               ik_pltp(k_count) = i
               jk_pltp(k_count) = j
               k_count = k_count + 1
          enddo
     enddo

     RETURN

END SUBROUTINE util_dist_index_pltp
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
