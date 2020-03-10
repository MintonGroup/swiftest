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
!    Arguments : num_comparisons     : length of the distance array
!              : ik_plpl           : i index to convert linear indices to matrix indices
!              : jk_plpl           : j index to convert linear indices to matrix indices
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_dist_index_plpl(npl, num_comparisons, ik, jk)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_dist_index_plpl(npl, l, ik, jk)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_index_plpl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     INTEGER(I4B), DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: ik, jk
     INTEGER(I4B), INTENT(OUT) :: l

! Internals
     INTEGER(I4B)              :: i,m,k_count, j
     INTEGER(I4B), DIMENSION(:),ALLOCATABLE :: k

! Executable code
     l = npl * (npl - 1) / 2 ! length of the distance matrix for a strict lower triangle, npl x npl
     l = l - (npl - 1)! however, swifter doesn't compare anything to the 
     m = ceiling(sqrt(2.*l))

     allocate(ik(l))
     allocate(jk(l))
     allocate(k(l))

     k = (/(i, i=1,l, 1)/)     

     k_count = 1

     do i = 2, npl
          do j = i + 1, npl
               ik(k_count) = i
               jk(k_count) = j
               k_count = k_count + 1
          enddo
     enddo

     RETURN

END SUBROUTINE util_dist_index_plpl
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
