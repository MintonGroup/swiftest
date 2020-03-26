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
SUBROUTINE util_dist_index_plpl(npl, num_comparisons, ik_plpl, jk_plpl)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_index_plpl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: npl
     INTEGER(I4B), DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: ik_plpl, jk_plpl
     INTEGER(I4B), INTENT(OUT) :: num_comparisons

! Internals
     INTEGER(I4B)              :: i,j,k,count,m
     ! INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: k

! Executable code
     num_comparisons = (npl - 1) * (npl - 2) / 2 ! number of entries in a strict lower triangle, npl x npl, minus first column
     m = ceiling(sqrt(2. * num_comparisons))
     
     allocate(ik_plpl(num_comparisons))
     allocate(jk_plpl(num_comparisons))

     ! this is a 'fancier' code, but so far i think it runs slower
     ! so leaving it in, but commenting it out
     ! i think it's because of the 'mod' call, but i haven't profiled it yet
     ! don't forget to uncomment the 'k' declaration up top!
     ! allocate(k(num_comparisons))


     ! k = (/(i, i=1,num_comparisons, 1)/)

     ! ik_plpl = m - nint( sqrt( dble(2) * (dble(1) + num_comparisons - k))) + 1
     ! jk_plpl = mod(k + (ik_plpl - 1) * ik_plpl / 2 - 1, m) + 2

     ! brute force the index creation

!$omp parallel do default(none) schedule(dynamic) &
!$omp shared (ik_plpl, jk_plpl, npl) &
!$omp private (i, j, count)
     do i = 2,npl
          count = (i - 2) * npl - i*(i-1)/2 + 2
          ik_plpl(count:count+(npl-(i+1))) = i
          jk_plpl(count:count+(npl-(i+1))) = (/(j, j=i+1,npl, 1)/)
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
