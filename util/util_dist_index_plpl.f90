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
!              : mass         : array of planet masses
!              : mtiny        : value of mtiny (we do not compare semi interacting particles with themselves)
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : num_comparisons     : length of the distance array
!              : k_plpl              : matrix of index conversions (linear index to i, j indices)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL util_dist_index_plpl(npl, num_comparisons, ik, jk)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE util_dist_index_plpl(npl, mass, mtiny, num_comparisons, k_plpl)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_index_plpl
     IMPLICIT NONE

! Arguments
     REAL(DP), DIMENSION(npl), INTENT(IN) :: mass
     INTEGER(I4B), INTENT(IN)  :: npl
     REAL (DP), INTENT(IN) :: mtiny
     INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE,INTENT(OUT) :: k_plpl
     INTEGER(I4B), INTENT(OUT) :: num_comparisons

! Internals
     INTEGER(I4B)              :: i,j,counter,nplm

! Executable code
     nplm = count(mass>mtiny)

     num_comparisons = ((npl - 1) * (npl - 2) / 2) - ( (npl-nplm-1) * ((npl-nplm-1)+1)/2 )! number of entries in a strict lower triangle, npl x npl, minus first column
     allocate(k_plpl(num_comparisons,2))

     ! this is a 'fancier' code, but so far i think it runs slower
     ! so leaving it in, but commenting it out
     ! i think it's because of the 'mod' call, but i haven't profiled it yet
     ! don't forget to uncomment the 'k' declaration up top!
     ! allocate(k(num_comparisons))

     ! m = ceiling(sqrt(2. * num_comparisons))

     ! k = (/(i, i=1,num_comparisons, 1)/)

     ! ik_plpl = m - nint( sqrt( dble(2) * (dble(1) + num_comparisons - k))) + 1
     ! jk_plpl = mod(k + (ik_plpl - 1) * ik_plpl / 2 - 1, m) + 2

     ! brute force the index creation

!$omp parallel do default(none) schedule(dynamic) &
!$omp shared (k_plpl, npl) &
!$omp private (i, j, counter)
     do i = 2,nplm
          counter = (i - 2) * npl - i*(i-1)/2 + 2
          k_plpl(counter:counter+(npl-(i+1)),1) = i
          k_plpl(counter:counter+(npl-(i+1)),2) = (/(j, j=i+1,npl, 1)/)
     enddo
!$omp end parallel do

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
