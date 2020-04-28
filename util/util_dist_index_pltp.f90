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
SUBROUTINE util_dist_index_pltp(nplm, ntp, num_comparisons, k_pltp)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_interfaces, EXCEPT_THIS_ONE => util_dist_index_pltp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)  :: nplm, ntp
     INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE,INTENT(OUT) :: k_pltp
     INTEGER(I4B), INTENT(OUT) :: num_comparisons

! Internals
     INTEGER(I4B)              :: i,j,ii,jj,nb,np,nt,counter, nplm1, ntp1

! Executable code
     num_comparisons = (nplm - 1) * ntp ! number of entries in our distance array

     allocate(k_pltp(num_comparisons,2))


!$omp parallel do schedule(static) default(none) &
!$omp shared(k_pltp, nplm, ntp) &
!$omp private(i, j, counter)
     do i = 2,nplm
          counter = (i-2) * ntp + 1
          do j = 1,ntp
               k_pltp(counter,1) = i
               k_pltp(counter,2) = j
               counter = counter + 1
          enddo
     enddo
!$omp end parallel do

    !  nb = 2 ! number of blocks
    !  np = (nplm1-1)/nb ! number of planets per block
    !  nt = ntp1/nb ! number of test particles per block
    !  counter = 1

    !  do i = 2,np*nb,np
    !       do j = 1,nt*nb,nt
    !            do ii = i, i+np-1
    !                 do jj = j, j+nt-1
    !                      print *,'i j ii jj: ',i,j,ii,jj
    !                      k_pltp(counter,1) = ii
    !                      k_pltp(counter,2) = jj
    !                      counter = counter + 1
    !                 enddo
    !            enddo
    !       enddo
    !  enddo

    !  print *,'break'
    !            do ii = i, nplm1
    !                 do jj = j, ntp1
    !                      print *,'i j ii jj: ',i,j,ii,jj
    !                      k_pltp(counter,1) = ii
    !                      k_pltp(counter,2) = jj
    !                      counter = counter + 1
    !                 enddo
    !  enddo

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
