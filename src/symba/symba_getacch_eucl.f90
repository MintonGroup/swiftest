!**********************************************************************************************************************************
!
!  Unit Name   : symba_getacch
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of planets
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplm         : number of planets with mass > mtiny
!                nplmax       : maximum allowed number of planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-test particle encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_pl1P, j2rp2, j4rp4, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_getacch.f
!
!                Accelerations in an encounter are not included here
!
!**********************************************************************************************************************************
SUBROUTINE symba_getacch_eucl(lextra_force, t, npl, nplm, nplmax, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list, &
     num_plpl_comparisons, k_plpl)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_getacch_eucl
     USE omp_lib
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                      :: lextra_force
     INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, nplplenc, num_plpl_comparisons
     REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
     TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
     TYPE(symba_plplenc), INTENT(INOUT)            :: plplenc_list
     INTEGER(I4B), DIMENSION(2,num_plpl_comparisons), INTENT(IN) :: k_plpl


! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j, index_i, index_j, k, counter
     REAL(DP)                                     :: rji2, irij3, faci, facj, r2, fac
     REAL(DP), DIMENSION(NDIM)                    :: dx
     REAL(DP), DIMENSION(NDIM, npl)               :: ah
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     ! REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: dist_plpl_array


!Executable code
 
     symba_plA%helio%ah(:,2:npl) = 0.0_DP
     ah(:,2:npl) = 0.0_DP
     
     ! CALL util_dist_eucl_plpl(npl,symba_plA%helio%swiftest%xh, num_plpl_comparisons, k_plpl, dist_plpl_array) ! does not care about mtiny

! There is floating point arithmetic round off error in this loop
! For now, we will keep it in the serial operation, so we can easily compare
! it to the older swifter versions

!$omp parallel do default(none) schedule(static) &
!$omp num_threads(min(omp_get_max_threads(),ceiling(num_plpl_comparisons/10000.))) &
!$omp shared (num_plpl_comparisons, k_plpl, symba_plA) &
!$omp private (i, j, k, dx, rji2, irij3, faci, facj) &
!$omp reduction(+:ah)
     DO k = 1, num_plpl_comparisons
          i = k_plpl(1,k)
          j = k_plpl(2,k)
          
          IF ((.NOT. symba_plA%lmerged(i) .OR. (.NOT. symba_plA%lmerged(j)) .OR. &
               (symba_plA%index_parent(i) /= symba_plA%index_parent(j)))) THEN
               
               dx(:) = symba_plA%helio%swiftest%xh(:,k_plpl(2,k)) - symba_plA%helio%swiftest%xh(:,k_plpl(1,k))
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP/(rji2*SQRT(rji2))
               faci = symba_plA%helio%swiftest%mass(i)*irij3
               facj = symba_plA%helio%swiftest%mass(j)*irij3
               ah(:,i) = ah(:,i) + facj*dx(:)
               ah(:,j) = ah(:,j) - faci*dx(:)

          ENDIF
     END DO
!$omp end parallel do

     symba_plA%helio%ah(:,2:npl) = ah(:,2:npl)

     DO i = 1, nplplenc
          index_i = plplenc_list%index1(i)
          index_j = plplenc_list%index2(i)
          IF ((.NOT. symba_plA%lmerged(index_i)) .OR. (.NOT. symba_plA%lmerged(index_j))  &
                .OR. (symba_plA%index_parent(index_i) /= symba_plA%index_parent(index_j))) THEN !need to update parent/children
               dx(:) = symba_plA%helio%swiftest%xh(:,index_j) - symba_plA%helio%swiftest%xh(:,index_i)
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP/(rji2*SQRT(rji2))
               faci = symba_plA%helio%swiftest%mass(index_i)*irij3
               facj = symba_plA%helio%swiftest%mass(index_j)*irij3
               symba_plA%helio%ah(:,index_i) = symba_plA%helio%ah(:,index_i) - facj*dx(:)
               symba_plA%helio%ah(:,index_j) = symba_plA%helio%ah(:,index_j) + faci*dx(:)
          END IF
     END DO

     IF (j2rp2 /= 0.0_DP) THEN
          IF (lmalloc) THEN
               ALLOCATE(xh(NDIM, npl), aobl(NDIM, npl), irh(npl))
               lmalloc = .FALSE.
          END IF
          DO i = 2, npl
               
               r2 = DOT_PRODUCT(symba_plA%helio%swiftest%xh(:,i), symba_plA%helio%swiftest%xh(:,i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_acc(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, symba_plA%helio%swiftest%xh(:,:), irh, aobl)
          DO i = 2, npl
               symba_plA%helio%ah(:,i) = symba_plA%helio%ah(:,i) + aobl(:, i) - aobl(:, 1)
          END DO
     END IF

     IF (lextra_force) CALL symba_user_getacch(t, npl, symba_plA)

     RETURN

END SUBROUTINE symba_getacch_eucl
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
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
