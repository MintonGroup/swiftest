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
SUBROUTINE symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_getacch
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                      :: lextra_force
     INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, nplplenc
     REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
     TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
     TYPE(symba_plplenc), INTENT(INOUT)            :: plplenc_list

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j, index_i, index_j
     REAL(DP)                                     :: rji2, irij3, faci, facj, r2, fac
     REAL(DP), DIMENSION(NDIM)                    :: dx
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
! Added by D. Minton
     REAL(DP), DIMENSION(NDIM) :: accsum

! Executable code
     
     !Removed by D. Minton
     !helio_pliP => symba_pl1P%helio
     !symba_pliP => symba_pl1P
     !^^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
     !$OMP PRIVATE(i,helio_pliP) &
     !$OMP SHARED(npl,symba_pl1P) 
     DO i = 2, npl
          !Removed by D. Minton
          !helio_pliP => helio_pliP%nextP
          !^^^^^^^^^^^^^^^^^^^^
          !Added by D. Minton
          !helio_pliP => symba_pl1P%symba_plPA(i)%thisP%helio
          !^^^^^^^^^^^^^^^^^^
          symba_plA%helio%ah(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO
     !$OMP END PARALLEL DO
     DO i = 2, nplm
          !Removed by D. Minton
          !symba_pljP => symba_pliP
          !^^^^^^^^^^^^^^^^^^^^
          !Added by D. Minton
          accsum(:)=0.0_DP
          !^^^^^^^^^^^^^^^^^^^
          ! OpenMP parallelization added by D. Minton
          !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
          !$OMP PRIVATE(j,symba_pljP,helio_pljP,dx,rji2,irij3,faci,facj) &
          !$OMP SHARED(i,npl,symba_pliP,helio_pliP,symba_pl1P) &
          !$OMP REDUCTION (+:accsum)
          DO j = i + 1, npl
               !Removed by D. Minton
               !symba_pljP => symba_pljP%nextP
               !^^^^^^^^^^^^^^^^^^^^
               !Added by D. Minton
               !symba_pljP=>symba_pl1P%symba_plPA(j)%thisP
               !^^^^^^^^^^^^^^^^^^
               IF ((.NOT. symba_plA%lmerged(i)) .OR. &
                    (.NOT. symba_plA%lmerged(j))) &
                    ! .OR.  &
                   !(.NOT. ASSOCIATED(symba_pliP%parentP, symba_pljP%parentP))) & !need to handle parents and children separately
                   THEN !need to handle parents and children separately
                    dx(:) = symba_plA%helio%swiftest%xh(:,i) - symba_plA%helio%swiftest%xh(:,j)
                    rji2 = DOT_PRODUCT(dx(:), dx(:))
                    irij3 = 1.0_DP/(rji2*SQRT(rji2))
                    faci = symba_plA%helio%swiftest%mass(i)*irij3
                    facj = symba_plA%helio%swiftest%mass(j)*irij3
                    !Removed by D. Minton
                    !helio_pliP%ah(:) = helio_pliP%ah(:) + facj*dx(:)
                    !^^^^^^^^^^^^^^^^^^^^
                    !Added by D. Minton
                    accsum(:) = accsum(:) + facj*dx(:)
                    !^^^^^^^^^^^^^^^^^^^^
                    symba_plA%helio%ah(:,j) = symba_plA%helio%ah(:,j) - faci*dx(:)
               END IF
          END DO
          !$OMP END PARALLEL DO
          !Added by D. Minton
          symba_plA%helio%ah(:,i)=symba_plA%helio%ah(:,i)+accsum(:)
          !^^^^^^^^^^^^^^^^^^^^
     END DO
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) &
     !$OMP PRIVATE(i,symba_pliP,symba_pljP,helio_pliP,helio_pljP,dx,rji2,irij3,faci,facj) &
     !$OMP SHARED(plplenc_list,nplplenc)
     DO i = 1, nplplenc
          index_i = plplenc_list%id1(i)
          index_j = plplenc_list%id2(i)
          IF ((.NOT. symba_plA%lmerged(index_i)) .OR. &
               (.NOT. symba_plA%lmerged(index_j)))  &
               ! .OR. &
              !(.NOT. ASSOCIATED(symba_pliP%parentP, symba_pljP%parentP))) & !need to update parent/children
               THEN !need to update parent/children
               dx(:) = symba_plA%helio%swiftest%xh(:,index_j) - symba_plA%helio%swiftest%xh(:,index_i)
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP/(rji2*SQRT(rji2))
               faci = symba_plA%helio%swiftest%mass(index_i)*irij3
               facj = symba_plA%helio%swiftest%mass(index_j)*irij3
               !$OMP CRITICAL
               symba_plA%helio%ah(:,index_i) = symba_plA%helio%ah(:,index_i) - facj*dx(:)
               symba_plA%helio%ah(:,index_j) = symba_plA%helio%ah(:,index_j) + faci*dx(:)
               !$OMP END CRITICAL
          END IF
     END DO
     !$OMP END PARALLEL DO
     IF (j2rp2 /= 0.0_DP) THEN
          IF (lmalloc) THEN
               ALLOCATE(aobl(NDIM, nplmax), irh(nplmax))
               lmalloc = .FALSE.
          END IF
          !Removed by D. Minton
          !swifter_plP => swifter_pl1P
          !^^^^^^^^^^^^^^^^^^^^^
          ! OpenMP parallelization added by D. Minton
          !$OMP PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) &
          !$OMP PRIVATE(i,swifter_plP,r2) &
          !$OMP SHARED(npl,symba_pl1P,xh,irh)
          DO i = 2, npl
               !Removed by D. Minton
               !swifter_plP => swifter_plP%nextP
               !^^^^^^^^^^^^^^^^^^
               !Added by D. Minton
               !swifter_plP=>symba_pl1P%symba_plPA(i)%thisP%helio%swifter
               !^^^^^^^^^^^^^^^^^^
               !xh(:, i) = swifter_plP%xh(:)
               r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          !$OMP END PARALLEL DO
          CALL obl_acc(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, symba_plA%helio%swiftest%xh(:,:), irh, aobl)
          !Removed by D. Minton
          !helio_pliP => symba_pl1P%helio
          !^^^^^^^^^^^^^^^^^^^^^
          ! OpenMP parallelization added by D. Minton
          !$OMP PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) &
          !$OMP PRIVATE(i,helio_pliP) &
          !$OMP SHARED(npl,symba_pl1P,aobl)
          DO i = 2, npl
               !Removed by D. Minton
               !helio_pliP => helio_pliP%nextP
               !^^^^^^^^^^^^^^^^^^^^^
               !Aded by D. Minton
               !helio_pliP => symba_pl1P%symba_plPA(i)%thisP%helio
               !^^^^^^^^^^^^^^^^^
               symba_plA%helio%ah(:,i) = symba_plA%helio%ah(:,i) + aobl(:, i) - aobl(:, 1)
          END DO
          !$OMP END PARALLEL DO
     END IF
     IF (lextra_force) CALL symba_user_getacch(t, npl, symba_plA)

     RETURN

END SUBROUTINE symba_getacch
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
