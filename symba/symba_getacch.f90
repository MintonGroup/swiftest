!**********************************************************************************************************************************
!
!  Unit Name   : symba_getacch
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_pl1P, j2rp2, j4rp4, nplplenc, plplenc_list)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_random_access, EXCEPT_THIS_ONE => symba_getacch
     USE module_interfaces, EXCEPT_THIS_ONE => symba_getacch
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                      :: lextra_force
     INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, nplplenc
     REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
     TYPE(symba_pl), POINTER                       :: symba_pl1P
     TYPE(symba_plplenc), DIMENSION(:), INTENT(IN) :: plplenc_list

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: rji2, irij3, faci, facj, r2, fac
     REAL(DP), DIMENSION(NDIM)                    :: dx
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_plP
     TYPE(helio_pl), POINTER                      :: helio_pliP, helio_pljP
     TYPE(symba_pl), POINTER                      :: symba_pliP, symba_pljP
! Added by D. Minton
     REAL(DP), DIMENSION(NDIM) :: accsum

! Executable code
     
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(AUTO) &
     !$OMP SHARED(npl) 
     DO i = 2, npl
          CALL get_point(i,helio_pliP)
          helio_pliP%ah(:) = 0.0_DP
     END DO
     !$OMP END PARALLEL DO

     symba_pliP => symba_pl1P
     DO i = 2, nplm
          symba_pliP => symba_pliP%nextP
          helio_pliP => symba_pliP%helio
          accsum(:) = 0.0_DP

          ! OpenMP parallelization added by D. Minton
          !$OMP PARALLEL DO DEFAULT (PRIVATE) SCHEDULE (AUTO) &
          !$OMP SHARED(i,npl,symba_pliP,helio_pliP) &
          !$OMP REDUCTION (+:accsum)
          DO j = i + 1, npl
               CALL get_point(j,symba_pljP)

               IF ((.NOT. symba_pliP%lmerged) .OR. (.NOT. symba_pljP%lmerged) .OR.                                                &
                   (.NOT. ASSOCIATED(symba_pliP%parentP, symba_pljP%parentP))) THEN
                    helio_pljP => symba_pljP%helio
                    dx(:) = helio_pljP%swifter%xh(:) - helio_pliP%swifter%xh(:)
                    rji2 = DOT_PRODUCT(dx(:), dx(:))
                    irij3 = 1.0_DP / (rji2 * SQRT(rji2))
                    faci = helio_pliP%swifter%mass * irij3
                    facj = helio_pljP%swifter%mass * irij3
                    accsum(:) = accsum(:) + facj*dx(:)
                    helio_pljP%ah(:) = helio_pljP%ah(:) - faci*dx(:)
               END IF
          END DO
          !$OMP END PARALLEL DO
          helio_pliP%ah(:) = helio_pliP%ah(:) + accsum(:)
     END DO

     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO DEFAULT (PRIVATE) SCHEDULE (AUTO) &
     !$OMP SHARED(nplplenc,plplenc_list)
     DO i = 1, nplplenc
          symba_pliP => plplenc_list(i)%pl1P
          symba_pljP => plplenc_list(i)%pl2P
          IF ((.NOT. symba_pliP%lmerged) .OR. (.NOT. symba_pljP%lmerged) .OR.                                                     &
              (.NOT. ASSOCIATED(symba_pliP%parentP, symba_pljP%parentP))) THEN
               helio_pliP => symba_pliP%helio
               helio_pljP => symba_pljP%helio
               dx(:) = helio_pljP%swifter%xh(:) - helio_pliP%swifter%xh(:)
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP / (rji2*SQRT(rji2))
               faci = helio_pliP%swifter%mass * irij3
               facj = helio_pljP%swifter%mass * irij3
               !$OMP CRITICAL
               helio_pliP%ah(:) = helio_pliP%ah(:) - facj*dx(:)
               helio_pljP%ah(:) = helio_pljP%ah(:) + faci*dx(:)
               !$OMP END CRITICAL
          END IF
     END DO
     !$OMP END PARALLEL DO
     IF (j2rp2 /= 0.0_DP) THEN
          swifter_pl1P => symba_pl1P%helio%swifter
          IF (lmalloc) THEN
               ALLOCATE(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax))
               lmalloc = .FALSE.
          END IF
          ! OpenMP parallelization added by D. Minton
          !$OMP PARALLEL DO DEFAULT (PRIVATE) SCHEDULE (AUTO) &
          !$OMP SHARED(npl,xh,irh)
          DO i = 2, npl
               CALL get_point(i,swifter_plP)
               xh(:, i) = swifter_plP%xh(:)
               r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
               irh(i) = 1.0_DP / SQRT(r2)
          END DO
          !$OMP END PARALLEL DO
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          ! OpenMP parallelization added by D. Minton
          !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE (AUTO) &
          !$OMP SHARED(npl,aobl)
          DO i = 2, npl
               CALL get_point(i,helio_pliP)
               helio_pliP%ah(:) = helio_pliP%ah(:) + aobl(:, i) - aobl(:, 1)
          END DO
          !$OMP END PARALLEL DO
     END IF
     IF (lextra_force) CALL symba_user_getacch(t, npl, symba_pl1P)

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
