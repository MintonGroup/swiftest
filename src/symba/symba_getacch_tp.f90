!**********************************************************************************************************************************
!
!  Unit Name   : symba_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of test particles
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplm         : number of planets with mass > mtiny
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P   : pointer to head of active SyMBA test particle structure linked-list
!                xh           : heliocentric positions of planets at time t
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                npltpenc     : number of planet-test particle encounters
!                pltpenc_list : array of planet-test particle encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_tp1P   : pointer to head of active SyMBA test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, xh, j2rp2, j4rp4,
!                                      npltpenc, pltpenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_getacch.f
!
!                Accelerations in an encounter are not included here
!
!**********************************************************************************************************************************
SUBROUTINE symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, xh, j2rp2, j4rp4, npltpenc,  &
     pltpenc_list)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_getacch_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                      :: lextra_force
     INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, ntp, ntpmax, npltpenc
     REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN)    :: xh
     TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                 :: symba_tpA
     TYPE(symba_pltpenc), INTENT(IN)               :: pltpenc_list

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j, index_pl, index_tp
     REAL(DP)                                     :: rji2, irij3, faci, facj, r2, fac, mu
     REAL(DP), DIMENSION(NDIM)                    :: dx
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh, irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: aobl, xht, aoblt

! Executable code
     !Removed by D. Minton
     !helio_tpP => symba_tp1P%helio
     !^^^^^^^^^^^^^^^^^^^^
     ! OpenMP parallelization added by D. Minton
     DO i = 1, ntp
          !Added by D. Minton
          !helio_tpP => symba_tp1P%symba_tpPA(i)%thisP%helio
          !^^^^^^^^^^^^^^^^^^
          symba_tpA%helio%ah(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          IF (symba_tpA%helio%swiftest%status(i) == ACTIVE) THEN
               !swifter_plP => swifter_pl1P
               !DO j = 2, nplm
               DO j = 2, nplm
                    !swifter_plP => swifter_plP%nextP
                    dx(:) = symba_tpA%helio%swiftest%xh(:,i) - xh(:, j)
                    r2 = DOT_PRODUCT(dx(:), dx(:))
                    fac = symba_PlA%helio%swiftest%mass(j)/(r2*SQRT(r2))
                    symba_tpA%helio%ah(:,i) = symba_tpA%helio%ah(:,i) - fac*dx(:)
               END DO
          END IF
          !Removed by D. Minton
          !helio_tpP => helio_tpP%nextP
          !^^^^^^^^^^^^^^^^^^^^
     END DO
     DO i = 1, npltpenc
          !swifter_plP => pltpenc_list(i)%plP%helio%swifter
          !helio_tpP => pltpenc_list(i)%tpP%helio
          index_pl = pltpenc_list%indexpl(i)
          index_tp = pltpenc_list%indextp(i)
          IF (symba_tpA%helio%swiftest%status(index_tp) == ACTIVE) THEN
               dx(:) = symba_tpA%helio%swiftest%xh(:,index_tp) - symba_plA%helio%swiftest%xh(:,index_pl)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               fac = symba_plA%helio%swiftest%mass(index_pl)/(r2*SQRT(r2))
               symba_tpA%helio%ah(:,index_tp) = symba_tpA%helio%ah(:,index_tp) + fac*dx(:)
          END IF
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          IF (lmalloc) THEN
               ALLOCATE(aobl(NDIM, nplmax), irh(nplmax), xht(NDIM, ntpmax), aoblt(NDIM, ntpmax), irht(ntpmax))
               lmalloc = .FALSE.
          END IF
          DO i = 2, npl
               r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_acc(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, symba_plA%helio%swiftest%xh(:,:), irh, aobl)
          mu = symba_plA%helio%swiftest%mass(1)
          DO i = 1, ntp
               xht(:, i) = symba_tpA%helio%swiftest%xh(:,i) !optimize
               r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
               irht(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
          DO i = 1, ntp
               IF (symba_tpA%helio%swiftest%status(i) == ACTIVE) &
               symba_tpA%helio%ah(:,i) = symba_tpA%helio%ah(:,i) + aoblt(:, i) - aobl(:, 1)
          END DO
     END IF
     IF (lextra_force) CALL symba_user_getacch_tp(t, ntp, symba_tpA)

     RETURN

END SUBROUTINE symba_getacch_tp
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
