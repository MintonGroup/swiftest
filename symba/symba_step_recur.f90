!**********************************************************************************************************************************
!
!  Unit Name   : symba_step_recur
!  Unit Type   : recursive subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step interacting planets and active test particles ahead in democratic heliocentric coordinates at the current
!                recursion level, if applicable, and descend to the next deeper level if necessary
!
!  Input
!    Arguments : lclose         : logical flag indicating whether to check for mergers
!                t              : time
!                ireci          : input recursion level
!                npl            : number of planets
!                nplm           : number of planets with mass > mtiny
!                ntp            : number of active test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                dt0            : time step (primary time step for overall integration)
!                eoffset        : energy offset (net energy lost in mergers)
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                eoffset        : energy offset (net energy lost in mergers)
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!    Terminal  : warning message
!    File      : none
!
!  Invocation  : CALL symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_pl1P, symba_tp1P, dt0, eoffset, nplplenc, npltpenc,
!                                      plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list,
!                                      encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_recur.F
!
!**********************************************************************************************************************************
RECURSIVE SUBROUTINE symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_pl1P, symba_tp1P, dt0, eoffset, nplplenc, npltpenc, &
     plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step_recur
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                         :: lclose
     INTEGER(I4B), INTENT(IN)                         :: ireci, npl, nplm, ntp, nplplenc, npltpenc
     INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt0
     REAL(DP), INTENT(INOUT)                          :: eoffset
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_pl), POINTER                          :: symba_pl1P
     TYPE(symba_tp), POINTER                          :: symba_tp1P
     TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
     TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list
     TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list

! Internals
     LOGICAL(LGT)              :: lencounter
     INTEGER(I4B)              :: i, j, irecp, icflg
     REAL(DP)                  :: dtl, dth, sgn
     REAL(DP), DIMENSION(NDIM) :: xr, vr, vbs
     TYPE(swifter_pl), POINTER :: swifter_pliP, swifter_pljP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP
     TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     dtl = dt0/(NTENC**ireci)
     dth = 0.5_DP*dtl
     IF (dtl/dt0 < TINY) THEN
          WRITE(*, *) "SWIFTER Warning:"
          WRITE(*, *) "   In symba_step_recur, local time step is too small"
          WRITE(*, *) "   Roundoff error will be important!"
     END IF
     irecp = ireci + 1
     IF (ireci == 0) THEN
          icflg = 0
          ! OpenMP parallelization added by D. Minton
          !$OMP PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) &
          !$OMP PRIVATE(i,symba_pliP,symba_pljP,swifter_pliP,swifter_pljP,xr,vr,lencounter) &
          !$OMP SHARED(plplenc_list,nplplenc,irecp,icflg,ireci,dtl)
          DO i = 1, nplplenc
               IF ((plplenc_list(i)%status == ACTIVE) .AND. (plplenc_list(i)%level == ireci)) THEN
                    symba_pliP => plplenc_list(i)%pl1P
                    symba_pljP => plplenc_list(i)%pl2P
                    swifter_pliP => symba_pliP%helio%swifter
                    swifter_pljP => symba_pljP%helio%swifter
                    xr(:) = swifter_pljP%xh(:) - swifter_pliP%xh(:)
                    vr(:) = swifter_pljP%vb(:) - swifter_pliP%vb(:)
                    CALL symba_chk(xr(:), vr(:), swifter_pliP%rhill, swifter_pljP%rhill, dtl, irecp, lencounter,                  &
                         plplenc_list(i)%lvdotr)
                    IF (lencounter) THEN
                         !Added by D. Minton
                         !$OMP CRITICAL
                         icflg = 1
                         symba_pliP%levelg = irecp
                         symba_pliP%levelm = MAX(irecp, symba_pliP%levelm)
                         symba_pljP%levelg = irecp
                         symba_pljP%levelm = MAX(irecp, symba_pljP%levelm)
                         plplenc_list(i)%level = irecp
                         !$OMP END CRITICAL
                    END IF
               END IF
          END DO
          !$OMP END PARALLEL DO
          DO i = 1, npltpenc
               IF ((pltpenc_list(i)%status == ACTIVE) .AND. (pltpenc_list(i)%level == ireci)) THEN
                    symba_pliP => pltpenc_list(i)%plP
                    symba_tpP => pltpenc_list(i)%tpP
                    swifter_pliP => symba_pliP%helio%swifter
                    swifter_tpP => symba_tpP%helio%swifter
                    xr(:) = swifter_tpP%xh(:) - swifter_pliP%xh(:)
                    vr(:) = swifter_tpP%vb(:) - swifter_pliP%vb(:)
                    CALL symba_chk(xr(:), vr(:), swifter_pliP%rhill, 0.0_DP, dtl, irecp, lencounter, pltpenc_list(i)%lvdotr)
                    IF (lencounter) THEN
                         icflg = 1
                         symba_pliP%levelg = irecp
                         symba_pliP%levelm = MAX(irecp, symba_pliP%levelm)
                         symba_tpP%levelg = irecp
                         symba_tpP%levelm = MAX(irecp, symba_tpP%levelm)
                         pltpenc_list(i)%level = irecp
                    END IF
               END IF
          END DO
          lencounter = (icflg == 1)
          sgn = 1.0_DP
          CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn)
          CALL symba_helio_drift(ireci, npl, symba_pl1P, dtl)
          IF (ntp > 0) CALL symba_helio_drift_tp(ireci, ntp, symba_tp1P, symba_pl1P%helio%swifter%mass, dtl)
          IF (lencounter) CALL symba_step_recur(lclose, t, irecp, npl, nplm, ntp, symba_pl1P, symba_tp1P, dt0, eoffset, nplplenc, &
               npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type)
          sgn = 1.0_DP
          CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn)
          IF (lclose) THEN
               vbs(:) = symba_pl1P%helio%swifter%vb(:)
               DO i = 1, nplplenc
                    IF (((plplenc_list(i)%status == ACTIVE) .AND.                                                                 &
                        (plplenc_list(i)%pl1P%levelg >= ireci) .AND.                                                              &
                        (plplenc_list(i)%pl2P%levelg >= ireci))) THEN
                        ! Create if statement to check for collisions (LS12) or merger depending on flag lfrag in param.in
                        ! Determines collisional regime if lfrag=.TRUE. for close encounter planets
                        ! CALL symba_frag_pl(...)
                        ! Determines if close encounter leads to merger if lfrag=.FALSE.   
                         IF (lfragmentation) THEN
                            CALL symba_fragmentation_pl(t, dtl, i, nplplenc, plplenc_list, nmergeadd, nmergesub,& 
                                   mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type)                                                       
                         ELSE
                            CALL symba_merge_pl(t, dtl, i, nplplenc, plplenc_list, nmergeadd, nmergesub,   & 
                                   mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type)
                         END IF
                     END IF
               END DO
               DO i = 1, npltpenc
                    IF ((pltpenc_list(i)%status == ACTIVE) .AND.                                                                  &
                        (pltpenc_list(i)%plP%levelg >= ireci) .AND.                                                               &
                        (pltpenc_list(i)%tpP%levelg >= ireci))                                                                    &
                         CALL symba_merge_tp(t, dtl, i, npltpenc, pltpenc_list, vbs, encounter_file, out_type)
               END DO
          END IF
          DO i = 1, nplplenc
               IF (plplenc_list(i)%pl1P%levelg == irecp) plplenc_list(i)%pl1P%levelg = ireci
               IF (plplenc_list(i)%pl2P%levelg == irecp) plplenc_list(i)%pl2P%levelg = ireci
               IF (plplenc_list(i)%level == irecp) plplenc_list(i)%level = ireci
          END DO
          DO i = 1, npltpenc
               IF (pltpenc_list(i)%plP%levelg == irecp) pltpenc_list(i)%plP%levelg = ireci
               IF (pltpenc_list(i)%tpP%levelg == irecp) pltpenc_list(i)%tpP%levelg = ireci
               IF (pltpenc_list(i)%level == irecp) pltpenc_list(i)%level = ireci
          END DO
     ELSE
          DO j = 1, NTENC
               icflg = 0
               ! OpenMP parallelization added by D. Minton
               !$OMP PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) &
               !$OMP PRIVATE(i,symba_pliP,symba_pljP,swifter_pliP,swifter_pljP,xr,vr,lencounter) &
               !$OMP SHARED(plplenc_list,nplplenc,irecp,icflg,ireci,dtl)
               DO i = 1, nplplenc
                    IF ((plplenc_list(i)%status == ACTIVE) .AND. (plplenc_list(i)%level == ireci)) THEN
                         symba_pliP => plplenc_list(i)%pl1P
                         symba_pljP => plplenc_list(i)%pl2P
                         swifter_pliP => symba_pliP%helio%swifter
                         swifter_pljP => symba_pljP%helio%swifter
                         xr(:) = swifter_pljP%xh(:) - swifter_pliP%xh(:)
                         vr(:) = swifter_pljP%vb(:) - swifter_pliP%vb(:)
                         CALL symba_chk(xr(:), vr(:), swifter_pliP%rhill, swifter_pljP%rhill, dtl, irecp, lencounter,             &
                              plplenc_list(i)%lvdotr)
                         IF (lencounter) THEN
                              !$OMP CRITICAL
                              icflg = 1
                              symba_pliP%levelg = irecp
                              symba_pliP%levelm = MAX(irecp, symba_pliP%levelm)
                              symba_pljP%levelg = irecp
                              symba_pljP%levelm = MAX(irecp, symba_pljP%levelm)
                              plplenc_list(i)%level = irecp
                              !$OMP END CRITICAL
                         END IF
                    END IF
               END DO
               DO i = 1, npltpenc
                    IF ((pltpenc_list(i)%status == ACTIVE) .AND. (pltpenc_list(i)%level == ireci)) THEN
                         symba_pliP => pltpenc_list(i)%plP
                         symba_tpP => pltpenc_list(i)%tpP
                         swifter_pliP => symba_pliP%helio%swifter
                         swifter_tpP => symba_tpP%helio%swifter
                         xr(:) = swifter_tpP%xh(:) - swifter_pliP%xh(:)
                         vr(:) = swifter_tpP%vb(:) - swifter_pliP%vb(:)
                         CALL symba_chk(xr(:), vr(:), swifter_pliP%rhill, 0.0_DP, dtl, irecp, lencounter, pltpenc_list(i)%lvdotr)
                         IF (lencounter) THEN
                              icflg = 1
                              symba_pliP%levelg = irecp
                              symba_pliP%levelm = MAX(irecp, symba_pliP%levelm)
                              symba_tpP%levelg = irecp
                              symba_tpP%levelm = MAX(irecp, symba_tpP%levelm)
                              pltpenc_list(i)%level = irecp
                         END IF
                    END IF
               END DO
               lencounter = (icflg == 1)
               sgn = 1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn)
               sgn = -1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn)
               CALL symba_helio_drift(ireci, npl, symba_pl1P, dtl)
               IF (ntp > 0) CALL symba_helio_drift_tp(ireci, ntp, symba_tp1P, symba_pl1P%helio%swifter%mass, dtl)
               IF (lencounter) CALL symba_step_recur(lclose, t, irecp, npl, nplm, ntp, symba_pl1P, symba_tp1P, dt0, eoffset,      &
                    nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list,           &
                    encounter_file, out_type)
               sgn = 1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn)
               sgn = -1.0_DP
               CALL symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn)
               IF (lclose) THEN
                    vbs(:) = symba_pl1P%helio%swifter%vb(:)
                    DO i = 1, nplplenc
                         IF ((plplenc_list(i)%status == ACTIVE) .AND.                                                             &
                             (plplenc_list(i)%pl1P%levelg >= ireci) .AND.                                                         &
                             (plplenc_list(i)%pl2P%levelg >= ireci))                                                              &
                              CALL symba_merge_pl(t, dtl, i, nplplenc, plplenc_list, nmergeadd, nmergesub, mergeadd_list,         &
                                   mergesub_list, eoffset, vbs, encounter_file, out_type)
                    END DO
                    DO i = 1, npltpenc
                         IF ((pltpenc_list(i)%status == ACTIVE) .AND.                                                             &
                             (pltpenc_list(i)%plP%levelg >= ireci) .AND.                                                          &
                             (pltpenc_list(i)%tpP%levelg >= ireci))                                                               &
                              CALL symba_merge_tp(t, dtl, i, npltpenc, pltpenc_list, vbs, encounter_file, out_type)
                    END DO
               END IF
               DO i = 1, nplplenc
                    IF (plplenc_list(i)%pl1P%levelg == irecp) plplenc_list(i)%pl1P%levelg = ireci
                    IF (plplenc_list(i)%pl2P%levelg == irecp) plplenc_list(i)%pl2P%levelg = ireci
                    IF (plplenc_list(i)%level == irecp) plplenc_list(i)%level = ireci
               END DO
               DO i = 1, npltpenc
                    IF (pltpenc_list(i)%plP%levelg == irecp) pltpenc_list(i)%plP%levelg = ireci
                    IF (pltpenc_list(i)%tpP%levelg == irecp) pltpenc_list(i)%tpP%levelg = ireci
                    IF (pltpenc_list(i)%level == irecp) pltpenc_list(i)%level = ireci
               END DO
          END DO
     END IF

     RETURN

END SUBROUTINE symba_step_recur
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
