!**********************************************************************************************************************************
!
!  Unit Name   : symba_merge_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check for merger between planets in SyMBA
!
!  Input
!    Arguments : t              : time
!                dt             : time step
!                index          : index of planet-planet encounter in array plplenc_list
!                nplplenc       : number of planet-planet encounters
!                plplenc_list   : array of planet-planet encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                eoffset        : energy offset (net energy lost in mergers)
!                vbs            : barycentric velocity of the Sun
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : plplenc_list   : array of planet-planet encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                eoffset        : energy offset (net energy lost in mergers)
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL symba_merge_pl(t, dt, index, nplplenc, plplenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list,
!                                    eoffset, vbs, encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_merge_pl(t, dt, index, nplplenc, plplenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, &
     encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_merge_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index, nplplenc
     INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
     TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list

! Internals
     LOGICAL(LGT)              :: lmerge
     INTEGER(I4B)              :: i, j, k, id1, id2, stat1, stat2
     REAL(DP)                  :: r2, rlim, rlim2, vdotr, tcr2, dt2, mtot, a, e, q, m1, m2, mtmp, mmax, eold, enew, rad1, rad2
     REAL(DP), DIMENSION(NDIM) :: xr, vr, x1, v1, x2, v2, xnew, vnew
     TYPE(swifter_pl), POINTER :: swifter_pliP, swifter_pljP, swifter_plP
     TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP, symba_plP

! Executable code
     lmerge = .FALSE.
     symba_pliP => plplenc_list(index)%pl1P
     symba_pljP => plplenc_list(index)%pl2P
     swifter_pliP => symba_pliP%helio%swifter
     swifter_pljP => symba_pljP%helio%swifter
     rlim = swifter_pliP%radius + swifter_pljP%radius
     xr(:) = swifter_pljP%xh(:) - swifter_pliP%xh(:)
     r2 = DOT_PRODUCT(xr(:), xr(:))
     rlim2 = rlim*rlim
     IF (rlim2 >= r2) THEN
          lmerge = .TRUE.
     ELSE
          vr(:) = swifter_pljP%vb(:) - swifter_pliP%vb(:)
          vdotr = DOT_PRODUCT(xr(:), vr(:))
          IF (plplenc_list(index)%lvdotr .AND. (vdotr > 0.0_DP)) THEN
               tcr2 = r2/DOT_PRODUCT(vr(:), vr(:))
               dt2 = dt*dt
               IF (tcr2 <= dt2) THEN
                    mtot = swifter_pliP%mass + swifter_pljP%mass
                    CALL orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
                    IF (q < rlim) lmerge = .TRUE.
               END IF
               IF (.NOT. lmerge) THEN
                    IF (encounter_file /= "") THEN
                         id1 = swifter_pliP%id
                         m1 = swifter_pliP%mass
                         rad1 = swifter_pliP%radius
                         x1(:) = swifter_pliP%xh(:)
                         v1(:) = swifter_pliP%vb(:) - vbs(:)
                         id2 = swifter_pljP%id
                         m2 = swifter_pljP%mass
                         rad2 = swifter_pljP%radius
                         x2(:) = swifter_pljP%xh(:)
                         v2(:) = swifter_pljP%vb(:) - vbs(:)
                         CALL io_write_encounter(t, id1, id2, m1, m2, rad1, rad2, x1(:), x2(:), v1(:), v2(:), encounter_file,     &
                              out_type)
                    END IF
               END IF
          END IF
     END IF
     IF (lmerge) THEN
          symba_pliP%lmerged = .TRUE.
          symba_pljP%lmerged = .TRUE.
          symba_pliP => symba_pliP%parentP
          m1 = symba_pliP%helio%swifter%mass
          x1(:) = m1*symba_pliP%helio%swifter%xh(:)
          v1(:) = m1*symba_pliP%helio%swifter%vb(:)
          mmax = m1
          id1 = symba_pliP%helio%swifter%id
          stat1 = symba_pliP%helio%swifter%status
          symba_plP => symba_pliP
          DO i = 1, symba_pliP%nchild
               symba_plP => symba_plP%childP
               mtmp = symba_plP%helio%swifter%mass
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    id1 = symba_plP%helio%swifter%id
                    stat1 = symba_plP%helio%swifter%status
               END IF
               m1 = m1 + mtmp
               x1(:) = x1(:) + mtmp*symba_plP%helio%swifter%xh(:)
               v1(:) = v1(:) + mtmp*symba_plP%helio%swifter%vb(:)
          END DO
          x1(:) = x1(:)/m1
          v1(:) = v1(:)/m1
          symba_pljP => symba_pljP%parentP
          m2 = symba_pljP%helio%swifter%mass
          x2(:) = m2*symba_pljP%helio%swifter%xh(:)
          v2(:) = m2*symba_pljP%helio%swifter%vb(:)
          mmax = m2
          id2 = symba_pljP%helio%swifter%id
          stat2 = symba_pljP%helio%swifter%status
          symba_plP => symba_pljP
          DO i = 1, symba_pljP%nchild
               symba_plP => symba_plP%childP
               mtmp = symba_plP%helio%swifter%mass
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    id2 = symba_plP%helio%swifter%id
                    stat2 = symba_plP%helio%swifter%status
               END IF
               m2 = m2 + mtmp
               x2(:) = x2(:) + mtmp*symba_plP%helio%swifter%xh(:)
               v2(:) = v2(:) + mtmp*symba_plP%helio%swifter%vb(:)
          END DO
          x2(:) = x2(:)/m2
          v2(:) = v2(:)/m2
          mtot = m1 + m2
          xnew(:) = (m1*x1(:) + m2*x2(:))/mtot
          vnew(:) = (m1*v1(:) + m2*v2(:))/mtot
          WRITE(*, *) "Merging particles ", id1, " and ", id2, " at time t = ",t
          nmergesub = nmergesub + 1
          mergesub_list(nmergesub)%id = id1
          mergesub_list(nmergesub)%status = MERGED
          mergesub_list(nmergesub)%xh(:) = x1(:)
          mergesub_list(nmergesub)%vh(:) = v1(:) - vbs(:)
          nmergesub = nmergesub + 1
          mergesub_list(nmergesub)%id = id2
          mergesub_list(nmergesub)%status = MERGED
          mergesub_list(nmergesub)%xh(:) = x2(:)
          mergesub_list(nmergesub)%vh(:) = v2(:) - vbs(:)
          nmergeadd = nmergeadd + 1
          IF (m2 > m1) THEN
               mergeadd_list(nmergeadd)%id = id2
               mergeadd_list(nmergeadd)%status = stat2
          ELSE
               mergeadd_list(nmergeadd)%id = id1
               mergeadd_list(nmergeadd)%status = stat1
          END IF
          mergeadd_list(nmergeadd)%ncomp = 2
          mergeadd_list(nmergeadd)%xh(:) = xnew(:)
          mergeadd_list(nmergeadd)%vh(:) = vnew(:) - vbs(:)
          eold = 0.5_DP*(m1*DOT_PRODUCT(v1(:), v1(:)) + m2*DOT_PRODUCT(v2(:), v2(:)))
          xr(:) = x2(:) - x1(:)
          eold = eold - m1*m2/SQRT(DOT_PRODUCT(xr(:), xr(:)))
          enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
          eoffset = eoffset + eold - enew
          DO k = 1, nplplenc
               IF (plplenc_list(k)%status == ACTIVE) THEN
                    symba_pliP => plplenc_list(index)%pl1P%parentP
                    DO i = 0, symba_pliP%nchild
                         symba_pljP => plplenc_list(index)%pl2P%parentP
                         DO j = 0, symba_pljP%nchild
                              IF (ASSOCIATED(plplenc_list(k)%pl1P, symba_pliP) .AND.                                              &
                                  ASSOCIATED(plplenc_list(k)%pl2P, symba_pljP)) THEN
                                   plplenc_list(k)%status = MERGED
                              ELSE IF (ASSOCIATED(plplenc_list(k)%pl1P, symba_pljP) .AND.                                         &
                                       ASSOCIATED(plplenc_list(k)%pl2P, symba_pliP)) THEN
                                   plplenc_list(k)%status = MERGED
                              END IF
                              symba_pljP => symba_pljP%childP
                         END DO
                         symba_pliP => symba_pliP%childP
                    END DO
               END IF
          END DO
          symba_pliP => plplenc_list(index)%pl1P%parentP
          symba_plP => symba_pliP
          swifter_plP => symba_plP%helio%swifter
          swifter_plP%xh(:) = xnew(:)
          swifter_plP%vb(:) = vnew(:)
          DO i = 1, symba_pliP%nchild
               symba_plP => symba_plP%childP
               swifter_plP => symba_plP%helio%swifter
               swifter_plP%xh(:) = xnew(:)
               swifter_plP%vb(:) = vnew(:)
          END DO
          symba_pljP => plplenc_list(index)%pl2P%parentP
          symba_plP%childP => symba_pljP
          DO i = 0, symba_pljP%nchild
               symba_plP => symba_plP%childP
               symba_plP%parentP => symba_pliP
               swifter_plP => symba_plP%helio%swifter
               swifter_plP%xh(:) = xnew(:)
               swifter_plP%vb(:) = vnew(:)
          END DO
          symba_pliP%nchild = symba_pliP%nchild + symba_pljP%nchild + 1
     END IF

     RETURN

END SUBROUTINE symba_merge_pl
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
