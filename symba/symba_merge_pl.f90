!**********************************************************************************************************************************
!
!  Unit Name   : symba_merge_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check whether or not bodies are colliding or on collision path
!                Set up the merger for symba_discard_merge_pl
!                
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
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list

! Internals
     LOGICAL(LGT)              :: lmerge
     INTEGER(I4B)              :: i, j, k, id1, id2, stat1, stat2, index1, index2, index_keep, index_rm, indexchild
     REAL(DP)                  :: r2, rlim, rlim2, vdotr, tcr2, dt2, mtot, a, e, q, m1, m2, mtmp, mmax, eold, enew, rad1, rad2, & 
     mass1, mass2
     REAL(DP), DIMENSION(NDIM) :: xr, vr, x1, v1, x2, v2, xnew, vnew
     TYPE(swiftest_pl)         :: swiftest_plA
     TYPE(swiftest_tp)         :: swiftest_tpA
     TYPE(symba_pl)            :: symba_plA
     TYPE(symba_tp)            :: symba_tpA
     !TYPE(swifter_pl), POINTER :: swifter_pliP, swifter_pljP, swifter_plP
     !TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP, symba_plP

! Executable code
     lmerge = .FALSE.

     index1 = plplenc_list%id1(index)
     index2 = plplenc_list%id2(index)
     rlim = symba_plA%helio%swiftest%radius(index1) + symba_plA%helio%swiftest%radius(index2)
     xr(:) = symba_plA%helio%swiftest%xh(:,index2) - symba_plA%helio%swiftest%xh(:,index1)
     r2 = DOT_PRODUCT(xr(:), xr(:))
     rlim2 = rlim*rlim
     ! checks if bodies are actively colliding in this time step
     IF (rlim2 >= r2) THEN 
          lmerge = .TRUE.
     ! if they are not actively colliding in  this time step, 
     !checks if they are going to collide next time step based on velocities and q
     ELSE 
          vr(:) = symba_plA%helio%swiftest%vb(:,index2) - symba_plA%helio%swiftest%vb(:,index1)
          vdotr = DOT_PRODUCT(xr(:), vr(:))
          IF (plplenc_list%lvdotr(index) .AND. (vdotr > 0.0_DP)) THEN
               tcr2 = r2/DOT_PRODUCT(vr(:), vr(:))
               dt2 = dt*dt
               IF (tcr2 <= dt2) THEN
                    mtot = symba_plA%helio%swiftest%mass(index1) + symba_plA%helio%swiftest%mass(index2)
                    CALL orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
                    IF (q < rlim) lmerge = .TRUE.
               END IF
               ! if no collision is going to happen, write as close encounter, not  merger
               IF (.NOT. lmerge) THEN
                    IF (encounter_file /= "") THEN
                         id1 = symba_plA%helio%swiftest%id(index1)
                         m1 = symba_plA%helio%swiftest%mass(index1)
                         rad1 = symba_plA%helio%swiftest%radius(index1)
                         x1(:) = symba_plA%helio%swiftest%xh(;,index1)
                         v1(:) = symba_plA%helio%swiftest%vb(;,index1) - vbs(:)
                         id2 = symba_plA%helio%swiftest%id(index2)
                         m2 = symba_plA%helio%swiftest%mass(index2)
                         rad2 = symba_plA%helio%swiftest%radius(index2)
                         x2(:) = symba_plA%helio%swiftest%xh(:,index2)
                         v2(:) = symba_plA%helio%swiftest%vb(:,index2) - vbs(:)
                         CALL io_write_encounter(t, id1, id2, m1, m2, rad1, rad2, x1(:), x2(:), v1(:), v2(:), encounter_file,     &
                              out_type)
                    END IF
               END IF
          END IF
     END IF
     !Set up the merger for symba_discard_merge_pl 
     IF (lmerge) THEN
          symba_plA%lmerged(index1) = .TRUE.
          symba_plA%lmerged(index2) = .TRUE.
          !symba_pliP => symba_pliP%parentP
          m1 = symba_plA%helio%swiftest%mass(index1)
          mass1 = m1 
          rad1 = symba_plA%helio%swiftest%radius(index1)
          x1(:) = m1*symba_plA%helio%swiftest%xh(:,index1)
          v1(:) = m1*symba_plA%helio%swiftest%vb(:,index1)
          mmax = m1
          id1 = symba_plA%helio%swiftest%id(index1)
          stat1 = symba_plA%helio%swiftest%status(index1)
          !symba_plP => symba_pliP
          DO i = 1, symba_plA%nchild(index1) ! initialize an array of children
               !symba_plP => symba_plP%childP
               indexchild = ????
               mtmp = symba_plA%helio%swiftest%mass(indexchild)
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    id1 = symba_plA%helio%swiftest%id(indexchild)
                    stat1 = symba_plA%helio%swiftest%status(indexchild)
               END IF
               m1 = m1 + mtmp
               x1(:) = x1(:) + mtmp*symba_plA%helio%swiftest%xh(indexchild)
               v1(:) = v1(:) + mtmp*symba_plA%helio%swiftest%vb(indexchild)
          END DO
          x1(:) = x1(:)/m1
          v1(:) = v1(:)/m1
          !symba_pljP => symba_pljP%parentP
          m2 = symba_plA%helio%swiftest%mass(index2)
          mass2 = m2
          x2(:) = m2*symba_plA%helio%swiftest%xh(:,index2)
          v2(:) = m2*symba_plA%helio%swiftest%vb(:,index2)
          mmax = m2
          id2 = symba_plA%helio%swiftest%id(index2)
          stat2 = symba_plA%helio%swiftest%status(index2)
          !symba_plP => symba_pljP
          DO i = 1, symba_plA%nchild(index2)
               !symba_plP => symba_plP%childP
               indexchild = ????
               mtmp = symba_plA%helio%swiftest%mass(indexchild)
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    id2 = symba_plA%helio%swiftest%id(indexchild)
                    stat2 = symba_plA%helio%swiftest%status(indexchild)
               END IF
               m2 = m2 + mtmp
               x2(:) = x2(:) + mtmp*symba_plA%helio%swiftest%xh(indexchild)
               v2(:) = v2(:) + mtmp*symba_plA%helio%swiftest%vb(indexchild)
          END DO
          x2(:) = x2(:)/m2
          v2(:) = v2(:)/m2
          mtot = m1 + m2
          xnew(:) = (m1*x1(:) + m2*x2(:))/mtot
          vnew(:) = (m1*v1(:) + m2*v2(:))/mtot
          WRITE(*, *) "Merging particles ", id1, " and ", id2, " at time t = ",t
          nmergesub = nmergesub + 1
          mergesub_list%id(nmergesub) = id1
          mergesub_list%status(nmergesub) = MERGED
          mergesub_list%xh(:,nmergesub) = x1(:)
          mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
          mergesub_list%mass(nmergesub) = mass1
          mergesub_list%radius(nmergesub) = rad1
          nmergesub = nmergesub + 1
          mergesub_list%id(nmergesub) = id2
          mergesub_list%status(nmergesub) = MERGED
          mergesub_list%xh(:,nmergesub) = x2(:)
          mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
          mergesub_list%mass(nmergesub) = mass2
          mergesub_list%radius(nmergesub) = rad2
          nmergeadd = nmergeadd + 1

          IF (m2 > m1) THEN
               index_keep = index2
               index_rm = index1
               mergeadd_list%id(nmergeadd) = id2
               mergeadd_list%status(nmergeadd) = stat2

          ELSE
               index_keep = index1
               index_rm = index2
               mergeadd_list%id(nmergeadd) = id1
               mergeadd_list%status(nmergeadd) = stat1

          END IF
          mergeadd_list%ncomp(nmergeadd) = 2
          mergeadd_list%xh(:,nmergeadd) = xnew(:)
          mergeadd_list%vh(:,nmergeadd) = vnew(:) - vbs(:)
          eold = 0.5_DP*(m1*DOT_PRODUCT(v1(:), v1(:)) + m2*DOT_PRODUCT(v2(:), v2(:)))
          xr(:) = x2(:) - x1(:)
          eold = eold - m1*m2/SQRT(DOT_PRODUCT(xr(:), xr(:)))
          enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
          eoffset = eoffset + eold - enew
          DO k = 1, nplplenc
               IF (plplenc_list%status(k) == ACTIVE) THEN
                    !symba_pliP => plplenc_list(index)%pl1P%parentP
                    index1 = plplenc_list%id1(index)
                    DO i = 0, symba_plA%nchild(index1)
                         !symba_pljP => plplenc_list(index)%pl2P%parentP
                         index2 = plplenc_list%id2(index)
                         DO j = 0, symba_plA%nchild(index2)
                              IF (index1 = plplenc_list%id1(k)) .AND. (index2 = plplenc_list%id2(k)) THEN
                                   plplenc_list%status(k) = MERGED
                              ELSE IF (index1 = plplenc_list%id2(k)) .AND. (index2 = plplenc_list%id1(k)) THEN
                                   plplenc_list%status(k) = MERGED
                              END IF
                              !symba_pljP => symba_pljP%childP
                              index1 = ???? 
                         END DO
                         !symba_pliP => symba_pliP%childP
                         index2 = ????
                    END DO
               END IF
          END DO
          symba_plA%helio%swiftest%xh(:,index_keep) = xnew(:)
          symba_plA%helio%swiftest%vb(:,index_keep) = vnew(:)

          DO i = 1, symba_plA%nchild(index_keep)
               !symba_plP => symba_plP%childP
               indexchild = ????
               !swifter_plP => symba_plP%helio%swifter
               symba_plA%helio%swiftest%xh(:,indexchild) = xnew(:)
               symba_plA%helio%swiftest%vb(:,indexchild) = vnew(:)
          END DO
          !symba_pljP => plplenc_list(index)%pl2P%parentP
          !symba_plP%childP => symba_pljP
          DO i = 0, symba_plA%nchild(index_rm)
               !symba_plP => symba_plP%childP
               !symba_plP%parentP => symba_pliP
               !swifter_plP => symba_plP%helio%swifter
               indexchild = ????
               symba_plA%helio%swiftest%xh(:,indexchild) = xnew(:)
               symba_plA%helio%swiftest%vb(:,indexchild) = vnew(:)
          END DO
          symba_plA%nchild(index_keep) = symba_plA%nchild(index_keep) + symba_plA%nchild(index_rm) + 1
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
