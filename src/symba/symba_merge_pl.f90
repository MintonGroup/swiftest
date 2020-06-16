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
!                index          : index of massive body-massive body encounter in array plplenc_list
!                nplplenc       : number of massive body-massive body encounters
!                plplenc_list   : array of massive body-massive body encounter structures
!                nmergeadd      : number of merged massive bodies to add
!                nmergesub      : number of merged massive bodies to subtract
!                mergeadd_list  : array of structures of merged massive bodies to add
!                mergesub_list  : array of structures of merged massive bodies to subtract
!                eoffset        : energy offset (net energy lost in mergers)
!                vbs            : barycentric velocity of the Sun
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : plplenc_list   : array of massive body-massive body encounter structures
!                nmergeadd      : number of merged massive bodies to add
!                nmergesub      : number of merged massive bodies to subtract
!                mergeadd_list  : array of structures of merged massive bodies to add
!                mergesub_list  : array of structures of merged massive bodies to subtract
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
SUBROUTINE symba_merge_pl(t, dt, index_enc, nplplenc, plplenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, &
     vbs, encounter_file, out_type, npl, symba_plA)

! Modules
     use swiftest, EXCEPT_THIS_ONE => symba_merge_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc, nplplenc
     INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, npl
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
     LOGICAL(LGT)                 :: lmerge
     INTEGER(I4B)                 :: i, j, k, stat1, stat2, index1, index2, index_keep, index_rm, indexchild
     INTEGER(I4B)                 :: index1_child, index2_child, index1_parent, index2_parent, index_big1, index_big2
     INTEGER(I4B)                 :: name1, name2
     REAL(DP)                     :: r2, rlim, rlim2, vdotr, tcr2, dt2, mtot, a, e, q, m1, m2, mtmp, mmax 
     REAL(DP)                     :: eold, enew, rad1, rad2, mass1, mass2
     REAL(DP), DIMENSION(NDIM)    :: xr, vr, x1, v1, x2, v2, xnew, vnew
     INTEGER(I4B), DIMENSION(npl) :: array_index1_child, array_index2_child, array_keep_child, array_rm_child

! Executable code
     lmerge = .FALSE.

     index1 = plplenc_list%index1(index_enc)
     index2 = plplenc_list%index2(index_enc)

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
          IF (plplenc_list%lvdotr(index_enc) .AND. (vdotr > 0.0_DP)) THEN 
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
                         name1 = symba_plA%helio%swiftest%name(index1)
                         m1 = symba_plA%helio%swiftest%mass(index1)
                         rad1 = symba_plA%helio%swiftest%radius(index1)
                         x1(:) = symba_plA%helio%swiftest%xh(:,index1)
                         v1(:) = symba_plA%helio%swiftest%vb(:,index1) - vbs(:)
                         name2 = symba_plA%helio%swiftest%name(index2)
                         m2 = symba_plA%helio%swiftest%mass(index2)
                         rad2 = symba_plA%helio%swiftest%radius(index2)
                         x2(:) = symba_plA%helio%swiftest%xh(:,index2)
                         v2(:) = symba_plA%helio%swiftest%vb(:,index2) - vbs(:)

                         CALL io_write_encounter(t, name1, name2, m1, m2, rad1, rad2, x1(:), x2(:), &
                              v1(:), v2(:), encounter_file, out_type)
                    END IF
               END IF
          END IF
     END IF
     !Set up the merger for symba_discard_merge_pl 
     IF (lmerge) THEN
          symba_plA%lmerged(index1) = .TRUE.
          symba_plA%lmerged(index2) = .TRUE.
          index1_parent = symba_plA%index_parent(index1)
          m1 = symba_plA%helio%swiftest%mass(index1_parent)
          mass1 = m1 
          rad1 = symba_plA%helio%swiftest%radius(index1_parent)
          x1(:) = m1*symba_plA%helio%swiftest%xh(:,index1_parent)
          v1(:) = m1*symba_plA%helio%swiftest%vb(:,index1_parent)
          mmax = m1
          name1 = symba_plA%helio%swiftest%name(index1_parent)
          index_big1 = index1_parent
          stat1 = symba_plA%helio%swiftest%status(index1_parent)
          array_index1_child(1:npl) = symba_plA%index_child(1:npl,index1_parent)
          DO i = 1, symba_plA%nchild(index1_parent) ! initialize an array of children
               index1_child = array_index1_child(i)
               mtmp = symba_plA%helio%swiftest%mass(index1_child)
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    name1 = symba_plA%helio%swiftest%name(index1_child)
                    index_big1 = index1_child
                    stat1 = symba_plA%helio%swiftest%status(index1_child)
               END IF
               m1 = m1 + mtmp
               x1(:) = x1(:) + mtmp*symba_plA%helio%swiftest%xh(:,index1_child)
               v1(:) = v1(:) + mtmp*symba_plA%helio%swiftest%vb(:,index1_child)
          END DO
          x1(:) = x1(:)/m1
          v1(:) = v1(:)/m1
          index2_parent = symba_plA%index_parent(index2)
          m2 = symba_plA%helio%swiftest%mass(index2_parent)
          mass2 = m2
          rad2 = symba_plA%helio%swiftest%radius(index2_parent)
          x2(:) = m2*symba_plA%helio%swiftest%xh(:,index2_parent)
          v2(:) = m2*symba_plA%helio%swiftest%vb(:,index2_parent)
          mmax = m2
          name2 = symba_plA%helio%swiftest%name(index2_parent)
          index_big2 = index2_parent
          stat2 = symba_plA%helio%swiftest%status(index2_parent)
          array_index2_child(1:npl) = symba_plA%index_child(1:npl,index2_parent)
          DO i = 1, symba_plA%nchild(index2_parent)
               index2_child = array_index2_child(i)
               mtmp = symba_plA%helio%swiftest%mass(index2_child)
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    name2 = symba_plA%helio%swiftest%name(index2_child)
                    index_big2 = index2_child
                    stat2 = symba_plA%helio%swiftest%status(index2_child)
               END IF
               m2 = m2 + mtmp
               x2(:) = x2(:) + mtmp*symba_plA%helio%swiftest%xh(:,index2_child)
               v2(:) = v2(:) + mtmp*symba_plA%helio%swiftest%vb(:,index2_child)
          END DO
          x2(:) = x2(:)/m2
          v2(:) = v2(:)/m2
          mtot = m1 + m2
          xnew(:) = (m1*x1(:) + m2*x2(:))/mtot
          vnew(:) = (m1*v1(:) + m2*v2(:))/mtot
          WRITE(*, *) "Merging particles ", name1, " and ", name2, " at time t = ",t
          nmergesub = nmergesub + 1
          mergesub_list%name(nmergesub) = name1
          mergesub_list%status(nmergesub) = MERGED
          mergesub_list%xh(:,nmergesub) = x1(:)
          mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
          mergesub_list%mass(nmergesub) = mass1
          mergesub_list%radius(nmergesub) = rad1
          nmergesub = nmergesub + 1
          mergesub_list%name(nmergesub) = name2
          mergesub_list%status(nmergesub) = MERGED
          mergesub_list%xh(:,nmergesub) = x2(:)
          mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
          mergesub_list%mass(nmergesub) = mass2
          mergesub_list%radius(nmergesub) = rad2
          nmergeadd = nmergeadd + 1
          IF (m2 > m1) THEN
               index_keep = index_big2
               index_rm = index_big1
               mergeadd_list%name(nmergeadd) = name2
               mergeadd_list%status(nmergeadd) = stat2

          ELSE
               index_keep = index_big1
               index_rm = index_big2
               mergeadd_list%name(nmergeadd) = name1
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

          !WRITE(*,*) "symba_merge_pl.f90 name", mergeadd_list%name(nmergeadd)
          !WRITE(*,*) "symba_merge_pl.f90 xh", mergeadd_list%xh(:,nmergeadd)
          !WRITE(*,*) "symba_merge_pl.f90 vh", mergeadd_list%vh(:,nmergeadd)
          !WRITE(*,*) "symba_merge_pl.f90 eoffset", eoffset
          DO k = 1, nplplenc                                          !go through the encounter list and for particles actively encoutering, get their children
               IF (plplenc_list%status(k) == ACTIVE) THEN
                    DO i = 0, symba_plA%nchild(index1_parent)
                         IF (i == 0) THEN 
                              index1_child = index1_parent
                         ELSE
                              index1_child = array_index1_child(i)
                         END IF 
                         DO j = 0, symba_plA%nchild(index2_parent)
                              IF (j == 0) THEN
                                   index2_child = index2_parent
                              ELSE
                                   index2_child = array_index2_child(j)
                              END IF
                              IF ((index1_child == plplenc_list%index1(k)) .AND. (index2_child == plplenc_list%index2(k))) THEN
                                   plplenc_list%status(k) = MERGED
                              ELSE IF ((index1_child == plplenc_list%index2(k)) .AND. (index2_child == plplenc_list%index1(k))) THEN
                                   plplenc_list%status(k) = MERGED
                              END IF
                         END DO
                    END DO
               END IF
          END DO
          symba_plA%helio%swiftest%xh(:,index1_parent) = xnew(:)
          symba_plA%helio%swiftest%vb(:,index1_parent) = vnew(:)
          symba_plA%helio%swiftest%xh(:,index2_parent) = xnew(:) 
          symba_plA%helio%swiftest%vb(:,index2_parent) = vnew(:) 
          array_keep_child(1:npl) = symba_plA%index_child(1:npl,index1_parent)
          DO i = 1, symba_plA%nchild(index1_parent)
               indexchild = array_keep_child(i)
               symba_plA%helio%swiftest%xh(:,indexchild) = xnew(:)
               symba_plA%helio%swiftest%vb(:,indexchild) = vnew(:)
          END DO

          symba_plA%index_child((symba_plA%nchild(index1_parent)+1),index1_parent) = index2_parent
          array_rm_child(1:npl) = symba_plA%index_child(1:npl,index2_parent)
          symba_plA%index_parent(index2) = index1_parent

          DO i = 1, symba_plA%nchild(index2_parent)
               symba_plA%index_parent(array_rm_child(i)) = index1_parent
               indexchild = array_rm_child(i)
               symba_plA%helio%swiftest%xh(:,indexchild) = xnew(:)
               symba_plA%helio%swiftest%vb(:,indexchild) = vnew(:)
          END DO
          DO i = 1, symba_plA%nchild(index2_parent)
               symba_plA%index_child(symba_plA%nchild(index1_parent)+i+1,index1_parent)= array_rm_child(i)
          END DO 
          symba_plA%nchild(index1_parent) = symba_plA%nchild(index1_parent) + symba_plA%nchild(index2_parent) + 1
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
