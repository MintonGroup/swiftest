!**********************************************************************************************************************************
!
!  Unit Name   : symba_casemerge
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge planets
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-planet encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_casemerge(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casemerge (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
     encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
     array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)

! Modules
     USE swiftest
     USE module_swiftest
     USE helio
     USE symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_casemerge
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc
     INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     INTEGER(I4B), DIMENSION(npl), INTENT(INOUT)      :: array_index1_child, array_index2_child

! Internals
 
     INTEGER(I4B)                 :: model, nres
     REAL(DP)                     :: mres, rres, pres, vres
     INTEGER(I4B)                 :: i, j, k, stat1, stat2, index1, index2, indexchild
     INTEGER(I4B)                 :: index1_child, index2_child, index1_parent, index2_parent, index_big1, index_big2
     INTEGER(I4B)                 :: name1, name2
     REAL(DP)                     :: r2, rlim, rlim2, vdotr, tcr2, dt2, mtot, a, e, q, mtmp, mmax 
     REAL(DP)                     :: eold, enew, mass1, mass2
     REAL(DP), DIMENSION(NDIM)    :: xr, vr, xnew, vnew
     INTEGER(I4B), DIMENSION(npl) :: array_keep_child, array_rm_child



! Executable code
               index1 = plplenc_list%index1(index_enc)
               index2 = plplenc_list%index2(index_enc)
               index1_parent = symba_plA%index_parent(index1)
               index2_parent = symba_plA%index_parent(index2)
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
                    mergeadd_list%name(nmergeadd) = name2
                    mergeadd_list%status(nmergeadd) = stat2

               ELSE
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

               DO k = 1, nplplenc
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
                                   ELSE IF ((index1_child == plplenc_list%index2(k)) .AND. &
                                        (index2_child == plplenc_list%index1(k))) THEN
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

               ! The children of parent one are the children we are keeping
               array_keep_child(1:npl) = symba_plA%index_child(1:npl,index1_parent)
               ! Go through the children of the kept parent and add those children to the array of kept children
               DO i = 1, symba_plA%nchild(index1_parent)
                    indexchild = array_keep_child(i)
                    symba_plA%helio%swiftest%xh(:,indexchild) = xnew(:)
                    symba_plA%helio%swiftest%vb(:,indexchild) = vnew(:)
               END DO
               ! the removed parent is assigned as a new child to the list of children of the kept parent
               ! gives kept parent a new child 
               symba_plA%index_child((symba_plA%nchild(index1_parent)+1),index1_parent) = index2_parent
               array_rm_child(1:npl) = symba_plA%index_child(1:npl,index2_parent)
               ! the parent of the removed parent is assigned to be the kept parent 
               ! gives removed parent a new parent
               symba_plA%index_parent(index2) = index1_parent
               ! go through the children of the removed parent and add those children to the array of removed children 
               DO i = 1, symba_plA%nchild(index2_parent)
                    symba_plA%index_parent(array_rm_child(i)) = index1_parent
                    indexchild = array_rm_child(i)
                    symba_plA%helio%swiftest%xh(:,indexchild) = xnew(:)
                    symba_plA%helio%swiftest%vb(:,indexchild) = vnew(:)
               END DO
               ! go through the children of the removed parent and add those children to the list of children of the kept parent
               DO i = 1, symba_plA%nchild(index2_parent)
                    symba_plA%index_child(symba_plA%nchild(index1_parent)+i+1,index1_parent)= array_rm_child(i)
               END DO 
               ! updates the number of children of the kept parent
               symba_plA%nchild(index1_parent) = symba_plA%nchild(index1_parent) + symba_plA%nchild(index2_parent) + 1

     RETURN 
END SUBROUTINE symba_casemerge

