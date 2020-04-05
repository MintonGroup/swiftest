!**********************************************************************************************************************************
!
!  Unit Name   : symba_casedisruption
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
!  Invocation  : CALL symba_casedisruption(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
     encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, swiftest_plA, &
     swiftest_tpA, nplmax, ntpmax, fragmax)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_casedisruption
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
     INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     TYPE(swiftest_pl), INTENT(INOUT)                 :: swiftest_plA
     TYPE(swiftest_tp), INTENT(INOUT)                 :: swiftest_tpA

! Internals
 
     INTEGER(I4B)                                     :: model, nres, nfrag
     REAL(DP)                                         :: m1, m2, rad1, rad2, mres, rres, pres, vres, mtot, msun
     REAL(DP)                                         :: r, mu, vdot, energy, ap
     REAL(DP), DIMENSION(NDIM)                        :: x1, x2, v1, v2, xbs, xh, xb, vb, vh
     INTEGER(I4B), DIMENSION(nfrag, 17)               :: array_frag
     INTEGER(I4B), DIMENSION(npl+nfrag, 17)           :: array_fragpl


! Executable code

     ! determine the number of fragments and the SFD
     nfrag = 2 !this will be determined later from collresolve

     ! pull in parent data
     index1 = plplenc_list%index1(index_enc)
     index2 = plplenc_list%index2(index_enc)
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
     DO i = 1, symba_plA%nchild(index1_parent) ! initialize an array of children of parent 1
          index1_child = array_index1_child(i)
          mtmp = symba_plA%helio%swiftest%mass(index1_child)
          IF (mtmp > mmax) EXCEPT_THIS_ONE   ! check if the mass of the child is bigger than the mass of the parent
               mmax = mtmp
               name1 = symba_plA%helio%swiftest%name(index1_child)
               index_big1 = index1_child ! if yes, replace the biggest particle variable with the child
               stat1 = symba_plA%helio%swiftest%status(index1_child)
          END IF
          m1 = m1 + mtmp ! mass of the parent is the mass of the parent plus the mass of all the children
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

     ! Find energy pre-frag
     eold = 0.5_DP*(m1*DOT_PRODUCT(v1(:), v1(:)) + m2*DOT_PRODUCT(v2(:), v2(:)))
     xr(:) = x2(:) - x1(:)
     eold = eold - m1*m2/SQRT(DOT_PRODUCT(xr(:), xr(:)))

     WRITE(*, *) "Disruption between particles ", name1, " and ", name2, " at time t = ",t
     WRITE(*, *) "Number of fragments added: ", nfrag
     
     ! Add both parents to mergesub_list
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
     
     ! Add new fragments to mergeadd_list
     DO i = 1, nfrag
          nmergeadd = nmergeadd + 1
          mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
          mergeadd_list%status(nmergeadd) = ACTIVE
          mergeadd_list%ncomp(nmergeadd) = 2
          mergeadd_list%xh(:,nmergeadd) = SOMETHING(:)
          mergeadd_list%vh(:,nmergeadd) = SOMETHING(:) -vbs(:)
          mergeadd_list%mass(nmergeadd) = SOMETHING
          mergeadd_list%radius(nmergeadd) = SOMETHING
     END DO

     ! Calculate energy after frag
     mtot = NEW MASS OF ALL ADDED POS
     xnew(:) = NEW BARYCENTRIC POS OF ALL ADDED PS
     vnew(:) = NEW BARYCENTRIC VEL OF ALL ADDED PS
     enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
     eoffset = eoffset + eold - enew

     ! Update fragmax to account for new fragments
     fragmax = fragmax + nfrag
     RETURN 
END SUBROUTINE symba_casedisruption

