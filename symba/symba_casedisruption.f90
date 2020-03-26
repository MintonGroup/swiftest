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
     swiftest_tpA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_casedisruption
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc
     INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc
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
 
     INTEGER(I4B)                                     :: model, nres
     REAL(DP)                                         :: m1, m2, rad1, rad2, mres, rres, pres, vres, mtot, msun
     REAL(DP), DIMENSION(NDIM)                        :: x1, x2, v1, v2, vbs, xbs
     INTEGER(I4B), DIMENSION(nfrag, 17)               :: array_frag
     INTEGER(I4B), DIMENSION(npl+nfrag, 17)           :: array_fragpl


! Executable code

     ! determine the number of fragments and the SFD
     nfrag = 1 !this will be determined later from collresolve

     ! pull in parent data
     msun = symba_plA%helio%swiftest%mass(1)
     vbs(:) = symba_plA%helio%swiftest%vb(:,1)
     xbs(:) = symba_plA%helio%swiftest%xb(:,1)
     index1 = plplenc_list%index1(index_enc)
     index2 = plplenc_list%index2(index_enc)
     symba_plA%lmerged(index1) = .TRUE.
     symba_plA%lmerged(index2) = .TRUE.
     index1_parent = symba_plA%index_parent(index1)
     index2_parent = symba_plA%index_parent(index2)
     name1 = symba_plA%helio%swiftest%name(index1_parent)
     m1 = symba_plA%helio%swiftest%mass(index1_parent)
     rad1 = symba_plA%helio%swiftest%radius(index1_parent)
     x1(:) = m1*symba_plA%helio%swiftest%xh(:,index1_parent)
     v1(:) = m1*symba_plA%helio%swiftest%vb(:,index1_parent)
     name2 = symba_plA%helio%swiftest%name(index2_parent)
     m2 = symba_plA%helio%swiftest%mass(index2_parent)
     rad2 = symba_plA%helio%swiftest%radius(index2_parent)
     x2(:) = m2*symba_plA%helio%swiftest%xh(:,index2_parent)
     v2(:) = m2*symba_plA%helio%swiftest%vb(:,index2_parent)

     mtot = m1 + m2

     ! initialize an array of fragments
     DO i = 1, nfrag 
          array_frag(i, 1) = npl + i !name
          array_frag(i, 2) = ACTIVE !status
          array_frag(i, 3) = mtot/nfrag !mass
          array_frag(i, 4) = ((mtot * rad1**3.0_DP) / m1)**(1.0_DP/3.0_DP) !radius    
          
          xh(:) = m1*x1(:) + m2*x2(:)
          array_frag(i, 6) =  xh(1,:) !xh
          array_frag(i, 7) =  xh(2,:) !xh
          array_frag(i, 8) =  xh(3,:) !xh

          xb(:) = xh(:) + xbs(:)
          array_frag(i, 9) =  xb(1,:)  !xb
          array_frag(i, 10) =  xb(2,:) !xb
          array_frag(i, 11) =  xb(3,:) !xb

          vb(:) = m1*v1(:) + m2*v2(:)
          array_frag(i, 12) =  vb(1,:) !vb
          array_frag(i, 13) =  vb(2,:) !vb
          array_frag(i, 14) =  vb(3,:) !vb

          vh(:) = vb(:) - vbs(:)
          array_frag(i, 15) =  vh(1,:) !vh
          array_frag(i, 16) =  vh(1,:) !vh
          array_frag(i, 17) =  vh(1,:) !vh

          r = SQRT(DOT_PRODUCT(xh(:), xh(:)))
          mu = msun*mtot/(msun + mtot)
          vdot = DOT_PRODUCT(vh(:), vh(:))
          energy = -1.0_DP*msun*mtot/r + 0.5_DP*mu*vdot
          ap = -1.0_DP*msun*mtot/(2.0_DP*energy)

          array_frag(i, 5) = ap*(((mu/msun)/3.0_DP)**(1.0_DP/3.0_DP)) !rhill
     END DO 

     ! concatenate array_frag and symba_plA into new array array_fragpl
     DO i = 1, npl
          array_fragpl(i, 1) = symba_plA%helio%swiftest%name(i)
          array_fragpl(i, 2) = symba_plA%helio%swiftest%status(i)
          array_fragpl(i, 3) = symba_plA%helio%swiftest%mass(i)
          array_fragpl(i, 4) = symba_plA%helio%swiftest%radius(i)
          array_fragpl(i, 5) = symba_plA%helio%swiftest%rhill(i)
          array_fragpl(i, 6) = symba_plA%helio%swiftest%xh(1, i)
          array_fragpl(i, 7) = symba_plA%helio%swiftest%xh(2, i)
          array_fragpl(i, 8) = symba_plA%helio%swiftest%xh(3, i)
          array_fragpl(i, 9) = symba_plA%helio%swiftest%xb(1, i)
          array_fragpl(i, 10) = symba_plA%helio%swiftest%xb(2, i)
          array_fragpl(i, 11) = symba_plA%helio%swiftest%xb(3, i)
          array_fragpl(i, 12) = symba_plA%helio%swiftest%vb(1, i)
          array_fragpl(i, 13) = symba_plA%helio%swiftest%vb(2, i)
          array_fragpl(i, 14) = symba_plA%helio%swiftest%vb(3, i)
          array_fragpl(i, 15) = symba_plA%helio%swiftest%vh(1, i)
          array_fragpl(i, 16) = symba_plA%helio%swiftest%vh(2, i)
          array_fragpl(i, 17) = symba_plA%helio%swiftest%vh(3, i)
     END DO

     DO i = npl, npl+nfrag
          array_fragpl(i, 1) = array_frag(i,1)
          array_fragpl(i, 2) = array_frag(i,2)
          array_fragpl(i, 3) = array_frag(i,3)
          array_fragpl(i, 4) = array_frag(i,4)
          array_fragpl(i, 5) = array_frag(i,5)
          array_fragpl(i, 6) = array_frag(i,6)
          array_fragpl(i, 7) = array_frag(i,7)
          array_fragpl(i, 8) = array_frag(i,8)
          array_fragpl(i, 9) = array_frag(i,9)
          array_fragpl(i, 10) = array_frag(i,10)
          array_fragpl(i, 11) = array_frag(i,11)
          array_fragpl(i, 12) = array_frag(i,12)
          array_fragpl(i, 13) = array_frag(i,13)
          array_fragpl(i, 14) = array_frag(i,14)
          array_fragpl(i, 15) = array_frag(i,15)
          array_fragpl(i, 16) = array_frag(i,16)
          array_fragpl(i, 17) = array_frag(i,17)
     END DO 

     ! set symba_plA == array_fragpl
     DO i = 1, npl+nfrag
          symba_plA%helio%swiftest%name(i) = array_fragpl(i, 1)
          symba_plA%helio%swiftest%status(i) = array_fragpl(i, 2)
          symba_plA%helio%swiftest%mass(i) = array_fragpl(i, 3)
          symba_plA%helio%swiftest%radius(i) = array_fragpl(i, 4) 
          symba_plA%helio%swiftest%rhill(i) = array_fragpl(i, 5)
          symba_plA%helio%swiftest%xh(1, i) = array_fragpl(i, 6) 
          symba_plA%helio%swiftest%xh(2, i) = array_fragpl(i, 7)
          symba_plA%helio%swiftest%xh(3, i) = array_fragpl(i, 8)
          symba_plA%helio%swiftest%xb(1, i) = array_fragpl(i, 9)
          symba_plA%helio%swiftest%xb(2, i) = array_fragpl(i, 10)
          symba_plA%helio%swiftest%xb(3, i) = array_fragpl(i, 11) 
          symba_plA%helio%swiftest%vb(1, i) = array_fragpl(i, 12)
          symba_plA%helio%swiftest%vb(2, i) = array_fragpl(i, 13)
          symba_plA%helio%swiftest%vb(3, i) = array_fragpl(i, 14)
          symba_plA%helio%swiftest%vh(1, i) = array_fragpl(i, 15)
          symba_plA%helio%swiftest%vh(2, i) = array_fragpl(i, 16)
          symba_plA%helio%swiftest%vh(3, i) = array_fragpl(i, 17)
     END DO
      
     ! reset npl
     npl = npl + nfrag

     ! print info to term.out
     WRITE(*, *) "Disruption between particles ", name1, " and ", name2, " at time t = ",t
     WRITE(*, *) "Number of fragments added: ", nfrag

     ! make sure that both parents get removed by symba_discard_merge_pl
     ! make sure that the ps+fragments get reordered by symba_rearray

     RETURN 
END SUBROUTINE symba_casedisruption