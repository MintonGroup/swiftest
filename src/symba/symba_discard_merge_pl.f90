!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_merge_pl
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
!  Invocation  : CALL symba_discard_merge_pl(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_merge_pl(t, npl, symba_plA, nplplenc, plplenc_list)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_merge_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                      :: nplplenc
     INTEGER(I4B), INTENT(INOUT)                   :: npl
     REAL(DP), INTENT(IN)                          :: t
     TYPE(symba_pl)                                :: symba_plA
     TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list

! Internals
     INTEGER(I4B)                  :: i, j, nchild, indexchild, enc_big, index1, index2, indexk 
     REAL(DP)                      :: m, mmax, mtot, r, r3, mu, energy, ap, v2, msun
     REAL(DP), DIMENSION(NDIM)     :: x, v, vbs
     INTEGER(I4B), DIMENSION(npl)  :: array_child

! Executable code
     msun = symba_plA%helio%swiftest%mass(1)
     vbs(:) = symba_plA%helio%swiftest%vb(:,1)
     DO i = 1, nplplenc
          IF (plplenc_list%status(i) == MERGED) THEN
               index1 = plplenc_list%index1(i)
               index2 = plplenc_list%index2(i)
               ! This IF statement is for if lfragmentation = FALSE
               IF ((symba_plA%helio%swiftest%status(index1) == ACTIVE) .AND.                                                    &
                   (symba_plA%helio%swiftest%status(index2) == ACTIVE)) THEN

                    enc_big = plplenc_list%index1(i)

                    m = symba_plA%helio%swiftest%mass(enc_big)
                    r = symba_plA%helio%swiftest%radius(enc_big)
                    r3 = r**3
                    mmax = m
                    mtot = m
                    x(:) = m*symba_plA%helio%swiftest%xh(:,enc_big)
                    v(:) = m*symba_plA%helio%swiftest%vb(:,enc_big)
                    indexk = enc_big

                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)

                    DO j = 1, nchild
                         indexchild = array_child(j)
                         m = symba_plA%helio%swiftest%mass(indexchild)
                         r = symba_plA%helio%swiftest%radius(indexchild)
                         r3 = r3 + r**3
                         mtot = mtot + m
                         x(:) = x(:) + m*symba_plA%helio%swiftest%xh(:,indexchild)
                         v(:) = v(:) + m*symba_plA%helio%swiftest%vb(:,indexchild)
                         IF (m > mmax) THEN
                              mmax = m
                              indexk = indexchild
                         END IF
                    END DO
                    x(:) = x(:)/mtot
                    v(:) = v(:)/mtot
                    r = r3**(1.0_DP/3.0_DP)
                    symba_plA%helio%swiftest%mass(indexk) = mtot
                    symba_plA%helio%swiftest%radius(indexk) = r
                    symba_plA%helio%swiftest%xh(:,indexk) = x(:)
                    symba_plA%helio%swiftest%vb(:,indexk) = v(:)
                    symba_plA%helio%swiftest%vh(:,indexk) = v(:) - vbs(:)
                    mu = msun*mtot/(msun + mtot)
                    r = SQRT(DOT_PRODUCT(x(:), x(:)))
                    v(:) = symba_plA%helio%swiftest%vh(:,indexk)
                    v2 = DOT_PRODUCT(v(:), v(:))
                    energy = -1.0_DP*msun*mtot/r + 0.5_DP*mu*v2
                    ap = -1.0_DP*msun*mtot/(2.0_DP*energy)
                    symba_plA%helio%swiftest%rhill(indexk) = ap*(((mu/msun)/3.0_DP)**(1.0_DP/3.0_DP))
                    array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
                    indexchild = enc_big
                    ldiscard = .TRUE.
                    DO j = 0, nchild
                         IF (indexchild /= indexk) THEN
                              symba_plA%helio%swiftest%status(indexchild) = MERGED
                         END IF
                         indexchild = array_child(j+1)
                    END DO

               ELSE IF ((symba_plA%helio%swiftest%status(index1) == DISRUPTION) .AND.    &                                                
                   (symba_plA%helio%swiftest%status(index2) == DISRUPTION)) THEN 

                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               ELSE IF ((symba_plA%helio%swiftest%status(index1) == SUPERCATASTROPHIC) .AND.   &                                                 
                   (symba_plA%helio%swiftest%status(index2) == SUPERCATASTROPHIC)) THEN 
                    
                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               ELSE IF ((symba_plA%helio%swiftest%status(index1) == HIT_AND_RUN) .AND.      &                                              
                   (symba_plA%helio%swiftest%status(index2) == HIT_AND_RUN)) THEN 

                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               ELSE IF ((symba_plA%helio%swiftest%status(index1) == GRAZE_AND_MERGE) .AND.  &                                                  
                   (symba_plA%helio%swiftest%status(index2) == GRAZE_AND_MERGE)) THEN 

                    enc_big = plplenc_list%index1(i)
                    nchild = symba_plA%nchild(enc_big)
                    array_child(1:npl) = symba_plA%index_child(1:npl,enc_big)
                    DO j = 1, nchild
                         symba_plA%helio%swiftest%status(array_child(j)) = INACTIVE
                    END DO
                    ldiscard = .TRUE.
               END IF
          END IF
     END DO

     RETURN
     
END SUBROUTINE symba_discard_merge_pl
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
