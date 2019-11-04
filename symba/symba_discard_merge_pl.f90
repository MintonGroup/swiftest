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
SUBROUTINE symba_discard_merge_pl(t, npl, nsppl, symba_plA, symba_pldA, nplplenc, plplenc_list)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_merge_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                      :: nplplenc
     INTEGER(I4B), INTENT(INOUT)                   :: npl, nsppl
     REAL(DP), INTENT(IN)                          :: t
     TYPE(symba_pl)                                :: symba_plA, symba_pldA
     TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list

! Internals
     INTEGER(I4B)              :: i, j, nchild
     REAL(DP)                  :: m, mmax, mtot, r, r3, mu, energy, ap, v2, msun
     REAL(DP), DIMENSION(NDIM) :: x, v, vbs
     TYPE(swiftest_pl)         :: swiftest_plA, swiftest_plsplA
     TYPE(symba_pl)            :: symba_plA, symba_plsplA

! Executable code
     !swifter_plP => symba_pl1P%helio%swifter
     msun = symba_plA%helio%swiftest%mass(1)
     vbs(:) = symba_plA%helio%swiftest%vb(:,1)
     DO i = 1, nplplenc
          IF (plplenc_list%status(i) == MERGED) THEN
               index1 = plplenc_list%id1(i)
               index2 = plplenc_list%id2(i)
               IF ((symba_plA%helio%swiftest%status(index1) == ACTIVE) .AND.                                                    &
                   (symba_plA%helio%swiftest%status(index2) == ACTIVE)) THEN
                    !symba_plkP => plplenc_list(i)%pl1P%parentP
                    !swifter_plP => symba_plkP%helio%swifter
                    m = symba_plA%helio%swiftest%mass(i)
                    r = symba_plA%helio%swiftest%radius(i)
                    r3 = r**3
                    mmax = m
                    mtot = m
                    x(:) = m*symba_plA%helio%swiftest%xh(:,i)
                    v(:) = m*symba_plA%helio%swiftest%vb(:,i)
                    !symba_plP => symba_plkP
                    indexk = index1
                    nchild = symba_plA%nchild(i)
                    DO j = 1, nchild
                         !symba_plP => symba_plP%childP
                         !swifter_plP => symba_plP%helio%swifter
                         indexchild = ????

                         m = symba_plA%helio%swiftest%mass(indexchild)
                         r = symba_plA%helio%swiftest%radius(indexchild)
                         r3 = r3 + r**3
                         mtot = mtot + m
                         x(:) = x(:) + m*symba_plA%helio%swiftest%xh(:,indexchild)
                         v(:) = v(:) + m*symba_plA%helio%swiftest%vb(:,indexchild)
                         IF (m > mmax) THEN
                              mmax = m
                              !symba_plkP => symba_plP
                              indexk = indexchild
                         END IF
                    END DO
                    x(:) = x(:)/mtot
                    v(:) = v(:)/mtot
                    r = r3**(1.0_DP/3.0_DP)
                    !symba_plP => plplenc_list(i)%pl1P%parentP
                    !DO j = 0, nchild
                         !indexchild = index1
                         !swifter_plP => symba_plP%helio%swifter
                         !IF (indexchiindexk) THEN
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
                         !ELSE
                    DO j = 0, nchild
                         indexchild = index1
                         IF (indexchild /= indexk) THEN
                              symba_plA%helio%swiftest%status(indexchild) = MERGED
                         END IF
                         indexchild = ????
                         !symba_plP => symba_plP%childP
                    END DO
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
