!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_merge_pl
!  Unit Type   : subroutine
!  Project     : Swifter
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
SUBROUTINE symba_discard_merge_pl(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)

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
     TYPE(symba_pl), POINTER                       :: symba_pl1P, symba_pld1P
     TYPE(symba_plplenc), DIMENSION(:), INTENT(IN) :: plplenc_list

! Internals
     INTEGER(I4B)              :: i, j, nchild
     REAL(DP)                  :: m, mmax, mtot, r, r3, mu, energy, ap, v2, msun
     REAL(DP), DIMENSION(NDIM) :: x, v, vbs
     TYPE(swifter_pl), POINTER :: swifter_plP, swifter_plspP
     TYPE(symba_pl), POINTER   :: symba_plP, symba_plkP, symba_plspP

! Executable code
     swifter_plP => symba_pl1P%helio%swifter
     msun = swifter_plP%mass
     vbs(:) = swifter_plP%vb(:)
     DO i = 1, nplplenc
          IF (plplenc_list(i)%status == MERGED) THEN
               IF ((plplenc_list(i)%pl1P%helio%swifter%status == ACTIVE) .AND.                                                    &
                   (plplenc_list(i)%pl2P%helio%swifter%status == ACTIVE)) THEN
                    symba_plkP => plplenc_list(i)%pl1P%parentP
                    swifter_plP => symba_plkP%helio%swifter
                    m = swifter_plP%mass
                    r = swifter_plP%radius
                    r3 = r**3
                    mmax = m
                    mtot = m
                    x(:) = m*swifter_plP%xh(:)
                    v(:) = m*swifter_plP%vb(:)
                    symba_plP => symba_plkP
                    nchild = symba_plP%nchild
                    DO j = 1, nchild
                         symba_plP => symba_plP%childP
                         swifter_plP => symba_plP%helio%swifter
                         m = swifter_plP%mass
                         r = swifter_plP%radius
                         r3 = r3 + r**3
                         mtot = mtot + m
                         x(:) = x(:) + m*swifter_plP%xh(:)
                         v(:) = v(:) + m*swifter_plP%vb(:)
                         IF (m > mmax) THEN
                              mmax = m
                              symba_plkP => symba_plP
                         END IF
                    END DO
                    x(:) = x(:)/mtot
                    v(:) = v(:)/mtot
                    r = r3**(1.0_DP/3.0_DP)
                    symba_plP => plplenc_list(i)%pl1P%parentP
                    DO j = 0, nchild
                         swifter_plP => symba_plP%helio%swifter
                         IF (ASSOCIATED(symba_plP, symba_plkP)) THEN
                              swifter_plP%mass = mtot
                              swifter_plP%radius = r
                              swifter_plP%xh(:) = x(:)
                              swifter_plP%vb(:) = v(:)
                              swifter_plP%vh(:) = v(:) - vbs(:)
                              mu = msun*mtot/(msun + mtot)
                              r = SQRT(DOT_PRODUCT(x(:), x(:)))
                              v(:) = swifter_plP%vh(:)
                              v2 = DOT_PRODUCT(v(:), v(:))
                              energy = -1.0_DP*msun*mtot/r + 0.5_DP*mu*v2
                              ap = -1.0_DP*msun*mtot/(2.0_DP*energy)
                              swifter_plP%rhill = ap*(((mu/msun)/3.0_DP)**(1.0_DP/3.0_DP))
                         ELSE
                              swifter_plP%status = MERGED
                         END IF
                         symba_plP => symba_plP%childP
                    END DO
               END IF
          END IF
     END DO
     symba_plP => symba_pl1P%nextP
     DO i = 2, npl
          symba_plspP => symba_plP
          symba_plP => symba_plP%nextP
          swifter_plspP => symba_plspP%helio%swifter
          IF (swifter_plspP%status /= ACTIVE) CALL symba_discard_spill_pl(npl, nsppl, symba_pld1P, symba_plspP)
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
