!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_getaccp_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Compute planetocentric accelerations of test particles closely encountering a planet
!
!  Input
!    Arguments : index        : inner substep number within current set
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                nenc         : number of test particles encountering current planet
!                ntpmax       : maximum allowed number of test particles
!                rmvs_pl1P    : pointer to head of RMVS planet structure linked-list
!                rmvs_pleP    : pointer to RMVS planet structure of planet being closely encountered
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_pleP    : pointer to RMVS planet structure of planet being closely encountered
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_getaccp_tp(index, lextra_force, t, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2, j4rp4)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_getaccp_tp(index, lextra_force, t, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_getaccp_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: index, npl, nplmax, nenc, ntpmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(rmvs_pl), POINTER   :: rmvs_pl1P, rmvs_pleP

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: r2, fac, mu
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: ir3p, irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xht, aoblt
     TYPE(rmvs_pl), POINTER                       :: rmvs_plP
     TYPE(rmvs_tp), POINTER                       :: rmvs_tpenc1P, rmvs_tpP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(ir3p(nplmax), xht(NDIM, ntpmax), aoblt(NDIM, ntpmax), irht(ntpmax))
          lmalloc = .FALSE.
     END IF
     rmvs_tpenc1P => rmvs_pleP%tpenc1P
     ah0(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     rmvs_plP => rmvs_pl1P
     DO i = 1, npl
          IF (.NOT. ASSOCIATED(rmvs_plP, rmvs_pleP)) THEN
               r2 = DOT_PRODUCT(rmvs_plP%xpc(:, index), rmvs_plP%xpc(:, index))
               ir3p(i) = 1.0_DP/(r2*SQRT(r2))
               fac = rmvs_plP%whm%swifter%mass*ir3p(i)
               ah0(:) = ah0(:) - fac*rmvs_plP%xpc(:, index)
          END IF
          rmvs_plP => rmvs_plP%nextP
     END DO
     CALL rmvs_getaccp_ah3_tp(index, npl, nenc, rmvs_pl1P, rmvs_pleP)
     rmvs_tpP => rmvs_tpenc1P
     DO i = 1, nenc
          rmvs_tpP%apc(:) = rmvs_tpP%apc(:) + ah0(:)
          rmvs_tpP => rmvs_tpP%tpencP
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          mu = rmvs_pl1P%whm%swifter%mass
          rmvs_tpP => rmvs_tpenc1P
          DO i = 1, nenc
               xht(:, i) = rmvs_tpP%whm%swifter%xh(:)
               r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
               irht(i) = 1.0_DP/SQRT(r2)
               rmvs_tpP => rmvs_tpP%tpencP
          END DO
          CALL obl_acc_tp(nenc, xht, j2rp2, j4rp4, irht, aoblt, mu)
          rmvs_tpP => rmvs_tpenc1P
          DO i = 1, nenc
               rmvs_tpP%apc(:) = rmvs_tpP%apc(:) + aoblt(:, i) - rmvs_pleP%aobl(:, index)
               rmvs_tpP => rmvs_tpP%tpencP
          END DO
     END IF
     IF (lextra_force) CALL rmvs_user_getacch_tp(t, nenc, rmvs_tpenc1P)

     RETURN

END SUBROUTINE rmvs_getaccp_tp
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
