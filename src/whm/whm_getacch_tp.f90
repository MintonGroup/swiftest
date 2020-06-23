!**********************************************************************************************************************************
!
!  Unit Name   : whm_getacch_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of test particles
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!                whm_tp1P     : pointer to head of active WHM test particle structure linked-list
!                xh           : heliocentric positions of planets at time t
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_tp1P     : pointer to head of active WHM test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_getacch_tp(lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xh, j2rp2, j4rp4)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_getacch_tp(lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xh, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_getacch_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                   :: lextra_force
     INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
     TYPE(whm_pl), POINTER                      :: whm_pl1P
     TYPE(whm_tp), POINTER                      :: whm_tp1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: r2, fac, mu
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh, ir3h
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: aobl
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xht, aoblt
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_plP
     TYPE(swifter_tp), POINTER                    :: swifter_tp1P, swifter_tpP
     TYPE(whm_pl), POINTER                        :: whm_plP
     TYPE(whm_tp), POINTER                        :: whm_tpP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(aobl(NDIM, nplmax), irh(nplmax), ir3h(nplmax), xht(NDIM, ntpmax), aoblt(NDIM, ntpmax), irht(ntpmax))
          lmalloc = .FALSE.
     END IF
     swifter_pl1P => whm_pl1P%swifter
     swifter_tp1P => whm_tp1P%swifter
     DO i = 2, npl
          r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
          irh(i) = 1.0_DP/SQRT(r2)
          ir3h(i) = irh(i)/r2
     END DO
     swifter_tpP => swifter_tp1P
     DO i = 1, ntp
          xht(:, i) = swifter_tpP%xh(:)
          r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
          irht(i) = 1.0_DP/SQRT(r2)
          swifter_tpP => swifter_tpP%nextP
     END DO
     ah0(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     swifter_plP => swifter_pl1P
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          fac = swifter_plP%mass*ir3h(i)
          ah0(:) = ah0(:) - fac*xh(:, i)
     END DO
     CALL whm_getacch_ah3_tp(npl, ntp, whm_pl1P, whm_tp1P, xh)
     whm_tpP => whm_tp1P
     DO i = 1, ntp
          whm_tpP%ah(:) = whm_tpP%ah(:) + ah0(:)
          whm_tpP => whm_tpP%nextP
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          mu = whm_pl1P%swifter%mass
          CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
          whm_tpP => whm_tp1P
          DO i = 1, ntp
               whm_tpP%ah(:) = whm_tpP%ah(:) + aoblt(:, i) - aobl(:, 1)
               whm_tpP => whm_tpP%nextP
          END DO
     END IF
     IF (lextra_force) CALL whm_user_getacch_tp(t, ntp, whm_tp1P)

     RETURN

END SUBROUTINE whm_getacch_tp
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
