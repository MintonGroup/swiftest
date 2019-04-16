!**********************************************************************************************************************************
!
!  Unit Name   : bs_getaccb_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Compute barycentric accelerations of test particles
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : bs_tp1P      : pointer to head of active BS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_getaccb_tp(lextra_force, t, npl, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Martin Duncan's Swift routine tu4_getaccb_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE bs_getaccb_tp(lextra_force, t, npl, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_getaccb_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, ntp, ntpmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(bs_pl), POINTER     :: bs_pl1P
     TYPE(bs_tp), POINTER     :: bs_tp1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: r2, fac, mu
     REAL(DP), DIMENSION(NDIM)                    :: dx, xsun
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xht, aoblt
     TYPE(swifter_pl), POINTER                    :: swifter_plP
     TYPE(swifter_tp), POINTER                    :: swifter_tpP
     TYPE(bs_pl), POINTER                         :: bs_plP
     TYPE(bs_tp), POINTER                         :: bs_tpP

! Executable code
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          bs_tpP%ab(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          bs_plP => bs_pl1P
          DO j = 1, npl
               dx(:) = bs_tpP%y(1:NDIM) - bs_plP%y(1:NDIM)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               fac = bs_plP%swifter%mass/(r2*SQRT(r2))
               bs_tpP%ab(:) = bs_tpP%ab(:) - fac*dx(:)
               bs_plP => bs_plP%nextP
          END DO
          bs_tpP => bs_tpP%nextP
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          IF (lmalloc) THEN
               ALLOCATE(xht(NDIM, ntpmax), aoblt(NDIM, ntpmax), irht(ntpmax))
               lmalloc = .FALSE.
          END IF
          xsun(:) = bs_pl1P%y(1:NDIM)
          mu = bs_pl1P%swifter%mass
          bs_tpP => bs_tp1P
          DO i = 1, ntp
               xht(:, i) = bs_tpP%y(1:NDIM) - xsun(:)
               r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
               irht(i) = 1.0_DP/SQRT(r2)
               bs_tpP => bs_tpP%nextP
          END DO
          CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
          bs_tpP => bs_tp1P
          DO i = 1, ntp
               bs_tpP%ab(:) = bs_tpP%ab(:) + aoblt(:, i)
               bs_tpP => bs_tpP%nextP
          END DO
     END IF
     IF (lextra_force) CALL bs_user_getaccb_tp(t, ntp, bs_tp1P)

     RETURN

END SUBROUTINE bs_getaccb_tp
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
