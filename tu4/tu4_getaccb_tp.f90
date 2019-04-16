!**********************************************************************************************************************************
!
!  Unit Name   : tu4_getaccb_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
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
!                tu4_pl1P     : pointer to head of TU4 planet structure linked-list
!                tu4_tp1P     : pointer to head of active TU4 test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : tu4_tp1P     : pointer to head of active TU4 test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL tu4_getaccb_tp(lextra_force, t, npl, ntp, ntpmax, tu4_pl1P, tu4_tp1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Martin Duncan's Swift routine tu4_getaccb_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE tu4_getaccb_tp(lextra_force, t, npl, ntp, ntpmax, tu4_pl1P, tu4_tp1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_getaccb_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, ntp, ntpmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(tu4_pl), POINTER    :: tu4_pl1P
     TYPE(tu4_tp), POINTER    :: tu4_tp1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: r2, fac, mu
     REAL(DP), DIMENSION(NDIM)                    :: dx, x1tmp, x2tmp, acc
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irht
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xht, aoblt
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_plP
     TYPE(swifter_tp), POINTER                    :: swifter_tpP
     TYPE(tu4_tp), POINTER                        :: tu4_tpP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xht(NDIM, ntpmax), aoblt(NDIM, ntpmax), irht(ntpmax))
          lmalloc = .FALSE.
     END IF
     swifter_pl1P => tu4_pl1P%swifter
     x1tmp(:) = swifter_pl1P%xb(:)
     tu4_tpP => tu4_tp1P
     DO i = 1, ntp
          swifter_tpP => tu4_tpP%swifter
          acc(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          xht(:, i) = swifter_tpP%xb(:) - x1tmp(:)
          r2 = DOT_PRODUCT(xht(:, i), xht(:, i))
          irht(i) = 1.0_DP/SQRT(r2)
          x2tmp(:) = swifter_tpP%xb(:)
          swifter_plP => swifter_pl1P
          DO j = 1, npl
               dx(:) = x2tmp(:) - swifter_plP%xb(:)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               fac = swifter_plP%mass/(r2*SQRT(r2))
               acc(:) = acc(:) - fac*dx(:)
               swifter_plP => swifter_plP%nextP
          END DO
          tu4_tpP%ab(:) = acc(:)
          tu4_tpP => tu4_tpP%nextP
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          mu = swifter_pl1P%mass
          CALL obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, mu)
          tu4_tpP => tu4_tp1P
          DO i = 1, ntp
               tu4_tpP%ab(:) = tu4_tpP%ab(:) + aoblt(:, i)
               tu4_tpP => tu4_tpP%nextP
          END DO
     END IF
     IF (lextra_force) CALL tu4_user_getaccb_tp(t, ntp, tu4_tp1P)

     RETURN

END SUBROUTINE tu4_getaccb_tp
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
