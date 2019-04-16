!**********************************************************************************************************************************
!
!  Unit Name   : tu4_getaccb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : tu4
!  Language    : Fortran 90/95
!
!  Description : Compute barycentric accelerations of planets
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                tu4_pl1P     : pointer to head of TU4 planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : tu4_pl1P     : pointer to head of TU4 planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL tu4_getaccb(lextra_force, t, npl, nplmax, tu4_pl1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Martin Duncan's Swift routine tu4_getaccb.f
!
!**********************************************************************************************************************************
SUBROUTINE tu4_getaccb(lextra_force, t, npl, nplmax, tu4_pl1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_tu4
     USE module_interfaces, EXCEPT_THIS_ONE => tu4_getaccb
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(tu4_pl), POINTER    :: tu4_pl1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: mu, massi, massj, r2, ir3
     REAL(DP), DIMENSION(NDIM)                    :: dx, acc
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_pliP, swifter_pljP
     TYPE(tu4_pl), POINTER                        :: tu4_pliP, tu4_pljP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax))
          lmalloc = .FALSE.
     END IF
     swifter_pl1P => tu4_pl1P%swifter
     mu = swifter_pl1P%mass
     tu4_pl1P%ab(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     tu4_pliP => tu4_pl1P
     DO i = 2, npl
          tu4_pliP => tu4_pliP%nextP
          swifter_pliP => tu4_pliP%swifter
          massi = swifter_pliP%mass
          dx(:) = swifter_pl1P%xb(:) - swifter_pliP%xb(:)
          xh(:, i) = -dx(:)
          r2 = DOT_PRODUCT(dx(:), dx(:))
          irh(i) = 1.0_DP/SQRT(r2)
          ir3 = irh(i)/r2
          acc(:) = ir3*dx(:)
          tu4_pl1P%ab(:) = tu4_pl1P%ab(:) - massi*acc(:)
          tu4_pliP%ab(:) = mu*acc(:)
     END DO
     tu4_pliP => tu4_pl1P
     DO i = 2, npl - 1
          tu4_pliP => tu4_pliP%nextP
          swifter_pliP => tu4_pliP%swifter
          massi = swifter_pliP%mass
          tu4_pljP => tu4_pliP
          DO j = i + 1, npl
               tu4_pljP => tu4_pljP%nextP
               swifter_pljP => tu4_pljP%swifter
               massj = swifter_pljP%mass
               dx(:) = swifter_pliP%xb(:) - swifter_pljP%xb(:)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               ir3 = 1.0_DP/(r2*SQRT(r2))
               acc(:) = ir3*dx(:)
               tu4_pliP%ab(:) = tu4_pliP%ab(:) - massj*acc(:)
               tu4_pljP%ab(:) = tu4_pljP%ab(:) + massi*acc(:)
          END DO
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          tu4_pliP => tu4_pl1P
          DO i = 1, npl
               tu4_pliP%ab(:) = tu4_pliP%ab(:) + aobl(:, i)
               tu4_pliP => tu4_pliP%nextP
          END DO
     END IF
     IF (lextra_force) CALL tu4_user_getaccb(t, npl, tu4_pl1P)

     RETURN

END SUBROUTINE tu4_getaccb
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
