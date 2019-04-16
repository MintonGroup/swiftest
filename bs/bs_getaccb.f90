!**********************************************************************************************************************************
!
!  Unit Name   : bs_getaccb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Compute barycentric accelerations of planets
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : bs_pl1P      : pointer to head of BS planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_getaccb(lextra_force, t, npl, nplmax, bs_pl1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Martin Duncan's Swift routine tu4_getaccb.f
!
!**********************************************************************************************************************************
SUBROUTINE bs_getaccb(lextra_force, t, npl, nplmax, bs_pl1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_getaccb
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(bs_pl), POINTER     :: bs_pl1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: mu, r2, ir3
     REAL(DP), DIMENSION(NDIM)                    :: dx, acc
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_pliP, swifter_pljP
     TYPE(bs_pl), POINTER                         :: bs_pliP, bs_pljP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax))
          lmalloc = .FALSE.
     END IF
     bs_pl1P%ab(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     mu = bs_pl1P%swifter%mass
     bs_pliP => bs_pl1P
     DO i = 2, npl
          bs_pliP => bs_pliP%nextP
          swifter_pliP => bs_pliP%swifter
          dx(:) = bs_pl1P%y(1:NDIM) - bs_pliP%y(1:NDIM)
          xh(:, i) = -dx(:)
          r2 = DOT_PRODUCT(dx(:), dx(:))
          irh(i) = 1.0_DP/SQRT(r2)
          ir3 = irh(i)/r2
          acc(:) = ir3*dx(:)
          bs_pl1P%ab(:) = bs_pl1P%ab(:) - swifter_pliP%mass*acc(:)
          bs_pliP%ab(:) = mu*acc(:)
     END DO
     bs_pliP => bs_pl1P
     DO i = 2, npl - 1
          bs_pliP => bs_pliP%nextP
          swifter_pliP => bs_pliP%swifter
          bs_pljP => bs_pliP
          DO j = i + 1, npl
               bs_pljP => bs_pljP%nextP
               swifter_pljP => bs_pljP%swifter
               dx(:) = bs_pliP%y(1:NDIM) - bs_pljP%y(1:NDIM)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               ir3 = 1.0_DP/(r2*SQRT(r2))
               acc(:) = ir3*dx(:)
               bs_pliP%ab(:) = bs_pliP%ab(:) - swifter_pljP%mass*acc(:)
               bs_pljP%ab(:) = bs_pljP%ab(:) + swifter_pliP%mass*acc(:)
          END DO
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          swifter_pl1P => bs_pl1P%swifter
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          bs_pliP => bs_pl1P
          DO i = 1, npl
               bs_pliP%ab(:) = bs_pliP%ab(:) + aobl(:, i)
               bs_pliP => bs_pliP%nextP
          END DO
     END IF
     IF (lextra_force) CALL bs_user_getaccb(t, npl, bs_pl1P)

     RETURN

END SUBROUTINE bs_getaccb
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
