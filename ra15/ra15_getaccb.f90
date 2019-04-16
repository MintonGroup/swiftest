!**********************************************************************************************************************************
!
!  Unit Name   : ra15_getaccb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
!  Language    : Fortran 90/95
!
!  Description : Compute barycentric accelerations of planets
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ra15_pl1P    : pointer to head of RA15 planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ra15_pl1P    : pointer to head of RA15 planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL ra15_getaccb(lextra_force, t, npl, nplmax, ra15_pl1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Martin Duncan's Swift routine tu4_getaccb.f
!
!**********************************************************************************************************************************
SUBROUTINE ra15_getaccb(lextra_force, t, npl, nplmax, ra15_pl1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => ra15_getaccb
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(ra15_pl), POINTER   :: ra15_pl1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: mu, massi, massj, r2, ir3
     REAL(DP), DIMENSION(NDIM)                    :: dx, acc
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_pliP, swifter_pljP
     TYPE(ra15_pl), POINTER                       :: ra15_pliP, ra15_pljP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax))
          lmalloc = .FALSE.
     END IF
     swifter_pl1P => ra15_pl1P%swifter
     mu = swifter_pl1P%mass
     ra15_pl1P%ab(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     ra15_pliP => ra15_pl1P
     DO i = 2, npl
          ra15_pliP => ra15_pliP%nextP
          swifter_pliP => ra15_pliP%swifter
          massi = swifter_pliP%mass
          dx(:) = swifter_pl1P%xb(:) - swifter_pliP%xb(:)
          xh(:, i) = -dx(:)
          r2 = DOT_PRODUCT(dx(:), dx(:))
          irh(i) = 1.0_DP/SQRT(r2)
          ir3 = irh(i)/r2
          acc(:) = ir3*dx(:)
          ra15_pl1P%ab(:) = ra15_pl1P%ab(:) - massi*acc(:)
          ra15_pliP%ab(:) = mu*acc(:)
     END DO
     ra15_pliP => ra15_pl1P
     DO i = 2, npl - 1
          ra15_pliP => ra15_pliP%nextP
          swifter_pliP => ra15_pliP%swifter
          massi = swifter_pliP%mass
          ra15_pljP => ra15_pliP
          DO j = i + 1, npl
               ra15_pljP => ra15_pljP%nextP
               swifter_pljP => ra15_pljP%swifter
               massj = swifter_pljP%mass
               dx(:) = swifter_pliP%xb(:) - swifter_pljP%xb(:)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               ir3 = 1.0_DP/(r2*SQRT(r2))
               acc(:) = ir3*dx(:)
               ra15_pliP%ab(:) = ra15_pliP%ab(:) - massj*acc(:)
               ra15_pljP%ab(:) = ra15_pljP%ab(:) + massi*acc(:)
          END DO
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          ra15_pliP => ra15_pl1P
          DO i = 1, npl
               ra15_pliP%ab(:) = ra15_pliP%ab(:) + aobl(:, i)
               ra15_pliP => ra15_pliP%nextP
          END DO
     END IF
     IF (lextra_force) CALL ra15_user_getaccb(t, npl, ra15_pl1P)

     RETURN

END SUBROUTINE ra15_getaccb
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
