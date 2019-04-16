!**********************************************************************************************************************************
!
!  Unit Name   : whm_getacch
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Compute heliocentric accelerations of planets
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P     : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_getacch(lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_getacch(lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_getacch
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(whm_pl), POINTER    :: whm_pl1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: r2, fac
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh, irj, ir3h, ir3j
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P, swifter_plP
     TYPE(whm_pl), POINTER                        :: whm_plP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xh(NDIM, nplmax), aobl(NDIM, nplmax), irh(nplmax), irj(nplmax), ir3h(nplmax), ir3j(nplmax))
          lmalloc = .FALSE.
     END IF
     swifter_pl1P => whm_pl1P%swifter
     whm_plP => whm_pl1P
     DO i = 2, npl
          whm_plP => whm_plP%nextP
          r2 = DOT_PRODUCT(whm_plP%xj(:), whm_plP%xj(:))
          irj(i) = 1.0_DP/SQRT(r2)
          ir3j(i) = irj(i)/r2
     END DO
     swifter_plP => swifter_pl1P
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          r2 = DOT_PRODUCT(swifter_plP%xh(:), swifter_plP%xh(:))
          irh(i) = 1.0_DP/SQRT(r2)
          ir3h(i) = irh(i)/r2
     END DO
     ah0(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     swifter_plP => swifter_pl1P%nextP
     DO i = 3, npl
          swifter_plP => swifter_plP%nextP
          fac = swifter_plP%mass*ir3h(i)
          ah0(:) = ah0(:) - fac*swifter_plP%xh(:)
     END DO
     CALL whm_getacch_ah1(npl, whm_pl1P, ir3h, ir3j)
     CALL whm_getacch_ah2(npl, whm_pl1P, ir3j)
     CALL whm_getacch_ah3(npl, whm_pl1P)
     whm_plP => whm_pl1P
     whm_plP%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     DO i = 2, npl
          whm_plP => whm_plP%nextP
          whm_plP%ah(:) = ah0(:) + whm_plP%ah1(:) + whm_plP%ah2(:) + whm_plP%ah3(:)
     END DO
     IF (j2rp2 /= 0.0_DP) THEN
          swifter_plP => swifter_pl1P
          DO i = 2, npl
               swifter_plP => swifter_plP%nextP
               xh(:, i) = swifter_plP%xh(:)
          END DO
          CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
          whm_plP => whm_pl1P
          DO i = 2, npl
               whm_plP => whm_plP%nextP
               whm_plP%ah(:) = whm_plP%ah(:) + aobl(:, i) - aobl(:, 1)
          END DO
     END IF
     IF (lextra_force) CALL whm_user_getacch(t, npl, whm_pl1P)

     RETURN

END SUBROUTINE whm_getacch
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
