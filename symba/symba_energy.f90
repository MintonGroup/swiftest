!**********************************************************************************************************************************
!
!  Unit Name   : symba_energy
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute total system angular momentum vector and kinetic, potential and total system energy
!
!  Input
!    Arguments : npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ke           : kinetic energy
!                pe           : potential energy
!                te           : total energy
!                htot         : angular momentum vector
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_energy(npl, nplmax, swifter_pl1P, j2rp2, j4rp4, ke, pe, te, htot)
!
!  Notes       : Adapted from Martin Duncan's Swift routine anal_energy.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_energy(npl, nplmax, swifter_pl1P, j2rp2, j4rp4, ke, pe, te, htot)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => symba_energy
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)               :: npl, nplmax
     REAL(DP), INTENT(IN)                   :: j2rp2, j4rp4
     REAL(DP), INTENT(OUT)                  :: ke, pe, te
     REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: htot
     TYPE(swifter_pl), POINTER              :: swifter_pl1P

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: mass, msys, r2, v2, oblpot
     REAL(DP), DIMENSION(NDIM)                    :: h, x, v, dx
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh
     TYPE(swifter_pl), POINTER                    :: swifter_pliP, swifter_pljP

! Executable code
     CALL coord_h2b(npl, swifter_pl1P, msys)
     htot = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     ke = 0.0_DP
     pe = 0.0_DP
     swifter_pliP => swifter_pl1P
     DO i = 1, npl - 1
          x(:) = swifter_pliP%xb(:)
          v(:) = swifter_pliP%vb(:)
          mass = swifter_pliP%mass
          h(1) = mass*(x(2)*v(3) - x(3)*v(2))
          h(2) = mass*(x(3)*v(1) - x(1)*v(3))
          h(3) = mass*(x(1)*v(2) - x(2)*v(1))
          htot(:) = htot(:) + h(:)
          v2 = DOT_PRODUCT(v(:), v(:))
          ke = ke + 0.5_DP*mass*v2
          swifter_pljP => swifter_pliP
          DO j = i + 1, npl
               swifter_pljP => swifter_pljP%nextP
               dx(:) = swifter_pljP%xb(:) - x(:)
               r2 = DOT_PRODUCT(dx(:), dx(:))
               pe = pe - mass*swifter_pljP%mass/SQRT(r2)
          END DO
          swifter_pliP => swifter_pliP%nextP
     END DO
     x(:) = swifter_pliP%xb(:)
     v(:) = swifter_pliP%vb(:)
     mass = swifter_pliP%mass
     h(1) = mass*(x(2)*v(3) - x(3)*v(2))
     h(2) = mass*(x(3)*v(1) - x(1)*v(3))
     h(3) = mass*(x(1)*v(2) - x(2)*v(1))
     htot(:) = htot(:) + h(:)
     v2 = DOT_PRODUCT(v(:), v(:))
     ke = ke + 0.5_DP*mass*v2
     IF (j2rp2 /= 0.0_DP) THEN
          IF (lmalloc) THEN
               ALLOCATE(xh(NDIM, nplmax), irh(nplmax))
               lmalloc = .FALSE.
          END IF
          swifter_pliP => swifter_pl1P
          DO i = 2, npl
               swifter_pliP => swifter_pliP%nextP
               xh(:, i) = swifter_pliP%xh(:)
               r2 = DOT_PRODUCT(xh(:, i), xh(:, i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_pot(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, oblpot)
          pe = pe + oblpot
     END IF
     te = ke + pe

     RETURN

END SUBROUTINE symba_energy
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
