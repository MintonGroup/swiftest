!**********************************************************************************************************************************
!
!  Unit Name   : symba_energy
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute total system angular momentum vector and kinetic, potential and total system energy
!
!  Input
!    Arguments : npl          : number of massive bodies
!                nplmax       : maximum allowed number of massive bodies
!                swiftest_pl1P : pointer to head of Swiftest massive body structure linked-list
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
!  Invocation  : CALL symba_energy(npl, nplmax, swiftest_pl1P, j2rp2, j4rp4, ke, pe, te, htot)
!
!  Notes       : Adapted from Martin Duncan's Swift routine anal_energy.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_energy(npl, nplmax, swiftest_plA, j2rp2, j4rp4, ke, pe, te, htot)

! Modules
     use swiftest, EXCEPT_THIS_ONE => symba_energy
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)               :: npl, nplmax
     REAL(DP), INTENT(IN)                   :: j2rp2, j4rp4
     REAL(DP), INTENT(OUT)                  :: ke, pe, te
     REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: htot
     TYPE(swiftest_pl), INTENT(INOUT)       :: swiftest_plA

! Internals
     LOGICAL(LGT)                                 :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j
     REAL(DP)                                     :: mass, msys, r2, v2, oblpot
     REAL(DP), DIMENSION(NDIM)                    :: h, x, v, dx
     REAL(DP), DIMENSION(npl)                     :: irh


! Executable code

     CALL coord_h2b(npl, swiftest_plA, msys)
     htot = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     ke = 0.0_DP
     pe = 0.0_DP

!$omp parallel do default(none) &
!$omp shared (swiftest_plA, npl) &
!$omp private (i, x, v, mass, h, htot, v2, dx, r2) &
!$omp reduction (+:ke, pe)
     DO i = 1, npl - 1
          x(:) = swiftest_plA%xb(:,i)
          v(:) = swiftest_plA%vb(:,i)
          mass = swiftest_plA%mass(i)
          h(1) = mass*(x(2)*v(3) - x(3)*v(2))
          h(2) = mass*(x(3)*v(1) - x(1)*v(3))
          h(3) = mass*(x(1)*v(2) - x(2)*v(1))
          htot(:) = htot(:) + h(:)
          v2 = DOT_PRODUCT(v(:), v(:))
          ke = ke + 0.5_DP*mass*v2
          DO j = i + 1, npl
               dx(:) = swiftest_plA%xb(:,j) - x(:) !this is 0 for the removed ps because swiftest_plA%xb(:,j) = swiftest_plA%xb(:,i)
               r2 = DOT_PRODUCT(dx(:), dx(:)) !this is 0 for the removed ps
               IF (r2 /= 0) THEN
                    pe = pe - mass*swiftest_plA%mass(j)/SQRT(r2) !DIVISION !!!!!!r2 is 0 for ps 12 which is the removed ps
               END IF
          END DO
     END DO
!$omp end parallel do
     i = npl ! needed to account for the parllelization above
     x(:) = swiftest_plA%xb(:,i)
     v(:) = swiftest_plA%vb(:,i)
     mass = swiftest_plA%mass(i)
     h(1) = mass*(x(2)*v(3) - x(3)*v(2))
     h(2) = mass*(x(3)*v(1) - x(1)*v(3))
     h(3) = mass*(x(1)*v(2) - x(2)*v(1))
     htot(:) = htot(:) + h(:)
     v2 = DOT_PRODUCT(v(:), v(:))
     ke = ke + 0.5_DP*mass*v2
     IF (j2rp2 /= 0.0_DP) THEN
          DO i = 2, npl
               r2 = DOT_PRODUCT(swiftest_plA%xh(:,i),swiftest_plA%xh(:,i))
               irh(i) = 1.0_DP/SQRT(r2)
          END DO
          CALL obl_pot(swiftest_plA, j2rp2, j4rp4, swiftest_plA%xh(:,:), irh, oblpot)
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
