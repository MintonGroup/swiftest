!**********************************************************************************************************************************
!
!  Unit Name   : bs_bsstep
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Take a Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust step size
!
!  Input
!    Arguments : npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!                x            : independent variable (time)
!                htry         : step size to be attempted
!                eps          : local truncation error tolerance
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                derivs       : name of subroutine used to compute planet and test particle velocities and accelerations
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!                x            : independent variable (time)
!                hdid         : step size actually accomplished
!                hnext        : estimated next step size
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL bs_bsstep(npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, x, htry, eps, hdid, hnext, j2rp2, j4rp4, lextra_force,
!                               derivs)
!
!  Notes       : Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
!                Vetterling, and Flannery, 2nd ed., pp. 1303-5
!
!**********************************************************************************************************************************
SUBROUTINE bs_bsstep(npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, x, htry, eps, hdid, hnext, j2rp2, j4rp4, lextra_force, derivs)

! Modules
     USE module_parameters
     USE module_bs
     USE module_nrutil
     USE module_interfaces, EXCEPT_THIS_ONE => bs_bsstep
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)     :: htry, eps, j2rp2, j4rp4
     REAL(DP), INTENT(INOUT)  :: x
     REAL(DP), INTENT(OUT)    :: hdid, hnext
     TYPE(bs_pl), POINTER     :: bs_pl1P
     TYPE(bs_tp), POINTER     :: bs_tp1P
     INTERFACE
          SUBROUTINE derivs(lextra_force, t, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
               USE module_parameters
               USE module_swifter
               USE module_bs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lextra_force
               INTEGER(I4B), INTENT(IN) :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
               TYPE(bs_pl), POINTER     :: bs_pl1P
               TYPE(bs_tp), POINTER     :: bs_tp1P
          END SUBROUTINE derivs
     END INTERFACE

! Internals
     LOGICAL(LGT)                            :: reduct
     LOGICAL(LGT), SAVE                      :: first = .TRUE.
     INTEGER(I4B), PARAMETER                 :: IMAX = 9, KMAXX = IMAX - 1
     INTEGER(I4B)                            :: i, k, km, ndum
     INTEGER(I4B), SAVE                      :: kopt, kmax
     INTEGER(I4B), DIMENSION(IMAX)           :: nseq = (/ 2, 4, 6, 8, 10, 12, 14, 16, 18 /)
     REAL(DP), PARAMETER                     :: SAFE1 = 0.25_DP, SAFE2 = 0.7_DP
     REAL(DP), PARAMETER                     :: REDMAX = 1.0E-5_DP, REDMIN = 0.7_DP, SCALMX = 0.1_DP
     REAL(DP)                                :: eps1, errmax, fact, h, red, scale, wrkmin, xest
     REAL(DP), SAVE                          :: epsold = -1.0_DP, xnew
     REAL(DP), DIMENSION(KMAXX)              :: err
     REAL(DP), DIMENSION(IMAX), SAVE         :: a
     REAL(DP), DIMENSION(KMAXX, KMAXX), SAVE :: alf
     TYPE(bs_pl), POINTER                    :: bs_plP
     TYPE(bs_tp), POINTER                    :: bs_tpP

! Executable code
     IF (eps /= epsold) THEN
          hnext = -1.0E29_DP
          xnew = -1.0E29_DP
          eps1 = SAFE1*eps
          a(:) = cumsum(nseq, 1)
          WHERE (upper_triangle(KMAXX, KMAXX))                                                                                    &
               alf = eps1**(outerdiff(a(2:), a(2:))/outerprod(arth(3.0_DP, 2.0_DP, KMAXX), (a(2:) - a(1) + 1.0_DP)))
          epsold = eps
          DO kopt = 2, KMAXX - 1
               IF (a(kopt+1) > a(kopt)*alf(kopt-1, kopt)) EXIT
          END DO
          kmax = kopt
     END IF
     h = htry
     bs_plP => bs_pl1P
     DO i = 1, npl
          bs_plP%ysav(:) = bs_plP%y(:)
          bs_plP%dydxsav(:) = bs_plP%dydx(:)
          bs_plP => bs_plP%nextP
     END DO
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          bs_tpP%ysav(:) = bs_tpP%y(:)
          bs_tpP%dydxsav(:) = bs_tpP%dydx(:)
          bs_tpP => bs_tpP%nextP
     END DO
     IF (h /= hnext .OR. x /= xnew) THEN
          first = .TRUE.
          kopt = kmax
     END IF
     reduct = .FALSE.
     MAIN_LOOP: DO
          DO k = 1, kmax
               xnew = x + h
               IF (xnew == x) THEN
                    WRITE(*, *) "Step size underflow in bs_bsstep"
                    CALL util_exit(FAILURE)
               END IF
               CALL bs_mmid(npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, x, h, nseq(k), j2rp2, j4rp4, lextra_force, derivs)
               xest = (h/nseq(k))**2
               bs_plP => bs_pl1P
               DO i = 1, npl
                    bs_plP%yseq(:) = bs_plP%y(:)
                    bs_plP => bs_plP%nextP
               END DO
               bs_tpP => bs_tp1P
               DO i = 1, ntp
                    bs_tpP%yseq(:) = bs_tpP%y(:)
                    bs_tpP => bs_tpP%nextP
               END DO
               CALL bs_pzextr(k, xest, npl, ntp, bs_pl1P, bs_tp1P)
               IF (k /= 1) THEN
                    errmax = 0.0_DP
                    bs_plP => bs_pl1P
                    DO i = 1, npl
                         errmax = MAX(errmax, MAXVAL(ABS(bs_plP%yerr(:)/bs_plP%yscal(:))))
                         bs_plP => bs_plP%nextP
                    END DO
                    bs_tpP => bs_tp1P
                    DO i = 1, ntp
                         errmax = MAX(errmax, MAXVAL(ABS(bs_tpP%yerr(:)/bs_tpP%yscal(:))))
                         bs_tpP => bs_tpP%nextP
                    END DO
                    errmax = MAX(TINYBS, errmax)/eps
                    km = k - 1
                    err(km) = (errmax/SAFE1)**(1.0_DP/(2*km + 1))
               END IF
               IF (k /= 1 .AND. (k >= kopt - 1 .OR. first)) THEN
                    IF (errmax < 1.0_DP) EXIT MAIN_LOOP
                    IF (k == kmax .OR. k == kopt + 1) THEN
                         red = SAFE2/err(km)
                         EXIT
                    ELSE IF (k == kopt) THEN
                         IF (alf(kopt-1, kopt) < err(km)) THEN
                              red = 1.0_DP/err(km)
                              EXIT
                         END IF
                    ELSE IF (kopt == kmax) THEN
                         IF (alf(km, kmax-1) < err(km)) THEN
                              red = alf(km, kmax-1)*SAFE2/err(km)
                              EXIT
                         END IF
                    ELSE IF (alf(km, kopt) < err(km)) THEN
                         red = alf(km, kopt-1)/err(km)
                         EXIT
                    END IF
               END IF
          END DO
          red = MAX(MIN(red, REDMIN), REDMAX)
          h = h*red
          reduct = .TRUE.
     END DO MAIN_LOOP
     x = xnew
     hdid = h
     first = .FALSE.
     kopt = 1 + iminloc(a(2:km+1)*MAX(err(1:km), SCALMX))
     scale = MAX(err(kopt-1), SCALMX)
     wrkmin = scale*a(kopt)
     hnext = h/scale
     IF (kopt >= k .AND. kopt /= kmax .AND. .NOT. reduct) THEN
          fact = MAX(scale/alf(kopt-1, kopt), SCALMX)
          IF (a(kopt+1)*fact <= wrkmin) THEN
               hnext = h/fact
               kopt = kopt + 1
          END IF
     END IF

     RETURN

END SUBROUTINE bs_bsstep
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
