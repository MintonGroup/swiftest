!**********************************************************************************************************************************
!
!  Unit Name   : ra15_sequence
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
!  Language    : Fortran 90/95
!
!  Description : Propagate the solution of the equations of motion by one time sequence using the 15th order RADAU method
!
!  Input
!    Arguments : niter        : number of iterations to perform to achieve desired accuracy
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                ra15_pl1P    : pointer to head of RA15 planet structure linked-list
!                ra15_tp1P    : pointer to head of active RA15 test particle structure linked-list
!                x            : independent variable (time)
!                htry         : time step to try
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ra15_pl1P    : pointer to head of RA15 planet structure linked-list
!                ra15_tp1P    : pointer to head of active RA15 test particle structure linked-list
!                x            : independent variable (time)
!                hdid         : time step actually accomplished
!                hnext        : recommended next time step
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL ra15_sequence(niter, npl, nplmax, ntp, ntpmax, ra15_pl1P, ra15_tp1P, x, htry, hdid, hnext, j2rp2, j4rp4,
!                                   lextra_force)
!
!  Notes       : Adapted from Edgar Everhart's RADAU15 routine ra15 (reference given in main program notes)
!
!                The length of the time sequence is self-adjusting to maintain the magnitude of the last term in the truncated
!                time series at or below the specified tolerance eps
!
!**********************************************************************************************************************************
SUBROUTINE ra15_sequence(niter, npl, nplmax, ntp, ntpmax, ra15_pl1P, ra15_tp1P, x, htry, hdid, hnext, j2rp2, j4rp4, lextra_force)

! Modules
     USE module_parameters
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => ra15_sequence
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: niter, npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)     :: htry, j2rp2, j4rp4
     REAL(DP), INTENT(INOUT)  :: x
     REAL(DP), INTENT(OUT)    :: hdid, hnext
     TYPE(ra15_pl), POINTER   :: ra15_pl1P
     TYPE(ra15_tp), POINTER   :: ra15_tp1P

! Internals
     INTEGER(I4B)           :: i, j, k, m, jd
     REAL(DP)               :: tpp, tm, t, t2, tval, s, y, z, a, temp, gk, hv, q
     REAL(DP), DIMENSION(7) :: b, g, e, bd
     TYPE(ra15_pl), POINTER :: ra15_plP
     TYPE(ra15_tp), POINTER :: ra15_tpP

! Executable code
     tpp = htry
     tm = x
     DO
          ra15_plP => ra15_pl1P
          DO i = 1, npl
               DO j = 1, NDIM
                    b(:) = ra15_plP%b(:, j)
                    g(1) = b(1) + d(1)*b(2) + d(2)*b(3) + d(4)*b(4) + d( 7)*b(5) + d(11)*b(6) + d(16)*b(7)
                    g(2) =             b(2) + d(3)*b(3) + d(5)*b(4) + d( 8)*b(5) + d(12)*b(6) + d(17)*b(7)
                    g(3) =                         b(3) + d(6)*b(4) + d( 9)*b(5) + d(13)*b(6) + d(18)*b(7)
                    g(4) =                                     b(4) + d(10)*b(5) + d(14)*b(6) + d(19)*b(7)
                    g(5) =                                                  b(5) + d(15)*b(6) + d(20)*b(7)
                    g(6) =                                                               b(6) + d(21)*b(7)
                    g(7) =                                                                            b(7)
                    ra15_plP%g(:, j) = g(:)
               END DO
               ra15_plP => ra15_plP%nextP
          END DO
          ra15_tpP => ra15_tp1P
          DO i = 1, ntp
               DO j = 1, NDIM
                    b(:) = ra15_tpP%b(:, j)
                    g(1) = b(1) + d(1)*b(2) + d(2)*b(3) + d(4)*b(4) + d( 7)*b(5) + d(11)*b(6) + d(16)*b(7)
                    g(2) =             b(2) + d(3)*b(3) + d(5)*b(4) + d( 8)*b(5) + d(12)*b(6) + d(17)*b(7)
                    g(3) =                         b(3) + d(6)*b(4) + d( 9)*b(5) + d(13)*b(6) + d(18)*b(7)
                    g(4) =                                     b(4) + d(10)*b(5) + d(14)*b(6) + d(19)*b(7)
                    g(5) =                                                  b(5) + d(15)*b(6) + d(20)*b(7)
                    g(6) =                                                               b(6) + d(21)*b(7)
                    g(7) =                                                                            b(7)
                    ra15_tpP%g(:, j) = g(:)
               END DO
               ra15_tpP => ra15_tpP%nextP
          END DO
          t = tpp
          t2 = t*t
          tval = ABS(t)
          DO m = 1, niter
               DO j = 2, 8
                    jd = j - 1
                    s = H(j)
                    ra15_plP => ra15_pl1P
                    DO i = 1, npl
                         DO k = 1, NDIM
                              b(:) = ra15_plP%b(:, k)
                              a = w(3)*b(3) + s*(w(4)*b(4) + s*(w(5)*b(5) + s*(w(6)*b(6) + s*w(7)*b(7))))
                              y = ra15_plP%xbsav(k) +                                                                             &
                                   s*(t*ra15_plP%vbsav(k) + t2*s*(w1*ra15_plP%absav(k) + s*(w(1)*b(1) + s*(w(2)*b(2) + s*a))))
                              ra15_plP%swifter%xb(k) = y
                              a = u(3)*b(3) + s*(u(4)*b(4) + s*(u(5)*b(5) + s*(u(6)*b(6) + s*u(7)*b(7))))
                              z = ra15_plP%vbsav(k) + s*t*(ra15_plP%absav(k) + s*(u(1)*b(1) + s*(u(2)*b(2) + s*a)))
                              ra15_plP%swifter%vb(k) = z
                         END DO
                         ra15_plP => ra15_plP%nextP
                    END DO
                    ra15_tpP => ra15_tp1P
                    DO i = 1, ntp
                         DO k = 1, NDIM
                              b(:) = ra15_tpP%b(:, k)
                              a = w(3)*b(3) + s*(w(4)*b(4) + s*(w(5)*b(5) + s*(w(6)*b(6) + s*w(7)*b(7))))
                              y = ra15_tpP%xbsav(k) +                                                                             &
                                   s*(t*ra15_tpP%vbsav(k) + t2*s*(w1*ra15_tpP%absav(k) + s*(w(1)*b(1) + s*(w(2)*b(2) + s*a))))
                              ra15_tpP%swifter%xb(k) = y
                              a = u(3)*b(3) + s*(u(4)*b(4) + s*(u(5)*b(5) + s*(u(6)*b(6) + s*u(7)*b(7))))
                              z = ra15_tpP%vbsav(k) + s*t*(ra15_tpP%absav(k) + s*(u(1)*b(1) + s*(u(2)*b(2) + s*a)))
                              ra15_tpP%swifter%vb(k) = z
                         END DO
                         ra15_tpP => ra15_tpP%nextP
                    END DO
                    CALL ra15_getaccb(lextra_force, tm+s*t, npl, nplmax, ra15_pl1P, j2rp2, j4rp4)
                    IF (ntp > 0) CALL ra15_getaccb_tp(lextra_force, tm+s*t, npl, ntp, ntpmax, ra15_pl1P, ra15_tp1P, j2rp2, j4rp4)
                    ra15_plP => ra15_pl1P
                    DO i = 1, npl
                         DO k = 1, NDIM
                              g(:) = ra15_plP%g(:, k)
                              gk = (ra15_plP%ab(k) - ra15_plP%absav(k))/s
                              SELECT CASE (j)
                                   CASE(2)
                                        ra15_plP%g(1, k) =       gk
                                   CASE(3)
                                        ra15_plP%g(2, k) =      (gk - g(1))*r( 1)
                                   CASE(4)
                                        ra15_plP%g(3, k) =     ((gk - g(1))*r( 2) - g(2))*r( 3)
                                   CASE(5)
                                        ra15_plP%g(4, k) =    (((gk - g(1))*r( 4) - g(2))*r( 5) - g(3))*r( 6)
                                   CASE(6)
                                        ra15_plP%g(5, k) =   ((((gk - g(1))*r( 7) - g(2))*r( 8) - g(3))*r( 9) - g(4))*r(10)
                                   CASE(7)
                                        ra15_plP%g(6, k) =  (((((gk - g(1))*r(11) - g(2))*r(12) - g(3))*r(13) - g(4))*r(14) -     &
                                             g(5))*r(15)
                                   CASE(8)
                                        ra15_plP%g(7, k) = ((((((gk - g(1))*r(16) - g(2))*r(17) - g(3))*r(18) - g(4))*r(19) -     &
                                             g(5))*r(20) - g(6))*r(21)
                              END SELECT
                              temp = ra15_plP%g(jd, k) - g(jd)
                              b(:) = ra15_plP%b(:, k)
                              b(jd) = b(jd) + temp
                              SELECT CASE (j)
                                   CASE(3)
                                        b(1) = b(1) + c( 1)*temp
                                   CASE(4)
                                        b(1) = b(1) + c( 2)*temp
                                        b(2) = b(2) + c( 3)*temp
                                   CASE(5)
                                        b(1) = b(1) + c( 4)*temp
                                        b(2) = b(2) + c( 5)*temp
                                        b(3) = b(3) + c( 6)*temp
                                   CASE(6)
                                        b(1) = b(1) + c( 7)*temp
                                        b(2) = b(2) + c( 8)*temp
                                        b(3) = b(3) + c( 9)*temp
                                        b(4) = b(4) + c(10)*temp
                                   CASE(7)
                                        b(1) = b(1) + c(11)*temp
                                        b(2) = b(2) + c(12)*temp
                                        b(3) = b(3) + c(13)*temp
                                        b(4) = b(4) + c(14)*temp
                                        b(5) = b(5) + c(15)*temp
                                   CASE(8)
                                        b(1) = b(1) + c(16)*temp
                                        b(2) = b(2) + c(17)*temp
                                        b(3) = b(3) + c(18)*temp
                                        b(4) = b(4) + c(19)*temp
                                        b(5) = b(5) + c(20)*temp
                                        b(6) = b(6) + c(21)*temp
                              END SELECT
                              ra15_plP%b(:, k) = b(:)
                         END DO
                         ra15_plP => ra15_plP%nextP
                    END DO
                    ra15_tpP => ra15_tp1P
                    DO i = 1, ntp
                         DO k = 1, NDIM
                              g(:) = ra15_tpP%g(:, k)
                              gk = (ra15_tpP%ab(k) - ra15_tpP%absav(k))/s
                              SELECT CASE (j)
                                   CASE(2)
                                        ra15_tpP%g(1, k) =       gk
                                   CASE(3)
                                        ra15_tpP%g(2, k) =      (gk - g(1))*r( 1)
                                   CASE(4)
                                        ra15_tpP%g(3, k) =     ((gk - g(1))*r( 2) - g(2))*r( 3)
                                   CASE(5)
                                        ra15_tpP%g(4, k) =    (((gk - g(1))*r( 4) - g(2))*r( 5) - g(3))*r( 6)
                                   CASE(6)
                                        ra15_tpP%g(5, k) =   ((((gk - g(1))*r( 7) - g(2))*r( 8) - g(3))*r( 9) - g(4))*r(10)
                                   CASE(7)
                                        ra15_tpP%g(6, k) =  (((((gk - g(1))*r(11) - g(2))*r(12) - g(3))*r(13) - g(4))*r(14) -     &
                                             g(5))*r(15)
                                   CASE(8)
                                        ra15_tpP%g(7, k) = ((((((gk - g(1))*r(16) - g(2))*r(17) - g(3))*r(18) - g(4))*r(19) -     &
                                             g(5))*r(20) - g(6))*r(21)
                              END SELECT
                              temp = ra15_tpP%g(jd, k) - g(jd)
                              b(:) = ra15_tpP%b(:, k)
                              b(jd) = b(jd) + temp
                              SELECT CASE (j)
                                   CASE(3)
                                        b(1) = b(1) + c( 1)*temp
                                   CASE(4)
                                        b(1) = b(1) + c( 2)*temp
                                        b(2) = b(2) + c( 3)*temp
                                   CASE(5)
                                        b(1) = b(1) + c( 4)*temp
                                        b(2) = b(2) + c( 5)*temp
                                        b(3) = b(3) + c( 6)*temp
                                   CASE(6)
                                        b(1) = b(1) + c( 7)*temp
                                        b(2) = b(2) + c( 8)*temp
                                        b(3) = b(3) + c( 9)*temp
                                        b(4) = b(4) + c(10)*temp
                                   CASE(7)
                                        b(1) = b(1) + c(11)*temp
                                        b(2) = b(2) + c(12)*temp
                                        b(3) = b(3) + c(13)*temp
                                        b(4) = b(4) + c(14)*temp
                                        b(5) = b(5) + c(15)*temp
                                   CASE(8)
                                        b(1) = b(1) + c(16)*temp
                                        b(2) = b(2) + c(17)*temp
                                        b(3) = b(3) + c(18)*temp
                                        b(4) = b(4) + c(19)*temp
                                        b(5) = b(5) + c(20)*temp
                                        b(6) = b(6) + c(21)*temp
                              END SELECT
                              ra15_tpP%b(:, k) = b(:)
                         END DO
                         ra15_tpP => ra15_tpP%nextP
                    END DO
               END DO
          END DO
          hv = ZERO
          ra15_plP => ra15_pl1P
          DO i = 1, npl
               DO j = 1, NDIM
                    hv = MAX(hv, ABS(ra15_plP%b(7, j)))
               END DO
               ra15_plP => ra15_plP%nextP
          END DO
          ra15_tpP => ra15_tp1P
          DO i = 1, ntp
               DO j = 1, NDIM
                    hv = MAX(hv, ABS(ra15_tpP%b(7, j)))
               END DO
               ra15_tpP => ra15_tpP%nextP
          END DO
          hv = hv*w(7)/tval**7
          tpp = (eps/hv)**PW
          IF ((niter /= 6) .OR. (tpp/t > ONE)) EXIT
          tpp = 0.8_DP*t
     END DO
     ra15_plP => ra15_pl1P
     DO i = 1, npl
          DO j = 1, NDIM
               b(:) = ra15_plP%b(:, j)
               ra15_plP%xbsav(j) = ra15_plP%xbsav(j) + t*ra15_plP%vbsav(j) + t2*(w1*ra15_plP%absav(j) + b(1)*w(1) + b(2)*w(2) +   &
                    b(3)*w(3) + b(4)*w(4) + b(5)*w(5) + b(6)*w(6) + b(7)*w(7))
               ra15_plP%vbsav(j) = ra15_plP%vbsav(j) + t*(ra15_plP%absav(j) + b(1)*u(1) + b(2)*u(2) + b(3)*u(3) + b(4)*u(4) +     &
                    b(5)*u(5) + b(6)*u(6) + b(7)*u(7))
          END DO
          ra15_plP%swifter%xb(:) = ra15_plP%xbsav(:)
          ra15_plP%swifter%vb(:) = ra15_plP%vbsav(:)
          ra15_plP => ra15_plP%nextP
     END DO
     ra15_tpP => ra15_tp1P
     DO i = 1, ntp
          DO j = 1, NDIM
               b(:) = ra15_tpP%b(:, j)
               ra15_tpP%xbsav(j) = ra15_tpP%xbsav(j) + t*ra15_tpP%vbsav(j) + t2*(w1*ra15_tpP%absav(j) + b(1)*w(1) + b(2)*w(2) +   &
                    b(3)*w(3) + b(4)*w(4) + b(5)*w(5) + b(6)*w(6) + b(7)*w(7))
               ra15_tpP%vbsav(j) = ra15_tpP%vbsav(j) + t*(ra15_tpP%absav(j) + b(1)*u(1) + b(2)*u(2) + b(3)*u(3) + b(4)*u(4) +     &
                    b(5)*u(5) + b(6)*u(6) + b(7)*u(7))
          END DO
          ra15_tpP%swifter%xb(:) = ra15_tpP%xbsav(:)
          ra15_tpP%swifter%vb(:) = ra15_tpP%vbsav(:)
          ra15_tpP => ra15_tpP%nextP
     END DO
     tm = tm + t
     x = tm
     hdid = t
     IF (tpp/t > SR) tpp = t*SR
     hnext = tpp
     q = tpp/t
     ra15_plP => ra15_pl1P
     DO i = 1, npl
          DO j = 1, NDIM
               b(:) = ra15_plP%b(:, j)
               bd(:) = ra15_plP%bd(:, j)
               e(:) = ra15_plP%e(:, j)
               IF (niter /= 6) bd(:) = b(:) - e(:)
               e(1) = q*(b(1) + 2.0_DP*b(2) + 3.0_DP*b(3) + 4.0_DP*b(4) +  5.0_DP*b(5) +  6.0_DP*b(6) +  7.0_DP*b(7))
               e(2) =            q**2*(b(2) + 3.0_DP*b(3) + 6.0_DP*b(4) + 10.0_DP*b(5) + 15.0_DP*b(6) + 21.0_DP*b(7))
               e(3) =                          q**3*(b(3) + 4.0_DP*b(4) + 10.0_DP*b(5) + 20.0_DP*b(6) + 35.0_DP*b(7))
               e(4) =                                        q**4*(b(4) +  5.0_DP*b(5) + 15.0_DP*b(6) + 35.0_DP*b(7))
               e(5) =                                                       q**5*(b(5) +  6.0_DP*b(6) + 21.0_DP*b(7))
               e(6) =                                                                      q**6*(b(6) +  7.0_DP*b(7))
               e(7) =                                                                                      q**7*b(7)
               b(:) = e(:) + bd(:)
               ra15_plP%b(:, j) = b(:)
               ra15_plP%bd(:, j) = bd(:)
               ra15_plP%e(:, j) = e(:)
          END DO
          ra15_plP => ra15_plP%nextP
     END DO
     ra15_tpP => ra15_tp1P
     DO i = 1, ntp
          DO j = 1, NDIM
               b(:) = ra15_tpP%b(:, j)
               bd(:) = ra15_tpP%bd(:, j)
               e(:) = ra15_tpP%e(:, j)
               IF (niter /= 6) bd(:) = b(:) - e(:)
               e(1) = q*(b(1) + 2.0_DP*b(2) + 3.0_DP*b(3) + 4.0_DP*b(4) +  5.0_DP*b(5) +  6.0_DP*b(6) +  7.0_DP*b(7))
               e(2) =            q**2*(b(2) + 3.0_DP*b(3) + 6.0_DP*b(4) + 10.0_DP*b(5) + 15.0_DP*b(6) + 21.0_DP*b(7))
               e(3) =                          q**3*(b(3) + 4.0_DP*b(4) + 10.0_DP*b(5) + 20.0_DP*b(6) + 35.0_DP*b(7))
               e(4) =                                        q**4*(b(4) +  5.0_DP*b(5) + 15.0_DP*b(6) + 35.0_DP*b(7))
               e(5) =                                                       q**5*(b(5) +  6.0_DP*b(6) + 21.0_DP*b(7))
               e(6) =                                                                      q**6*(b(6) +  7.0_DP*b(7))
               e(7) =                                                                                      q**7*b(7)
               b(:) = e(:) + bd(:)
               ra15_tpP%b(:, j) = b(:)
               ra15_tpP%bd(:, j) = bd(:)
               ra15_tpP%e(:, j) = e(:)
          END DO
          ra15_tpP => ra15_tpP%nextP
     END DO

     RETURN

END SUBROUTINE ra15_sequence
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
