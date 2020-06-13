!**********************************************************************************************************************************
!
!  Unit Name   : drift_dan
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : drift
!  Language    : Fortran 90/95
!
!  Description : Perform Kepler drift, solving Kepler's equation in appropriate variables
!
!  Input
!    Arguments : mu    : G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
!                x0    : position of body to drift
!                v0    : velocity of body to drift
!                dt0   : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : x0    : position of body to drift
!                v0    : velocity of body to drift
!                iflag : error status flag for Kepler drift (0 = OK, nonzero = NO CONVERGENCE)
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL drift_dan(mu, x0, v0, dt0, iflag)
!
!  Notes       : Adapted from Hal Levison and Martin Duncan's Swift routine drift_dan.f
!
!**********************************************************************************************************************************
pure subroutine drift_dan(mu, x0, v0, dt0, iflag)

! modules
     use swiftest
     use drift, except_this_one => drift_dan
     implicit none

! arguments
     real(DP), intent(in)                     :: mu, dt0
     real(DP), dimension(:), intent(inout) :: x0, v0
     integer(I4B), intent(out) :: iflag

! internals
     real(DP)                  :: dt, f, g, fdot, gdot, c1, c2, c3, u, alpha, fp, r0
     real(DP)                  :: v0s, a, asq, en, dm, ec, es, esq, xkep, fchk, s, c
     real(DP), dimension(NDIM) :: x, v

! executable code
     iflag = 0
     dt = dt0
     r0 = sqrt(dot_product(x0(:), x0(:)))
     v0s = dot_product(v0(:), v0(:))
     u = dot_product(x0(:), v0(:))
     alpha = 2.0_DP*mu/r0 - v0s
     if (alpha > 0.0_DP) then
          a = mu/alpha
          asq = a*a
          en = sqrt(mu/(a*asq))
          ec = 1.0_DP - r0/a
          es = u/(en*asq)
          esq = ec*ec + es*es
          dm = dt*en - int(dt*en/twopi)*twopi
          dt = dm/en
          if ((esq < e2max) .and. (dm*dm < dm2max) .and. (esq*dm*dm < e2dm2max)) then
               call drift_kepmd(dm, es, ec, xkep, s, c)
               fchk = (xkep - ec*s + es*(1.0_DP - c) - dm)
! dek - original code compared fchk*fchk with danbyb, but i think it should
! dek - be compared with danbyb*danbyb, and i changed it accordingly - please
! dek - check with hal and/or martin about this
               if (fchk*fchk > danbyb*danbyb) then
                    iflag = 1
                    return
               end if
               fp = 1.0_DP - ec*c + es*s
               f = a/r0*(c - 1.0_DP) + 1.0_DP
               g = dt + (s - xkep)/en
               fdot = -(a/(r0*fp))*en*s
               gdot = (c - 1.0_DP)/fp + 1.0_DP
               x(:) = x0(:)*f + v0(:)*g
               v(:) = x0(:)*fdot + v0(:)*gdot
               x0(:) = x(:)
               v0(:) = v(:)
               iflag = 0
               return
          end if
     end if
     call drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
     if (iflag == 0) then
          f = 1.0_DP - mu/r0*c2
          g = dt - mu*c3
          fdot = -mu/(fp*r0)*c1
          gdot = 1.0_DP - mu/fp*c2
          x(:) = x0(:)*f + v0(:)*g
          v(:) = x0(:)*fdot + v0(:)*gdot
          x0(:) = x(:)
          v0(:) = v(:)
     end if

     return

end subroutine drift_dan
!**********************************************************************************************************************************
!
!  author(s)   : david e. kaufmann
!
!  revision control system (rcs) information
!
!  source file : $rcsfile$
!  full path   : $source$
!  revision    : $revision$
!  date        : $date$
!  programmer  : $author$
!  locked by   : $locker$
!  state       : $state$
!
!  modification history:
!
!  $log$
!**********************************************************************************************************************************
