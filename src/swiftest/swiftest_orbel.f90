! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_orbel
   real(DP), parameter :: TINYVALUE = 4.0e-15_DP 
      !! Tiny value used to prevent floating point errors. Value set based on the Swifter TINY parameter.
contains

   module subroutine swiftest_orbel_el2xv_vec(self, cb)
      !! author: David A. Minton
      !!
      !! A wrapper method that converts all of the orbital element vectors into cartesian position and velocity vectors for a 
      !! Swiftest body object. This method deallocates all of the orbital elements after it is finished.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self !! Swiftest body object
      class(swiftest_cb),  intent(inout) :: cb !! Swiftest central body objec
      ! Internals
      integer(I4B) :: i, n
   
      if (self%nbody == 0) return
      n = self%nbody

      call self%set_mu(cb)
      associate(mu => self%mu, a => self%a, e => self%e, inc => self%inc, capom => self%capom, omega => self%omega, &
                capm => self%capm, rh => self%rh, vh => self%vh)
#ifdef DOCONLOC
         do concurrent (i = 1:n) shared(mu, a, e, inc, capom, omega, capm, rh, vh)
#else
         do concurrent (i = 1:n)
#endif
            call swiftest_orbel_el2xv(mu(i), a(i), e(i), inc(i), capom(i), omega(i), capm(i), & 
                                      rh(1,i), rh(2,i), rh(3,i), vh(1,i), vh(2,i), vh(3,i))
         end do
      end associate
      return
   end subroutine swiftest_orbel_el2xv_vec


   pure elemental module subroutine swiftest_orbel_el2xv(mu, a, ie, inc, capom, omega, capm, rx, ry, rz, vx, vy, vz)
      !! author: David A. Minton
      !!
      !! Compute osculating orbital elements from relative C)rtesian position and velocity
      !!  All angular measures are returned in radians
      !!      If inclination < TINY, longitude of the ascending node is arbitrarily set to 0
      !!
      !!      If eccentricity < sqrt(TINY), argument of pericenter is arbitrarily set to 0
      !!
      !!      ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
      !!
      !! Adapted from Martin Duncan's el2xv.f
      !! DATE WRITTEN:  May 11, 1992.
      !! REVISIONS: May 26 - now use better Kepler solver for ellipses
      !!  and hyperbolae called EHYBRID.F and FHYBRID.F
      implicit none
      real(DP), intent(in)  :: mu
      real(DP), intent(in)  :: a, ie, inc, capom, omega, capm
      real(DP), intent(out) :: rx, ry, rz, vx, vy, vz

      integer(I4B) :: iorbit_type
      real(DP) :: e, cape, capf, zpara, em1
      real(DP) :: sip, cip, so, co, si, ci
      real(DP) :: d11, d12, d13, d21, d22, d23
      real(DP) :: scap, ccap, shcap, chcap
      real(DP) :: sqe, sqgma, xfac1, xfac2, ri, vfac1, vfac2

      if(ie < 0.0_DP) then
         !write(*,*) ' ERROR in swiftest_orbel_el2xv: e<0, setting e=0!!1'
         e = 0.0_DP
         iorbit_type = ELLIPSE
      else
         e = ie
         em1 = e - 1._DP
         if (abs(em1) < TINYVALUE) then
            iorbit_type = PARABOLA
         else if (e > 1.0_DP)  then
            iorbit_type = HYPERBOLA
         else
            iorbit_type = ELLIPSE
         end if
      endif

      call swiftest_orbel_scget(omega,sip,cip)
      call swiftest_orbel_scget(capom,so,co)
      call swiftest_orbel_scget(inc,si,ci)
      d11 = cip * co - sip * so * ci
      d12 = cip * so + sip * co * ci
      d13 = sip * si
      d21 = -sip * co - cip * so * ci
      d22 = -sip * so + cip * co * ci
      d23 = cip * si

      !--
      ! Get the other quantities depending on orbit type 
      !
      if (iorbit_type == ELLIPSE) then
         cape = swiftest_orbel_ehybrid(e,capm)
         call swiftest_orbel_scget(cape,scap,ccap)
         sqe = sqrt(1._DP - e**2)
         sqgma = sqrt(mu* a)
         xfac1 = a * (ccap - e)
         xfac2 = a * sqe * scap
         ri = 1._DP / (a * (1._DP - e* ccap))
         vfac1 = -ri *  sqgma *  scap
         vfac2 = ri *  sqgma *  sqe *  ccap
      endif
      !--
      if (iorbit_type == HYPERBOLA) then
         capf = swiftest_orbel_fhybrid(e,capm)
         call swiftest_orbel_schget(capf,shcap,chcap)
         sqe = sqrt(e**2 - 1._DP )
         sqgma = sqrt(mu * a)
         xfac1 = a * (e - chcap)
         xfac2 = a * sqe * shcap
         ri = 1._DP / (a * (e * chcap - 1._DP))
         vfac1 = -ri * sqgma * shcap
         vfac2 = ri * sqgma * sqe * chcap
      endif
      !--
      if (iorbit_type == PARABOLA) then
         zpara = swiftest_orbel_zget(capm)
         sqgma = sqrt(2 * mu * a)
         xfac1 = a * (1._DP - zpara * zpara)
         xfac2 = 2 * a * zpara
         ri = 1._DP / (a * (1._DP + zpara * zpara))
         vfac1 = -ri  *  sqgma  *  zpara
         vfac2 = ri  *  sqgma 
      endif
      !--
      rx = d11 * xfac1 + d21 * xfac2
      ry = d12 * xfac1 + d22 * xfac2
      rz = d13 * xfac1 + d23 * xfac2
      vx = d11 * vfac1 + d21 * vfac2
      vy = d12 * vfac1 + d22 * vfac2
      vz = d13 * vfac1 + d23 * vfac2

      return
   end subroutine swiftest_orbel_el2xv


   pure module subroutine swiftest_orbel_scget(angle, sx, cx)
      !! author: David A. Minton
      !!
      !! Efficiently compute the sine and cosine of an input angle
      !!      Input angle must be in radians
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: orbel_scget.f90
      !! Adapted from Hal Levison's Swift routine orbel_scget.f
      implicit none
      ! Arguments
      real(DP), intent(in)  :: angle
      real(DP), intent(out) :: sx, cx
      ! Internals
      integer(I4B) :: nper
      real(DP)   :: x

      nper = angle / TWOPI
      x = angle - nper * TWOPI
      if (x < 0.0_DP) x = x + TWOPI
      sx = sin(x)
      cx = sqrt(1.0_DP - sx**2)
      if ((x > PIBY2) .and. (x < PI3BY2)) cx = -cx

      return

   end subroutine swiftest_orbel_scget


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !**********************************************************************
   !                   swiftest_orbel_SCHGET.F
   !**********************************************************************
   !     PURPOSE:  Given an angle, efficiently compute sinh and cosh.
   ! 
   !        Input:
   !             angle ==> angle in radians (real scalar)
   !
   !        Output:
   !             shx    ==>  sinh(angle)  (real scalar)
   !             chx    ==>  cosh(angle)  (real scalar)
   !
   !     ALGORITHM: Obvious from the code
   !     REMARKS: Based on the routine SCGET for sine's and cosine's.
   !       We use the sqrt rather than cosh (it's faster)
   !       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
   !       OR OVERFLOWS WILL OCCUR!
   !     AUTHOR:  M. Duncan.
   !     DATE WRITTEN:  May 6, 1992.
   !     REVISIONS:
   !**********************************************************************
   pure subroutine swiftest_orbel_schget(angle,shx,chx)
      real(DP), intent(in)  ::  angle
      real(DP), intent(out) :: shx,chx
      
      shx = sinh(angle)
      chx= sqrt(1._DP + shx * shx)
      
      return
   end subroutine swiftest_orbel_schget


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !        !                    swiftest_orbel_FLON.F
   !**********************************************************************
   !     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
   !
   !             Input:
   !                           e ==> eccentricity anomaly. (real scalar)
   !                        capn ==> hyperbola mean anomaly. (real scalar)
   !             Returns:
   !                  swiftest_orbel_flon ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: Uses power series for N in terms of F and Newton,s method
   !     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 26, 1992.
   !     REVISIONS:
   !**********************************************************************
   real(DP) pure function swiftest_orbel_flon(e,icapn)
      implicit none
      real(DP), intent(in) ::  e, icapn
      integer(I4B) :: iflag,i
      real(DP) ::  a,b,sq,biga,bigb, capn
      real(DP) ::  x,x2
      real(DP) ::  f,fp,dx
      real(DP) ::  a0,a1
      real(DP) ::  b1
      integer(I4B), parameter :: IMAX = 10
      real(DP), parameter :: a11 = 156._DP, a9 = 17160._DP, a7 = 1235520._DP
      real(DP), parameter :: a5 = 51891840._DP,  a3 = 1037836800._DP
      real(DP), parameter :: b11 = 11 * a11, b9 = 9 * a9, b7 = 7 * a7
      real(DP), parameter :: b5 = 5 * a5, b3 = 3 * a3

      ! Function to solve "Kepler's eqn" for F (here called
      ! x) for given e and CAPN. Only good for smallish CAPN

      iflag = 0
      if (icapn < 0._DP) then
         iflag = 1
         capn = -icapn
      else
         capn = icapn
      end if

      a1 = 6227020800._DP * (1._DP - 1._DP / e)
      a0 = -6227020800._DP * capn / e
      b1 = a1

      !  set iflag nonzero if capn < 0., in which case solve for -capn
      !  and change the sign of the final answer for f.
      !  Begin with a reasonable guess based on solving the cubic for small F
      a = 6 * ( e - 1.d0) / e
      b = -6 * capn / e
      sq = SQRT(0.25_DP * b**2 + a**3 / 27._DP)
      biga =  (-0.5_DP * b + sq)**(1.0_DP / 3.0_DP)
      bigb = -(+0.5_DP * b + sq)**(1.0_DP / 3.0_DP) 
      x = biga + bigb
      swiftest_orbel_flon = x
      ! If capn is TINYVALUE (or zero) no need to go further than cubic even for
      ! e =1.
      if( capn >= TINYVALUE) then
         do i = 1,IMAX
            x2 = x**2
            f = a0 + x * (a1 + x2 * (a3 + x2 * (a5 + x2 * (a7 + x2 * (a9 + x2 * (a11 + x2))))))
            fp = b1 + x2 * (b3 + x2 * (b5 + x2 * (b7 + x2 * (b9 + x2 * (b11 + 13 * x2)))))
            dx = -f / fp
            swiftest_orbel_flon = x + dx
            !   if we have converged here there's no point in going on
            if(abs(dx) <= TINYVALUE) exit
            x = swiftest_orbel_flon
         end do
      end if

      !  normal return here, but check if capn was originally negative
      if(iflag == 1) then
         swiftest_orbel_flon = -swiftest_orbel_flon
         capn = -capn
      end if

      return
   end function swiftest_orbel_flon


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
            !                    swiftest_orbel_FGET.F
   !**********************************************************************
   !     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
   !
   !             Input:
   !                           e ==> eccentricity anomaly. (real scalar)
   !                        capn ==> hyperbola mean anomaly. (real scalar)
   !             Returns:
   !                  swiftest_orbel_fget ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
   !           Cel. Mech. ".  Quartic convergence from Danby's book.
   !     REMARKS:
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 11, 1992.
   !     REVISIONS: 2/26/93 hfl
   !**********************************************************************
   real(DP) pure function swiftest_orbel_fget(e,capn)
      implicit none

      real(DP), intent(in) ::  e,capn

      integer :: i
      real(DP) ::  tmp,x,shx,chx
      real(DP) ::  esh,ech,f,fp,fpp,fppp,dx
      integer(I4B), parameter :: IMAX = 10

      !----
      !...  executable code

      ! function to solve "kepler's eqn" for f (here called
      ! x) for given e and capn.

      !  begin with a guess proposed by danby
      if( capn < 0.d0) then
         tmp = -2 * capn / e + 1.8_DP
         x = -log(tmp)
      else
         tmp = +2 * capn / e + 1.8_DP
         x = log(tmp)
      end if

      swiftest_orbel_fget = x

      do i = 1, IMAX
         call swiftest_orbel_schget(x,shx,chx)
         esh = e * shx
         ech = e * chx
         f = esh - x - capn
      !   write(6,*) 'i,x,f : ',i,x,f
         fp = ech - 1.d0
         fpp = esh
         fppp = ech
         dx = -f / fp
         dx = -f / (fp + dx * fpp / 2._DP)
         dx = -f / (fp + dx * fpp / 2._DP + dx**2 * fppp / 6._DP)
         swiftest_orbel_fget = x + dx
      !   if we have converged here there's no point in going on
         if(abs(dx) <= TINYVALUE) return
         x = swiftest_orbel_fget
      end do

      !write(*,*) 'fget : returning without complete convergence'
      return
   end function swiftest_orbel_fget


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !                    swiftest_orbel_ZGET.F
   !**********************************************************************
   !     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola
   !          given Q (Fitz. notation.)
   !  
   !             Input:
   !                           q ==>  parabola mean anomaly. (real scalar)
   !             Returns:
   !                  swiftest_orbel_zget ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
   !     REMARKS: For a parabola we can solve analytically.
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 11, 1992.
   !     REVISIONS: May 27 - corrected it for negative Q and use power
   !       series for small Q.
   !**********************************************************************
   real(DP) pure function swiftest_orbel_zget(iq)
      implicit none

      real(DP), intent(in)  :: iq

      integer(I4B) :: iflag
      real(DP) ::  x,tmp,q
      
      iflag = 0
      if (iq < 0.0_DP) then
         iflag = 1
         q = -iq
      else
         q = iq
      end if

      if (q < 1.e-3_DP) then
         swiftest_orbel_zget = q * (1._DP - (q**2 / 3._DP) * (1._DP - q**2))
      else
         x = 0.5_DP * (3 * q + sqrt(9 * q**2 + 4._DP))
         tmp = x**(1._DP / 3._DP)
         swiftest_orbel_zget = tmp - 1._DP / tmp
      end if

      if(iflag == 1) then
         swiftest_orbel_zget = -swiftest_orbel_zget
         q = -q
      end if

      return
   end function swiftest_orbel_zget


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !                    swiftest_orbel_ESOLMD.F
   !**********************************************************************
   !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
   !
   !             Input:
   !                           e ==> eccentricity anomaly. (real scalar)
   !                           m ==> mean anomaly. (real scalar)
   !             Returns:
   !                swiftest_orbel_esolmd ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: Some sort of quartic convergence from Wisdom.
   !     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
   !         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
   !         ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI
   !     INCLUDES: needs SCGET.F
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 7, 1992.
   !     REVISIONS: 2/26/93 hfl
   !**********************************************************************
   real(DP) pure function swiftest_orbel_esolmd(e,m)
      implicit none

      real(DP), intent(in)  :: e
      real(DP), intent(in)  :: m

      real(DP) :: x,sm,cm,sx,cx
      real(DP) :: es,ec,f,fp,fpp,fppp,dx

      !...    function to solve kepler's eqn for e (here called
      !...    x) for given e and m. returns value of x.

      call swiftest_orbel_scget(m,sm,cm)
      x = m + e * sm * (1._DP + e * ( cm + e * (1._DP - 1.5_DP * sm**2)))

      call swiftest_orbel_scget(x,sx,cx)
      es = e * sx
      ec = e * cx
      f = x - es  - m
      fp = 1._DP - ec
      fpp = es
      fppp = ec
      dx = -f / fp
      dx = -f / (fp + dx * fpp / 2._DP)
      dx = -f / (fp + dx * fpp / 2._DP + dx**2 * fppp / 6._DP)

      swiftest_orbel_esolmd = x + dx

      return   
   end function swiftest_orbel_esolmd


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !                    swiftest_orbel_EHIE.F
   !**********************************************************************
   !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
   !
   !             Input:
   !                           e ==> eccentricity anomaly. (real scalar)
   !                           m ==> mean anomaly. (real scalar)
   !             Returns:
   !              swiftest_orbel_ehybrid ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: Use Danby's quartic for 3 iterations.
   !                Eqn. is f(x) = x - e*sin(x+M). Note  that
   !          E = x + M. First guess is very good for e near 1.
   !          Need to first get M between 0. and PI and use
   !   symmetry to return right answer if M between PI and 2PI
   !     REMARKS: Modifies M so that both E and M are in range (0,TWOPI)
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 25,1992.
   !     REVISIONS:
   !**********************************************************************
   real(DP) pure function swiftest_orbel_ehie(e,im)
      implicit none

      real(DP), intent(in) :: e,im

      integer(I4B) :: iflag,nper,niter
      real(DP) :: dx,x,sa,ca,esa,eca,f,fp,m

      integer(I4B), parameter  :: NMAX = 3

      ! in this section, bring m into the range (0,TWOPI) and if
      ! the result is greater than pi, solve for (TWOPI - m).
      iflag = 0
      nper = im / TWOPI
      m = im - nper * TWOPI
      if (m < 0._DP) m = m + TWOPI

      if (m > PI) then
         m = TWOPI - m
         iflag = 1
      end if

      ! make a first guess that works well for e near 1.
      x = (6 * m)**(1._DP / 3._DP) - m
      niter =0

      ! iteration loop
      do niter =1,NMAX
         call swiftest_orbel_scget(x + m,sa,ca)
         esa = e * sa
         eca = e * ca
         f = x - esa
         fp = 1._DP -eca
         dx = -f / fp
         dx = -f / (fp + 0.5_DP * dx * esa)
         dx = -f / (fp + 0.5_DP * dx * (esa + eca * dx / 3.0_DP))
         x = x + dx
      end do

      swiftest_orbel_ehie = m + x

      if (iflag == 1) then
         swiftest_orbel_ehie = TWOPI - swiftest_orbel_ehie
         m = TWOPI - m
      end if

      return
   end function swiftest_orbel_ehie


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !                             swiftest_orbel_EGET.F
   !**********************************************************************
   !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
   !
   !             Input:
   !                           e ==> eccentricity anomaly. (real scalar)
   !                           m ==> mean anomaly. (real scalar)
   !             Returns:
   !                  swiftest_orbel_eget ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: Quartic convergence from Danby
   !     REMARKS: For results very near roundoff, give it M between
   !           0 and 2*pi. One can condition M before calling EGET
   !           by calling my double precision function MOD2PI(M).
   !           This is not done within the routine to speed it up
   !           and because it works fine even for large M.
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 7, 1992.
   !     REVISIONS: May 21, 1992.  Now have it go through EXACTLY two iterations
   !                with the premise that it will only be called if
   !          we have an ellipse with e between 0.15 and 0.8
   !**********************************************************************
   real(DP) pure function swiftest_orbel_eget(e,m)
      implicit none
      
      real(DP), intent(in) ::  e,m
      real(DP) ::  x,sm,cm,sx,cx
      real(DP) ::  es,ec,f,fp,fpp,fppp,dx

      ! function to solve kepler's eqn for e (here called
      ! x) for given e and m. returns value of x.
      ! may 21 : for e < 0.18 use esolmd for speed and sufficient accuracy
      ! may 21 : for e > 0.8 use ehie - this one may not converge fast enough.

      call swiftest_orbel_scget(m,sm,cm)

      !  begin with a guess accurate to order ecc**3
      x = m + e * sm * ( 1._DP + e * (cm + e * (1._DP - 1.5_DP * sm * sm)))

      !  go through one iteration for improved estimate
      call swiftest_orbel_scget(x,sx,cx)
      es = e * sx
      ec = e * cx
      f = x - es  - m
      fp = 1._DP - ec
      fpp = es
      fppp = ec
      dx = -f / fp
      dx = -f / (fp + dx * fpp / 2._DP)
      dx = -f / (fp + dx * fpp / 2._DP + dx*2 * fppp / 6._DP)
      swiftest_orbel_eget = x + dx

      ! do another iteration.
      ! for m between 0 and 2*pi this seems to be enough to
      ! get near roundoff error for eccentricities between 0 and 0.8

      x = swiftest_orbel_eget
      call swiftest_orbel_scget(x,sx,cx)
      es = e * sx
      ec = e * cx
      f = x - es  - m
      fp = 1._DP - ec
      fpp = es
      fppp = ec
      dx = -f / fp
      dx = -f / (fp + dx * fpp / 2._DP)
      dx = -f / (fp + dx * fpp / 2._DP + dx**2 * fppp / 6._DP)

      swiftest_orbel_eget = x + dx

      return
   end function swiftest_orbel_eget


   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !                    swiftest_orbel_EHYBRID.F
   !**********************************************************************
   !     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
   !
   !             Input:
   !                           e ==> eccentricity anomaly. (real scalar)
   !                           m ==> mean anomaly. (real scalar)
   !             Returns:
   !              swiftest_orbel_ehybrid ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: For e < 0.18 uses fast routine ESOLMD
   !          For larger e but less than 0.8, uses EGET
   !          For e > 0.8 uses EHIE
   !     REMARKS: Only EHIE brings M and E into range (0,TWOPI)
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 25,1992.
   !     REVISIONS: 2/26/93 hfl
   !**********************************************************************
   real(DP) pure function swiftest_orbel_ehybrid(e,m)
      implicit none

      real(DP), intent(in) :: e,m
      !real(DP) :: swiftest_orbel_esolmd,orbel_eget,orbel_ehie

      if (e < 0.18_DP) then
         swiftest_orbel_ehybrid = swiftest_orbel_esolmd(e,m)
      else
         if( e <= 0.8_DP) then
            swiftest_orbel_ehybrid = swiftest_orbel_eget(e,m)
         else
            swiftest_orbel_ehybrid = swiftest_orbel_ehie(e,m)
         end if
      end if
      return
   end function swiftest_orbel_ehybrid

   !**********************************************************************
   ! Code converted to Modern Fortran by David A. Minton
   ! Date: 2020-06-29  
   !                    swiftest_orbel_FHYBRID.F
   !**********************************************************************
   !     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
   !
   !             Input:
   !                           e ==> eccentricity anomaly. (real scalar)
   !                           n ==> hyperbola mean anomaly. (real scalar)
   !             Returns:
   !               swiftest_orbel_fhybrid ==>  eccentric anomaly. (real scalar)
   !
   !     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON
   !          For larger N, uses FGET
   !     REMARKS:
   !     AUTHOR: M. Duncan
   !     DATE WRITTEN: May 26,1992.
   !     REVISIONS::
   !     REVISIONS: 2/26/93 hfl
   !**********************************************************************
   real(DP) pure function swiftest_orbel_fhybrid(e,n)
      implicit none
      real(DP), intent(in) :: e,n

      real(DP) :: abn

      abn = n
      if(n < 0._DP) abn = -abn

      if(abn < 0.636_DP * e - 0.6_DP) then
         swiftest_orbel_fhybrid = swiftest_orbel_flon(e,n)
      else
         swiftest_orbel_fhybrid = swiftest_orbel_fget(e,n)
      end if

      return
   end function swiftest_orbel_fhybrid
   

   pure elemental module subroutine swiftest_orbel_xv2aeq(mu, rx, ry, rz, vx, vy, vz, a, e, q)
      !! author: David A. Minton
      !!
      !! Compute semimajor axis, eccentricity, and pericentric distance from relative Cartesian position and velocity
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: orbel_xv2aeq.f90
      !! Adapted from Luke Dones' Swift routine orbel_xv2aeq.f
      implicit none
      !! Arguments
      real(DP), intent(in)  :: mu !! Gravitational constant
      real(DP), intent(in)  :: rx,ry,rz  !! Position vector
      real(DP), intent(in)  :: vx,vy,vz  !! Velocity vector
      real(DP), intent(out) :: a  !! semimajor axis
      real(DP), intent(out) :: e  !! eccentricity
      real(DP), intent(out) :: q  !! periapsis
      ! Internals
      integer(I4B) :: iorbit_type
      real(DP)   :: hx, hy, hz, r, v2, h2, energy, fac

      a = 0.0_DP
      e = 0.0_DP
      q = 0.0_DP

      r = sqrt(rx*rx + ry*ry + rz*rz)
      v2 = vx*vx + vy*vy + vz*vz

      hx = ry*vz - rz*vy
      hy = rz*vx - rx*vz
      hz = rx*vy - ry*vx
      h2 = hx*hx + hy*hy + hz*hz

      if (h2 < tiny(h2)) return
      energy = 0.5_DP * v2 - mu / r
      if (abs(energy * r / mu) < sqrt(TINYVALUE)) then
         iorbit_type = PARABOLA
      else
         a = -0.5_DP * mu / energy
         if (a < 0.0_DP) then
            fac = -h2 / (mu * a)
            if (fac > TINYVALUE) then
               iorbit_type = HYPERBOLA
            else
               iorbit_type = PARABOLA
            end if
         else
            iorbit_type = ELLIPSE
         end if
      end if
      select case (iorbit_type)
      case (ELLIPSE)
         fac = 1.0_DP - h2 / (mu * a)
         if (fac > TINYVALUE) e = sqrt(fac)
         q = a * (1.0_DP - e)
      case (PARABOLA)
         a = 0.5_DP * h2 / mu
         e = 1.0_DP
         q = a
      case (HYPERBOLA)
         e = sqrt(1.0_DP + fac)
         q = a * (1.0_DP - e)
      end select

      return

   end subroutine swiftest_orbel_xv2aeq


   pure module subroutine swiftest_orbel_xv2aqt(mu, rx, ry, rz, vx, vy, vz, a, q, capm, tperi)
      !! author: David A. Minton
      !!
      !! Compute semimajor axis, pericentric distance, mean anomaly, and time to nearest pericenter passage from
      !! relative Cartesian position and velocity
      !!      tperi > 0 means nearest pericenter passage is in the future
      !!      tperi < 0 means nearest pericenter passage is in the past
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: orbel_xv2aqt.f90
      implicit none
      ! Arguments
      real(DP), intent(in)  :: mu    !! Gravitational constant
      real(DP), intent(in)  :: rx,ry,rz !! Position vector
      real(DP), intent(in)  :: vx,vy,vz !! Velocity vector
      real(DP), intent(out) :: a     !! semimajor axis
      real(DP), intent(out) :: q     !! periapsis
      real(DP), intent(out) :: capm  !! mean anomaly
      real(DP), intent(out) :: tperi !! time of pericenter passage
      ! Internals
      integer(I4B) :: iorbit_type
      real(DP)   :: hx, hy, hz, r, v2, h2, rdotv, energy, fac, w, face, cape, e, tmpf, capf, mm

      a = 0.0_DP
      q = 0.0_DP
      capm = 0.0_DP
      tperi = 0.0_DP
      hx = ry*vz - rz*vy
      hy = rz*vx - rx*vz
      hz = rx*vy - ry*vx
      h2 = hx*hx + hy*hy + hz*hz
      if (h2 < tiny(h2)) return

      r = sqrt(rx*rx + ry*ry + rz*rz)
      v2 = vx*vx + vy*vy + vz*vz
      rdotv = rx*vx + ry*vy + rz*vz 
      energy = 0.5_DP * v2 - mu / r
      if (abs(energy * r / mu) < sqrt(TINYVALUE)) then
         iorbit_type = PARABOLA
      else
         a = -0.5_DP * mu / energy
         if (a < 0.0_DP) then
            fac = -h2 / (mu * a)
            if (fac > TINYVALUE) then
               iorbit_type = HYPERBOLA
            else
               iorbit_type = PARABOLA
            end if
         else
            iorbit_type = ELLIPSE
         end if
      end if
      select case (iorbit_type)
      case (ELLIPSE)
         fac = 1.0_DP - h2 / (mu * a)
         if (fac > TINYVALUE) then
            e = sqrt(fac)
            cape = 0.0_DP
            face = (a - r) / (a * e)
            if (face < -1.0_DP) then
               cape = PI
            else if (face < 1.0_DP) then
               cape = acos(face)
            end if
            if (rdotv < 0.0_DP) cape = TWOPI - cape
         else
            e = 0.0_DP
            cape = 0.0_DP
         end if
         capm = cape - e * sin(cape)
         q = a * (1.0_DP - e)
         mm = sqrt(mu / a**3)
         if (capm < PI) then
            tperi = -1.0_DP * capm / mm
         else
            tperi = -1.0_DP * (capm - TWOPI) / mm
         end if
      case (PARABOLA)
         a = 0.5_DP * h2 / mu
         e = 1.0_DP
         w = 0.0_DP
         fac = 2 * a / r - 1.0_DP
         if (fac < -1.0_DP) then
            w = PI
         else if (fac < 1.0_DP) then
            w = acos(fac)
         end if
         if (rdotv < 0.0_DP) w = TWOPI - w
         tmpf = tan(0.5_DP * w)
         capm = tmpf*(1.0_DP + tmpf * tmpf / 3.0_DP)
         q = a
         mm = sqrt(0.5_DP * mu / q**3)
         tperi = -1.0_DP * capm / mm
      case (HYPERBOLA)
         e = sqrt(1.0_DP + fac)
         tmpf = (a - r) / (a * e)
         if (tmpf < 1.0_DP) tmpf = 1.0_DP
         capf = log(tmpf + sqrt(tmpf * tmpf - 1.0_DP))
         if (rdotv < 0.0_DP) capf = -capf
         capm = e * sinh(capf) - capf
         q = a * (1.0_DP - e)
         mm = sqrt(-mu / a**3)
         tperi = -1.0_DP * capm / mm
      end select

      return

   end subroutine swiftest_orbel_xv2aqt


   module subroutine swiftest_orbel_xv2el_vec(self, cb)
      !! author: David A. Minton
      !!
      !! A wrapper method that converts all of the cartesian position and velocity vectors of a Swiftest body object to orbital 
      !! elements.
      implicit none
      ! Arguments
      class(swiftest_body), intent(inout) :: self !! Swiftest body object
      class(swiftest_cb),   intent(inout) :: cb   !! Swiftest central body object
      ! internals
      integer(I4B) :: i
      real(DP) :: varpi, lam, f, cape, capf
    
      if (self%nbody == 0) return

      call self%set_mu(cb)
      if (allocated(self%a)) deallocate(self%a); allocate(self%a(self%nbody))
      if (allocated(self%e)) deallocate(self%e); allocate(self%e(self%nbody))
      if (allocated(self%inc)) deallocate(self%inc); allocate(self%inc(self%nbody))
      if (allocated(self%capom)) deallocate(self%capom); allocate(self%capom(self%nbody))
      if (allocated(self%omega)) deallocate(self%omega); allocate(self%omega(self%nbody))
      if (allocated(self%capm)) deallocate(self%capm);  allocate(self%capm(self%nbody))
#ifdef DOCONLOC
      do concurrent (i = 1:self%nbody) shared(self) local(varpi,lam,f,cape,capf)
#else
      do concurrent (i = 1:self%nbody)
#endif
         call swiftest_orbel_xv2el(self%mu(i), self%rh(1,i), self%rh(2,i), self%rh(3,i), &
                                      self%vh(1,i), self%vh(2,i), self%vh(3,i), &
                                      self%a(i), self%e(i), self%inc(i),  &
                                      self%capom(i), self%omega(i), self%capm(i), &
                                      varpi, lam, f, cape, capf)
      end do

      return
   end subroutine swiftest_orbel_xv2el_vec 


   pure elemental module subroutine swiftest_orbel_xv2el(mu, rx, ry, rz, vx, vy, vz, &
                                                         a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf)
      !! author: David A. Minton
      !!
      !! Compute osculating orbital elements from relative Cartesian position and velocity
      !!  All angular measures are returned in radians
      !!      If inclination < TINY, longitude of the ascending node is arbitrarily set to 0
      !!
      !!      If eccentricity < sqrt(TINY), argument of pericenter is arbitrarily set to 0
      !!
      !!      References: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 201 - 206.
      !!              Fitzpatrick, P. M. 1970. Principles of Celestial Mechanics, (Academic Press), 69 - 73.
      !!              Roy, A. E. 1982. Orbital Motion, (Adam Hilger, Ltd.), 75 - 95
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: orbel_xv2el.f90
      !! Adapted from Martin Duncan's Swift routine orbel_xv2el.f
      implicit none
      ! Arguments
      real(DP), intent(in)  :: mu    !! Gravitational constant
      real(DP), intent(in)  :: rx,ry,rz    !! Position vector
      real(DP), intent(in)  :: vx,vy,vz !! Velocity vector
      real(DP), intent(out) :: a     !! semimajor axis
      real(DP), intent(out) :: e     !! eccentricity
      real(DP), intent(out) :: inc   !! inclination
      real(DP), intent(out) :: capom !! longitude of ascending node
      real(DP), intent(out) :: omega !! argument of periapsis
      real(DP), intent(out) :: capm  !! mean anomaly
      real(DP), intent(out) :: varpi !! longitude of periapsis
      real(DP), intent(out) :: lam   !! mean longitude
      real(DP), intent(out) :: f     !! true anomaly
      real(DP), intent(out) :: cape  !! eccentric anomaly (eccentric orbits)
      real(DP), intent(out) :: capf  !! hyperbolic anomaly (hyperbolic orbits)
      ! Internals
      integer(I4B) :: iorbit_type
      real(DP)   :: hx, hy, hz, r, v2, h2, h, rdotv, energy, fac, u, w, cw, sw, face, tmpf, sf, cf, rdot

      a = 0.0_DP
      e = 0.0_DP
      inc = 0.0_DP
      capom = 0.0_DP
      omega = 0.0_DP
      capm = 0.0_DP
      varpi = 0.0_DP
      lam = 0.0_DP
      f = 0.0_DP
      cape = 0.0_DP
      capf = 0.0_DP

      hx = ry*vz - rz*vy
      hy = rz*vx - rx*vz
      hz = rx*vy - ry*vx
      h2 = hx*hx + hy*hy +hz*hz
      h  = sqrt(h2)
      if(hz>h) then                 ! Hal's fix
         hz = h
         hx = 0.0_DP
         hy = 0.0_DP
      endif
      if (h2 <= 10 * tiny(0.0_DP)) return
      h = SQRT(h2)

      r = sqrt(rx*rx + ry*ry + rz*rz)
      v2 = vx*vx + vy*vy + vz*vz
      rdotv = rx*vx + ry*vy + rz*vz
      energy = 0.5_DP * v2 - mu / r
      fac = hz / h
      if (fac < -1.0_DP) then
         inc = PI
      else if (fac < 1.0_DP) then
         inc = acos(fac)
      end if
      fac = sqrt(hx**2 + hy**2) / h
      if (fac**2 < TINYVALUE) then
         u = atan2(ry, rx)
         if (hz < 0.0_DP) u = -u
      else
         capom = atan2(hx, -hy)
         u = atan2(rz / sin(inc), rx * cos(capom) + ry * sin(capom))
      end if
      if (capom < 0.0_DP) capom = capom + TWOPI
      if (u < 0.0_DP) u = u + TWOPI
      if (abs(energy * r / mu) < sqrt(TINYVALUE)) then
         iorbit_type = parabola
      else
         a = -0.5_DP * mu / energy
         if (a < 0.0_DP) then
            fac = -h2 / (mu * a)
            if (fac > TINYVALUE) then
               iorbit_type = HYPERBOLA
            else
               iorbit_type = PARABOLA
            end if
         else
            iorbit_type = ELLIPSE
         end if
      end if
      select case (iorbit_type)
         case (ELLIPSE)
            fac = 1.0_DP - h2 / (mu * a)
            if (fac > TINYVALUE) then
               e = sqrt(fac)
               cape = 0.0_DP
               face = (a - r) / (a * e)
               if (face < -1.0_DP) then
                  cape = PI
               else if (face < 1.0_DP) then
                  cape = acos(face)
               end if
               if (rdotv < 0.0_DP) cape = TWOPI - cape
               fac = 1.0_DP - e * cos(cape)
               cw = (cos(cape) - e) / fac
               sw = sqrt(1.0_DP - e**2) * sin(cape) / fac
               w = atan2(sw, cw)
               if (w < 0.0_DP) w = w + TWOPI
            else
               cape = u
               w = u
            end if
            capm = cape - e * sin(cape)
         case (PARABOLA)
            a = 0.5_DP * h2 / mu
            e = 1.0_DP
            w = 0.0_DP
            fac = 2 * a / r - 1.0_DP
            if (fac < -1.0_DP) then
               w = PI
            else if (fac < 1.0_DP) then
               w = acos(fac)
            end if
            if (rdotv < 0.0_DP) w = TWOPI - w
            tmpf = tan(0.5_DP * w)
            capm = tmpf * (1.0_DP + tmpf * tmpf / 3.0_DP)
         case (HYPERBOLA)
            e = sqrt(1.0_DP + fac)
            tmpf = max((a - r) / (a * e), 1.0_DP)
            capf = log(tmpf + sqrt(tmpf**2 - 1.0_DP))
            if (rdotv < 0.0_DP) capf = -capf
            fac = e * cosh(capf) - 1.0_DP
            cw = (e - cosh(capf)) / fac
            sw = sqrt(e * e - 1.0_DP) * sinh(capf) / fac
            w = atan2(sw, cw)
            if (w < 0.0_DP) w = w + TWOPI
            capm = e * sinh(capf) - capf
            a = -a
      end select
      omega = u - w
      if (omega < 0.0_DP) omega = omega + TWOPI
      varpi = mod(omega + capom, TWOPI)
      lam = mod(capm + varpi, TWOPI)
      if (e > TINYVALUE) then
         cf = 1.0_DP / e * (a * (1.0_DP - e**2)/r - 1.0_DP)
         rdot = sign(sqrt(max(v2 - (h / r)**2,0.0_DP)),rdotv)
         sf = a * (1.0_DP - e**2) / (h * e) * rdot
         f = atan2(sf,cf)
         if (f < 0.0_DP) f = f + TWOPI
      else
         f = u
      end if


      return
   end subroutine swiftest_orbel_xv2el

   
end submodule s_swiftest_orbel
