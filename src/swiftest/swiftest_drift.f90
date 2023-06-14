!! Coryright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a cory of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_drift
   !> Integration control parameters:
   real(DP), parameter :: E2MAX    = 0.36_DP      
   real(DP), parameter :: DM2MAX   = 0.16_DP
   real(DP), parameter :: E2DM2MAX = 0.0016_DP
   real(DP),     parameter :: DANBYB   = 1.0E-13_DP
   integer(I2B), parameter :: NLAG1    = 50
   integer(I2B), parameter :: NLAG2    = 40

contains

   module subroutine swiftest_drift_body(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Loop bodies and call Danby drift routine on the heliocentric position and velocities.
      !!
      !! Adapted from Hal Levison's Swift routine drift_tp.f 
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f90
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! Swiftest test particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B)                              :: i   
      integer(I4B), dimension(:), allocatable   :: iflag
      character(len=STRMAX) :: message

      associate(n => self%nbody)
         allocate(iflag(n))
         iflag(:) = 0
         call swiftest_drift_all(self%mu, self%rh, self%vh, self%nbody, param, dt, self%lmask, iflag)
         if (any(iflag(1:n) /= 0)) then
            where(iflag(1:n) /= 0) self%status(1:n) = DISCARDED_DRIFTERR
            do i = 1, n
               if (iflag(i) /= 0) then
                  write(message, *) " Body ", self%id(i), " lost due to error in Danby drift"
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
               end if
            end do
         end if
      end associate

      deallocate(iflag)

      return
   end subroutine swiftest_drift_body


   module subroutine swiftest_drift_all(mu, x, v, n, param, dt, lmask, iflag)
      !! author: David A. Minton
      !!
      !! Loop bodies and call Danby drift routine on all bodies for the given position and velocity vector.
      !!
      !! Adapted from Hal Levison's Swift routine drift_tp.f 
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f9
      implicit none
      ! Arguments
      real(DP), dimension(:),     intent(in)    :: mu    !! Vector of gravitational constants
      real(DP), dimension(:,:),   intent(inout) :: x, v  !! Position and velocity vectors
      integer(I4B),               intent(in)    :: n     !! number of bodies
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      real(DP),                   intent(in)    :: dt    !! Stepsize
      logical, dimension(:),      intent(in)    :: lmask !! Logical mask of size self%nbody that determines which bodies to drift.
      integer(I4B), dimension(:), intent(out)   :: iflag !! Vector of error flags. 0 means no problem
      ! Internals
      integer(I4B)                              :: i   
      real(DP)                                  :: energy, vmag2, rmag  !! Variables used in GR calculation
      real(DP), dimension(:), allocatable       :: dtp

      if (n == 0) return

      allocate(dtp(n))
      if (param%lgr) then
#ifdef DOCONLOC
         do concurrent(i = 1:n, lmask(i)) shared(param,lmask,x,v,mu,dtp,dt) local(rmag,vmag2,energy)
#else
         do concurrent(i = 1:n, lmask(i))
#endif
            rmag = norm2(x(:, i))
            vmag2 = dot_product(v(:, i), v(:, i))
            energy = 0.5_DP * vmag2 - mu(i) / rmag
            dtp(i) = dt * (1.0_DP + 3 * param%inv_c2 * energy)
         end do
      else
         where(lmask(1:n)) dtp(1:n) = dt
      end if 

      !!$omp simd ! SIMD does not yet work
      do i = 1, n
         if (lmask(i)) call swiftest_drift_one(mu(i), x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i), dtp(i), iflag(i))
      end do
      !!$omp end simd

      deallocate(dtp)

      return
   end subroutine swiftest_drift_all


   pure elemental module subroutine swiftest_drift_one(mu, rx, ry, rz, vx, vy, vz, dt, iflag) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Perform Danby drift for one body, redoing drift with smaller substeps if original accuracy is insufficient
      !!
      !! Adapted from David E. Kaufmann's Swifter routine routine drift_one.f90
      !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_one.f
      implicit none
      ! Arguments
      real(DP), intent(in)      :: mu                     !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body to drift
      real(DP), intent(inout)   :: rx, ry, rz, vx, vy, vz !! Position and velocity of body to drift
      real(DP), intent(in)      :: dt                     !! Step size
      integer(I4B), intent(out) :: iflag                  !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      ! Internals
      integer(I4B) :: i
      real(DP)   :: dttmp

      call swiftest_drift_dan(mu, rx, ry, rz, vx, vy, vz, dt, iflag)
      if (iflag /= 0) then
         dttmp = 0.1_DP * dt
         do i = 1, 10
            call swiftest_drift_dan(mu, rx, ry, rz, vx, vy, vz, dttmp, iflag)
            if (iflag /= 0) exit
         end do
      end if

      return
   end subroutine swiftest_drift_one


   pure subroutine swiftest_drift_dan(mu, rx0, ry0, rz0, vx0, vy0, vz0, dt0, iflag)
      !! author: David A. Minton
      !!
      !! Perform Kepler drift, solving Kepler's equation in appropriate variables
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_dan.f90
      !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_dan.f
      implicit none
      ! Arguments
      real(DP),     intent(in)    :: mu            !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
      real(DP),     intent(inout) :: rx0, ry0, rz0 !! position of body to drift
      real(DP),     intent(inout) :: vx0, vy0, vz0 !! velocity of body to drift     
      real(DP),     intent(in)    :: dt0           !! time step
      integer(I4B), intent(out)   :: iflag         !! error status flag for Kepler drift (0 = OK, nonzero = NO CONVERGENCE)
      ! Internals
      real(DP)        :: rx, ry, rz, vx, vy, vz, dt
      real(DP)        :: f, g, fdot, gdot, c1, c2, c3, u, alpha, fp, r0
      real(DP)        :: v0s, a, asq, en, dm, ec, es, esq, xkep, fchk, s, c

      ! Executable code
      iflag = 0
      dt = dt0
      r0 = sqrt(rx0*rx0 + ry0*ry0 + rz0*rz0)
      v0s = vx0*vx0 + vy0*vy0 + vz0*vz0
      u = rx0*vx0 + ry0*vy0 + rz0*vz0
      alpha = 2 * mu / r0 - v0s
      if (alpha > 0.0_DP) then
         a = mu / alpha
         asq = a**2
         en = sqrt(mu / (a * asq))
         ec = 1.0_DP - r0 / a
         es = u / (en * asq)
         esq = ec**2 + es**2
         dm = dt * en - int(dt * en / TWOPI, kind = I4B) * TWOPI
         dt = dm / en
         if ((esq < E2MAX) .and. (dm**2 < DM2MAX) .and. (esq * dm**2 < E2DM2MAX)) then
            call swiftest_drift_kepmd(dm, es, ec, xkep, s, c)
            fchk = (xkep - ec * s + es * (1.0_DP - c) - dm)
            ! DEK - original code compared fchk*fchk with DANBYB, but i think it should
            ! DEK - be compared with DANBYB*DANBYB, and i changed it accordingly - please
            ! DEK - check with hal and/or martin about this
            if (fchk**2 > DANBYB**2) then
               iflag = 1
               return
            end if
            fp = 1.0_DP - ec * c + es * s
            f = a / r0 * (c - 1.0_DP) + 1.0_DP
            g = dt + (s - xkep) / en
            fdot = -(a / (r0 * fp)) * en * s
            gdot = (c - 1.0_DP) / fp + 1.0_DP
            rx = rx0 * f + vx0 * g
            ry = ry0 * f + vy0 * g
            rz = rz0 * f + vz0 * g
            vx = rx0 * fdot + vx0 * gdot
            vy = ry0 * fdot + vy0 * gdot
            vz = rz0 * fdot + vz0 * gdot

            rx0 = rx
            ry0 = ry
            rz0 = rz
            vx0 = vx
            vy0 = vy
            vz0 = vz
            iflag = 0
            return
         end if
      end if

      call swiftest_drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      if (iflag == 0) then
         f = 1.0_DP - mu / r0 * c2
         g = dt - mu * c3
         fdot = -mu / (fp * r0) * c1
         gdot = 1.0_DP - mu / fp * c2
         rx = rx0 * f + vx0 * g
         ry = ry0 * f + vy0 * g
         rz = rz0 * f + vz0 * g
         vx = rx0 * fdot + vx0 * gdot
         vy = ry0 * fdot + vy0 * gdot
         vz = rz0 * fdot + vz0 * gdot

         rx0 = rx
         ry0 = ry
         rz0 = rz
         vx0 = vx
         vy0 = vy
         vz0 = vz
      end if

      return
   end subroutine swiftest_drift_dan


   pure subroutine swiftest_drift_kepmd(dm, es, ec, x, s, c)
      !! author: David A. Minton
      !!
      !! Solve Kepler's equation in difference form for an ellipse for small input dm and eccentricity
      !!    Original disclaimer: built for speed, does not check how well the original equation is solved
      !!    Can do that in calling routine by checking how close (x - ec*s + es*(1.0 - c) - dm) is to zero
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepmd.f90
      !! Adapted from Martin Duncan's Swift routine drift_kepmd.f
      implicit none
      ! Arguments
      real(DP), intent(in)  :: dm !! increment in mean anomaly
      real(DP), intent(in)  :: es !! eccentricity times the sine of eccentric anomaly
      real(DP), intent(in)  :: ec !! eccentricity times the cosine of eccentric anomaly
      real(DP), intent(out) :: x  !! solution to Kepler's equation in difference form (x = dE)
      real(DP), intent(out) :: s  !! sine of x
      real(DP), intent(out) :: c  !! cosine of x
      ! Internals
      real(DP), parameter :: a0 = 39916800.0_DP, a1 = 6652800.0_DP, a2 = 332640.0_DP, a3 = 7920.0_DP, a4 = 110.0_DP
      real(DP)      :: dx, fac1, fac2, q, y, f, fp, fpp, fppp
      
      ! executable code
      fac1 = 1.0_DP / (1.0_DP - ec)
      q = fac1 * dm
      fac2 = es*es*fac1 - ec / 3.0_DP
      x = q * (1.0_DP - 0.5_DP * fac1 * q * (es - q * fac2))
      y = x*x
      s = x * (a0 - y * (a1 - y * (a2 - y * (a3 - y * (a4 - y))))) / a0
      c = sqrt(1.0_DP - s*s)
      f = x - ec * s + es * (1.0_DP - c) - dm
      fp = 1.0_DP - ec * c + es * s
      fpp = ec * s + es * c
      fppp = ec * c - es * s
      dx = -f / fp
      dx = -f / (fp + dx * fpp/2.0_DP)
      dx = -f / (fp + dx * fpp/2.0_DP + dx*dx * fppp * SIXTH)
      x = x + dx
      y = x*x
      s = x * (a0 - y * (a1 - y * (a2 - y * (a3 - y * (a4 - y))))) / a0
      c = sqrt(1.0_DP - s*s)
      
      return
   end subroutine swiftest_drift_kepmd


   pure subroutine swiftest_drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
      !! author: David A. Minton
      !!
      !! Solve Kepler's equation in universal variables
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu.f90
      !! Adapted from Hal Levison's Swift routine drift_kepu.f
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(out)     :: fp, c1, c2, c3      
      real(DP) :: s, st, fo, fn
      
      ! executable code
      call swiftest_drift_kepu_guess(dt, r0, mu, alpha, u, s)
      st = s
      call swiftest_drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      if (iflag /= 0) then
         call swiftest_drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo)
         call swiftest_drift_kepu_fchk(dt, r0, mu, alpha, u, s, fn)
         if (abs(fo) < abs(fn)) s = st
         call swiftest_drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      end if
      
      return
   end subroutine swiftest_drift_kepu


   pure subroutine swiftest_drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
      !! author: David A. Minton
      !!
      !! Computes the value of f, the function whose root we are trying to find in universal variables
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_fchk.f90
      !! Adapted from Martin Duncan's Swift routine drift_kepu_fchk.f
      implicit none
      ! Internals
      real(DP), intent(in)  :: dt    !! time step
      real(DP), intent(in)  :: r0    !! distance between two bodies
      real(DP), intent(in)  :: mu    !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
      real(DP), intent(in)  :: alpha !! twice the binding energy
      real(DP), intent(in)  :: u     !! dot product of position and velocity vectors
      real(DP), intent(in)  :: s     !! universal variable (approximate root of f)
      real(DP), intent(out) :: f     !! function value
      ! Arguments
      real(DP) :: x, c0, c1, c2, c3

      x = s**2 * alpha
      call swiftest_drift_kepu_stumpff(x, c0, c1, c2, c3)
      c1 = c1 * s
      c2 = c2 * s**2
      c3 = c3 * s**3
      f = r0 * c1 + u * c2 + mu * c3 - dt

      return
   end subroutine swiftest_drift_kepu_fchk


   pure subroutine swiftest_drift_kepu_guess(dt, r0, mu, alpha, u, s)
      !! author: David A. Minton
      !!
      !! Compute initial guess for solving Kepler's equation using universal variables
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_guess.f90
      !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_kepu_guess.f
      implicit none
      ! Arguments
      real(DP), intent(in)  :: dt    !! time ste4p
      real(DP), intent(in)  :: r0    !! distance between two bodies
      real(DP), intent(in)  :: mu    !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
      real(DP), intent(in)  :: alpha !! twice the binding energy
      real(DP), intent(in)  :: u     !! dot product of position and velocity vectors
      real(DP), intent(out) :: s     !! initial guess for the value of the universal variable
      ! Internals
      integer(I4B) :: iflag
      real(DP), parameter :: thresh = 0.4_DP, danbyk = 0.85_DP
      real(DP)        :: y, sy, cy, sigma, es, x, a, en, ec, e

      if (alpha > 0.0_DP) then
         if (dt / r0 <= thresh) then
            s = dt/r0 - (dt*dt*u)/(2.0_DP*r0*r0*r0)
         else
            a = mu/alpha
            en = sqrt(mu/(a*a*a))
            ec = 1.0_DP - r0/a
            es = u/(en*a*a)
            e = sqrt(ec*ec + es*es)
            y = en*dt - es
            call swiftest_orbel_scget(y, sy, cy)
            sigma = sign(1.0_DP, es*cy + ec*sy)
            x = y + sigma*DANBYK*e
            s = x/sqrt(alpha)
         end if
      else
         call swiftest_drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
         if (iflag /= 0) s = dt / r0
      end if

      return
   end subroutine swiftest_drift_kepu_guess


   pure subroutine swiftest_drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      !! author: David A. Minton
      !!
      !! Solve Kepler's equation in universal variables using Laguerre's method
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 178 - 180.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_lag.f90
      !! Adapted from Hal Levison's Swift routine drift_kepu_lag.f
      implicit none
      ! Arguments
      real(DP),     intent(inout) :: s     !! universal variable 
      real(DP),     intent(in)    :: dt    !! time step
      real(DP),     intent(in)    :: r0    !! distance between two bodies 
      real(DP),     intent(in)    :: mu    !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift 
      real(DP),     intent(in)    :: alpha !! twice the binding energy 
      real(DP),     intent(in)    :: u     !! dot product of position and velocity vectors
      real(DP),     intent(out)   :: fp    !! first derivative of Kepler's equation in universal variables with respect to s (see Danby, p. 175)
      real(DP),     intent(out)   :: c1    !! Stumpff function c1 times s
      real(DP),     intent(out)   :: c2    !! Stumpff function c2 times s**2
      real(DP),     intent(out)   :: c3    !! Stumpff function c3 times s**3
      integer(I4B), intent(out)   :: iflag !! error status flag for convergence (0 = CONVERGED, nonzero = NOT CONVERGED)
      ! Internals
      integer(I4B) :: nc, ncmax
      real(DP)   :: x, fpp, ds, c0, f, fdt
      integer(I4B), parameter :: ln = 5

      if (alpha < 0.0_DP) then
         ncmax = NLAG2
      else
         ncmax = NLAG1
      end if
      do nc = 0, ncmax
         x = s * s * alpha
         call swiftest_drift_kepu_stumpff(x, c0, c1, c2, c3)
         c1 = c1*s
         c2 = c2*s*s
         c3 = c3*s*s*s
         f = r0 * c1 + u * c2 + mu * c3 - dt
         fp = r0 * c0 + u * c1 + mu * c2
         fpp = (-r0 * alpha + mu) * c1 + u * c0
         ds = -ln*f/(fp + sign(1.0_DP, fp)*sqrt(abs((ln - 1.0_DP)*(ln - 1.0_DP)*fp*fp - (ln - 1.0_DP)*ln*f*fpp)))
         s = s + ds
         fdt = f / dt
         if (fdt*fdt < DANBYB*DANBYB) then
            iflag = 0
            return
         end if
      end do
      iflag = 2

      return
   end subroutine swiftest_drift_kepu_lag


   pure subroutine swiftest_drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      !! author: David A. Minton
      !!
      !! Solve Kepler's equation in universal variables using Newton's method
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 174 - 175.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_new.f90
      !! Adapted from Hal Levison's Swift routine drift_kepu_new.f
      implicit none
      ! Arguments
      real(DP),     intent(inout) :: s     !! universal variable
      real(DP),     intent(in)    :: dt    !! time step
      real(DP),     intent(in)    :: r0    !! distance between two bodies
      real(DP),     intent(in)    :: mu    !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
      real(DP),     intent(in)    :: alpha !! twice the binding energy
      real(DP),     intent(in)    :: u     !! dot product of position and velocity vectors
      real(DP),     intent(out)   :: fp    !! first derivative of Kepler's equation in universal variables with respect to s (see Danby, p. 175)
      real(DP),     intent(out)   :: c1    !! Stumpff function c1 times s
      real(DP),     intent(out)   :: c2    !! Stumpff function c2 times s**2
      real(DP),     intent(out)   :: c3    !! Stumpff function c3 times s**3
      integer(I4B), intent(out)   :: iflag !! error status flag for convergence (0 = CONVERGED, nonzero = NOT CONVERGED)
      ! Internals
      integer( I4B) :: nc
      real(DP)   :: x, c0, ds, f, fpp, fppp, fdt

      do nc = 0, 6
         x = s*s*alpha
         call swiftest_drift_kepu_stumpff(x, c0, c1, c2, c3)
         c1 = c1*s
         c2 = c2*s*s
         c3 = c3*s*s*s
         f = r0*c1 + u*c2 + mu*c3 - dt
         fp = r0*c0 + u*c1 + mu*c2
         fpp = (-r0*alpha + mu)*c1 + u*c0
         fppp = (-r0*alpha + mu)*c0 - u*alpha*c1
         ds = -f/fp
         ds = -f/(fp + ds*fpp/2.0_DP)
         ds = -f/(fp + ds*fpp/2.0_DP + ds*ds*fppp/6.0_DP)
         s = s + ds
         fdt = f/dt
         if (fdt*fdt < DANBYB*DANBYB) then
              iflag = 0
              return
         end if
      end do
      iflag = 1

      return
   end subroutine swiftest_drift_kepu_new


   pure subroutine swiftest_drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
      !! author: David A. Minton
      !!
      !! Computes real root of cubic involved in setting initial guess for solving Kepler's equation in universal variables
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 177 - 178.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_p3solve.f90
      !! Adapted from Martin Duncan's Swift routine drift_kepu_p3solve.f
      implicit none
      ! Arguments
      real(DP), intent(in)      :: dt    !! time step
      real(DP), intent(in)      :: r0    !! distance between two bodies
      real(DP), intent(in)      :: mu    !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
      real(DP), intent(in)      :: alpha !! twice the binding energy
      real(DP), intent(in)      :: u     !! dot product of position and velocity vectors
      real(DP), intent(out)     :: s     !! s     : real solution of cubic equation
      integer(I4B), intent(out) :: iflag !! error status flag for solution (0 = OK, nonzero = ERROR)
      ! Internals
      real(DP) :: denom, a0, a1, a2, q, r, sq2, sq, p1, p2

      denom = (mu - alpha * r0) * SIXTH
      a2 = 0.5_DP * u / denom
      a1 = r0 / denom
      a0 = -dt / denom
      q = (a1 - a2*a2 * THIRD) * THIRD
      r = (a1*a2 - 3 * a0) * SIXTH - (a2*a2*a2)/27.0_DP
      sq2 = q*q*q + r*r
      if (sq2 >= 0.0_DP) then
         sq = sqrt(sq2)
         if ((r + sq) <= 0.0_DP) then
            p1 = -(-(r + sq))**(THIRD)
         else
            p1 = (r + sq)**(THIRD)
         end if
         if ((r - sq) <= 0.0_DP) then
            p2 = -(-(r - sq))**(THIRD)
         else
            p2 = (r - sq)**(THIRD)
         end if
         iflag = 0
         s = p1 + p2 - a2 * THIRD
      else
         iflag = 1
         s = 0.0_DP
      end if

      return
   end subroutine swiftest_drift_kepu_p3solve


   pure subroutine swiftest_drift_kepu_stumpff(x, c0, c1, c2, c3)
      !! author: David A. Minton
      !!
      !! Compute Stumpff functions needed for Kepler drift in universal variables
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 171 - 172.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_stumpff.f90
      !! Adapted from Hal Levison's Swift routine drift_kepu_stumpff.f
      implicit none
      ! Arguments
      real(DP), intent(inout) :: x  !! argument of Stumpff functions
      real(DP), intent(out)   :: c0 !! zeroth Stumpff function
      real(DP), intent(out)   :: c1 !! first Stumpff function
      real(DP), intent(out)   :: c2 !! second Stumpff function
      real(DP), intent(out)   :: c3 !! third Stumpff function
      ! Internals
      integer(I4B) :: i, n
      real(DP)   :: xm

      n = 0
      xm = 0.1_DP
      do while (abs(x) >= xm)
         n = n + 1
         x = x / 4.0_DP
      end do
      c2 = (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * &
            (1.0_DP - x / 182.0_DP) / 132.0_DP) / 90.0_DP) / 56.0_DP) /        &
            30.0_DP) / 12.0_DP) / 2.0_DP
      c3 = (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * (1.0_DP - x * &
            (1.0_DP - x / 210.0_DP) / 156.0_DP) / 110.0_DP) / 72.0_DP) /       &
            42.0_DP) / 20.0_DP ) / 6.0_DP
      c1 = 1.0_DP - x * c3
      c0 = 1.0_DP - x * c2
      if (n /= 0) then
         do i = n, 1, -1
            c3 = (c2 + c0 * c3) / 4.0_DP
            c2 = c1*c1 / 2.0_DP
            c1 = c0 * c1
            c0 = 2 * c0*c0 - 1.0_DP
            x = x * 4
         end do
      end if

      return
   end subroutine swiftest_drift_kepu_stumpff


end submodule s_swiftest_drift
