submodule (swiftest_classes) drift_implementation
   use swiftest
   !> Integration control parameters:
   real(DP), parameter :: E2MAX    = 0.36_DP      
   real(DP), parameter :: DM2MAX   = 0.16_DP
   real(DP), parameter :: E2DM2MAX = 0.0016_DP
   real(DP),     parameter :: DANBYB   = 1.0E-13_DP
   integer(I2B), parameter :: NLAG1    = 50
   integer(I2B), parameter :: NLAG2    = 40

contains

   module subroutine drift_body(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Loop bodies and call Danby drift routine on the heliocentric position and velocities.
      !!
      !! Adapted from Hal Levison's Swift routine drift_tp.f 
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f90
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! Swiftest test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B)                              :: i   
      integer(I4B), dimension(:), allocatable   :: iflag

      associate(n => self%nbody)
         allocate(iflag(n))
         iflag(:) = 0
         call drift_all(self%mu, self%xh, self%vh, self%nbody, param, dt, self%lmask, iflag)
         if (any(iflag(1:n) /= 0)) then
            where(iflag(1:n) /= 0) self%status(1:n) = DISCARDED_DRIFTERR
            do i = 1, n
               if (iflag(i) /= 0) write(*, *) " Body ", self%id(i), " lost due to error in Danby drift"
            end do
         end if
      end associate

      return
   end subroutine drift_body


   module subroutine drift_all(mu, x, v, n, param, dt, lmask, iflag)
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
         do concurrent(i = 1:n, lmask(i))
            rmag = norm2(x(:, i))
            vmag2 = dot_product(v(:, i), v(:, i))
            energy = 0.5_DP * vmag2 - mu(i) / rmag
            dtp(i) = dt * (1.0_DP + 3 * param%inv_c2 * energy)
         end do
      else
         where(lmask(1:n)) dtp(1:n) = dt
      end if 
      !$omp simd
      do i = 1, n
         if (lmask(i)) call drift_one(mu(i), x(1,i), x(2,i), x(3,i), v(1,i), v(2,i), v(3,i), dtp(i), iflag(i))
      end do
      !$omp end simd

      return
   end subroutine drift_all


   module pure elemental subroutine drift_one(mu, px, py, pz, vx, vy, vz, dt, iflag) 
      !$omp declare simd(drift_one)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Perform Danby drift for one body, redoing drift with smaller substeps if original accuracy is insufficient
      !!
      !! Adapted from David E. Kaufmann's Swifter routine routine drift_one.f90
      !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_one.f
      implicit none
      ! Arguments
      real(DP), intent(in)      :: mu                     !! G * (Mcb + m), G = gravitational constant, Mcb = mass of central body, m = mass of body to drift
      real(DP), intent(inout)   :: px, py, pz, vx, vy, vz !! Position and velocity of body to drift
      real(DP), intent(in)      :: dt                     !! Step size
      integer(I4B), intent(out) :: iflag                  !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      ! Internals
      integer(I4B) :: i
      real(DP)   :: dttmp

      call drift_dan(mu, px, py, pz, vx, vy, vz, dt, iflag)
      if (iflag /= 0) then
         dttmp = 0.1_DP * dt
         do i = 1, 10
            call drift_dan(mu, px, py, pz, vx, vy, vz, dttmp, iflag)
            if (iflag /= 0) exit
         end do
      end if
   
      return
   end subroutine drift_one


   pure subroutine drift_dan(mu, px0, py0, pz0, vx0, vy0, vz0, dt0, iflag)
      !$omp declare simd(drift_dan)
      !! author: David A. Minton
      !!
      !! Perform Kepler drift, solving Kepler's equation in appropriate variables
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_dan.f90
      !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_dan.f
      implicit none
      integer(I4B), intent(out)                :: iflag
      real(DP), intent(in)                     :: mu, dt0
      real(DP), intent(inout)    :: px0, py0, pz0, vx0, vy0, vz0      
      real(DP)        :: dt, f, g, fdot, gdot, c1, c2, c3, u, alpha, fp, r0
      real(DP)        :: v0s, a, asq, en, dm, ec, es, esq, xkep, fchk, s, c
      real(DP), dimension(NDIM) :: x, v, x0, v0

      ! Executable code
      iflag = 0
      dt = dt0
      x0 = [px0, py0, pz0]
      v0 = [vx0, vy0, vz0]
      r0 = sqrt(dot_product(x0(:), x0(:))) 
      v0s = dot_product(v0(:), v0(:))
      u = dot_product(x0(:),  v0(:))
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
            call drift_kepmd(dm, es, ec, xkep, s, c)
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
            x(:) = x0(:) * f + v0(:) * g
            v(:) = x0(:) * fdot + v0(:) * gdot
            px0 = x(1); py0 = x(2); pz0 = x(3)
            vx0 = v(1); vy0 = v(2); vz0 = v(3)
            iflag = 0
            return
         end if
      end if

      call drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      if (iflag == 0) then
         f = 1.0_DP - mu / r0 * c2
         g = dt - mu * c3
         fdot = -mu / (fp * r0) * c1
         gdot = 1.0_DP - mu / fp * c2
         x(:) = x0(:) * f + v0(:) * g
         v(:) = x0(:) * fdot + v0(:) * gdot
         px0 = x(1); py0 = x(2); pz0 = x(3)
         vx0 = v(1); vy0 = v(2); vz0 = v(3)
      end if

      return
   end subroutine drift_dan


   pure subroutine drift_kepmd(dm, es, ec, x, s, c)
      !$omp declare simd(drift_kepmd)
      !! author: David A. Minton
      !!
      !! Solve Kepler's equation in difference form for an ellipse for small input dm and eccentricity
      !!    Original disclaimer: built for speed, does not check how well the original equation is solved
      !!    Can do that in calling routine by checking how close (x - ec*s + es*(1.0 - c) - dm) is to zero
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepmd.f90
      !! Adapted from Martin Duncan's Swift routine drift_kepmd.f
      implicit none
      real(DP), intent(in)  :: dm, es, ec
      real(DP), intent(out) :: x, s, c      
      real(DP), parameter :: a0 = 39916800.0_DP, a1 = 6652800.0_DP, a2 = 332640.0_DP, a3 = 7920.0_DP, a4 = 110.0_DP
      real(DP)      :: dx, fac1, fac2, q, y, f, fp, fpp, fppp
      
      ! executable code
      fac1 = 1.0_DP / (1.0_DP - ec)
      q = fac1 * dm
      fac2 = es**2 * fac1 - ec / 3.0_DP
      x = q * (1.0_DP - 0.5_DP * fac1 * q * (es - q * fac2))
      y = x**2
      s = x * (a0 - y * (a1 - y * (a2 - y * (a3 - y * (a4 - y))))) / a0
      c = sqrt(1.0_DP - s**2)
      f = x - ec * s + es * (1.0_DP - c) - dm
      fp = 1.0_DP - ec * c + es * s
      fpp = ec * s + es * c
      fppp = ec * c - es * s
      dx = -f / fp
      dx = -f / (fp + dx * fpp / 2.0_DP)
      dx = -f / (fp + dx * fpp / 2.0_DP + dx**2* fppp / 6.0_DP)
      x = x + dx
      y = x**2
      s = x * (a0 - y * (a1 - y * (a2 - y * (a3 - y * (a4 - y))))) / a0
      c = sqrt(1.0_DP - s**2)
      
      return
   end subroutine drift_kepmd


   pure subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
      !$omp declare simd(drift_kepu)
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
      call drift_kepu_guess(dt, r0, mu, alpha, u, s)
      st = s
      call drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      if (iflag /= 0) then
         call drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo)
         call drift_kepu_fchk(dt, r0, mu, alpha, u, s, fn)
         if (abs(fo) < abs(fn)) s = st
         call drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      end if
      
      return
   end subroutine drift_kepu


   pure subroutine drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
      !$omp declare simd(drift_kepu_fchk)
      !! author: David A. Minton
      !!
      !! Computes the value of f, the function whose root we are trying to find in universal variables
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_fchk.f90
      !! Adapted from Martin Duncan's Swift routine drift_kepu_fchk.f
      implicit none
      real(DP), intent(in)  :: dt, r0, mu, alpha, u, s
      real(DP), intent(out) :: f
      real(DP) :: x, c0, c1, c2, c3

      x = s**2 * alpha
      call drift_kepu_stumpff(x, c0, c1, c2, c3)
      c1 = c1 * s
      c2 = c2 * s**2
      c3 = c3 * s**3
      f = r0 * c1 + u * c2 + mu * c3 - dt

      return
   end subroutine drift_kepu_fchk


   pure subroutine drift_kepu_guess(dt, r0, mu, alpha, u, s)
      !$omp declare simd(drift_kepu_guess)
      !! author: David A. Minton
      !!
      !! Compute initial guess for solving Kepler's equation using universal variables
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_guess.f90
      !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_kepu_guess.f
      implicit none
      real(DP), intent(in)  :: dt, r0, mu, alpha, u
      real(DP), intent(out) :: s      
      integer(I4B)      :: iflag
      real(DP), parameter :: thresh = 0.4_DP, danbyk = 0.85_DP
      real(DP)        :: y, sy, cy, sigma, es, x, a, en, ec, e

      if (alpha > 0.0_DP) then
         if (dt / r0 <= thresh) then
            s = dt / r0 - (dt**2 * u) / (2 * r0**3)
         else
            a = mu / alpha
            en = sqrt(mu / a**3)
            ec = 1.0_DP - r0 / a
            es = u / (en * a**2)
            e = sqrt(ec**2 + es**2)
            y = en * dt - es
            call orbel_scget(y, sy, cy)
            sigma = sign(1.0_DP, es * cy + ec * sy)
            x = y + sigma * danbyk * e
            s = x / sqrt(alpha)
         end if
      else
         call drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
         if (iflag /= 0) s = dt / r0
      end if

      return
   end subroutine drift_kepu_guess


   pure subroutine drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      !$omp declare simd(drift_kepu_lag)
      !! author: David A. Minton
      !!
      !! Solve Kepler's equation in universal variables using Laguerre's method
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 178 - 180.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_lag.f90
      !! Adapted from Hal Levison's Swift routine drift_kepu_lag.f
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(inout)   :: s
      real(DP), intent(out)     :: fp, c1, c2, c3      
      integer( I4B) :: nc, ncmax
      real(DP)   :: ln, x, fpp, ds, c0, f, fdt
   
      if (alpha < 0.0_DP) then
         ncmax = NLAG2
      else
         ncmax = NLAG1
      end if
      ln = 5.0_DP
      do nc = 0, ncmax
         x = s * s * alpha
         call drift_kepu_stumpff(x, c0, c1, c2, c3)
         c1 = c1 * s
         c2 = c2 * s**2
         c3 = c3 * s**3
         f = r0 * c1 + u * c2 + mu * c3 - dt
         fp = r0 * c0 + u * c1 + mu * c2
         fpp = (-r0 * alpha + mu) * c1 + u * c0
         ds = -ln * f / (fp + sign(1.0_DP, fp) * sqrt(abs((ln - 1.0_DP)**2 * fp**2 - (ln - 1.0_DP) * ln * f * fpp)))
         s = s + ds
         fdt = f / dt
         if (fdt**2 < DANBYB**2) then
            iflag = 0
            return
         end if
      end do
      iflag = 2
   
      return
   end subroutine drift_kepu_lag


   pure subroutine drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      !$omp declare simd(drift_kepu_new)
      !! author: David A. Minton
      !!
      !! Solve Kepler's equation in universal variables using Newton's method
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 174 - 175.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_new.f90
      !! Adapted from Hal Levison's Swift routine drift_kepu_new.f
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(inout)   :: s
      real(DP), intent(out)     :: fp, c1, c2, c3      
      integer( I4B) :: nc
      real(DP)   :: x, c0, ds, f, fpp, fppp, fdt
   
      do nc = 0, 6
         x = s**2 * alpha
         call drift_kepu_stumpff(x, c0, c1, c2, c3)
         c1 = c1 * s
         c2 = c2 * s**2
         c3 = c3 * s**3
         f = r0 * c1 + u * c2 + mu * c3 - dt
         fp = r0 * c0 + u * c1 + mu * c2
         fpp = (-r0 * alpha + mu) * c1 + u * c0
         fppp = (-r0 * alpha + mu) * c0 - u * alpha * c1
         ds = -f / fp
         ds = -f / (fp + ds * fpp / 2.0_DP)
         ds = -f / (fp + ds * fpp / 2.0_DP + ds**2 * fppp / 6.0_DP)
         s = s + ds
         fdt = f / dt
         if (fdt**2 < DANBYB**2) then
            iflag = 0
            return
         end if
      end do
      iflag = 1
   
      return
   end subroutine drift_kepu_new


   pure subroutine drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
      !$omp declare simd(drift_kepu_p3solve)
      !! author: David A. Minton
      !!
      !! Computes real root of cubic involved in setting initial guess for solving Kepler's equation in universal variables
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 177 - 178.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_p3solve.f90
      !! Adapted from Martin Duncan's Swift routine drift_kepu_p3solve.f
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(out)     :: s
      real(DP) :: denom, a0, a1, a2, q, r, sq2, sq, p1, p2

      denom = (mu - alpha * r0) / 6.0_DP
      a2 = 0.5_DP * u / denom
      a1 = r0 / denom
      a0 = -dt / denom
      q = (a1 - a2**2 / 3.0_DP) / 3.0_DP
      r = (a1 * a2 - 3 * a0) / 6.0_DP - a2**3 / 27.0_DP
      sq2 = q**3 + r**2
      if (sq2 >= 0.0_DP) then
         sq = sqrt(sq2)
         if ((r + sq) <= 0.0_DP) then
            p1 = -(-(r + sq))**(1.0_DP / 3.0_DP)
         else
            p1 = (r + sq)**(1.0_DP / 3.0_DP)
         end if
         if ((r - sq) <= 0.0_DP) then
            p2 = -(-(r - sq))**(1.0_DP / 3.0_DP)
         else
            p2 = (r - sq)**(1.0_DP / 3.0_DP)
         end if
         iflag = 0
         s = p1 + p2 - a2 / 3.0_DP
      else
         iflag = 1
         s = 0.0_DP
      end if

      return
   end subroutine drift_kepu_p3solve
   

   pure subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)
      !$omp declare simd(drift_kepu_stumpff)
      !! author: David A. Minton
      !!
      !! Compute Stumpff functions needed for Kepler drift in universal variables
      !!      Reference: Danby, J. M. A. 1988. Fundamentals of Celestial Mechanics, (Willmann-Bell, Inc.), 171 - 172.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: drift_kepu_stumpff.f90
      !! Adapted from Hal Levison's Swift routine drift_kepu_stumpff.f
      implicit none
      real(DP), intent(inout) :: x
      real(DP), intent(out)   :: c0, c1, c2, c3
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
            c2 = c1**2 / 2.0_DP
            c1 = c0 * c1
            c0 = 2 * c0**2 - 1.0_DP
            x = x * 4
         end do
      end if

      return
   end subroutine drift_kepu_stumpff


end submodule drift_implementation
