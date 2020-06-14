module drift_vec
use swiftest_globals
implicit none

interface

   pure subroutine drift_vec_dan(mu, x0, v0, dt0, iflag)
      use swiftest_globals
      implicit none
      integer(I4B), intent(out)                :: iflag
      real(DP), intent(in)                     :: mu, dt0
      real(DP), dimension(:), intent(inout)    :: x0, v0
   end subroutine drift_vec_dan

   pure subroutine drift_vec_kepmd(dm, es, ec, x, s, c)
      use swiftest_globals
      implicit none
      real(DP), intent(in)  :: dm, es, ec
      real(DP), intent(out) :: x, s, c
   end subroutine drift_vec_kepmd

   pure subroutine drift_vec_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
      use swiftest_globals
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(out)     :: fp, c1, c2, c3
   end subroutine drift_vec_kepu

   pure subroutine drift_vec_kepu_fchk(dt, r0, mu, alpha, u, s, f)
      use swiftest_globals
      implicit none
      real(DP), intent(in)  :: dt, r0, mu, alpha, u, s
      real(DP), intent(out) :: f
   end subroutine drift_vec_kepu_fchk

   pure subroutine drift_vec_kepu_guess(dt, r0, mu, alpha, u, s)
      use swiftest_globals
      implicit none
      real(DP), intent(in)  :: dt, r0, mu, alpha, u
      real(DP), intent(out) :: s
   end subroutine drift_vec_kepu_guess

   pure subroutine drift_vec_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      use swiftest_globals
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(inout)   :: s
      real(DP), intent(out)     :: fp, c1, c2, c3
   end subroutine drift_vec_kepu_lag

   pure subroutine drift_vec_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
      use swiftest_globals
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(inout)   :: s
       real(DP), intent(out)     :: fp, c1, c2, c3
   end subroutine drift_vec_kepu_new

   pure subroutine drift_vec_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
      use swiftest_globals
      implicit none
      integer(I4B), intent(out) :: iflag
      real(DP), intent(in)      :: dt, r0, mu, alpha, u
      real(DP), intent(out)     :: s
   end subroutine drift_vec_kepu_p3solve

   pure subroutine drift_vec_kepu_stumpff(x, c0, c1, c2, c3)
      use swiftest_globals
      implicit none
      real(DP), intent(inout) :: x
      real(DP), intent(out)   :: c0, c1, c2, c3
   end subroutine drift_vec_kepu_stumpff

   module elemental subroutine drift_vec_one(mu, posx, posy, posz, vx, vy, vz, dt, iflag)
      use swiftest_globals
      implicit none
      real(DP), intent(in)      :: mu                !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
      real(DP), intent(inout)   :: posx, posy, posz  !! Position of body to drift
      real(DP), intent(inout)   :: vx, vy, vz        !! Velocity of body to drift
      real(DP), intent(in)      :: dt                !! Step size
      integer(I4B), intent(out) :: iflag             !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
   end subroutine drift_vec_one

end interface

end module drift_vec


module drift

   interface
      subroutine drift_dan(mu, x0, v0, dt0, iflag)
         use swiftest_globals
         implicit none
         integer(I4B), intent(out)                :: iflag
         real(DP), intent(in)                     :: mu, dt0
         real(DP), dimension(:), intent(inout)    :: x0, v0
      end subroutine drift_dan

      subroutine drift_kepmd(dm, es, ec, x, s, c)
         use swiftest_globals
         implicit none
         real(DP), intent(in)  :: dm, es, ec
         real(DP), intent(out) :: x, s, c
      end subroutine drift_kepmd

      subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
         use swiftest_globals
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(out)     :: fp, c1, c2, c3
      end subroutine drift_kepu

      subroutine drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
         use swiftest_globals
         implicit none
         real(DP), intent(in)  :: dt, r0, mu, alpha, u, s
         real(DP), intent(out) :: f
      end subroutine drift_kepu_fchk

      subroutine drift_kepu_guess(dt, r0, mu, alpha, u, s)
         use swiftest_globals
         implicit none
         real(DP), intent(in)  :: dt, r0, mu, alpha, u
         real(DP), intent(out) :: s
      end subroutine drift_kepu_guess

      subroutine drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
         use swiftest_globals
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(inout)   :: s
         real(DP), intent(out)     :: fp, c1, c2, c3
      end subroutine drift_kepu_lag

      subroutine drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
         use swiftest_globals
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(inout)   :: s
         real(DP), intent(out)     :: fp, c1, c2, c3
      end subroutine drift_kepu_new

      subroutine drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
         use swiftest_globals
         implicit none
         integer(I4B), intent(out) :: iflag
         real(DP), intent(in)      :: dt, r0, mu, alpha, u
         real(DP), intent(out)     :: s
      end subroutine drift_kepu_p3solve

      subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)
         use swiftest_globals
         implicit none
         real(DP), intent(inout) :: x
         real(DP), intent(out)   :: c0, c1, c2, c3
      end subroutine drift_kepu_stumpff

      subroutine drift_one(mu, x, v, dt, iflag)
         use swiftest_globals
         implicit none
         real(DP), intent(in)                    :: mu      !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
         real(DP), dimension(:), intent(inout)   :: x       !! Position of body to drift
         real(DP), dimension(:), intent(inout)   :: v       !! Velocity of body to drift
         real(DP), intent(in)                    :: dt      !! Step size
         integer(I4B), intent(out)               :: iflag   !! iflag : error status flag for Danby drift (0 = OK, nonzero = ERROR)
      end subroutine drift_one

   end interface

end module drift
