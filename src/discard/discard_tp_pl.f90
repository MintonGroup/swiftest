submodule (swiftest_classes) s_discard_tp_pl
contains
   module procedure discard_tp_pl
   !! author: David A. Minton
   !!
   !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: discard_pl.f90
   !! Adapted from Hal Levison's Swift routine discard_pl.f
   use swiftest
   implicit none

   integer(I4B)              :: i, j, isp, ntp, npl
   real(DP)                  :: r2min, radius
   real(DP), dimension(NDIM) :: dx, dv

   ntp = swiftest_tpA%nbody
   npl = swiftest_plA%nbody
   do i = 1, ntp
      if (swiftest_tpA%status(i) == ACTIVE) then
         do j = 2, npl
            dx(:) = swiftest_tpA%xh(:,i) - swiftest_plA%xh(:,i)
            dv(:) = swiftest_tpA%vh(:,i) - swiftest_plA%vh(:,i)
            radius = swiftest_plA%radius(i)
            call discard_pl_close(dx(:), dv(:), dt, radius * radius, isp, r2min)
            if (isp /= 0) then
               swiftest_tpA%status(i) = DISCARDED_PLR
               swiftest_plA%ldiscard = .true.
               write(*, *) "Particle ", swiftest_tpA%name(i), " too close to massive body ", swiftest_plA%name(i), " at t = ", t
               swiftest_tpA%ldiscard = .true.
               exit
            end if
         end do
      end if
   end do

   return

   end procedure discard_tp_pl
end submodule s_discard_tp_pl
