submodule (symba) s_symba_discard_sun_pl
contains
   module procedure symba_discard_sun_pl
   !! author: David A. Minton
   !!
   !! Check to see if planets should be discarded based on their positions relative to the Sun
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_discard_sun_pl.f90
   !! Adapted from Hal Levison's Swift routine discard_massive5.f
use swiftest
implicit none
   integer(I4B)          :: i
   real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2


! executable code
   rmin2 = rmin*rmin
   rmax2 = rmax*rmax
   rmaxu2 = rmaxu*rmaxu
   do i = 2, npl
      if (swiftest_pla%status(i) == active) then
         rh2 = dot_product(swiftest_pla%xh(:,i), swiftest_pla%xh(:,i))
         if ((rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
            ldiscards = .true.
            swiftest_pla%status(i) = discarded_rmax
            write(*, *) "particle ",  swiftest_pla%name(i), " too far from sun at t = ", t
            print *,'rmax: ',rmax
            print *,'rh2: ',rh2
         else if ((rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
            ldiscards = .true.
            swiftest_pla%status(i) = discarded_rmin
            write(*, *) "particle ", swiftest_pla%name(i), " too close to sun at t = ", t
         else if (rmaxu >= 0.0_DP) then
            rb2 = dot_product(swiftest_pla%xb(:,i), swiftest_pla%xb(:,i))
            vb2 = dot_product(swiftest_pla%vb(:,i), swiftest_pla%vb(:,i))
            energy = 0.5_DP*vb2 - msys/sqrt(rb2)
            if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
               ldiscards = .true.
               swiftest_pla%status(i) = discarded_rmaxu
               write(*, *) "particle ", swiftest_pla%name(i), " is unbound and too far from barycenter at t = ", t
            end if
         end if
      end if
   end do

   return

   end procedure symba_discard_sun_pl
end submodule s_symba_discard_sun_pl
