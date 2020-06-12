submodule (helio) s_helio_drift_pl
contains
   module procedure helio_drift_pl
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Read in parameters for the integration
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_init_param.
   implicit none

   integer(I4B)                   :: i, iflag, i
   real(DP)                       :: mu

   iflag = helio_plA%drift_one(dt) 
   do i = 2, npl
      call drift_one(mu, swiftest_plA%xh(:,i), swiftest_plA%vb(:,i), dt, iflag)
      if (iflag /= 0) then
          write(*, *) " plAnet ", swiftest_plA%name(i), " is lost!!!!!!!!!!"
          write(*, *) mu, dt
          write(*, *) swiftest_plA%xh(:,i)
          write(*, *) swiftest_plA%vb(:,i)
          write(*, *) " stopping "
          call util_exit(failure)
      end if
   end do

   return

   end procedure helio_drift_pl
end submodule s_helio_drift_pl
