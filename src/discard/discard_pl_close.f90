submodule (swiftest_data_structures) s_discard_pl_close
contains
   module procedure discard_pl_close
   !! author: David A. Minton
   !!
   !!  Check to see if a test particle and massive body are having, or will have within the next time step, an encounter such
   !!          that the separation distance r is less than some critical radius rcrit (or r**2 < rcrit**2 = r2crit)
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: discard_pl_close.f90
   !! Adapted from Hal Levison's Swift routine discard_pl_close.f
   use swiftest
   real(DP) :: r2, v2, vdotr, tmin

! executable code
   r2 = dot_product(dx(:), dx(:))
   if (r2 <= r2crit) then
      iflag = 1
   else
      vdotr = dot_product(dx(:), dv(:))
      if (vdotr > 0.0_DP) then
         iflag = 0
      else
         v2 = dot_product(dv(:), dv(:))
         tmin = -vdotr / v2
         if (tmin < dt) then
            r2min = r2 - vdotr * vdotr / v2
         else
            r2min = r2 + 2 * vdotr * dt + v2 * dt * dt
         end if
         r2min = min(r2min, r2)
         if (r2min <= r2crit) then
            iflag = 1
         else
            iflag = 0
         end if
      end if
   end if

   return

   end procedure discard_pl_close
end submodule s_discard_pl_close
