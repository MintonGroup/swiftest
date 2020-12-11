submodule (rmvs_classes) s_rmvs_chk_ind
contains
   module procedure rmvs_chk_ind
   !! author: David A. Minton
   !!
   !! Determine whether a test particle and planet are having or will have an encounter within the next time step
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk_ind.f90
   !! Adapted from Hal Levison's Swift routine rmvs_chk_ind.f
   use swiftest
   implicit none
   real(DP) :: r2, v2, vdotr, tmin, r2min

   iflag = 0
   r2 = dot_product(xr(:), xr(:))
   if (r2 < r2crit) then
      iflag = 1
   else
      vdotr = dot_product(vr(:), xr(:))
      if (vdotr < 0.0_DP) then
         v2 = dot_product(vr(:), vr(:))
         tmin = -vdotr / v2
         if (tmin < dt) then
            r2min = r2 - vdotr**2 / v2
         else
            r2min = r2 + 2 * vdotr * dt + v2 * dt**2
         end if
         r2min = min(r2min, r2)
         if (r2min <= r2crit) iflag = 1
      end if
   end if

   return

   end procedure rmvs_chk_ind
end submodule s_rmvs_chk_ind
