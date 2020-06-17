submodule (symba) s_symba_chk
contains
   module procedure symba_chk
   !! author: David A. Minton
   !!
   !! Check for an encounter
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_chk.f90
   !! Adapted from Hal Levison's Swift routine symba5_chk.f
use swiftest
implicit none
   integer(I4B) :: iflag
   real(DP)   :: rcrit, r2crit, vdotr

   lencounter = .false.
   rcrit = (rhill1 + rhill2)*rhscale*(rshell**(irec))
   r2crit = rcrit*rcrit
   call rmvs_chk_ind(xr(:), vr(:), dt, r2crit, iflag)
   if (iflag /= 0) lencounter = .true.
   vdotr = dot_product(vr(:), xr(:))
   lvdotr = (vdotr < 0.0_DP)

   return

   end procedure symba_chk
end submodule s_symba_chk
