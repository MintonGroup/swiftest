submodule (symba) s_symba_chk_eucl_pltp
contains
   module procedure symba_chk_eucl_pltp
   !! author: Jacob R. Elliott
   !!
   !! Check for an encounter between planets and test particles using the single-loop blocking for the Euclidean distance 
   !!    matrix evaluation
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_chk.f90
   !! Adapted from Hal Levison's Swift routine symba5_chk.f 
use swiftest
implicit none
   ! logical(lgt) :: iflag lvdotr_flag
   real(DP)   :: rcrit, r2crit, vdotr, r2, v2, tmin, r2min, term2, rcritmax, r2critmax
   integer(I4B) :: k
   real(DP), dimension(ndim) :: xr, vr

! executable code

   npltpenc = 0
   
   term2 = rhscale*(rshell**0)

   rcritmax = symba_pla%helio%swiftest%rhill(2) * term2
   r2critmax = rcritmax * rcritmax

!$omp parallel do default(none) schedule(static) &
!$omp private(k, rcrit, r2crit, r2, vdotr, v2, tmin, r2min, xr, vr) &
!$omp shared(num_encounters, lvdotr, lencounter, k_pltp, dt, term2, r2critmax, symba_pla, symba_tpa) &
!$omp reduction(+:npltpenc)

   do k = 1,num_encounters
      xr(:) = symba_tpa%helio%swiftest%xh(:,k_pltp(2,k)) - symba_pla%helio%swiftest%xh(:,k_pltp(1,k))
      r2 = dot_product(xr(:), xr(:)) 
      if (r2<r2critmax) then

         rcrit = symba_pla%helio%swiftest%rhill(k_pltp(1,k))*term2
         r2crit = rcrit*rcrit 
         vr(:) = symba_tpa%helio%swiftest%vh(:,k_pltp(2,k)) - symba_pla%helio%swiftest%vh(:,k_pltp(1,k))

         vdotr = dot_product(vr(:), xr(:))

         if (vdotr < 0.0_DP) lvdotr(k) = k

         if (r2 < r2crit) then
            lencounter(k) = k
            npltpenc = npltpenc + 1
         else
            if (vdotr < 0.0_DP) then
               v2 = dot_product(vr(:), vr(:))
               tmin = -vdotr/v2
               if (tmin < dt) then
                  r2min = r2 - vdotr*vdotr/v2
               else
                  r2min = r2 + 2.0_DP*vdotr*dt + v2*dt*dt
               end if
               r2min = min(r2min, r2)
               if (r2min <= r2crit)then
                  lencounter(k) = k
                  npltpenc = npltpenc + 1
               endif
            end if
         end if
      endif
   enddo

!$omp end parallel do

   return

   end procedure symba_chk_eucl_pltp
end submodule s_symba_chk_eucl_pltp
