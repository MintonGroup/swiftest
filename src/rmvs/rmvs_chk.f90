submodule (rmvs_classes) s_rmvs_chk
contains
   module procedure rmvs_chk
   !! author: David A. Minton
   !!
   !! Determine whether a test particle and planet are having or will have an encounter within the next time step
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk.f90
   !! Adapted from Hal Levison's Swift routine rmvs3_chk.f
   use swiftest
   implicit none
   ! Internals
   integer(I4B)           :: i, j, k, iflag, nenc
   real(DP)              :: r2crit
   real(DP), dimension(NDIM) :: xht, vht, xr, vr

! executable code
   associate(ntp => tp%nbody, npl => pl%nbody)
      lencounter = .false.
      do i = 1, npl
          pl%nenc(i) = 0
          plp%tpenc1p(i)) = 0
      end do
      do i = 1, ntp
         if (tp%status == ACTIVE) then
            iflag = 0
            xht(:) = tp%xh(:, i)
            vht(:) = tp%vh(:, i)
            do j = 1, npl
               r2crit = (rts * plrhill(j))**2
               xr(:) = xht(:) - xh(:, j)
               vr(:) = vht(:) - vh(:, j)
               call rmvs_chk_ind(xr(:), vr(:), dt, r2crit, iflag)
               if (iflag /= 0) then
                     lencounter = .true.
                     pl%nenc(j) = pl%nenc(j) + 1
                     nenc = pl%nenc(j)
                     if (nenc == 1) then
                        pl%tpenc1p(j) = i
                     else
                        rmvs_tpencp => rmvs_plp%tpenc1p
                        do k = 2, nenc - 1
                           rmvs_tpencp => rmvs_tpencp%tpencp
                        end do
                        rmvs_tpencp%tpencp => rmvs_tpp
                     end if
                     rmvs_tpp%plencp => rmvs_plp
                     exit
               end if
            end do
         end if
         rmvs_tpp => rmvs_tpp%nextp
      end do

   end associate
   return


   end procedure rmvs_chk

   module procedure rmvs_chk_ind
   !! author: david a. minton
   !!
   !! determine whether a test particle and planet are having or will have an encounter within the next time step
   !!
   !! adapted from david e. kaufmann's swifter routine: rmvs_chk_ind.f90
   !! adapted from hal levison's swift routine rmvs_chk_ind.f
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
end submodule s_rmvs_chk
