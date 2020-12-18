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
      integer(I4B)              :: i, j, k, iflag, nenc
      real(DP)                  :: r2crit
      real(DP), dimension(NDIM) :: xht, vht, xr, vr
      integer(I4B)              :: tpencPindex

      associate(ntp => tp%nbody, npl => pl%nbody)
         lencounter = .false.
         pl%nenc(:) = 0
         pl%tpenc1p(:) = 0
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               tp%plencP(i) = 0
               tp%tpencP(i) = 0
               iflag = 0
               xht(:) = tp%xh(:, i)
               vht(:) = tp%vh(:, i)
               do j = 1, npl
                  r2crit = (rts * pl%rhill(j))**2
                  xr(:) = xht(:) - xh(:, j)
                  vr(:) = vht(:) - vh(:, j)
                  iflag = rmvs_chk_ind(xr(:), vr(:), dt, r2crit)
                  if (iflag /= 0) then
                        lencounter = .true.
                        pl%nenc(j) = pl%nenc(j) + 1
                        nenc = pl%nenc(j)
                        if (nenc == 1) then
                           pl%tpenc1P(j) = i
                        else
                           tpencPindex = pl%tpenc1P(j)
                           do k = 2, nenc - 1
                              tpencPindex = tp%tpencP(tpencPindex)
                           end do
                           tp%tpencP(tpencPindex) = i
                        end if
                        tp%plencP(i) = j
                        exit
                  end if
               end do
            end if
         end do
      end associate
      return
   end procedure rmvs_chk

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
end submodule s_rmvs_chk
