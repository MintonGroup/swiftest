submodule (rmvs_classes) s_rmvs_chk
   use swiftest
contains

   module function rmvs_encounter_check_tp(self, system, dt) result(lencounter)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk.f90
      !! Adapted from Hal Levison's Swift routine rmvs3_chk.f
      implicit none
      ! Arguments
      class(rmvs_tp),           intent(inout) :: self   !! RMVS test particle object  
      class(rmvs_nbody_system), intent(inout) :: system !! RMVS nbody system object
      real(DP),                 intent(in)    :: dt     !! step size
      ! Result
      logical                                 :: lencounter  !! Returns true if there is at least one close encounter
      ! Internals
      integer(I4B)                            :: i, j
      real(DP)                                :: r2, v2, vdotr
      real(DP), dimension(NDIM)               :: xr, vr
      real(DP), dimension(system%pl%nbody)    :: r2crit
      logical                                 :: lflag

      select type(pl => system%pl)
      class is (rmvs_pl)
         associate(tp => self, ntp => self%nbody, npl => pl%nbody, rts => system%rts)
            r2crit(:) = (rts * pl%rhill(:))**2
            tp%plencP(:) = 0
            do j = 1, npl
               do i = 1, ntp
                  if ((.not.tp%lmask(i)).or.(tp%plencP(i) /= 0)) cycle
                  xr(:) = tp%xh(:, i) - pl%xbeg(:, j)
                  vr(:) = tp%vh(:, i) - pl%vbeg(:, j)
                  r2 = dot_product(xr(:), xr(:))
                  v2 = dot_product(vr(:), vr(:))
                  vdotr = dot_product(vr(:), xr(:))
                  lflag = rmvs_chk_ind(r2, v2, vdotr, dt, r2crit(j))
                  if (lflag) tp%plencP(i) = j
               end do
               pl%nenc(j) = count(tp%plencP(:) == j)
            end do
            lencounter = any(pl%nenc(:) > 0)
         end associate
      end select
      return
   end function rmvs_encounter_check_tp


   module elemental function rmvs_chk_ind(r2, v2, vdotr, dt, r2crit) result(lflag)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk_ind.f90
      !! Adapted from Hal Levison's Swift routine rmvs_chk_ind.f
      implicit none
      ! Arguments
      real(DP), intent(in)       :: r2, v2, vdotr, dt, r2crit
      logical                    :: lflag
      ! Internals
      real(DP) ::  tmin, r2min

      lflag = .false.
      if (r2 < r2crit) then
         lflag = .true.
      else
         if (vdotr < 0.0_DP) then
            tmin = -vdotr / v2
            if (tmin < dt) then
               r2min = r2 - vdotr**2 / v2
            else
               r2min = r2 + 2 * vdotr * dt + v2 * dt**2
            end if
            r2min = min(r2min, r2)
            lflag = (r2min <= r2crit)
         end if
      end if

      return
   end function rmvs_chk_ind

end submodule s_rmvs_chk
