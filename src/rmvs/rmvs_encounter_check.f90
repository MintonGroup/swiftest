submodule (rmvs_classes) s_rmvs_chk
contains
   module function rmvs_encounter_check_tp(self, cb, pl, dt, rts) result(lencounter)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk.f90
      !! Adapted from Hal Levison's Swift routine rmvs3_chk.f
      use swiftest
      implicit none
      ! Arguments
      class(rmvs_tp),            intent(inout) :: self        !! RMVS test particle object  
      class(rmvs_cb),            intent(inout) :: cb          !! RMVS central body object  
      class(rmvs_pl),            intent(inout) :: pl          !! RMVS massive body object  
      real(DP),                  intent(in)    :: dt          !! step size
      real(DP),                  intent(in)    :: rts         !! fraction of Hill's sphere radius to use as radius of encounter regio
      logical                                  :: lencounter  !! Returns true if there is at least one close encounter
      ! Internals
      integer(I4B)                             :: i, j, k, nenc
      real(DP)                                 :: r2crit
      real(DP), dimension(NDIM)                :: xht, vht, xr, vr
      integer(I4B)                             :: tpencPindex
      logical                                  :: lflag
      logical, save                            :: lfirst = .true.

      associate(tp => self, ntp => self%nbody, npl => pl%nbody, rhill => pl%rhill)
         if (.not.allocated(pl%encmask)) allocate(pl%encmask(ntp, npl))
         pl%encmask(:,:) = .false.
         lencounter = .false.
         pl%nenc(:) = 0
         pl%tpenc1P(:) = 0
         ! if first time through, calc hill's sphere for the planets
         !if (lfirst) then
           
         !   lfirst = .false.
         !end if
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               tp%plencP(i) = 0
               tp%tpencP(i) = 0
               lflag = .false. 
               xht(:) = tp%xh(:, i)
               vht(:) = tp%vh(:, i)

               do j = 1, npl
                  r2crit = (rts * pl%rhill(j))**2
                  xr(:) = xht(:) - tp%xbeg(:, j)
                  vr(:) = vht(:) - tp%vbeg(:, j)
                  lflag = rmvs_chk_ind(xr(:), vr(:), dt, r2crit)
                  if (lflag) then
                     lencounter = .true.
                     pl%encmask(i,j) = .true.
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
   end function rmvs_encounter_check_tp

   module function rmvs_chk_ind(xr, vr, dt, r2crit) result(lflag)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk_ind.f90
      !! Adapted from Hal Levison's Swift routine rmvs_chk_ind.f
      use swiftest
      implicit none
      ! Arguments
      real(DP), intent(in)                     :: dt, r2crit
      real(DP), dimension(:), intent(in)       :: xr, vr
      logical                                  :: lflag
      ! Internals
      real(DP) :: r2, v2, vdotr, tmin, r2min

      lflag = .false.
      r2 = dot_product(xr(:), xr(:))
      if (r2 < r2crit) then
         lflag = .true.
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
            if (r2min <= r2crit) lflag = .true.
         end if
      end if

      return

   end function rmvs_chk_ind
end submodule s_rmvs_chk
