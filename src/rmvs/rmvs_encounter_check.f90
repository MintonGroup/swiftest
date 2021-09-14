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
      real(DP)                                :: xr, yr, zr, vxr, vyr, vzr
      real(DP), dimension(system%pl%nbody)    :: r2crit
      logical                                 :: lflag, lvdotr

      ! The minimization and linear solvers can sometimes lead to floating point exceptions. Rather than halting the code entirely if this occurs, we
      ! can simply fail the attempt and try again. So we need to turn off any floating point exception halting modes temporarily 

      lencounter = .false.
      if (self%nbody == 0) return

      select type(pl => system%pl)
      class is (rmvs_pl)
         associate(tp => self, ntp => self%nbody, npl => pl%nbody, rts => system%rts)
            r2crit(1:npl) = (rts * pl%rhill(1:npl))**2
            tp%plencP(1:ntp) = 0
            !$omp parallel do default(private)&
            !$omp shared(npl, ntp, tp, pl, dt, r2crit)
            do j = 1, npl
               do i = 1, ntp
                  if ((.not.tp%lmask(i)).or.(tp%plencP(i) /= 0)) cycle
                  xr = tp%xh(1, i) - pl%xbeg(1, j)
                  yr = tp%xh(2, i) - pl%xbeg(2, j)
                  zr = tp%xh(3, i) - pl%xbeg(3, j)
                  vxr = tp%vh(1, i) - pl%vbeg(1, j)
                  vyr = tp%vh(2, i) - pl%vbeg(2, j)
                  vzr = tp%vh(3, i) - pl%vbeg(3, j)
                  call rmvs_chk_ind(xr, yr, zr, vxr, vyr, vzr, dt, r2crit(j), lflag, lvdotr)
                  if (lflag) tp%plencP(i) = j
               end do
               pl%nenc(j) = count(tp%plencP(1:ntp) == j)
            end do
            !$omp end parallel do
            lencounter = any(pl%nenc(1:npl) > 0)
         end associate
      end select

      return
   end function rmvs_encounter_check_tp


   module pure subroutine rmvs_chk_ind(xr, yr, zr, vxr, vyr, vzr, dt, r2crit, lencounter, lvdotr)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk_ind.f90
      !! Adapted from Hal Levison's Swift routine rmvs_chk_ind.f
      implicit none
      ! Arguments
      real(DP), intent(in)  :: xr, yr, zr    !! Relative distance vector components
      real(DP), intent(in)  :: vxr, vyr, vzr !! Relative velocity vector components
      real(DP), intent(in)  :: dt            !! Step size
      real(DP), intent(in)  :: r2crit        !! Square of the critical encounter distance
      logical,  intent(out) :: lencounter    !! Flag indicating that an encounter has occurred
      logical,  intent(out) :: lvdotr        !! Logical flag indicating the direction of the v .dot. r vector
      ! Internals
      real(DP) :: r2min, r2, v2, vdotr

      r2 = xr**2 + yr**2 + zr**2
      lencounter = (r2 < r2crit) 
      if (lencounter) return 

      vdotr = vxr * xr + vyr * yr + vzr * zr
      lvdotr = (vdotr < 0.0_DP)
      if (.not.lvdotr) return
     
      v2 = vxr**2 + vyr**2 + vzr**2

      if (-vdotr < v2 * dt) then
         r2min = r2 - vdotr**2 / v2
      else
         r2min = r2 + 2 * vdotr * dt + v2 * dt**2
      end if
      lencounter = (r2min <= r2crit)

      return
   end subroutine rmvs_chk_ind

end submodule s_rmvs_chk
