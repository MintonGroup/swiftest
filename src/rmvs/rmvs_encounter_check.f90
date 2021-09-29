submodule (rmvs_classes) s_rmvs_chk
   use swiftest
contains

   module function rmvs_encounter_check_tp(self, param, system, dt) result(lencounter)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk.f90
      !! Adapted from Hal Levison's Swift routine rmvs3_chk.f
      implicit none
      ! Arguments
      class(rmvs_tp),             intent(inout) :: self   !! RMVS test particle object  
      class(swiftest_parameters), intent(in)    :: param  !! Current swiftest run configuration parameters
      class(rmvs_nbody_system),   intent(inout) :: system !! RMVS nbody system object
      real(DP),                   intent(in)    :: dt     !! step size
      ! Result
      logical                                 :: lencounter  !! Returns true if there is at least one close encounter
      ! Internals
      integer(I4B)                            :: i, j, nenc
      real(DP)                                :: xr, yr, zr, vxr, vyr, vzr
      real(DP), dimension(system%pl%nbody)    :: rcrit
      logical                                 :: lflag
      logical, dimension(:),      allocatable :: lvdotr
      integer(I4B), dimension(:), allocatable :: index1, index2

      ! The minimization and linear solvers can sometimes lead to floating point exceptions. Rather than halting the code entirely if this occurs, we
      ! can simply fail the attempt and try again. So we need to turn off any floating point exception halting modes temporarily 

      lencounter = .false.
      if (self%nbody == 0) return

      select type(pl => system%pl)
      class is (rmvs_pl)
         associate(tp => self, ntp => self%nbody, npl => pl%nbody)
            tp%plencP(1:ntp) = 0
            call encounter_check_all_sort_and_sweep_pltp(npl, ntp, pl%xbeg, pl%vbeg, tp%xh, tp%vh, pl%renc, dt, lvdotr, index1, index2, nenc)

            lencounter = (nenc > 0)
            if (lencounter) then
               tp%plencP(index2(1:nenc)) = index1(1:nenc)
               do j = 1, npl
                  pl%nenc(j) = count(tp%plencP(1:ntp) == j)
               end do
            end if
         end associate
      end select

      return
   end function rmvs_encounter_check_tp


end submodule s_rmvs_chk
