submodule (symba_classes) s_symba_encounter_check
   use swiftest
contains
   module function symba_encounter_check_pl(self, system, dt, irec) result(lencounter)
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: self       !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      integer(I4B),              intent(in)    :: irec       !! Current recursion level
      ! Result
      logical                                  :: lencounter !! Returns true if there is at least one close encounter      
      ! Internals
      integer(I4B)               :: i, j
      real(DP)                   :: r2, v2, vdotr
      real(DP), dimension(NDIM)  :: xr, vr
      real(DP)                   :: rcrit, r2crit
      logical                    :: lflag

      associate(pl => self, npl => self%nbody)
         lencounter = .false.
         do concurrent(j = 1:npl, pl%lminty(j))
            do i = 1, npl
               rcrit = (pl%rhill(i) + pl%rhill(j)) * RHSCALE * (RSHELL**(irec))
               r2crit = r2crit**2
               xr(:) = pl%xh(:, j) - pl%xh(:, i)
               vr(:) = pl%vh(:, j) - pl%vh(:, i)
               v2 = dot_product(vr(:), vr(:))
               vdotr = dot_product(vr(:), xr(:))
               lflag = rmvs_chk_ind(r2, v2, vdotr, dt, r2crit)
               if (lflag) lencounter = .true.
            end do
         end do
      end associate
      return
   end function symba_encounter_check_pl

   module function symba_encounter_check_tp(self, system, dt) result(lencounter)
      implicit none
      ! Arguments
      class(symba_tp),           intent(inout) :: self       !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      ! Result
      logical                                  :: lencounter !! Returns true if there is at least one close encounter      

      lencounter = .false.
      return
   end function symba_encounter_check_tp

end submodule s_symba_encounter_check