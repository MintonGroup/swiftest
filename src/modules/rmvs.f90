module rmvs
   implicit none
   interface

      subroutine rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: dt, r2crit
         real(DP), dimension(ndim), intent(in) :: xr, vr
         integer(I4B), intent(out)      :: iflag
      end subroutine rmvs_chk_ind

   end interface
end module rmvs
