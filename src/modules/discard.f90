module discard
   implicit none
   interface


      subroutine discard_tp(t, dt, npl, ntp, swiftest_plA, swiftest_tpA, rmin, rmax, rmaxu, qmin,  &
         qmin_alo, qmin_ahi, qmin_coord, lclose, lrhill_present)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         logical(lgt), intent(in)  :: lclose, lrhill_present
         integer(I4B), intent(in)  :: npl, ntp
         real(DP), intent(in)   :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
         character(*), intent(in)  :: qmin_coord
         type(swiftest_pl), intent(inout) :: swiftest_plA
         type(swiftest_tp), intent(inout) :: swiftest_tpA
      end subroutine discard_tp

      subroutine discard_peri(t, npl, ntp, swiftest_plA, swiftest_tpA, msys, qmin, qmin_alo, & 
         qmin_ahi, qmin_coord, lrhill_present)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         logical(lgt), intent(in)  :: lrhill_present
         integer(I4B), intent(in)  :: npl, ntp
         real(DP), intent(in)   :: t, msys, qmin, qmin_alo, qmin_ahi
         character(*), intent(in)  :: qmin_coord
         type(swiftest_pl), intent(inout) :: swiftest_plA
         type(swiftest_tp), intent(inout) :: swiftest_tpA
      end subroutine discard_peri

      subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
         use swiftest_globals
         implicit none
         integer(I4B), intent(out)      :: iflag
         real(DP), intent(in)      :: dt, r2crit
         real(DP), dimension(ndim), intent(in) :: dx, dv
         real(DP), intent(out)      :: r2min
      end subroutine discard_pl_close

      subroutine discard_pl(t, dt, npl, ntp, swiftest_plA, swiftest_tpA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: npl, ntp
         real(DP), intent(in)   :: t, dt
         type(swiftest_pl), intent(inout) :: swiftest_plA
         type(swiftest_tp), intent(inout) :: swiftest_tpA
      end subroutine discard_pl

      subroutine discard_sun(t, ntp, msys, swifter_tpA, rmin, rmax, rmaxu)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: ntp
         real(DP), intent(in)   :: t, msys, rmin, rmax, rmaxu
         type(swiftest_tp), intent(inout) :: swifter_tpA
      end subroutine discard_sun

   end interface
end module discard
