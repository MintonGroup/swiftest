module orbel
   implicit none
   interface
      pure subroutine orbel_scget(angle, sx, cx)
         use swiftest_globals
         implicit none
         real(DP), intent(in)  :: angle
         real(DP), intent(out) :: sx, cx
      end subroutine orbel_scget

      subroutine orbel_xv2aeq(x, v, mu, a, e, q)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: mu
         real(DP), dimension(ndim), intent(in) :: x, v
         real(DP), intent(out)      :: a, e, q
      end subroutine orbel_xv2aeq

      subroutine orbel_xv2aqt(x, v, mu, a, q, capm, tperi)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: mu
         real(DP), dimension(ndim), intent(in) :: x, v
         real(DP), intent(out)      :: a, q, capm, tperi
      end subroutine orbel_xv2aqt

      subroutine orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: mu
         real(DP), dimension(ndim), intent(in) :: x, v
         real(DP), intent(out)      :: a, e, inc, capom, omega, capm
      end subroutine orbel_xv2el

   end interface
end module orbel
