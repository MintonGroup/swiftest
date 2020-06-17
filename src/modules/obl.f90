module obl
   implicit none
   interface

      subroutine obl_acc(swiftest_plA, j2rp2, j4rp4, xh, irh, aobl)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         class(swiftest_pl), intent(in)     :: swiftest_plA
         real(DP), intent(in)        :: j2rp2, j4rp4
         real(DP), dimension(:), intent(in)   :: irh
         real(DP), dimension(:, :), intent(in)  :: xh
         real(DP), dimension(:, :), intent(out) :: aobl
      end subroutine obl_acc

      subroutine obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)      :: ntp
         real(DP), intent(in)        :: j2rp2, j4rp4, msun
         real(DP), dimension(:), intent(in)   :: irht
         real(DP), dimension(:, :), intent(in)  :: xht
         real(DP), dimension(:, :), intent(out) :: aoblt
      end subroutine obl_acc_tp

      subroutine obl_pot(swiftest_plA, j2rp2, j4rp4, xh, irh, oblpot)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         class(swiftest_pl), intent(in)    :: swiftest_plA
         real(DP), intent(in)       :: j2rp2, j4rp4
         real(DP), intent(out)        :: oblpot
         real(DP), dimension(:), intent(in)   :: irh
         real(DP), dimension(:, :), intent(in) :: xh
      end subroutine obl_pot

   end interface
end module obl
