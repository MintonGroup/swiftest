module obl
   implicit none
   interface

      subroutine obl_acc(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, aobl)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)      :: npl
         real(DP), intent(in)        :: j2rp2, j4rp4
         real(DP), dimension(npl), intent(in)   :: irh
         real(DP), dimension(ndim, npl), intent(in)  :: xh
         real(DP), dimension(ndim, npl), intent(out) :: aobl
         type(swiftest_pl), intent(inout)     :: swiftest_plA
      end subroutine obl_acc

      subroutine obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)      :: ntp
         real(DP), intent(in)        :: j2rp2, j4rp4, msun
         real(DP), dimension(ntp), intent(in)   :: irht
         real(DP), dimension(ndim, ntp), intent(in)  :: xht
         real(DP), dimension(ndim, ntp), intent(out) :: aoblt
      end subroutine obl_acc_tp

      subroutine obl_pot(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, oblpot)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)       :: npl
         real(DP), intent(in)       :: j2rp2, j4rp4
         real(DP), intent(out)        :: oblpot
         real(DP), dimension(npl), intent(in)   :: irh
         real(DP), dimension(ndim, npl), intent(in) :: xh
         type(swiftest_pl), intent(inout)    :: swiftest_plA
      end subroutine obl_pot

   end interface
end module obl
