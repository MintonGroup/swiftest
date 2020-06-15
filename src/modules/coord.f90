module coord
   implicit none

   interface 
      subroutine coord_h2b(npl, swiftest_plA, msys)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl
         real(DP), intent(out)     :: msys
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine coord_h2b

      subroutine coord_h2b_tp(ntp, swiftest_tpA, swiftest_plA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: ntp
         type(swiftest_tp), intent(inout) :: swiftest_tpA
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine coord_h2b_tp

      subroutine coord_vb2vh(npl, swiftest_plA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine coord_vb2vh

      subroutine coord_vb2vh_tp(ntp, swiftest_tpA, vs)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)     :: ntp
         real(DP), dimension(ndim), intent(in) :: vs
         type(swiftest_tp), intent(inout)   :: swiftest_tpA
      end subroutine coord_vb2vh_tp

      subroutine coord_vh2vb(npl, swiftest_plA, msys)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)   :: npl
         real(DP), intent(out)    :: msys
         type(swiftest_pl),intent(inout) :: swiftest_plA
      end subroutine coord_vh2vb

      subroutine coord_vh2vb_tp(ntp, swiftest_tpA, vs)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)     :: ntp
         real(DP), dimension(ndim), intent(in) :: vs
         type(swiftest_tp), intent(inout)   :: swiftest_tpA
      end subroutine coord_vh2vb_tp
   end interface
end module coord
