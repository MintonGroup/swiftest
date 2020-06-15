module module_interfaces
  implicit none





      subroutine python_io_write_frame_pl(t, symba_plA, npl, out_stat)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         real(DP), intent(in)   :: t
         type(symba_pl),intent(in) :: symba_plA
         integer, intent(in)   :: npl
         character(*), intent(in)  :: out_stat
      end subroutine python_io_write_frame_pl

      subroutine python_io_write_frame_tp(t, symba_tpA, ntp, out_stat)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         real(DP), intent(in)   :: t
         type(symba_tp),intent(in) :: symba_tpA
         integer, intent(in)   :: ntp
         character(*), intent(in)  :: out_stat
      end subroutine python_io_write_frame_tp




   end interface
end module module_interfaces
