module pydriver_module
   use swiftest, only : swiftest_driver
   implicit none

   contains 
   subroutine driver(integrator, param_file_name, display_style)

      !! author: David A. Minton
      !!
      !! Driver program for the Swiftest integrators. Unlike the earlier Swift and Swifter drivers, in Swiftest all integrators 
      !!    are run from this single program. 
      !!
      !! Adapted from Swifter by David E. Kaufmann's Swifter driver programs swifter_[bs,helio,ra15,rmvs,symba,tu4,whm].f90
      !! Adapted from Hal Levison and Martin Duncan's Swift driver programs
      implicit none

      ! Arguments
      character(len=:), intent(in), allocatable :: integrator      !! Symbolic code of the requested integrator  
      character(len=:), intent(in), allocatable :: param_file_name !! Name of the input parameters file
      character(len=:), intent(in), allocatable :: display_style   !! Style of the output display {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD")

      call swiftest_driver(integrator, param_file_name, display_style)

      return
   end subroutine driver
end module pydriver_module