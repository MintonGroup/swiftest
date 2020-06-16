module util
   implicit none
   interface

      subroutine util_exit(code)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in) :: code
      end subroutine util_exit

      subroutine util_dist_index_plpl(npl, nplm, num_comparisons, k_plpl)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: npl, nplm
         integer(I4B), dimension(:,:),allocatable,intent(out) :: k_plpl
         integer(I4B), intent(out) :: num_comparisons
      end subroutine

      subroutine util_dist_index_pltp(nplm, ntp, num_comparisons, k_pltp)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: nplm, ntp
         integer(I4B), dimension(:,:),allocatable,intent(out) :: k_pltp
         integer(I4B), intent(out) :: num_comparisons
      end subroutine util_dist_index_pltp

      subroutine util_dist_eucl_plpl(npl, invar, num_comparisons, k_plpl, outvar)
         use swiftest_globals
         use swiftest_data_structures
         use symba
         implicit none
         integer(I4B), intent(in)  :: npl
         integer(I4B), dimension(2,num_comparisons),intent(in) :: k_plpl
         integer(I4B), intent(in) :: num_comparisons
         real(DP),dimension(ndim,npl),intent(in) :: invar
         real(DP), dimension(ndim,num_comparisons),intent(inout) :: outvar
      end subroutine util_dist_eucl_plpl

      subroutine util_dist_eucl_pltp(npl, ntp, massive bodies, test_particles, num_pltp_comparisons, k_pltp, outvar)
         use swiftest_globals
         use swiftest_data_structures
         use symba
         implicit none
         integer(I4B), intent(in) :: npl, ntp
         integer(I4B), dimension(num_pltp_comparisons,2),intent(in) :: k_pltp
         integer(I4B), intent(in) :: num_pltp_comparisons
         real(DP),dimension(ndim,npl),intent(in) :: massive bodies
         real(DP),dimension(ndim,ntp),intent(in) :: test_particles
         real(DP), dimension(ndim,num_pltp_comparisons),intent(inout) :: outvar
      end subroutine

      subroutine util_hills(npl, swiftest_plA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine util_hills

      subroutine util_index(arr, index)
         use swiftest_globals
         use module_nrutil
         implicit none
         integer(I4B), dimension(:), intent(out) :: index
         real(DP), dimension(:), intent(in)   :: arr
      end subroutine util_index

      subroutine util_resize_pl(symba_plA, npl_new, npl_old)
         use swiftest_globals
         use symba
         use swiftest_data_structures
         use helio
         use module_nrutil
         implicit none
         type(symba_pl), intent(inout) :: symba_plA
         integer(I4B), intent(in)   :: npl_old, npl_new
      end subroutine util_resize_pl

      subroutine util_sort_I4B(arr)
         use swiftest_globals
         implicit none
         integer(I4B), dimension(:), intent(inout) :: arr
      end subroutine util_sort_I4B
      subroutine util_sort_sp(arr)
         use swiftest_globals
         implicit none
         real(sp), dimension(:), intent(inout) :: arr
      end subroutine util_sort_sp
      subroutine util_sort_DP(arr)
         use swiftest_globals
         implicit none
         real(DP), dimension(:), intent(inout) :: arr
      end subroutine util_sort_DP

      subroutine util_toupper(string)
         use swiftest_globals
         implicit none
         character(*), intent(inout) :: string
      end subroutine util_toupper

      subroutine util_valid(npl, ntp, swiftest_plA, swiftest_tpA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl, ntp
         type(swiftest_pl), intent(in) :: swiftest_plA
         type(swiftest_tp), intent(in) :: swiftest_tpA
      end subroutine util_valid

      subroutine util_version
         use swiftest_globals
         implicit none
      end subroutine util_version

      function util_kahan_sum(xsum_current, xi, xerror) 
         use swiftest_globals
         implicit none
         real(DP)     :: util_kahan_sum
         real(DP), intent(in)   :: xsum_current, xi
         real(DP), intent(inout) :: xerror
      end function

   end interface
end module util
