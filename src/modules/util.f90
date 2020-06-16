module util
   implicit none
   interface

      module subroutine util_exit(code)
         implicit none
         integer(I4B), intent(in) :: code
      end subroutine util_exit

     module subroutine util_dist_index_plpl(npl, nplm, num_comparisons, k_plpl)
         implicit none
         integer(I4B), intent(in)  :: npl, nplm
         integer(I4B), dimension(:,:),allocatable,intent(out) :: k_plpl
         integer(I4B), intent(out) :: num_comparisons
      end subroutine

      module subroutine util_dist_index_pltp(nplm, ntp, num_comparisons, k_pltp)
         implicit none
         integer(I4B), intent(in)  :: nplm, ntp
         integer(I4B), dimension(:,:),allocatable,intent(out) :: k_pltp
         integer(I4B), intent(out) :: num_comparisons
      end subroutine util_dist_index_pltp

      module subroutine util_dist_eucl_plpl(npl, invar, num_comparisons, k_plpl, outvar)
         implicit none
         integer(I4B), intent(in)  :: npl
         integer(I4B), dimension(:,:),intent(in) :: k_plpl
         integer(I4B), intent(in) :: num_comparisons
         real(DP),dimension(:,:),intent(in) :: invar
         real(DP), dimension(:,:),intent(inout) :: outvar
      end subroutine util_dist_eucl_plpl

      module subroutine util_dist_eucl_pltp(npl, ntp, planets, test_particles, num_pltp_comparisons, k_pltp, outvar)
         implicit none
         integer(I4B), intent(in) :: npl, ntp
         integer(I4B), dimension(:,:),intent(in) :: k_pltp
         integer(I4B), intent(in) :: num_pltp_comparisons
         real(DP), dimension(:,:),intent(in) :: planets
         real(DP), dimension(:,:),intent(in) :: test_particles
         real(DP), dimension(:, :),intent(inout) :: outvar
      end subroutine

      module subroutine util_hills(npl, swiftest_plA)
         implicit none
         integer(I4B), intent(in)    :: npl
         class(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine util_hills

      module subroutine util_index(arr, index)
         implicit none
         integer(I4B), dimension(:), intent(out) :: index
         real(DP), dimension(:), intent(in)   :: arr
      end subroutine util_index

      module subroutine util_peri(lfirst, ntp, swiftest_tpa, mu, msys, qmin_coord)
         logical, intent(in)               :: lfirst       !! Logical flag indicating whether current invocation is the first
         integer(I4B), intent(in)          :: ntp          !! Number of active test particles
         class(swiftest_tp), intent(inout) :: swiftest_tpA !! Swiftest test particle class
         real(DP), intent(in)              :: mu           !! G * (m1 + m2) = mass of the Sun in this routine
         real(DP), intent(in)              :: msys         !! Total system masse
         character(len=*), intent(in)      :: qmin_coord   !! Coordinate frame for qmin (see swiftest_globals for symbolic definitions)
      end subroutine util_peri

      module subroutine util_sort_i4b(arr)
         implicit none
         integer(I4B), dimension(:), intent(inout) :: arr
      end subroutine util_sort_i4b

      module subroutine util_sort_sp(arr)
         implicit none
         real(SP), dimension(:), intent(inout) :: arr
      end subroutine util_sort_sp

      module subroutine util_sort_dp(arr)
         implicit none
         real(DP), dimension(:), intent(inout) :: arr
      end subroutine util_sort_dp

      module subroutine util_toupper(string)
         implicit none
         character(*), intent(inout) :: string
      end subroutine util_toupper

     module subroutine util_valid(swiftest_plA, swiftest_tpA)
         implicit none
         class(swiftest_pl), intent(in) :: swiftest_plA
         class(swiftest_tp), intent(in) :: swiftest_tpA
      end subroutine util_valid

      module subroutine util_version
         implicit none
      end subroutine util_version

      module function util_kahan_sum(xsum_current, xi, xerror) 
         implicit none
         real(DP)     :: util_kahan_sum
         real(DP), intent(in)   :: xsum_current, xi
         real(DP), intent(inout) :: xerror
      end function

   end interface
end module util
