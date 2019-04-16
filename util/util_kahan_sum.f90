  FUNCTION util_kahan_sum(xsum_current, xi, xerror) RESULT(xsum_new)
  !-------------------------------------------------------------------------------------------
  !					UTIL_KAHAN_SUM.F90
  !-------------------------------------------------------------------------------------------
  ! Sums two floating point scalars more accurately utilitizing the Kahan summation formula
  ! This function is designed to be used inside a do or while loop, where the initial value of
  ! of xsum_current is initialized appropriately and the initial value of xerror is 0.0_dp
  !
  ! N.B. Use this function if the summation is being performed for more than *three* terms
  !
  ! Input:  xsum_current - Current value of the sum
  !         xi           - i-th term to be added to the sum
  !         xerror       - Error term from the previous term of the sum
  !
  ! Output: xsum_new     - The updated value of the sum
  !         xerror       - The error term for this term of the sum
  !
  ! By: Chris Capobianco
  ! Date: 05/04/09
  !       01/26/10 D. Minton - added to SWIFTER and modified types
  !-------------------------------------------------------------------------------------------
! Modules
   USE module_parameters
   USE module_interfaces, EXCEPT_THIS_ONE => util_kahan_sum
  implicit none

  ! Input/Output variables
  real(dp), intent(in) :: xsum_current, xi
  real(dp), intent(inout) :: xerror
  real(dp) :: xsum_new

  ! Internal variables
  real(dp) :: low_bits

  low_bits = xi - xerror
  xsum_new = xsum_current + low_bits
  xerror = (xsum_new - xsum_current) - low_bits

  return
  end function util_kahan_sum
