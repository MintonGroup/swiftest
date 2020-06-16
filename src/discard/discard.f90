submodule (swiftest_data_structures) s_discard
contains
   module procedure discard
   !! author: David A. Minton
   !!
   !! Check to see if test particles should be discarded based on their positions or because they are unbound from the systemj
   use swiftest
!subroutine discard(swiftest_pla, swiftes_tpa, config, t, dt)

! modules
!   use module_parameters
!   use module_swifter
!   use module_interfaces, except_this_one => discard
!   implicit none
!
!! arguments
!   logical, intent(in)  :: lclose, lrhill_present
!   integer(I4B), intent(in)  :: npl, ntp
!   real(DP), intent(in)    :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
!   character(*), intent(in)  :: qmin_coord
!   type(swifter_pl), pointer :: swifter_pl1p
!   type(swifter_tp), pointer :: swifter_tp1p
!
! internals
   real(DP) :: msys

! executable code
   if ((config%rmin  >= 0.0_DP) .or. &
       (config%rmax  >= 0.0_DP) .or. &
       (config%rmaxu >= 0.0_DP) .or. &
       ((config%qmin >= 0.0_DP) .and. (config%qmin_coord == "BARY"))) then
      call swiftest_plA%h2b()
      call swiftest_tpA%h2b(swifter_plA)
   end if
   if ((config%rmin >= 0.0_DP) .or. &
       (config%rmax >= 0.0_DP) .or. &
       (config%rmaxu >= 0.0_DP)) call discard_sun(swiftest_tpA, config, swiftest_plA%msys ntp, msys, swifter_tp1p, rmin, rmax,  &
      rmaxu)
   if (qmin >= 0.0_DP) call discard_peri(t, npl, ntp, swifter_pl1p, swifter_tp1p, msys, qmin, qmin_alo, qmin_ahi, qmin_coord,   &
      lrhill_present)
   if (lclose) call discard_pl(t, dt, npl, ntp, swifter_pl1p, swifter_tp1p)

   return

   end procedure discard
end submodule s_discard