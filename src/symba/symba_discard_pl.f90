submodule (symba) s_symba_discard_pl
contains
   module procedure symba_discard_pl
   !! author: David A. Minton
   !!
   !! Check to see if planets should be discarded based on their positions or because they are unbound
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_discard_pl.f90
   !! Adapted from Hal Levison's Swift routine discard_massive5.f
use swiftest
implicit none
   logical(lgt)          :: ldiscards
   integer(I4B)          :: i
   real(DP)            :: msys, ke, pe, tei, tef
   real(DP), dimension(ndim) :: htot

! executable code
   ldiscards = .false.
   if ((rmin >= 0.0_DP) .or. (rmax >= 0.0_DP) .or. (config%rmaxu >= 0.0_DP) .or. ((qmin >= 0.0_DP) .and. (qmin_coord == "bary")))    &
      call coord_h2b(npl, symba_plA, msys)
   if ((rmin >= 0.0_DP) .or. (rmax >= 0.0_DP) .or. (config%rmaxu >= 0.0_DP))                                     &
      call symba_discard_sun_pl(t, npl, msys, symba_plA, rmin, rmax, config%rmaxu, ldiscards)
   if (qmin >= 0.0_DP) call symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)

   return

   end procedure symba_discard_pl
end submodule s_symba_discard_pl
