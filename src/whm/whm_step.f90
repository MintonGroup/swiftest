submodule(whm_classes) s_whm_step
contains
   module procedure whm_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, j2rp2, j4rp4, dt, c2)
   !! author: David A. Minton
   !!
   !! Step planets and active test particles ahead in heliocentric coordinates
   !!
   !! Adapted from Hal Levison's Swift routine step_kdk.f
   !! Adapted from David E. Kaufmann's Swifter routine whm_step.f90
   use swiftest
   implicit none
   logical(lgt)                     :: lfirsttp
   logical(lgt), save                 :: lmalloc = .true.
   integer(I4B)                     :: i
   real(DP), dimension(:, :), allocatable, save :: xbeg, xend
   type(whm_pl), pointer                :: whm_plp

! executable code
   lfirsttp = lfirst
   if (ntp > 0) then
      if (lmalloc) then
         allocate(xbeg(ndim, nplmax), xend(ndim, nplmax))
         lmalloc = .false.
      end if
      whm_plp => whm_pl1p
      do i = 2, npl
         whm_plp => whm_plp%nextp
         xbeg(:, i) = whm_plp%swifter%xh(:)
      end do
   end if
   call whm_step_pl(lfirst, lextra_force, t, npl, nplmax, whm_pl1p, j2rp2, j4rp4, dt, c2)
   if (ntp > 0) then
      whm_plp => whm_pl1p
      do i = 2, npl
         whm_plp => whm_plp%nextp
         xend(:, i) = whm_plp%swifter%xh(:)
      end do
      call whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1p, whm_tp1p, xbeg, xend, j2rp2, j4rp4, dt,c2)
   end if

   return

   end procedure whm_step
end submodule s_whm_step
