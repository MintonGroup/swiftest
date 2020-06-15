submodule (helio) s_helio_getacch_pl
contains
   module procedure helio_getacch_pl   
   !! author: David A. Minton
   !!
   !! Compute heliocentric accelerations of plAnets
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_getacch.f90
   !! Adapted from Hal Levison's Swift routine helio_getacch.f
   use swiftest
   logical, save                    :: lmalloc = .true.
   integer(I4B)                     :: i, npl
   real(DP)                       :: r2
   real(DP), dimension(:), allocatable, save    :: irh
   real(DP), dimension(:, :), allocatable, save :: xh, aobl

   npl = helio_plA%nbody
   if (lflag) then
      do i = 2, npl
         helio_plA%ahi(:,i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      end do
      call helio_getacch_int_pl(helio_plA)
   end if
   if (config%j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(xh(NDIM, config%nplmax), aobl(NDIM, config%nplmax), irh(config%nplmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         xh(:, i) = helio_plA%xh(:,i)
         r2 = dot_product(xh(:, i), xh(:, i))
         irh(i) = 1.0_DP / sqrt(r2)
      end do
      call obl_acc(npl, helio_plA, config%j2rp2, config%j4rp4, xh, irh, aobl)
      do i = 2, npl
         helio_plA%ah(:,i) = helio_plA%ahi(:,i) + aobl(:, i) - aobl(:, 1)
      end do
   else
      do i = 2, npl
         helio_plA%ah(:,i) = helio_plA%ahi(:,i)
      end do
   end if
   if (config%lextra_force) call helio_user_getacch(helio_plA, t)

   return

   end procedure helio_getacch_pl
end submodule s_helio_getacch_pl