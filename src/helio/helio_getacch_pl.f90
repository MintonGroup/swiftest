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
   real(DP), dimension(:, :), allocatable, save :: xh_loc, aobl

   npl = self%nbody
   if (lflag) then
      self%ahi(:,2:npl) = 0.0_DP
      call helio_getacch_int_pl(self)
   end if
   if (config%j2rp2 /= 0.0_DP) then
      if (lmalloc) then
         allocate(xh_loc(NDIM, config%nplmax), aobl(NDIM, config%nplmax), irh(config%nplmax))
         lmalloc = .false.
      end if
      do i = 2, npl
         xh_loc(:, i) = self%xh(:,i)
         r2 = dot_product(xh_loc(:, i), xh_loc(:, i))
         irh(i) = 1.0_DP / sqrt(r2)
      end do
      call obl_acc(self, config%j2rp2, config%j4rp4, xh_loc, irh, aobl)
      self%ah(:,2:npl) = self%ahi(:,2:npl) + aobl(:, 2:npl) - aobl(:, 1)
   else
      self%ah(:,:) = self%ahi(:,:)
   end if
   if (config%lextra_force) call helio_user_getacch(self, t)

   return

   end procedure helio_getacch_pl
end submodule s_helio_getacch_pl