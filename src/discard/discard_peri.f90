submodule (swiftest_data_structures) s_discard_peri
contains
   module procedure discard_peri
   !! author: David A. Minton
   !!
   !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: discard_peri.f90
   !! Adapted from Hal Levison's Swift routine discard_peri.f
   use swiftest
   logical, save             :: lfirst = .true.
   integer(I4B)              :: i, j, ih
   real(DP)                  :: r2
   real(DP), dimension(NDIM) :: dx


   if (lfirst) then
      if (.not. config%lrhill_present) call util_hills(npl, swiftest_plA)
      call util_peri(lfirst, swiftest_tpA%nbody, swiftest_tpA, swiftest_plA%mass(1), swiftest_plA%msys, config%qmin_coord)
      lfirst = .false.
   else
      call util_peri(swiftest_tpA, swiftest_plA, config, lfirst)
      do i = 1, ntp
         if (swiftest_tpA%status(i) == ACTIVE) then
            if (swiftest_tpA%isperi(i) == 0) then
               ih = 1
               do j = 2, npl
                  dx(:) = swiftest_tpA%xh(:,i) - swiftest_plA%xh(:,j)
                  r2 = dot_product(dx(:), dx(:))
                  if (r2 <= swiftest_plA%rhill(j) * swiftest_plA%rhill(j)) ih = 0
               end do
               if (ih == 1) then
                  if ((swiftest_tpA%atp(i) >= qmin_alo) .and.    &
                     (swiftest_tpA%atp(i) <= qmin_ahi) .and.    &           
                      (swiftest_tpA%peri(i) <= qmin)) then
                     swiftest_tpA%status(i) = DISCARDED_PERI
                     write(*, *) "Particle ", swiftest_tpA%name(i), " perihelion distance too small at t = ", t
                     swiftest_tpA%ldiscard = .true.
                  end if
               end if
            end if
         end if
      end do
   end if

   return

   end procedure discard_peri
end submodule s_discard_peri