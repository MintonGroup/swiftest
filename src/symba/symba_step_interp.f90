submodule (symba) s_symba_step_interp
contains
   module procedure symba_step_interp
   !! author: David A. Minton
   !!
   !! Step planets and active test particles ahead in democratic heliocentric coordinates, calling the recursive
   !! subroutine to descend to the appropriate level to handle close encounters
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_step_interp.f90
   !! Adapted from Hal Levison's Swift routine symba5_step_interp.f
   use swiftest
   implicit none
   logical , save                 :: lmalloc = .true.
   integer( I4B)                     :: i, irec
   real(DP)                       :: dth, msys
   real(DP), dimension(NDIM)            :: ptb, pte
   real(DP), dimension(:, :), allocatable, save :: xbeg, xend

! executable code

   if (lmalloc) then
      allocate(xbeg(NDIM, config%nplmax), xend(NDIM, config%nplmax))
      lmalloc = .false.
   end if
   dth = 0.5_DP*dt

   call coord_vh2vb(npl, symba_plA, msys)

   call helio_lindrift(npl, symba_plA, dth, ptb)
   if (ntp > 0) then
      call coord_vh2vb_tp(ntp, symba_tpA, -ptb)
      call helio_lindrift_tp(ntp, symba_tpA, dth, ptb) 
      do i = 2, npl
         xbeg(:, i) = symba_plA%xh(:,i)
      end do
   end if

   call symba_getacch(lextra_force, t, npl, nplm, symba_plA, config%j2rp2, config%j4rp4, nplplenc, plplenc_list)
   if (ntp > 0) call symba_getacch_tp(lextra_force, t, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, xbeg, config%j2rp2,   &
      config%j4rp4, npltpenc, pltpenc_list)

   call helio_kickvb(npl, symba_plA, dth)
   if (ntp > 0) call helio_kickvb_tp(ntp, symba_tpA, dth)
   irec = -1

   call symba_helio_drift(irec, npl, symba_plA, dt)
   if (ntp > 0) call symba_helio_drift_tp(irec, ntp, symba_tpA, symba_plA%mass(1), dt)
   irec = 0

   call symba_step_recur(lclose, t, irec, npl, nplm, ntp, symba_plA, symba_tpA, dt, eoffset, nplplenc, npltpenc,          &
      plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type, &
      config%nplmax, config%ntpmax, fragmax, param)
   if (ntp > 0) then
      do i = 2, npl
         xend(:, i) = symba_plA%xh(:,i)
      end do
   end if
   call symba_getacch(lextra_force, t+dt, npl, nplm, symba_plA, config%j2rp2, config%j4rp4, nplplenc, plplenc_list)
   if (ntp > 0) call symba_getacch_tp(lextra_force, t+dt, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, xend, config%j2rp2,  &
      config%j4rp4, npltpenc, pltpenc_list)
   call helio_kickvb(npl, symba_plA, dth)
   if (ntp > 0) call helio_kickvb_tp(ntp, symba_tpA, dth)
   call coord_vb2vh(npl, symba_plA)
   call helio_lindrift(npl, symba_plA, dth, pte)
   if (ntp > 0) then
      call coord_vb2vh_tp(ntp, symba_tpA, -pte)
      call helio_lindrift_tp(ntp, symba_tpA, dth, pte)
   end if
   return

   end procedure symba_step_interp
end submodule s_symba_step_interp
