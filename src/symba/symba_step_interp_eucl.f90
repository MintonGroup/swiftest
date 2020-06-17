submodule (symba) s_symba_step_interp_eucl
contains
   module procedure symba_step_interp_eucl
   !! author: Jacob R. Elliott
   !!
   !! Same as symba_step_interp, but with th new single loop-blocking Euclidean distance matrix evaluation
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_step_interp.f90
   !! Adapted from Hal Levison's Swift routine symba5_step_interp.f
   use swiftest
   implicit none
   logical(lgt), save                 :: lmalloc = .true.
   integer(I4B)                     :: i, irec
   real(DP)                       :: dth, msys
   real(DP), dimension(ndim)            :: ptb, pte
   real(DP), dimension(:, :), allocatable, save :: xbeg, xend

! executable code

   if (lmalloc) then
      allocate(xbeg(ndim, nplmax), xend(ndim, nplmax))
      lmalloc = .false.
   end if
   dth = 0.5_DP*dt

   call coord_vh2vb(npl, symba_pla%helio%swiftest, msys)

   call helio_lindrift(npl, symba_pla%helio%swiftest, dth, ptb)
   if (ntp > 0) then
      call coord_vh2vb_tp(ntp, symba_tpa%helio%swiftest, -ptb)
      call helio_lindrift_tp(ntp, symba_tpa%helio%swiftest, dth, ptb) 
      do i = 2, npl
         xbeg(:, i) = symba_pla%helio%swiftest%xh(:,i)
      end do
   end if

   call symba_getacch_eucl(lextra_force, t, npl, nplm, nplmax, symba_pla, j2rp2, j4rp4, nplplenc, plplenc_list, &
      num_plpl_comparisons, k_plpl)
   if (ntp > 0) call symba_getacch_tp_eucl(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_pla, symba_tpa, xbeg, j2rp2,&
      j4rp4, npltpenc, pltpenc_list, num_pltp_comparisons, k_pltp)

   call helio_kickvb(npl, symba_pla%helio, dth)
   if (ntp > 0) call helio_kickvb_tp(ntp, symba_tpa%helio, dth)
   irec = -1

   call symba_helio_drift(irec, npl, symba_pla, dt)
   if (ntp > 0) call symba_helio_drift_tp(irec, ntp, symba_tpa, symba_pla%helio%swiftest%mass(1), dt)
   irec = 0

   call symba_step_recur(lclose, t, irec, npl, nplm, ntp, symba_pla, symba_tpa, dt, eoffset, nplplenc, npltpenc,&
      plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type)
   if (ntp > 0) then
      do i = 2, npl
         xend(:, i) = symba_pla%helio%swiftest%xh(:,i)
      end do
   end if
   call symba_getacch_eucl(lextra_force, t+dt, npl, nplm, nplmax, symba_pla, j2rp2, j4rp4, nplplenc, plplenc_list, &
      num_plpl_comparisons, k_plpl)
   if (ntp > 0) call symba_getacch_tp_eucl(lextra_force, t+dt, npl, nplm, nplmax, ntp, ntpmax, symba_pla, symba_tpa, xend, &
      j2rp2,j4rp4, npltpenc, pltpenc_list, num_pltp_comparisons, k_pltp)
   call helio_kickvb(npl, symba_pla%helio, dth)
   if (ntp > 0) call helio_kickvb_tp(ntp, symba_tpa%helio, dth)
   call coord_vb2vh(npl, symba_pla%helio%swiftest)
   call helio_lindrift(npl, symba_pla%helio%swiftest, dth, pte)
   if (ntp > 0) then
      call coord_vb2vh_tp(ntp, symba_tpa%helio%swiftest, -pte)
      call helio_lindrift_tp(ntp, symba_tpa%helio%swiftest, dth, pte)
   end if

   return

   end procedure symba_step_interp_eucl
end submodule s_symba_step_interp_eucl
