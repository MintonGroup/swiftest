submodule (symba) s_symba_step_recur
contains
   module procedure symba_step_recur
   !! author: David A. Minton
   !!
   !! Step interacting planets and active test particles ahead in democratic heliocentric coordinates at the current
   !!    recursion level, if applicable, and descend to the next deeper level if necessary
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_step_recur.f90
   !! Adapted from Hal Levison's Swift routine symba5_step_recur.F
use swiftest
implicit none
   logical(lgt)          :: lencounter
   integer(I4B)          :: i, j, irecp, icflg, index_i, index_j, index_pl, index_tp
   real(DP)            :: dtl, dth, sgn
   real(DP), dimension(ndim) :: xr, vr, vbs

! executable code
   dtl = dt0/(ntenc**ireci)
   dth = 0.5_DP*dtl
   if (dtl/dt0 < tiny) then
      write(*, *) "swiftest warning:"
      write(*, *) "   in symba_step_recur, local time step is too small"
      write(*, *) "   roundoff error will be important!"
      call util_exit(failure)
   end if
   irecp = ireci + 1

   if (ireci == 0) then
      icflg = 0
      do i = 1, nplplenc
         if ((plplenc_list%status(i) == active) .and. (plplenc_list%level(i) == ireci)) then
            index_i  = plplenc_list%index1(i)
            index_j  = plplenc_list%index2(i)
            xr(:) = symba_pla%helio%swiftest%xh(:,index_j) - symba_pla%helio%swiftest%xh(:,index_i)
            vr(:) = symba_pla%helio%swiftest%vb(:,index_j) - symba_pla%helio%swiftest%vb(:,index_i)
            call symba_chk(xr(:), vr(:), symba_pla%helio%swiftest%rhill(index_i),   &  
               symba_pla%helio%swiftest%rhill(index_j), dtl, irecp, lencounter,            &
               plplenc_list%lvdotr(i))
            if (lencounter) then
               icflg = 1
               symba_pla%levelg(index_i) = irecp
               symba_pla%levelm(index_i) = max(irecp, symba_pla%levelm(index_i))
               symba_pla%levelg(index_j) = irecp
               symba_pla%levelm(index_j) = max(irecp, symba_pla%levelm(index_j))
               plplenc_list%level(i) = irecp
            end if
         end if
      end do
      do i = 1, npltpenc
         if ((pltpenc_list%status(i) == active) .and. (pltpenc_list%level(i) == ireci)) then
            index_pl  = pltpenc_list%indexpl(i)
            index_tp  = pltpenc_list%indextp(i)
            
            xr(:) = symba_tpa%helio%swiftest%xh(:,index_tp) - symba_pla%helio%swiftest%xh(:,index_pl)
            vr(:) = symba_tpa%helio%swiftest%vb(:,index_tp) - symba_pla%helio%swiftest%vb(:,index_pl)
            call symba_chk(xr(:), vr(:), symba_pla%helio%swiftest%rhill(index_pl), 0.0_DP,   &
               dtl, irecp, lencounter, pltpenc_list%lvdotr(i))
            if (lencounter) then
               icflg = 1
               symba_pla%levelg(index_pl) = irecp
               symba_pla%levelm(index_pl) = max(irecp, symba_pla%levelm(index_pl))
               symba_tpa%levelg(index_tp) = irecp
               symba_tpa%levelm(index_tp) = max(irecp, symba_tpa%levelm(index_tp))
               pltpenc_list%level(i) = irecp
            end if
         end if
      end do
      lencounter = (icflg == 1)
      sgn = 1.0_DP
      call symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_pla, symba_tpa)
      call symba_helio_drift(ireci, npl, symba_pla, dtl)
      if (ntp > 0) call symba_helio_drift_tp(ireci, ntp, symba_tpa, symba_pla%helio%swiftest%mass(1), dtl)
      if (lencounter) call symba_step_recur(lclose, t, irecp, npl, nplm, ntp, symba_pla, symba_tpa, dt0, eoffset, nplplenc, &
         npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type, &
         nplmax, ntpmax, fragmax, param)
      sgn = 1.0_DP
      call symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_pla, symba_tpa) 
      if (lclose) then
         vbs(:) = symba_pla%helio%swiftest%vb(:,1)
         do i = 1, nplplenc
            index_i  = plplenc_list%index1(i) 
            index_j  = plplenc_list%index2(i)
            if (((plplenc_list%status(i) == active) .and.                                       &
                (symba_pla%levelg(index_i) >= ireci) .and.                                      &
                (symba_pla%levelg(index_j) >= ireci))) then
                ! create if statement to check for collisions (ls12) or merger depending on flag lfrag in param.in
                ! determines collisional regime if lfrag=.true. for close encounter massive bodies
                ! call symba_frag_pl(...)
                ! determines if close encounter leads to merger if lfrag=.false.   
               if (config%lfragmentation) then
                  call symba_fragmentation (t, dtl, i, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
                        eoffset, vbs, encounter_file, out_type, npl, symba_pla, nplplenc, plplenc_list, nplmax, &
                        ntpmax, fragmax)
               else
                  call symba_merge_pl(t, dtl, i, nplplenc, plplenc_list, nmergeadd, nmergesub, mergeadd_list, &
                  mergesub_list, eoffset, vbs, encounter_file, out_type, npl, symba_pla)
               end if
             end if
         end do
         do i = 1, npltpenc
            index_pl  = pltpenc_list%indexpl(i) 
            index_tp  = pltpenc_list%indextp(i) 
            if ((pltpenc_list%status(i) == active) .and.                          &
                (symba_pla%levelg(index_pl) >= ireci) .and.                         &
                (symba_tpa%levelg(index_tp) >= ireci)) then                          
               call symba_merge_tp(t, dtl, i, pltpenc_list, vbs, encounter_file, out_type, symba_pla, symba_tpa)            !check later 
            end if
         end do
      end if
      do i = 1, nplplenc
         index_i  = plplenc_list%index1(i) 
         index_j  = plplenc_list%index2(i) 
         if (symba_pla%levelg(index_i) == irecp) symba_pla%levelg(index_i) = ireci
         if (symba_pla%levelg(index_j) == irecp) symba_pla%levelg(index_j) = ireci
         if (plplenc_list%level(i) == irecp) plplenc_list%level(i) = ireci
      end do
      do i = 1, npltpenc
         index_pl  = pltpenc_list%indexpl(i) 
         index_tp  = pltpenc_list%indextp(i) 
         if (symba_pla%levelg(index_pl) == irecp) symba_pla%levelg(index_pl) = ireci
         if (symba_tpa%levelg(index_tp) == irecp) symba_tpa%levelg(index_tp) = ireci
         if (pltpenc_list%level(i) == irecp) pltpenc_list%level(i) = ireci
      end do
   else
      do j = 1, ntenc
         icflg = 0
         do i = 1, nplplenc
            if ((plplenc_list%status(i) == active) .and. (plplenc_list%level(i) == ireci)) then
               index_i  = plplenc_list%index1(i) 
               index_j  = plplenc_list%index2(i) 
               xr(:) = symba_pla%helio%swiftest%xh(:,index_j) - symba_pla%helio%swiftest%xh(:,index_i)
               vr(:) = symba_pla%helio%swiftest%vb(:,index_j) - symba_pla%helio%swiftest%vb(:,index_i)
               call symba_chk(xr(:), vr(:), symba_pla%helio%swiftest%rhill(index_i),    &
                  symba_pla%helio%swiftest%rhill(index_j), dtl, irecp, lencounter,         &
                  plplenc_list%lvdotr(i))
               if (lencounter) then
                  icflg = 1
                  symba_pla%levelg(index_i) = irecp
                  symba_pla%levelm(index_i) = max(irecp, symba_pla%levelm(index_i))
                  symba_pla%levelg(index_j) = irecp
                  symba_pla%levelm(index_j) = max(irecp, symba_pla%levelm(index_j))
                  plplenc_list%level(i) = irecp
               end if
            end if
         end do
         do i = 1, npltpenc
            if ((pltpenc_list%status(i) == active) .and. (pltpenc_list%level(i) == ireci)) then
               index_pl  = pltpenc_list%indexpl(i) 
               index_tp  = pltpenc_list%indextp(i) 
               xr(:) = symba_tpa%helio%swiftest%xh(:,index_tp) - symba_pla%helio%swiftest%xh(:,index_pl)
               vr(:) = symba_tpa%helio%swiftest%vb(:,index_tp)  - symba_pla%helio%swiftest%vb(:,index_pl) 
               call symba_chk(xr(:), vr(:), symba_pla%helio%swiftest%rhill(index_pl), 0.0_DP, &   
                  dtl, irecp, lencounter, pltpenc_list%lvdotr(i))
               if (lencounter) then
                  icflg = 1
                  symba_pla%levelg(index_pl) = irecp
                  symba_pla%levelm(index_pl) = max(irecp, symba_pla%levelm(index_pl))
                  symba_tpa%levelg(index_tp) = irecp
                  symba_tpa%levelm(index_tp) = max(irecp, symba_tpa%levelm(index_tp))
                  pltpenc_list%level(i) = irecp
               end if
            end if
         end do
         lencounter = (icflg == 1)
         sgn = 1.0_DP
         call symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_pla, symba_tpa) 
         sgn = -1.0_DP
         call symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_pla, symba_tpa)
         call symba_helio_drift(ireci, npl, symba_pla, dtl)
         if (ntp > 0) call symba_helio_drift_tp(ireci, ntp, symba_tpa, symba_pla%helio%swiftest%mass(1), dtl)
         if (lencounter) call symba_step_recur(lclose, t, irecp, npl, nplm, ntp, symba_pla, symba_tpa, dt0, eoffset,    &
            nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list,       &
            encounter_file, out_type, nplmax, ntpmax, fragmax, param)
         sgn = 1.0_DP
         call symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_pla, symba_tpa) 
         sgn = -1.0_DP
         call symba_kick(irecp, nplplenc, npltpenc, plplenc_list, pltpenc_list, dth, sgn,symba_pla, symba_tpa)
         if (lclose) then
            vbs(:) = symba_pla%helio%swiftest%vb(:,1)
            do i = 1, nplplenc
               index_i  = plplenc_list%index1(i) 
               index_j  = plplenc_list%index2(i) 
               if ((plplenc_list%status(i) == active) .and.                                     &
                   (symba_pla%levelg(index_i) >= ireci) .and.                                   &
                   (symba_pla%levelg(index_j) >= ireci))  then    
                  if (config%lfragmentation) then
                     call symba_fragmentation (t, dtl, i, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
                        eoffset, vbs, encounter_file, out_type, npl, symba_pla, nplplenc, plplenc_list, nplmax, &
                        ntpmax, fragmax)
                  else
                     call symba_merge_pl(t, dtl, i, nplplenc, plplenc_list, nmergeadd, nmergesub, mergeadd_list, &
                        mergesub_list, eoffset, vbs, encounter_file, out_type, npl, symba_pla)
                  end if
               end if
            end do
            do i = 1, npltpenc
               index_pl  = pltpenc_list%indexpl(i) 
               index_tp  = pltpenc_list%indextp(i) 
               if ((pltpenc_list%status(i) == active) .and.                                     &
                   (symba_pla%levelg(index_pl) >= ireci) .and.                                    &
                   (symba_tpa%levelg(index_tp) >= ireci))                                       &
                  call symba_merge_tp(t, dtl, i, pltpenc_list, vbs, encounter_file, out_type, symba_pla, symba_tpa)          !check that later
            end do
         end if
         do i = 1, nplplenc
            index_i  = plplenc_list%index1(i) 
            index_j  = plplenc_list%index2(i) 
            if (symba_pla%levelg(index_i) == irecp) symba_pla%levelg(index_i) = ireci
            if (symba_pla%levelg(index_j) == irecp) symba_pla%levelg(index_j) = ireci
            if (plplenc_list%level(i) == irecp) plplenc_list%level(i) = ireci
         end do
         do i = 1, npltpenc
            index_pl  = pltpenc_list%indexpl(i) 
            index_tp  = pltpenc_list%indextp(i) 
            if (symba_pla%levelg(index_pl) == irecp) symba_pla%levelg(index_pl) = ireci
            if (symba_tpa%levelg(index_tp) == irecp) symba_tpa%levelg(index_tp) = ireci
            if (pltpenc_list%level(i) == irecp) pltpenc_list%level(i) = ireci
         end do
      end do
   end if

   return

   end procedure symba_step_recur
end submodule s_symba_step_recur
