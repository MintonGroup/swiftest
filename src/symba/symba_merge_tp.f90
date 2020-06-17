submodule (symba) s_symba_merge_tp
contains
   module procedure symba_merge_tp
   !! author: David A. Minton
   !!
   !! Check for merger between planet and test particle in SyMBAs
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_merge_tp.f90
   !! Adapted from Hal Levison's Swift routine symba5_merge.f
use swiftest
implicit none
   logical(lgt)          :: lmerge
   integer(I4B)          :: name1, name2, indexpl, indextp
   real(DP)            :: r2, rlim, rlim2, vdotr, tcr2, dt2, mu, a, e, q, rad1
   real(DP), dimension(ndim) :: xr, vr, xh1, vh1, xh2, vh2

! executable code
   lmerge = .false.
   
   indexpl = pltpenc_list%indexpl(index_enc)
   indextp = pltpenc_list%indextp(index_enc)

   rlim = symba_pla%helio%swiftest%radius(indexpl)
   xr(:) = symba_tpa%helio%swiftest%xh(:,indextp) - symba_pla%helio%swiftest%xh(:,indexpl)
   r2 = dot_product(xr(:), xr(:))
   rlim2 = rlim*rlim
   if (rlim2 >= r2) then
      lmerge = .true.
   else
      vr(:) = symba_tpa%helio%swiftest%vb(:,indextp) - symba_pla%helio%swiftest%vb(:,indexpl)
      vdotr = dot_product(xr(:), vr(:))
      if (pltpenc_list%lvdotr(index_enc) .and. (vdotr > 0.0_DP)) then
         mu = symba_pla%helio%swiftest%mass(indexpl)
         tcr2 = r2/dot_product(vr(:), vr(:))
         dt2 = dt*dt
         if (tcr2 <= dt2) then
            call orbel_xv2aeq(xr(:), vr(:), mu, a, e, q)
            if (q < rlim) lmerge = .true.
         end if
         if (.not. lmerge) then
            if (encounter_file /= "") then
               name1 = symba_pla%helio%swiftest%name(indexpl)
               rad1 = symba_pla%helio%swiftest%radius(indexpl)
               xh1(:) = symba_pla%helio%swiftest%xh(:,indexpl)
               vh1(:) = symba_pla%helio%swiftest%vb(:,indexpl) - vbs(:)
               name2 = symba_tpa%helio%swiftest%name(indextp)
               xh2(:) = symba_tpa%helio%swiftest%xh(:,indextp)
               vh2(:) = symba_tpa%helio%swiftest%vb(:,indextp) - vbs(:)
               call io_write_encounter(t, name1, name2, mu, 0.0_DP, rad1, 0.0_DP, &
                  xh1(:), xh2(:), vh1(:), vh2(:), encounter_file, out_type)
            end if
         end if
      end if
   end if
   if (lmerge) then
      pltpenc_list%status(index_enc) = merged
      symba_tpa%helio%swiftest%status = discarded_plr
      write(*, *) "particle ", symba_tpa%helio%swiftest%name, " too close to massive body ", &
      symba_pla%helio%swiftest%name, " at t = ", t
   end if

   return

   end procedure symba_merge_tp
end submodule s_symba_merge_tp
