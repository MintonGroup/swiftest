submodule (symba) s_symba_kick
contains
   module procedure symba_kick
   !! author: David A. Minton
   !!
   !! Kick barycentric velocities of planets and active test particles within SyMBA recursion
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_kick.f90
   !! Adapted from Hal Levison's Swift routine symba5_kick.f
use swiftest
implicit none
   integer(I4B)          :: i, irm1, irecl, index_i,index_j,index_tp,index_pl
   real(DP)            :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
   real(DP), dimension(ndim) :: dx

! executable code
   irm1 = irec - 1
   if (sgn < 0.0_DP) then
      irecl = irec - 1
   else
      irecl = irec
   end if
   do i = 1, nplplenc
      index_i  = plplenc_list%index1(i)
      index_j  = plplenc_list%index2(i)
      symba_pla%helio%ah(:,index_i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      symba_pla%helio%ah(:,index_j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do
   do i = 1, npltpenc
      index_tp  = pltpenc_list%indextp(i)
      symba_tpa%helio%ah(:,index_tp) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do
   do i = 1, nplplenc
      if (plplenc_list%status(i) == active) then
         index_i  = plplenc_list%index1(i) 
         index_j  = plplenc_list%index2(i) 
         if ((symba_pla%levelg(index_i) >= irm1) .and. (symba_pla%levelg(index_j) >= irm1)) then
            ri = ((symba_pla%helio%swiftest%rhill(index_i) &
               + symba_pla%helio%swiftest%rhill(index_j))**2)*(rhscale**2)*(rshell**(2*irecl))
            rim1 = ri*(rshell**2)
            dx(:) = symba_pla%helio%swiftest%xh(:,index_j) - symba_pla%helio%swiftest%xh(:,index_i)
            r2 = dot_product(dx(:), dx(:))
            if (r2 < rim1) then
               fac = 0.0_DP
            else if (r2 < ri) then
               ris = sqrt(ri)
               r = sqrt(r2)
               rr = (ris - r)/(ris*(1.0_DP - rshell))
               fac = (r2**(-1.5_DP))*(1.0_DP - 3.0_DP*(rr**2) + 2.0_DP*(rr**3))
            else
               ir3 = 1.0_DP/(r2*sqrt(r2))
               fac = ir3
            end if
            faci = fac*symba_pla%helio%swiftest%mass(index_i)
            facj = fac*symba_pla%helio%swiftest%mass(index_j)
            symba_pla%helio%ah(:,index_i) = symba_pla%helio%ah(:,index_i) + facj*dx(:)
            symba_pla%helio%ah(:,index_j) = symba_pla%helio%ah(:,index_j) - faci*dx(:)
         end if
      end if
   end do
   do i = 1, npltpenc
      if (pltpenc_list%status(i) == active) then
         index_pl  = pltpenc_list%indexpl(i) 
         index_tp  = pltpenc_list%indextp(i) 
         if ((symba_pla%levelg(index_pl) >= irm1) .and. (symba_tpa%levelg(index_tp) >= irm1)) then
            ri = ((symba_pla%helio%swiftest%rhill(index_pl))**2)*(rhscale**2)*(rshell**(2*irecl))
            rim1 = ri*(rshell**2)
            dx(:) = symba_tpa%helio%swiftest%xh(:,index_tp) - symba_pla%helio%swiftest%xh(:,index_pl)
            r2 = dot_product(dx(:), dx(:))
            if (r2 < rim1) then
               fac = 0.0_DP
            else if (r2 < ri) then
               ris = sqrt(ri)
               r = sqrt(r2)
               rr = (ris - r)/(ris*(1.0_DP - rshell))
               fac = (r2**(-1.5_DP))*(1.0_DP - 3.0_DP*(rr**2) + 2.0_DP*(rr**3))
            else
               ir3 = 1.0_DP/(r2*sqrt(r2))
               fac = ir3
            end if
            faci = fac*symba_pla%helio%swiftest%mass(index_pl)
            symba_tpa%helio%ah(:,index_tp) = symba_tpa%helio%ah(:,index_tp) - faci*dx(:)
         end if
      end if
   end do
   do i = 1, nplplenc
      index_i  = plplenc_list%index1(i) 
      index_j  = plplenc_list%index2(i) 
      symba_pla%helio%swiftest%vb(:,index_i) = symba_pla%helio%swiftest%vb(:,index_i) + sgn*dt*symba_pla%helio%ah(:,index_i)
      symba_pla%helio%swiftest%vb(:,index_j) = symba_pla%helio%swiftest%vb(:,index_j) + sgn*dt*symba_pla%helio%ah(:,index_j)
      symba_pla%helio%ah(:,index_i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      symba_pla%helio%ah(:,index_j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do
   do i = 1, npltpenc
      index_tp  = pltpenc_list%indextp(i)
      if (symba_tpa%helio%swiftest%status(index_tp) == active)   &
      symba_tpa%helio%swiftest%vb(:,index_tp) = symba_tpa%helio%swiftest%vb(:,index_tp) + sgn*dt*symba_tpa%helio%ah(:,index_tp)
      symba_tpa%helio%ah(:,index_tp) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do

   return

   end procedure symba_kick
end submodule s_symba_kick
