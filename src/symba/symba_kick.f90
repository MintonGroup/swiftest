submodule (symba) s_symba_kick
contains
   module procedure symba_kick
   !! author: David A. Minton
   !!
   !! Kick barycentric velocities of planets and ACTIVE test particles within SyMBA recursion
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_kick.f90
   !! Adapted from Hal Levison's Swift routine symba5_kick.f
use swiftest
implicit none
   integer(I4B)          :: i, irm1, irecl, index_i,index_j,index_tp,index_pl
   real(DP)            :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
   real(DP), dimension(NDIM) :: dx

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
      symba_plA%ah(:,index_i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      symba_plA%ah(:,index_j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do
   do i = 1, npltpenc
      index_tp  = pltpenc_list%indextp(i)
      symba_tpA%ah(:,index_tp) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do
   do i = 1, nplplenc
      if (plplenc_list%status(i) == ACTIVE) then
         index_i  = plplenc_list%index1(i) 
         index_j  = plplenc_list%index2(i) 
         if ((symba_plA%levelg(index_i) >= irm1) .and. (symba_plA%levelg(index_j) >= irm1)) then
            ri = ((symba_plA%rhill(index_i) &
               + symba_plA%rhill(index_j))**2)*(rhscale**2)*(rshell**(2*irecl))
            rim1 = ri*(rshell**2)
            dx(:) = symba_plA%xh(:,index_j) - symba_plA%xh(:,index_i)
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
            faci = fac*symba_plA%mass(index_i)
            facj = fac*symba_plA%mass(index_j)
            symba_plA%ah(:,index_i) = symba_plA%ah(:,index_i) + facj*dx(:)
            symba_plA%ah(:,index_j) = symba_plA%ah(:,index_j) - faci*dx(:)
         end if
      end if
   end do
   do i = 1, npltpenc
      if (pltpenc_list%status(i) == ACTIVE) then
         index_pl  = pltpenc_list%indexpl(i) 
         index_tp  = pltpenc_list%indextp(i) 
         if ((symba_plA%levelg(index_pl) >= irm1) .and. (symba_tpA%levelg(index_tp) >= irm1)) then
            ri = ((symba_plA%rhill(index_pl))**2)*(rhscale**2)*(rshell**(2*irecl))
            rim1 = ri*(rshell**2)
            dx(:) = symba_tpA%xh(:,index_tp) - symba_plA%xh(:,index_pl)
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
            faci = fac*symba_plA%mass(index_pl)
            symba_tpA%ah(:,index_tp) = symba_tpA%ah(:,index_tp) - faci*dx(:)
         end if
      end if
   end do
   do i = 1, nplplenc
      index_i  = plplenc_list%index1(i) 
      index_j  = plplenc_list%index2(i) 
      symba_plA%vb(:,index_i) = symba_plA%vb(:,index_i) + sgn*dt*symba_plA%ah(:,index_i)
      symba_plA%vb(:,index_j) = symba_plA%vb(:,index_j) + sgn*dt*symba_plA%ah(:,index_j)
      symba_plA%ah(:,index_i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
      symba_plA%ah(:,index_j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do
   do i = 1, npltpenc
      index_tp  = pltpenc_list%indextp(i)
      if (symba_tpA%status(index_tp) == ACTIVE)   &
      symba_tpA%vb(:,index_tp) = symba_tpA%vb(:,index_tp) + sgn*dt*symba_tpA%ah(:,index_tp)
      symba_tpA%ah(:,index_tp) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   end do

   return

   end procedure symba_kick
end submodule s_symba_kick
