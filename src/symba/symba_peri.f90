submodule (symba) s_symba_peri
contains
   module procedure symba_peri
   !! author: David A. Minton
   !!
   !! Determine system pericenter passages for planets in SyMBA
   !!      If the coordinate system used is barycentric, then this routine assumes that the barycentric coordinates in the
   !!      massive body structures are up-to-date and are not recomputed
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_peri.f90
   !! Adapted from Hal Levison's Swift routine util_mass_peri.f
use swiftest
implicit none
   integer(I4B)          :: i
   real(DP)            :: vdotr, e, mu, msun

! executable code
   msun = symba_plA%mass(1)
   if (lfirst) then
      if (qmin_coord == "helio") then
         do i = 2, npl
            if (symba_plA%status(i) == ACTIVE) then
               vdotr = dot_product(symba_plA%xh(:,i), symba_plA%vh(:,i))
               if (vdotr > 0.0_DP) then
                  symba_plA%isperi(i) = 1
               else
                  symba_plA%isperi(i) = -1
               end if
            end if
         end do
      else
         do i = 2, npl
            if (symba_plA%status(i) == ACTIVE) then
               vdotr = dot_product(symba_plA%xb(:,i), symba_plA%vb(:,i))
               if (vdotr > 0.0_DP) then
                  symba_plA%isperi(i) = 1
               else
                  symba_plA%isperi(i) = -1
               end if
            end if
         end do
      end if
   else
      if (qmin_coord == "helio") then
         do i = 2, npl
            if (symba_plA%status(i) == ACTIVE) then
               vdotr = dot_product(symba_plA%xh(:,i), symba_plA%vh(:,i))
               if (symba_plA%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     symba_plA%isperi(i) = 0
                     mu = msun + symba_plA%mass(i)
                     call orbel_xv2aeq(symba_plA%xh(:,i), &
                        symba_plA%vh(:,i), mu, symba_plA%atp(i), e, symba_plA%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     symba_plA%isperi(i) = 1
                  else
                     symba_plA%isperi(i) = -1
                  end if
               end if
            end if
         end do
      else
         do i = 2, npl
            if (symba_plA%status(i) == ACTIVE) then
               vdotr = dot_product(symba_plA%xb(:,i), symba_plA%vb(:,i))
               if (symba_plA%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     symba_plA%isperi(i) = 0
                     call orbel_xv2aeq(symba_plA%xb(:,i), & 
                        symba_plA%vb(:,i), msys, symba_plA%atp(i), e, symba_plA%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     symba_plA%isperi(i) = 1
                  else
                     symba_plA%isperi(i) = -1
                  end if
               end if
            end if
         end do
      end if
   end if

   return

   end procedure symba_peri
end submodule s_symba_peri
