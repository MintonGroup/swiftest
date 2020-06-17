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
   msun = symba_pla%helio%swiftest%mass(1)
   if (lfirst) then
      if (qmin_coord == "helio") then
         do i = 2, npl
            if (symba_pla%helio%swiftest%status(i) == active) then
               vdotr = dot_product(symba_pla%helio%swiftest%xh(:,i), symba_pla%helio%swiftest%vh(:,i))
               if (vdotr > 0.0_DP) then
                  symba_pla%isperi(i) = 1
               else
                  symba_pla%isperi(i) = -1
               end if
            end if
         end do
      else
         do i = 2, npl
            if (symba_pla%helio%swiftest%status(i) == active) then
               vdotr = dot_product(symba_pla%helio%swiftest%xb(:,i), symba_pla%helio%swiftest%vb(:,i))
               if (vdotr > 0.0_DP) then
                  symba_pla%isperi(i) = 1
               else
                  symba_pla%isperi(i) = -1
               end if
            end if
         end do
      end if
   else
      if (qmin_coord == "helio") then
         do i = 2, npl
            if (symba_pla%helio%swiftest%status(i) == active) then
               vdotr = dot_product(symba_pla%helio%swiftest%xh(:,i), symba_pla%helio%swiftest%vh(:,i))
               if (symba_pla%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     symba_pla%isperi(i) = 0
                     mu = msun + symba_pla%helio%swiftest%mass(i)
                     call orbel_xv2aeq(symba_pla%helio%swiftest%xh(:,i), &
                        symba_pla%helio%swiftest%vh(:,i), mu, symba_pla%atp(i), e, symba_pla%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     symba_pla%isperi(i) = 1
                  else
                     symba_pla%isperi(i) = -1
                  end if
               end if
            end if
         end do
      else
         do i = 2, npl
            if (symba_pla%helio%swiftest%status(i) == active) then
               vdotr = dot_product(symba_pla%helio%swiftest%xb(:,i), symba_pla%helio%swiftest%vb(:,i))
               if (symba_pla%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     symba_pla%isperi(i) = 0
                     call orbel_xv2aeq(symba_pla%helio%swiftest%xb(:,i), & 
                        symba_pla%helio%swiftest%vb(:,i), msys, symba_pla%atp(i), e, symba_pla%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     symba_pla%isperi(i) = 1
                  else
                     symba_pla%isperi(i) = -1
                  end if
               end if
            end if
         end do
      end if
   end if

   return

   end procedure symba_peri
end submodule s_symba_peri
