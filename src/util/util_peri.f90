submodule (util) s_util_peri
contains
   module procedure util_peri
   !! author: David A. Minton
   !!
   !! Determine system pericenter passages for test particles
   !! Note:  If the coordinate system used is barycentric, then this routine assumes that the barycentric coordinates in the
   !!        test particle structures are up-to-date and are not recomputed
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: util_peri.f90
   !! Adapted from Hal Levison's Swift routine util_peri.f
   use swiftest
   implicit none
   
   integer(I4B) :: i
   real(DP)     :: e
   real(DP), dimension(:), allocatable :: vdotr

   associate(ntp => tp%nbody, xht => tp%xh, vht => tp%vh, status => tp%status, isperi => tp%isperi, &
             xbt => tp%xb, vbt => tp%vb, atp => tp%atp, peri => tp%peri)
      if (lfirst) then
         if (qmin_coord == "HELIO") then
            where(status(1:ntp) == ACTIVE)
               vdotr(1:ntp) = xht(:, 1:ntp) .dot. vht(:, 1:ntp)
            elsewhere
               vdotr(1:ntp) = 0.0_DP
            end where
         else
            where(status(1:ntp) == ACTIVE)
               vdotr(1:ntp) = xbt(:, 1:ntp) .dot. vbt(:, 1:ntp)
            elsewhere
               vdotr(1:ntp) = 0.0_DP
            end where 
         end if
         where(vdotr(1:ntp) > 0.0_DP)
            isperi(1:ntp) = 1
         elsewhere
            isperi = -1
         end where
      else
         if (qmin_coord == "HELIO") then
            do concurrent (i = 1:ntp, status(i) == ACTIVE)
               vdotr(i) = xht(:, i) .dot. vht(:, i)
               if (isperi(i) == -1) then
                  if (vdotr(i) >= 0.0_DP) then
                     isperi(i) = 0
                     call orbel_xv2aeq(mu, xht(1, i), xht(2, i), xht(3, i), &
                                           vht(1, i), vht(2, i), vht(3, i), & 
                                           atp(i), e, peri(i))
                  end if
               else
                  if (vdotr(i) > 0.0_DP) then
                     isperi(i) = 1
                  else
                     isperi(i) = -1
                  end if
               end if
            end do
         else
            do concurrent (i = 1:ntp, status(i) == ACTIVE)
               vdotr(i) = xbt(:, i) .dot. vbt(:, i)
               if (isperi(i) == -1) then
                  if (vdotr(i) >= 0.0_DP) then
                     isperi(i) = 0
                     call orbel_xv2aeq(msys, xbt(1, i), xbt(2, i), xbt(3, i), &
                                           vbt(1, i), vbt(2, i), vbt(3, i), & 
                                           atp(i), e, peri(i))
                  end if
               else
                  if (vdotr(i) > 0.0_DP) then
                     isperi(i) = 1
                  else
                     isperi(i) = -1
                  end if
               end if
            end do
         end if
      end if
   end associate
   return

   end procedure util_peri
end submodule s_util_peri
