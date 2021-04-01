submodule (util) s_util_peri
   use swiftest
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
   implicit none
   
   integer(I4B) :: i
   real(DP)     :: e
   real(DP), dimension(:), allocatable, save :: vdotr

   associate(ntp => tp%nbody, xht => tp%xh, vht => tp%vh, status => tp%status, isperi => tp%isperi, &
             xbt => tp%xb, vbt => tp%vb, atp => tp%atp, peri => tp%peri)
      if (lfirst) then
         if (.not. allocated(vdotr)) allocate(vdotr(ntp))
         if (qmin_coord == "HELIO") then
            !do concurrent(i = 1:ntp, status(i) == ACTIVE)
            do i = 1, ntp
               vdotr(i) = dot_product(xht(:, i), vht(:, i))
            end do
         else
            !do concurrent(i = 1:ntp, status(i) == ACTIVE)
            do i = 1, ntp
               vdotr(i) = dot_product(xbt(:, i), vbt(:, i))
            end do
         end if
         where(vdotr(1:ntp) > 0.0_DP)
            isperi(1:ntp) = 1
         elsewhere
            isperi = -1
         end where
      else
         if (qmin_coord == "HELIO") then
            !do concurrent (i = 1:ntp, status(i) == ACTIVE)
            do i = 1, ntp
               vdotr(i) = dot_product(xht(:, i), vht(:, i))
               if (isperi(i) == -1) then
                  if (vdotr(i) >= 0.0_DP) then
                     isperi(i) = 0
                     call orbel_xv2aeq(mu, xht(:, i), vht(:, i), atp(i), e, peri(i))
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
            !do concurrent (i = 1:ntp, status(i) == ACTIVE)
            do i = 1, ntp
               vdotr(i) = dot_product(xbt(:, i), vbt(:, i))
               if (isperi(i) == -1) then
                  if (vdotr(i) >= 0.0_DP) then
                     isperi(i) = 0
                     call orbel_xv2aeq(msys, xbt(:, i), vbt(:, i), atp(i), e, peri(i))
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
