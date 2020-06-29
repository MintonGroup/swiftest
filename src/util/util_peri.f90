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
   real(DP)     :: vdotr, e

   if (lfirst) then
      if (qmin_coord == "HELIO") then
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               vdotr = dot_product(tp%xh(:,i), tp%vh(:,i))
               if (vdotr > 0.0_DP) then
                  tp%isperi(i) = 1
               else
                  tp%isperi(i) = -1
               end if
            end if
         end do
      else
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               vdotr = dot_product(tp%xb(:,i), tp%vb(:,i))
               if (vdotr > 0.0_DP) then
                  tp%isperi(i) = 1
               else
                  tp%isperi(i) = -1
               end if
            end if
         end do
      end if
   else
      if (qmin_coord == "HELIO") then
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               vdotr = dot_product(tp%xh(:,i), tp%vh(:,i))
               if (tp%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     tp%isperi(i) = 0
                     call orbel_xv2aeq(mu, tp%xh(i,1), tp%xh(i,2), tp%xh(i,3), &
                                           tp%vh(i,1), tp%vh(i,2), tp%vh(i,3), & 
                                           tp%atp(i), e, tp%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     tp%isperi(i) = 1
                  else
                     tp%isperi(i) = -1
                  end if
               end if
            end if
         end do
      else
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               vdotr = dot_product(tp%xb(:,i), tp%vb(:,i))
               if (tp%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     tp%isperi(i) = 0
                     call orbel_xv2aeq(msys, tp%xb(i,1), tp%xh(i,2), tp%xb(i,3), &
                                           tp%vb(i,1), tp%vb(i,2), tp%vb(i,3), & 
                                           tp%atp(i), e, tp%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     tp%isperi(i) = 1
                  else
                     tp%isperi(i) = -1
                  end if
               end if
            end if
         end do
      end if
   end if

   return

   end procedure util_peri
end submodule s_util_peri
