submodule (util) s_util_peri
contains
   module procedure util_peri
   !! author: David A. Minton
   !!
   !! Determine system pericenter passages for test particles
   !! Note:  If the coordinate system used is barycentric, then this routine assumes that the barycentric coordinates in the
   !!        test particle structures are up-to-date and are not recomputed
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: util_peri.f90
   !! Adapted from Hal Levison's Swift routine util_peri.f
   use swiftest
   implicit none
   
   integer(I4B) :: i
   real(DP)     :: vdotr, e

   if (lfirst) then
      if (qmin_coord == "HELIO") then
         do i = 1, ntp
            if (swiftest_tpA%status(i) == ACTIVE) then
               vdotr = dot_product(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i))
               if (vdotr > 0.0_DP) then
                  swiftest_tpA%isperi(i) = 1
               else
                  swiftest_tpA%isperi(i) = -1
               end if
            end if
         end do
      else
         do i = 1, ntp
            if (swiftest_tpA%status(i) == ACTIVE) then
               vdotr = dot_product(swiftest_tpA%xb(:,i), swiftest_tpA%vb(:,i))
               if (vdotr > 0.0_DP) then
                  swiftest_tpA%isperi(i) = 1
               else
                  swiftest_tpA%isperi(i) = -1
               end if
            end if
         end do
      end if
   else
      if (qmin_coord == "HELIO") then
         do i = 1, ntp
            if (swiftest_tpA%status(i) == ACTIVE) then
               vdotr = dot_product(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i))
               if (swiftest_tpA%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     swiftest_tpA%isperi(i) = 0
                     call orbel_xv2aeq(swiftest_tpA%xh(:,i), swiftest_tpA%vh(:,i), mu, & 
                        swiftest_tpA%atp(i), e, swiftest_tpA%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     swiftest_tpA%isperi(i) = 1
                  else
                     swiftest_tpA%isperi(i) = -1
                  end if
               end if
            end if
         end do
      else
         do i = 1, ntp
            if (swiftest_tpA%status(i) == ACTIVE) then
               vdotr = dot_product(swiftest_tpA%xb(:,i), swiftest_tpA%vb(:,i))
               if (swiftest_tpA%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     swiftest_tpA%isperi(i) = 0
                     call orbel_xv2aeq(swiftest_tpA%xb(:,i), swiftest_tpA%vb(:,i), msys, & 
                      swiftest_tpA%atp(i), e, swiftest_tpA%peri(i))
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     swiftest_tpA%isperi(i) = 1
                  else
                     swiftest_tpA%isperi(i) = -1
                  end if
               end if
            end if
         end do
      end if
   end if

   return

   end procedure util_peri
end submodule s_util_peri
