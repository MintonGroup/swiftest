submodule(rmvs_classes) s_rmvs_peri
contains   
   module procedure rmvs_peri_tp 
      !! author: David A. Minton
      !!
      !! Determine planetocentric pericenter passages for test particles in close encounters with a planet
      !! 
      !! Adapted from Hal Levison's Swift routine Adapted from Hal Levison's Swift routine util_peri.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_peri.f90
      use swiftest
      implicit none

      integer(I4B)              :: i, id1, id2
      real(DP)                  :: r2, mu, rhill2, vdotr, a, peri, capm, tperi, rpl
      real(DP), dimension(NDIM) :: xh1, xh2, vh1, vh2

      rhill2 = pl%Rhill(ipleP)**2
      link = pl%tpenc1P(ipleP)
      mu = pl%Gmass(ipleP)
      associate(nenc => self%nbody)
      if (lfirst) then
         do i = 1, nenc
            if (self%status(i) == ACTIVE) then
               vdotr = dot_product(self%xh(:, i), self%vh(:, i))
               if (vdotr > 0.0_DP) then
                  self%isperi(i) = 1
               else
                  self%isperi(i) = -1
               end if
            end if
         end do
      else
         do i = 1, nenc
            if (tp%status(i) == ACTIVE) then
               vdotr = dot_product(tp%xh(:, i),tp%vh(:, i))
               if (self%isperi(i) == -1) then
                  if (vdotr >= 0.0_DP) then
                     self%isperi(i) = 0
                     call orbel_xv2aqt(mu, self%xh(:, i), self%vh(:, i), a, q, capm, tperi)
                     r2 = dot_product(self%xh(:, i), self%xh(:, i))
                     if ((abs(tperi) > FACQDT * dt) .or. (r2 > rhill2)) peri = sqrt(r2)
                     if (config%encounter_file /= "") then
                        id1 = pl%name(pleP)
                        rpl = pl%radius(pleP)
                        xh1(:) = pl%xin(:, pleP, index)
                        vh1(:) = pl%vin(:, pleP, index)
                        id2 = self%name(link)
                        xh2(:) = self%xh(:, link) + pl%xin(:, pleP, index) 
                        vh2(:) = self%vh(:, link) + pl%vin(:, pleP, index)
                        call io_write_encounter(t, id1, id2, mu, 0.0_DP, rpl, 0.0_DP, xh1(:), xh2(:), vh1(:), vh2(:),  &
                           config%encounter_file, config%out_type)
                     end if
                     if (self%lperi(i)) then
                        if (peri < self%peri(i)) then
                           self%peri(i) = peri
                           self%plperP(i) = ipleP
                        end if
                     else
                        self%lperi(i) = .true.
                        self%peri(i) = peri
                        self%plperP(i) = ipleP
                     end if
                  end if
               else
                  if (vdotr > 0.0_DP) then
                     self%isperi(i) = 1
                  else
                     self%isperi(i) = -1
                  end if
               end if
            end if
         end do                   
      end if

      return

   end procedure rmvs_peri_tp
end submodule s_rmvs_peri