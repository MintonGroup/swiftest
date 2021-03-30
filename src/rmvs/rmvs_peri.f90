submodule(rmvs_classes) s_rmvs_peri
contains   
   module subroutine rmvs_peri_tp(self, cb, pl, t, dt, lfirst, index, nenc, ipleP, config)
      !! author: David A. Minton
      !!
      !! Determine planetocentric pericenter passages for test particles in close encounters with a planet
      !! 
      !! Adapted from Hal Levison's Swift routine Adapted from Hal Levison's Swift routine util_peri.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_peri.f90
      use swiftest
      implicit none
      ! Arguments
      class(rmvs_tp),                intent(inout) :: self   !! RMVS test particle object  
      class(rmvs_cb),                intent(inout) :: cb     !! RMVS central body object  
      class(rmvs_pl),                intent(inout) :: pl     !! RMVS massive body object  
      real(DP),                      intent(in)    :: t      !! current time
      real(DP),                      intent(in)    :: dt     !! step size
      logical,                       intent(in)    :: lfirst !! Logical flag indicating whether current invocation is the first
      integer(I4B),                  intent(in)    :: index !! outer substep number within current set
      integer(I4B),                  intent(in)    :: nenc  !! number of test particles encountering current planet 
      integer(I4B),                  intent(in)    :: ipleP !!index of RMVS planet being closely encountered
      class(swiftest_configuration), intent(in)    :: config  !! Input collection of  configuration parameters
      ! Internals
      integer(I4B)                                 :: i, id1, id2
      real(DP)                                     :: r2, mu, rhill2, vdotr, a, peri, capm, tperi, rpl
      real(DP), dimension(NDIM)                    :: xh1, xh2, vh1, vh2

      rhill2 = pl%Rhill(ipleP)**2
      mu = pl%Gmass(ipleP)
      associate(tp => self, nenc => self%nbody)
         if (lfirst) then
            do i = 1, nenc
               if (tp%status(i) == ACTIVE) then
                  vdotr = dot_product(tp%xh(:, i), tp%vh(:, i))
                  if (vdotr > 0.0_DP) then
                     tp%isperi(i) = 1
                  else
                     tp%isperi(i) = -1
                  end if
               end if
            end do
         else
            do i = 1, nenc
               if (tp%status(i) == ACTIVE) then
                  vdotr = dot_product(tp%xh(:, i), tp%vh(:, i))
                  if (tp%isperi(i) == -1) then
                     if (vdotr >= 0.0_DP) then
                        tp%isperi(i) = 0
                        call orbel_xv2aqt(mu, tp%xh(:, i), tp%vh(:, i), a, peri, capm, tperi)
                        r2 = dot_product(tp%xh(:, i), tp%xh(:, i))
                        if ((abs(tperi) > FACQDT * dt) .or. (r2 > rhill2)) peri = sqrt(r2)
                        if (config%encounter_file /= "") then
                           id1 = pl%name(ipleP)
                           rpl = pl%radius(ipleP)
                           xh1(:) = pl%xin(:, ipleP, index)
                           vh1(:) = pl%vin(:, ipleP, index)
                           id2 = tp%name(i)
                           xh2(:) = tp%xh(:, i) + pl%xin(:, ipleP, index) 
                           vh2(:) = tp%vh(:, i) + pl%vin(:, ipleP, index)
                           call io_write_encounter(t, id1, id2, mu, 0.0_DP, rpl, 0.0_DP, xh1(:), xh2(:), vh1(:), vh2(:),  &
                              config%encounter_file, config%out_type)
                        end if
                        if (tp%lperi(i)) then
                           if (peri < tp%peri(i)) then
                              tp%peri(i) = peri
                              tp%plperP(i) = ipleP
                           end if
                        else
                           tp%lperi(i) = .true.
                           tp%peri(i) = peri
                           tp%plperP(i) = ipleP
                        end if
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
         end associate
      return

   end subroutine rmvs_peri_tp
end submodule s_rmvs_peri