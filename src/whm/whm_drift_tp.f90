submodule(whm_classes) s_whm_drift_tp
contains
   module procedure whm_drift_tp
      !! author: David A. Minton
      !!
      !! Loop through test particles and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift_tp.f 
      !! Includes 
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f90
      use swiftest
      implicit none
      integer(I4B) :: i
      real(DP),     dimension(:), allocatable :: dtp, energy, vmag2, rmag  !! Variables used in GR calculation
      integer(I4B), dimension(:), allocatable :: iflag

      associate( ntp => self%nbody)
         if (config%lgr) then
            where(self%status(1:ntp) == ACTIVE)
               rmag(:) = .mag. self%xh(1:ntp, :) 
               vmag2(:) = ().mag. self%vh(1:ntp, :))**2! dot_product(self%vh(:,i), self%vh(:,i))
               energy(:) = 0.5_DP * vmag2(:) - cb%Gmass / rmag(:)
               dtp(:) = dt * (1.0_DP + 3 * config%inv_c2 * energy(:))
            end where
         else
            dtp = dt
         end if
         !do i = 1, ntp
         where(self%status(1:ntp) == ACTIVE)

            iflag = self%drift(cb, config, dtp) !drift_one(mu, swifter_tpp%xh(:), swifter_tpp%vh(:), dtp, iflag)
            if (iflag /= 0) then
               self%status(i) = DISCARDED_DRIFTERR
               write(*, *) "Particle ", self%id, " lost due to error in danby drift"
            end if
         end if
         whm_tpp => whm_tpp%nextp
      end do
      end associate

      return

      end procedure whm_drift_tp
end submodule s_whm_drift_tp
