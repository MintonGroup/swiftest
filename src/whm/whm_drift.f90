submodule(whm_classes) whm_drift
   use swiftest
contains
   module subroutine whm_drift_pl(self, cb, config, dt)
      !! author: David A. Minton
      !!
      !! Loop through planets and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift.f90
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structur
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B)                                 :: i
      real(DP)                                     :: energy, vmag2, rmag  !! Variables used in GR calculation
      integer(I4B), dimension(:), allocatable      :: iflag
      real(DP), dimension(:), allocatable          :: dtp

      associate(npl    => self%nbody, &
         xj     => self%xj, &
         vj     => self%vj, &
         status => self%status,&
         mu     => self%muj)

         if (npl == 0) return

         allocate(iflag(npl))
         iflag(:) = 0
         allocate(dtp(npl))

         if (config%lgr) then
            do i = 1,npl
               rmag = norm2(xj(:, i))
               vmag2 = dot_product(vj(:, i),  vj(:, i))
               energy = 0.5_DP * vmag2 - mu(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * config%inv_c2 * energy)
            end do
         else
            dtp(:) = dt
         end if 

         !!$omp simd
         call drift_one(mu(1:npl), xj(1,1:npl), xj(2,1:npl), xj(3,1:npl), &
                                   vj(1,1:npl), vj(2,1:npl), vj(3,1:npl), &
                                   dtp(1:npl), iflag(1:npl))
         if (any(iflag(1:npl) /= 0)) then
            do i = 1, npl
               if (iflag(i) /= 0) then
                  write(*, *) " Planet ", self%name(i), " is lost!!!!!!!!!!"
                  write(*, *) xj(:,i)
                  write(*, *) vj(:,i)
                  write(*, *) " stopping "
                  call util_exit(FAILURE)
               end if
            end do
         end if
      end associate

      return

   end subroutine whm_drift_pl

   module subroutine whm_drift_tp(self, cb, config, dt)
      !! author: David A. Minton
      !!
      !! Loop through test particles and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift_tp.f 
      !! Includes 
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                 intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structuree
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B)                                 :: i   
      real(DP)                                     :: energy, vmag2, rmag  !! Variables used in GR calculation
      integer(I4B), dimension(:), allocatable      :: iflag
      real(DP), dimension(:), allocatable          :: dtp


      associate(ntp    => self%nbody, &
         xh     => self%xh, &
         vh     => self%vh, &
         status => self%status,&
         mu     => self%mu)
         
         if (ntp == 0) return
         allocate(iflag(ntp))
         iflag(:) = 0
         allocate(dtp(ntp))
         if (config%lgr) then
            do i = 1,ntp
               rmag = norm2(xh(:, i))
               vmag2 = dot_product(vh(:, i), vh(:, i))
               energy = 0.5_DP * vmag2 - cb%Gmass / rmag
               dtp(i) = dt * (1.0_DP + 3 * config%inv_c2 * energy)
            end do
         else
            dtp(:) = dt
         end if 
         do concurrent(i = 1:ntp, status(i) == ACTIVE)
            call drift_one(mu(i), xh(1,i), xh(2,i), xh(3,i), &
                                  vh(1,i), vh(2,i), vh(3,i), &
                                  dtp(i), iflag(i))
         end do
         if (any(iflag(1:ntp) /= 0)) then
            where(iflag(1:ntp) /= 0) status(1:ntp) = DISCARDED_DRIFTERR
            do i = 1, ntp
               if (iflag(i) /= 0) write(*, *) "Particle ", self%name(i), " lost due to error in Danby drift"
            end do
         end if
      end associate

      return

      end subroutine whm_drift_tp
end submodule whm_drift
