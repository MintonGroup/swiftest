submodule(whm_classes) whm_drift
   use swiftest
contains
   module subroutine whm_drift_pl(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Loop through planets and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift.f90
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_nbody_system),  intent(inout) :: system !! WHM nbody system object
      class(swiftest_parameters),    intent(in)    :: param  !! Current run configuration parameters of 
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B)                              :: i
      real(DP)                                  :: energy, vmag2, rmag  !! Variables used in GR calculation
      integer(I4B), dimension(:), allocatable   :: iflag
      real(DP), dimension(:), allocatable       :: dtp

      associate(pl => self, npl => self%nbody)

         if (npl == 0) return

         allocate(iflag(npl))
         iflag(:) = 0
         allocate(dtp(npl))

         if (param%lgr) then
            do i = 1,npl
               rmag = norm2(pl%xj(:, i))
               vmag2 = dot_product(pl%vj(:, i), pl%vj(:, i))
               energy = 0.5_DP * vmag2 - pl%muj(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * param%inv_c2 * energy)
            end do
         else
            dtp(:) = dt
         end if 

         call drift_one(pl%muj(1:npl), pl%xj(1,1:npl), pl%xj(2,1:npl), pl%xj(3,1:npl), &
                                   pl%vj(1,1:npl), pl%vj(2,1:npl), pl%vj(3,1:npl), &
                                   dtp(1:npl), iflag(1:npl))
         if (any(iflag(1:npl) /= 0)) then
            do i = 1, npl
               if (iflag(i) /= 0) then
                  write(*, *) " Planet ", self%name(i), " is lost!!!!!!!!!!"
                  write(*, *) pl%xj(:,i)
                  write(*, *) pl%vj(:,i)
                  write(*, *) " stopping "
                  call util_exit(FAILURE)
               end if
            end do
         end if
      end associate

      return

   end subroutine whm_drift_pl

   module subroutine whm_drift_tp(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Loop through test particles and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift_tp.f 
      !! Includes 
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters of 
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B)                              :: i   
      real(DP)                                  :: energy, vmag2, rmag  !! Variables used in GR calculation
      integer(I4B), dimension(:), allocatable   :: iflag
      real(DP), dimension(:), allocatable       :: dtp

      associate(tp => self, ntp => self%nbody)
         if (ntp == 0) return
         allocate(iflag(ntp))
         iflag(:) = 0
         allocate(dtp(ntp))
         if (param%lgr) then
            do i = 1, ntp
               rmag = norm2(tp%xh(:, i))
               vmag2 = dot_product(tp%vh(:, i), tp%vh(:, i))
               energy = 0.5_DP * vmag2 - tp%mu(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * param%inv_c2 * energy)
            end do
         else
            dtp(:) = dt
         end if 
         do concurrent(i = 1:ntp, tp%status(i) == ACTIVE)
            call drift_one(tp%mu(i), tp%xh(1,i), tp%xh(2,i), tp%xh(3,i), &
                                     tp%vh(1,i), tp%vh(2,i), tp%vh(3,i), &
                                     dtp(i), iflag(i))
         end do
         if (any(iflag(1:ntp) /= 0)) then
            where(iflag(1:ntp) /= 0) tp%status(1:ntp) = DISCARDED_DRIFTERR
            do i = 1, ntp
               if (iflag(i) /= 0) write(*, *) "Particle ", self%name(i), " lost due to error in Danby drift"
            end do
         end if
      end associate

      return
      end subroutine whm_drift_tp
end submodule whm_drift
