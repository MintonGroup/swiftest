submodule(whm_classes) whm_drift
   use swiftest
contains

   module subroutine whm_drift_pl(self, system, param, dt, mask)
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
      class(swiftest_parameters),    intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                      intent(in)    :: dt     !! Stepsize
      logical, dimension(:),         intent(in)    :: mask   !! Logical mask of size self%nbody that determines which bodies to drift
      ! Internals
      integer(I4B)                              :: i
      integer(I4B), dimension(:), allocatable   :: iflag

      associate(pl => self, npl => self%nbody)
         if (npl == 0) return

         allocate(iflag(npl))
         iflag(:) = 0
         call drift_all(pl%muj, pl%xj, pl%vj, npl, param, dt, mask, iflag)
         if (any(iflag(1:npl) /= 0)) then
            where(iflag(1:npl) /= 0) pl%status(1:npl) = DISCARDED_DRIFTERR
            do i = 1, npl
               if (iflag(i) /= 0) then 
                  write(*, *) " Planet ", pl%id(i), " is lost!!!!!!!!!!!!"
                  WRITE(*, *) pl%muj(i), dt
                  WRITE(*, *) pl%xj(:,i)
                  WRITE(*, *) pl%vj(:,i)
                  WRITE(*, *) " STOPPING "
               end if
            end do
            call util_exit(FAILURE)
         end if
      end associate

      return
   end subroutine whm_drift_pl

end submodule whm_drift
