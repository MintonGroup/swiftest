submodule (rmvs_classes) s_rmvs_interp
   use swiftest
contains
   module subroutine rmvs_interp_in(self, cb, dt)
      !! author: David A. Minton
      !!
      !! Interpolate planet positions between two Keplerian orbits in inner encounter regio
      !!
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_interp_in.f90 
      !!
      !! Adapted from Hal Levison's Swift routine rmvs3_interp.f
      implicit none
      ! Arguments
      class(rmvs_pl), intent(inout)   :: self !! RMVS test particle object
      class(rmvs_cb), intent(in)      :: cb   !! RMVS central body particle type
      real(DP), intent(in)            :: dt   !! Step size
      ! Internals
      integer(I4B)                    :: i, j, iflag
      real(DP)                        :: dti, frac, dntphenc
      real(DP), dimension(NDIM)       :: xtmp, vtmp

      dntphenc = real(NTPHENC, kind=DP)
      associate (msun => cb%Gmass, npl => self%nbody)
         dti = dt / dntphenc
         do i = 1, npl
            xtmp(:) = self%xin(:, i, 0)
            vtmp(:) = self%vin(:, i, 0)
            do j = 1, NTPHENC - 1
               call drift_one(msun, xtmp(:), vtmp(:), dti, iflag)
               if (iflag /= 0) then
                  write(*, *) " Planet ", self%name(i), " is lost!!!!!!!!!!"
                  write(*, *) msun, dti
                  write(*, *) xtmp(:)
                  write(*, *) vtmp(:)
                  write(*, *) " STOPPING "
                  call util_exit(failure)
               end if
               frac = 1.0_DP - j / dntphenc
               self%xin(:, i, j) = frac * xtmp(:)
               self%vin(:, i, j) = frac * vtmp(:)
            end do
            xtmp(:) = self%xin(:, i, NTPHENC)
            vtmp(:) = self%vin(:, i, NTPHENC)
            do j = NTPHENC - 1, 1, -1
               call drift_one(msun, xtmp(:), vtmp(:), -dti, iflag)
               if (iflag /= 0) then
                  write(*, *) " Planet ", self%name(i), " is lost!!!!!!!!!!"
                  write(*, *) msun, -dti
                  write(*, *) xtmp(:)
                  write(*, *) vtmp(:)
                  write(*, *) " STOPPING "
                  call util_exit(failure)
               end if
               frac = j / dntphenc
               self%xin(:, i, j) = self%xin(:, i, j) + frac * xtmp(:)
               self%vin(:, i, j) = self%vin(:, i, j) + frac * vtmp(:)
            end do
         end do
      end associate
      return

   end subroutine rmvs_interp_in

   module subroutine rmvs_interp_out(self, cb, dt)
      !! author: David A. Minton
      !!
      !! Interpolate planet positions between two Keplerian orbits in outer encounter region
      !!
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_interp_out.f90 
      !!
      !! Adapted from Hal Levison's Swift routine rmvs3_interp.f
      implicit none
      ! Arguments
      class(rmvs_pl), intent(inout)   :: self !! RMVS test particle object
      class(rmvs_cb), intent(in)      :: cb   !! RMVS central body particle type
      real(DP), intent(in)            :: dt   !! Step size
      ! Internals
      integer(I4B)                    :: i, j, iflag
      real(DP)                        :: dto, frac, dntenc
      real(DP), dimension(NDIM)       :: xtmp, vtmp

   ! executable code
      dntenc = real(NTENC, DP)
      associate (msun => cb%Gmass, npl => self%nbody)
         dto = dt / dntenc
         do i = 1, npl
            xtmp(:) = self%xout(:, i, 0)
            vtmp(:) = self%vout(:, i, 0)
            do j = 1, NTENC - 1
               call drift_one(msun, xtmp(:), vtmp(:), dto, iflag)
               if (iflag /= 0) then
                  write(*, *) " Planet ", self%name(i), " is lost!!!!!!!!!!"
                  write(*, *) msun, dto
                  write(*, *) xtmp(:)
                  write(*, *) vtmp(:)
                  write(*, *) " STOPPING "
                  call util_exit(FAILURE)
               end if
               frac = 1.0_DP - j / dntenc
               self%xout(:, i, j) = frac*xtmp(:)
               self%vout(:, i, j) = frac*vtmp(:)
            end do
            xtmp(:) = self%xout(:, i, NTENC)
            vtmp(:) = self%vout(:, i, NTENC)
            do j = NTENC - 1, 1, -1
               call drift_one(msun, xtmp(:), vtmp(:), -dto, iflag)
               if (iflag /= 0) then
                  write(*, *) " Planet ", self%name(i), " is lost!!!!!!!!!!"
                  write(*, *) msun, -dto
                  write(*, *) xtmp(:)
                  write(*, *) vtmp(:)
                  write(*, *) " STOPPING "
                  call util_exit(FAILURE)
               end if
               frac = j / dntenc
               self%xout(:, i, j) = self%xout(:, i, j) + frac * xtmp(:)
               self%vout(:, i, j) = self%vout(:, i, j) + frac * vtmp(:)
            end do
         end do
      end associate

      return

   end subroutine rmvs_interp_out   
end submodule s_rmvs_interp