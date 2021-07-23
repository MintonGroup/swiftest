submodule (helio_classes) s_helio_drift
   use swiftest
contains
   module subroutine helio_drift_body(self, system, param, dt, mask)

      !! author: David A. Minton
      !!
      !! Loop through bodies and call Danby drift routine on democratic heliocentric coordinates
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_drift.f90
      !! Adapted from Hal Levison's Swift routine drift.f
      implicit none
      ! Arguments
      class(swiftest_body),         intent(inout) :: self   !! Swiftest body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical, dimension(:),        intent(in)    :: mask   !! Logical mask of size self%nbody that determines which bodies to drift.
      ! Internals
      integer(I4B) :: i !! Loop counter
      real(DP) :: rmag, vmag2, energy
      integer(I4B), dimension(:),allocatable :: iflag !! Vectorized error code flag
      real(DP), dimension(:), allocatable    :: dtp, mu

      associate(pl => self, npl => self%nbody, cb => system%cb)
         if (npl == 0) return

         allocate(iflag(npl))
         iflag(:) = 0
         allocate(dtp(npl))
         allocate(mu(npl))
         mu(:) =  cb%Gmass

         if (param%lgr) then
            do concurrent(i = 1:npl, mask(i))
               rmag = norm2(pl%xh(:, i))
               vmag2 = dot_product(pl%vb(:, i), pl%vb(:, i))
               energy = 0.5_DP * vmag2 - mu(i) / rmag
               dtp(i) = dt * (1.0_DP + 3 * param%inv_c2 * energy)
            end do
         else
            where(mask(1:npl)) dtp(1:npl) = dt
         end if 

         do concurrent(i = 1:npl, mask(i))
            call drift_one(mu(i), pl%xh(1,i), pl%xh(2,i), pl%xh(3,i), &
                                  pl%vb(1,i), pl%vb(2,i), pl%vb(3,i), &
                                  dtp(i), iflag(i))
         end do
         if (any(iflag(1:npl) /= 0)) then
            do i = 1, npl
               if (iflag(i) /= 0) then
                  write(*, *) " Body", self%id(i), " is lost!!!!!!!!!!"
                  write(*, *) pl%xh(:,i)
                  write(*, *) pl%vb(:,i)
                  write(*, *) " stopping "
                  call util_exit(FAILURE)
               end if
            end do
         end if
      end associate

      return

   end subroutine helio_drift_body

   module subroutine helio_drift_pl(self, system, param, dt, mask)
      !! author: David A. Minton
      !!
      !! Wrapper function used to call the body drift routine from a helio_pl structure
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical, dimension(:),        intent(in)    :: mask   !! Logical mask of size self%nbody that determines which bodies to drift.

      call helio_drift_body(self, system, param, dt, mask)
      return
   end subroutine helio_drift_pl

   module subroutine helio_drift_tp(self, system, param, dt, mask)
      !! author: David A. Minton
      !!
      !! Wrapper function used to call the body drift routine from a helio_pl structure
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self   !! Helio massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize
      logical, dimension(:),        intent(in)    :: mask   !! Logical mask of size self%nbody that determines which bodies to drift.

      call helio_drift_body(self, system, param, dt, mask)
      return
   end subroutine helio_drift_tp
   
   module subroutine helio_drift_linear_pl(self, cb, dt, lbeg)
      !! author: David A. Minton
      !!
      !! Perform linear drift of massive bodies due to barycentric momentum of Sun
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift.f90
      !! Adapted from Hal Levison's Swift routine helio_lindrift.f
      implicit none
      ! Arguments
      class(helio_pl), intent(inout) :: self   !! Helio massive body object
      class(helio_cb), intent(inout) :: cb   !! Helio central bod
      real(DP),        intent(in)    :: dt     !! Stepsize
      logical,         intent(in)    :: lbeg   !! Argument that determines whether or not this is the beginning or end of the step
      ! Internals
      real(DP), dimension(NDIM)      :: pt     !! negative barycentric velocity of the central body

      associate(pl => self, npl => self%nbody)
         pt(1) = sum(pl%Gmass(1:npl) * pl%vb(1,1:npl))
         pt(2) = sum(pl%Gmass(1:npl) * pl%vb(2,1:npl))
         pt(3) = sum(pl%Gmass(1:npl) * pl%vb(3,1:npl))
         pt(:) = pt(:) / cb%Gmass
         pl%xh(1,1:npl) = pl%xh(1,1:npl) + pt(1) * dt
         pl%xh(2,1:npl) = pl%xh(2,1:npl) + pt(2) * dt
         pl%xh(3,1:npl) = pl%xh(3,1:npl) + pt(3) * dt

         if (lbeg) then
            cb%ptbeg = pt(:)
         else
            cb%ptend = pt(:) 
         end if
      end associate
   
      return
   end subroutine helio_drift_linear_pl

   module subroutine helio_drift_linear_tp(self, cb, dt, lbeg)
      !! author: David A. Minton
      !!
      !! Perform linear drift of test particles due to barycentric momentum of Sun
      !! New vectorized version included
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_lindrift_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_lindrift_tp.f
      implicit none
      ! Arguments
      class(helio_tp), intent(inout) :: self !! Helio test particleb object
      class(helio_cb), intent(in)    :: cb   !! Helio central body
      real(DP),        intent(in)    :: dt   !! Stepsize
      logical,         intent(in)    :: lbeg !! Argument that determines whether or not this is the beginning or end of the step
      ! Internals
      real(DP), dimension(NDIM)      :: pt     !! negative barycentric velocity of the central body
        
      associate(tp => self, ntp => self%nbody)
         if (lbeg) then
            pt(:) = cb%ptbeg
         else
            pt(:) = cb%ptend
         end if
         where (tp%status(1:ntp) == ACTIVE)
            tp%xh(1, 1:ntp) = tp%xh(1, 1:ntp) + pt(1) * dt
            tp%xh(2, 1:ntp) = tp%xh(2, 1:ntp) + pt(2) * dt
            tp%xh(3, 1:ntp) = tp%xh(3, 1:ntp) + pt(3) * dt
         end where
      end associate
   
      return
   end subroutine helio_drift_linear_tp

end submodule s_helio_drift
