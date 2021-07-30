submodule (helio_classes) s_helio_drift
   use swiftest
contains

   module subroutine helio_drift_body(self, system, param, dt)
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
      ! Internals
      integer(I4B) :: i !! Loop counter
      real(DP) :: rmag, vmag2, energy
      integer(I4B), dimension(:),allocatable :: iflag !! Vectorized error code flag
      real(DP), dimension(:), allocatable    :: dtp, mu

      if (self%nbody == 0) return

      associate(n => self%nbody)
         allocate(iflag(n))
         iflag(:) = 0
         allocate(mu(n))
         mu(:) = system%cb%Gmass
         call drift_all(mu, self%xh, self%vb, self%nbody, param, dt, self%lmask, iflag)
         if (any(iflag(1:n) /= 0)) then
            where(iflag(1:n) /= 0) self%status(1:n) = DISCARDED_DRIFTERR
            do i = 1, n
               if (iflag(i) /= 0) write(*, *) " Body ", self%id(i), " lost due to error in Danby drift"
            end do
         end if
      end associate

      return
   end subroutine helio_drift_body


   module subroutine helio_drift_pl(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Wrapper function used to call the body drift routine from a helio_pl structure
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize

      call helio_drift_body(self, system, param, dt)

      return
   end subroutine helio_drift_pl


   module subroutine helio_drift_tp(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Wrapper function used to call the body drift routine from a helio_pl structure
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self   !! Helio massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize

      call helio_drift_body(self, system, param, dt)

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
      class(helio_pl),               intent(inout) :: self !! Helio massive body object
      class(helio_cb),               intent(inout) :: cb   !! Helio central body
      real(DP),                      intent(in)    :: dt   !! Stepsize
      logical,                       intent(in)    :: lbeg !! Argument that determines whether or not this is the beginning or end of the step
      ! Internals
      real(DP), dimension(NDIM) :: pt     !! negative barycentric velocity of the central body
      integer(I4B)              :: i    

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         if (npl == 0) return
         pt(1) = sum(pl%Gmass(1:npl) * pl%vb(1,1:npl), self%lmask(1:npl))
         pt(2) = sum(pl%Gmass(1:npl) * pl%vb(2,1:npl), self%lmask(1:npl))
         pt(3) = sum(pl%Gmass(1:npl) * pl%vb(3,1:npl), self%lmask(1:npl))
         pt(:) = pt(:) / cb%Gmass
         do concurrent(i = 1:npl, self%lmask(i))
            pl%xh(:,i) = pl%xh(:,i) + pt(:) * dt
         end do

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
      class(helio_tp),               intent(inout) :: self !! Helio test particleb object
      class(helio_cb),               intent(in)    :: cb   !! Helio central body
      real(DP),                      intent(in)    :: dt   !! Stepsize
      logical,                       intent(in)    :: lbeg !! Argument that determines whether or not this is the beginning or end of the step
      ! Internals
      real(DP), dimension(NDIM)      :: pt     !! negative barycentric velocity of the central body
        
      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         if (ntp == 0) return
         if (lbeg) then
            pt(:) = cb%ptbeg
         else
            pt(:) = cb%ptend
         end if
         where (self%lmask(1:ntp))
            tp%xh(1, 1:ntp) = tp%xh(1, 1:ntp) + pt(1) * dt
            tp%xh(2, 1:ntp) = tp%xh(2, 1:ntp) + pt(2) * dt
            tp%xh(3, 1:ntp) = tp%xh(3, 1:ntp) + pt(3) * dt
         end where
      end associate
   
      return
   end subroutine helio_drift_linear_tp

end submodule s_helio_drift
