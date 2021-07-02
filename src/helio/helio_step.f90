submodule(helio_classes) s_helio_step
   use swiftest
contains
   module subroutine helio_step_system(self, param)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
      implicit none
      ! Arguments
      class(helio_nbody_system),     intent(inout) :: self    !! Helio nbody system object
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters of on parameters 

      select type(cb => self%cb)
      class is (helio_cb)
      select type(pl => self%pl)
      class is (helio_pl)
      select type(tp => self%tp)
      class is (helio_tp)
      associate(ntp => tp%nbody, npl => pl%nbody, t => param%t, dt => param%dt)
         call pl%set_rhill(cb)
         call tp%set_beg_end(xbeg = pl%xh)
         call pl%step(cb, param, t, dt)
         if (ntp > 0) then
            call tp%set_beg_end(xend = pl%xh)
            call tp%step(cb, pl, param, t, dt)
         end if
      end associate
      end select
      end select
      end select
      return
   end subroutine helio_step_system 

   module subroutine helio_step_pl(self, cb, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies ahead Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
      !! Adapted from Hal Levison's Swift routine helio_step_pl.f
      implicit none
      ! Arguments
      class(helio_pl),               intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Helio central body particle data structure
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of 
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals 
      integer(I4B)     :: i,npl
      real(DP)         :: dth, msys
      real(DP), dimension(NDIM) :: ptbeg, ptend !! TODO: Incorporate these into the tp structure
      logical, save :: lfirst = .true.
   
      npl = self%nbody
      dth = 0.5_DP * dt
      if (lfirst) then
         call self%vh2vb(cb)
         lfirst = .false.
      end if
      call self%lindrift(cb, dth, ptbeg)
      call self%getacch(cb, param, t)
      call self%kickvb(dth)

      call self%drift(cb, param, dt)
      call self%getacch(cb, param, t + dt)
      call self%kickvb(dth)
      call self%lindrift(cb, dth, ptend)
      call self%vb2vh(cb)
   
      return
   
   end subroutine helio_step_pl

   module subroutine helio_step_tp(self, cb, pl, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_step_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_step_tp.f
      implicit none
      ! Arguments
      class(helio_tp),                 intent(inout) :: self !! Helio test particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structure
      class(whm_pl),                 intent(inout) :: pl     !! WHM massive body data structure
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of 
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals
      logical, save  :: lfirst = .true. !! Flag to indicate that this is the first call
      real(DP) :: dth   !! Half step size
      real(DP) :: mu    !! Central mass term
   
   ! executable code
      dth = 0.5_DP * dt
      mu = cb%Gmass
      if (lfirst) then
         call self%vh2vb(vbcb = -self%ptbeg)
         lfirst = .false.
      end if
      call self%lindrift(dth, self%ptbeg)
      call self%getacch(cb, pl, param, t, self%xbeg)
      call self%kickvb(dth)
      call self%drift(cb, param, dt)
      call self%getacch(cb, pl, param, t + dt, self%xend)
      call self%kickvb(dth)
      call self%lindrift(dth, self%ptend)
      call self%vb2vh(vbcb = -self%ptend)
   
      return
   
   end subroutine helio_step_tp

end submodule s_helio_step
