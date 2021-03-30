submodule(helio_classes) s_helio_step
contains
   module subroutine helio_step_system(cb, pl, tp, config)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
      use swiftest
      implicit none
      !! Arguments
      class(helio_cb),               intent(inout) :: cb      !! WHM central body object  
      class(helio_pl),               intent(inout) :: pl      !! WHM central body object  
      class(helio_tp),               intent(inout) :: tp      !! WHM central body object  
      class(swiftest_configuration), intent(in)    :: config  !! Input collection of on parameters 

      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt)
         call pl%set_rhill(cb)
         call tp%set_beg_end(xbeg = pl%xh)
         call pl%step(cb, config, t, dt)
         if (ntp > 0) then
            call tp%set_beg_end(xend = pl%xh)
            call tp%step(cb, pl, config, t, dt)
         end if
      end associate
   end subroutine helio_step_system 

   module subroutine helio_step_pl(self, cb, config, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies ahead Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
      !! Adapted from Hal Levison's Swift routine helio_step_pl.f
      use swiftest
      implicit none
      ! Arguments
      class(helio_pl),               intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! WHM central body particle data structure
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals 
      integer(I4B)     :: i,npl
      real(DP)         :: dth, msys
      real(DP), dimension(NDIM) :: ptbeg, ptend !! TODO: Incorporate these into the tp structure
   
      npl = self%nbody
      dth = 0.5_DP * dt
      if (lfirst) then
         call self%vh2vb()
         lfirst = .false.
      end if
      call self%lindrift(dth, ptbeg)
      call self%getacch(config, t)
      call self%kickvb(dth)
      call self%drift(dt)
      call self%getacch(config, t + dt)
      call self%kickvb(dth)
      call self%lindrift(dth, ptend)
      call self%vb2vh()
   
      return
   
   end subroutine helio_step_pl

   module subroutine helio_step_tp(self, cb, pl, config, t, dt)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_step_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_step_tp.f
      use swiftest
      implicit none
      ! Arguments
      class(helio_tp),                 intent(inout) :: self !! Helio test particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structure
      class(whm_pl),                 intent(inout) :: pl     !! WHM massive body data structure
      class(swiftest_configuration), intent(in)    :: config !! Input collection of 
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
         call self%vh2vb(vs = -self%ptbeg)
         lfirst = .false.
      end if
      call self%lindrift(dth, self%ptbeg)
      call self%getacch(config, t, pl, self%xbeg)
      call self%kickvb(dth)
      call self%drift(mu, dt)
      call self%getacch(config, t + dt, pl, self%xend)
      call self%kickvb(dth)
      call self%lindrift(dth, self%ptend)
      call self%vb2vh(vs = -self%ptend)
   
      return
   
   end subroutine helio_step_tp

end submodule s_helio_step
