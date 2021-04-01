submodule(rmvs_classes) s_rmvs_step
   use swiftest
contains
   module subroutine rmvs_step_system(cb, pl, tp, config)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step.f90
      implicit none
      ! Arguments
      class(rmvs_cb),                    intent(inout)  :: cb      !! RMVS central body object  
      class(rmvs_pl),                    intent(inout)  :: pl      !! RMVS central body object  
      class(rmvs_tp),                    intent(inout)  :: tp      !! RMVS central body object  
      class(swiftest_configuration),     intent(in)     :: config  !! Input collection of  configuration parameters 
      ! Internals
      logical :: lencounter, lfirstpl, lfirsttp 
      real(DP) :: rts
      real(DP), dimension(:,:), allocatable :: xbeg, xend, vbeg
      integer(I4B) :: i
 
      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt, &
         xh => pl%xh, vh => pl%vh, xj => pl%xj, vj => pl%vj, ah => pl%ah,  eta => pl%eta, & ! These two lines of associations aid in debugging with gdb
         xht => tp%xh, vht => tp%vh)
         allocate(xbeg, source=pl%xh)
         allocate(vbeg, source=pl%vh)
         call pl%set_rhill(cb)
         call tp%set_beg_end(xbeg = xbeg, vbeg = vbeg)
         ! ****** Check for close encounters ***** !
         rts = RHSCALE
         lencounter = tp%encounter_check(cb, pl, dt, rts)
         if (lencounter) then
            lfirstpl = pl%lfirst
            lfirsttp = tp%lfirst
            pl%xout(:,:,0) = tp%xbeg(:,:)
            pl%vout(:,:,0) = tp%vbeg(:,:)
            call pl%step(cb, config, t, dt) 
            pl%xout(:,:,NTENC) = pl%xh(:,:)
            pl%vout(:,:,NTENC) = pl%vh(:,:)
            call tp%set_beg_end(xend = pl%xh)
            call pl%interp_out(cb, dt)
            call pl%step_out(cb, tp, dt, config)
            call tp%reverse_status()
            call tp%set_beg_end(xbeg = xbeg, xend = xend)
            tp%lfirst = .true.
            call tp%step(cb, pl, config, t, dt)
            where (tp%status(:) == INACTIVE) tp%status(:) = ACTIVE
            pl%lfirst = lfirstpl
            tp%lfirst = lfirsttp
         else
            call whm_step_system(cb, pl, tp, config)
         end if
      end associate

   end subroutine rmvs_step_system 

   module subroutine rmvs_step_out(self, cb, tp, dt, config)
      !! author: David A. Minton
      !!
      !! Step ACTIVE test particles ahead in the outer encounter region
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step_out.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_out.f90
      implicit none
      ! Arguments
      class(rmvs_pl),                    intent(inout)  :: self !! RMVS massive body object
      class(rmvs_cb),                    intent(inout)  :: cb   !! RMVS central body object
      class(rmvs_tp),                    intent(inout)  :: tp   !! RMVS test particle object
      real(DP),                          intent(in)     :: dt   !! Step size
      class(swiftest_configuration),     intent(in)     :: config  !! Input collection of  configuration parameters
      ! Internals
      integer(I4B)                                      :: i, j, k, nenc, itpp
      real(DP)                                          :: dto, time

   ! executable code
      associate(pl => self, npl => self%nbody, ntp => tp%nbody, t => config%t, xht => tp%xh, vht => tp%vh, &
         status => tp%status)
         dto = dt / NTENC
         where(tp%plencp(:) == 0)
            tp%status(:) = INACTIVE
         elsewhere  
            tp%lperi(:) = .false.
         end where
         do i = 1, NTENC
            time = t + (i - 1) * dto
            call rmvs_step_out2(cb, pl, tp, time, dto, i, config)
            do j = 1, npl
               nenc = pl%nenc(j) 
               if (nenc > 0) then
                  itpp = pl%tpenc1P(j)
                  do k = 1, nenc
                     if (tp%status(itpp) == INACTIVE) tp%status(itpp) = ACTIVE
                     itpp = tp%tpencp(itpp)
                  end do
               end if
            end do
         end do
      end associate
      return

   end subroutine rmvs_step_out   

   module subroutine rmvs_step_out2(cb, pl, tp, t, dt, index, config)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead in the outer encounter region, setting up and calling the inner region
      !!    integration if necessary
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step_out.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_out.f90 
      implicit none
      ! Arguments
      class(rmvs_cb),                 intent(inout)  :: cb   !! RMVS central body object
      class(rmvs_pl),                 intent(inout)  :: pl   !! RMVS massive body object
      class(rmvs_tp),                 intent(inout)  :: tp   !! RMVS test particle object
      real(DP),                       intent(in)     :: t     !! Simulation time
      real(DP),                       intent(in)     :: dt    !! Step size
      integer(I4B),                   intent(in)     :: index !! outer substep number within current set
      class(swiftest_configuration),  intent(in)     :: config  !! Input collection of  configuration parameters
      ! Internals
      logical                                        :: lfirsttp
      real(DP)                                       :: rts
      logical                                        :: lencounter

      associate(xht => tp%xh, vht => tp%vh, xbeg => tp%xbeg, xend => tp%xend, aht => tp%ah,&
         xout1 => pl%xout(:, :, index-1), vout1 => pl%vout(:,:,index-1), xout2 => pl%xout(:,:,index))

         call tp%set_beg_end(xbeg = pl%xout(:, :, index-1), &
                              vbeg = pl%vout(:, :, index-1), &
                              xend = pl%xout(:, :, index))
         rts = RHPSCALE
         lencounter = tp%encounter_check(cb, pl, dt, rts) 
         if (lencounter) then
            ! Interpolate planets in inner encounter region
            pl%xin(:,:,0) = tp%xbeg(:,:)
            pl%vin(:,:,0) = tp%vbeg(:,:)
            pl%xin(:,:,NTPHENC) = tp%xend(:, :)
            pl%vin(:,:,NTPHENC) = pl%vout(:, :, index)
            call pl%interp_in(cb, dt)
            call pl%step_in(cb, tp, config, t, dt)
            lfirsttp = tp%lfirst
            tp%lfirst = .true.
            call tp%step(cb, pl, config, t, dt)
            tp%lfirst = lfirsttp
         else
            call tp%step(cb, pl, config, t, dt)
         end if
      end associate
      return

   end subroutine rmvs_step_out2

   module subroutine rmvs_step_in_pl(self, cb, tp, config, t, dt)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead in the inner encounter region
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step_in.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_in.f90
      implicit none
      ! Arguments
      class(rmvs_pl),                intent(inout)  :: self !! RMVS massive body object
      class(rmvs_cb),                intent(inout)  :: cb   !! RMVS central body object
      class(rmvs_tp),                intent(inout)  :: tp   !! RMVS test particle object
      class(swiftest_configuration), intent(in)     :: config  !! Input collection of configuration parameters 
      real(DP),                      intent(in)     :: t    !! Current time
      real(DP),                      intent(in)     :: dt   !! Step size
      ! Internals
      logical                                        :: lfirsttp
      integer(I4B)                                   :: i, j, nenc
      real(DP)                                       :: dti, time

      dti = dt / NTPHENC
      associate(pl => self, npl => self%nbody, xht => tp%xh, vht => tp%vh)
         if (config%loblatecb) call pl%obl_acc_in(cb)
         call pl%make_planetocentric(cb, tp, config)
         do i = 1, npl
            nenc = pl%nenc(i) 
            if (nenc > 0) then
            ! There are inner encounters with this planet...switch to planetocentric coordinates to proceed
               time = t
               call pl%tpenc(i)%peri_pass(cb, pl, time, dti, .true., 0, nenc, i, config) 
               ! now step the encountering test particles fully through the inner encounter
               lfirsttp = .true.
               associate(index => pl%tpenc(i)%index, &
                  xpc => self%tpenc(i)%xh, vpc => self%tpenc(i)%vh, apc => self%tpenc(i)%ah)
                  pl%tpenc(i)%lfirst = .true.
                  do index = 1, NTPHENC ! Integrate over the encounter region, using the "substitute" planetocentric systems at each level
                     call pl%tpenc(i)%set_beg_end(xbeg = pl%plenc(i,index - 1)%xh, xend = pl%plenc(i, index)%xh)
                     call pl%tpenc(i)%step(pl%cbenc(i), pl%plenc(i, index), config, time, dti)
                     do j = 1, NDIM
                        pl%tpenc(i)%xheliocen(j, :) = pl%tpenc(i)%xh(j, :) + pl%xin(j, i, index)
                     end do
                     time = config%t + j * dti
                     call pl%tpenc(i)%peri_pass(cb, pl, time, dti, .false., index, nenc, i, config) 
                  end do
                  where(pl%tpenc(i)%status(:) == ACTIVE) pl%tpenc(i)%status(:) = INACTIVE
               end associate
            end if
         end do
         call pl%end_planetocentric(cb,tp)

      end associate

      return
   end subroutine rmvs_step_in_pl

   module subroutine rmvs_step_make_planetocentric(self, cb, tp, config)
      !! author: David A. Minton
      !!
      !! When encounters are detected, this method will call the interpolation methods for the planets and 
      !! creates a Swiftest test particle structure for each planet's encountering test particles to simplify the 
      !! planetocentric calculations. This subroutine is not based on an existing one from Swift and Swifter
      !!
      implicit none
      ! Arguments
      class(rmvs_pl),                 intent(inout)  :: self !! RMVS test particle object
      class(rmvs_cb),                 intent(inout)  :: cb   !! RMVS central body particle type
      class(rmvs_tp),                 intent(inout)  :: tp   !! RMVS test particle object
      class(swiftest_configuration),  intent(in)     :: config !! Input collection of configuration parameters 
      ! Internals
      integer(I4B)                                   :: i, j, k
      type(rmvs_pl)                                  :: cb_as_pl, tmp
      logical, dimension(:), allocatable             :: copyflag

      allocate(self%tpenc(self%nbody))
      allocate(self%plenc(self%nbody, 0:NTPHENC))
      allocate(self%cbenc(self%nbody))
      allocate(copyflag(self%nbody))
      associate(pl => self, npl => self%nbody, nenc => self%nenc, tpenc => self%tpenc, cbenc => self%cbenc, &
         plenc => self%plenc)
         ! Indicate that this is a planetocentric close encounter structor

         ! Save the original central body object so that it can be passed down as needed through the planetocentric structures
         call cb_as_pl%setup(1)
         call tmp%setup(1)
         cb_as_pl%name(1) = 0
         cb_as_pl%Gmass(1)  = cb%Gmass
         cb_as_pl%mass(1)   = cb%mass
         cb_as_pl%radius(1) = cb%radius

         do i = 1, npl
            if (nenc(i) > 0) then ! There are inner encounters with this planet
               ! Make the planet a central body for the encounter
               cbenc(i)        = cb
               cbenc(i)%Gmass  = pl%Gmass(i)
               cbenc(i)%mass   = pl%mass(i)
               cbenc(i)%radius = pl%radius(i)

               ! Create space for the heliocentric position values for those acceleration calculations that need them
               allocate(tpenc(i)%xheliocen(NDIM, nenc(i)))

               ! Save the index value of the planet corresponding to this encounter 
               tpenc(i)%ipleP = i
               tpenc(i)%lplanetocentric = .true.
               tpenc(i)%cb = cb 
               call pl%tpenc(i)%setup(nenc(i))
               tpenc(i)%status(:) = ACTIVE
               call tp%spill(pl%tpenc(i), pl%encmask(:,i))
               ! Grab all the encountering test particles and convert them to a planetocentric frame
               do j = 1, nenc(i)
                  tpenc(i)%xheliocen(:, j) = tpenc(i)%xh(:, j)
                  tpenc(i)%xh(:, j)        = tpenc(i)%xh(:, j) - pl%xin(:, i, 0)
                  tpenc(i)%vh(:, j)        = tpenc(i)%vh(:, j) - pl%vin(:, i, 0)
               end do

               ! Make sure that the test particles get the planetocentric value of mu 
               call tpenc(i)%set_mu(pl%cbenc(i)) 

               ! Save the heliocentric position values of the encountering planet
               allocate(tpenc(i)%xh_pl(NDIM,0:NTPHENC))
               tpenc(i)%xh_pl(:,:) = pl%xin(:,i,0:NTPHENC)
               ! Save the encountering planet's values of oblateness acceleration 
               if (config%loblatecb) then
                  allocate(tpenc(i)%aoblin_pl(NDIM,0:NTPHENC))
                  tpenc(i)%aoblin_pl(:,:) = pl%aoblin(:,i,0:NTPHENC)
               end if
               ! Now create a planetocentric "planet" structure containing the *other* planets (plus the Sun) in it at each point along
               ! the innner encounter trajectory of the planet
               do k = 0, NTPHENC
                  call plenc(i, k)%setup(npl)
                  ! Copy all the basic planet parameters and positions
                  copyflag(:) =  .true.
                  call plenc(i, k)%copy(pl,copyflag)
                  ! Give each planet a position and velocity vector that is planetocentric wrt planet i (the encountering planet)
                  do j = 1, npl
                     plenc(i, k)%xh(:, j)  = pl%xin(:, j, k) - pl%xin(:, i, k)
                     plenc(i, k)%vh(:, j)  = pl%vin(:, j, k) - pl%vin(:, i, k)
                  end do
                  cb_as_pl%xh(:, 1)  = cb%xin(:) - pl%xin(:, i, k)
                  cb_as_pl%vh(:, 1)  = cb%vin(:) - pl%vin(:, i, k)
                  ! Pull the encountering body out of the massive body list
                  copyflag(:) = .false.
                  copyflag(i) = .true.
                  call plenc(i, k)%spill(tmp, copyflag)
                  copyflag(:) = .false.
                  copyflag(1) = .true. ! Put the central body as planet 1 like in the original Swifter version
                  call plenc(i, k)%fill(cb_as_pl, copyflag)
               end do
            end if
         end do
      end associate
   end subroutine rmvs_step_make_planetocentric

   module subroutine rmvs_step_end_planetocentric(self, cb, tp)
      !! author: David A. Minton
      !!
      !! Deallocates all of the encountering particle data structures for next time
      !!
      implicit none
      ! Arguments
      class(rmvs_pl),                 intent(inout)  :: self !! RMVS test particle object
      class(rmvs_cb),                 intent(inout)  :: cb      !!  RMVS central body object
      class(rmvs_tp),                 intent(inout)  :: tp   !! RMVS test particle object
      ! Internals
      integer(I4B) :: i, j

      associate(pl => self, nenc => self%nenc, npl => self%nbody, name => tp%name, &
         tpenc => self%tpenc, plenc => self%plenc, cbenc => self%cbenc, encmask => self%encmask)

         do i = 1, npl
            if (nenc(i) == 0) cycle
            do j = 1, nenc(i)
               ! Copy the results of the integration back over
               tpenc(i)%xh(:, j) = tpenc(i)%xh(:, j) + pl%xin(:, i, NTPHENC) 
               tpenc(i)%vh(:, j) = tpenc(i)%vh(:, j) + pl%vin(:, i, NTPHENC) 
               ! Put back heliocentric value of mu
               call tpenc(i)%set_mu(cb)
            end do
            ! Replace the old test particle with the new one
            call tp%fill(tpenc(i), encmask(:,i))
         end do
      end associate
      deallocate(self%tpenc)
      deallocate(self%cbenc)
      deallocate(self%plenc)
      deallocate(self%encmask)

      return

   end subroutine rmvs_step_end_planetocentric
   
end submodule s_rmvs_step