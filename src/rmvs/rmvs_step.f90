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
      integer(I4B)                                      :: i, j, k
      real(DP)                                          :: dto, time

   ! executable code
      associate(pl => self, npl => self%nbody, ntp => tp%nbody, t => config%t, xht => tp%xh, vht => tp%vh)
         dto = dt / NTENC
         where(tp%plencP(:) == 0)
            tp%status(:) = INACTIVE
         elsewhere  
            tp%lperi(:) = .false.
         end where
         do i = 1, NTENC
            time = t + (i - 1) * dto
            call rmvs_step_out2(cb, pl, tp, time, dto, i, config)
            do j = 1, npl
               if (pl%nenc(j) == 0) cycle
               where((tp%plencP(:) == j) .and. (tp%status(:) == INACTIVE)) tp%status(:) = ACTIVE
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
      integer(I4B)                                   :: i, j
      real(DP)                                       :: dti, time
      real(DP), dimension(NDIM, self%nbody)          :: xbeg, xend, vbeg

      dti = dt / NTPHENC

      associate(pl => self, npl => self%nbody, xht => tp%xh, vht => tp%vh, plind => self%plind, nenc => self%nenc)
         if (config%loblatecb) call pl%obl_acc_in(cb)
         call pl%make_planetocentric(cb, tp, config)
         do i = 1, npl
            if (nenc(i) == 0) cycle
            ! There are inner encounters with this planet...switch to planetocentric coordinates to proceed
            time = t
            call pl%tpenc(i)%peri_pass(cb, pl, time, dti, .true., 0, i, config) 
            ! now step the encountering test particles fully through the inner encounter
            lfirsttp = .true.
            associate(index => pl%tpenc(i)%index, &
               xpc => self%tpenc(i)%xh, vpc => self%tpenc(i)%vh, apc => self%tpenc(i)%ah)
               pl%tpenc(i)%lfirst = .true.
               do index = 1, NTPHENC ! Integrate over the encounter region, using the "substitute" planetocentric systems at each level
                  xbeg(:,1) = cb%xin(:) - pl%xin(:, i, index - 1) 
                  xend(:,1) = cb%xin(:) - pl%xin(:, i, index)  
                  vbeg(:,1) = cb%vin(:) - pl%vin(:, i, index - 1)
                  do j = 1, NDIM
                     xbeg(j,2:npl) = pl%xin(j, plind(i,2:npl), index - 1) - pl%xin(j, i, index - 1)  
                     xend(j,2:npl) = pl%xin(j, plind(i,2:npl), index) - pl%xin(j, i, index)  
                     vbeg(j,2:npl) = pl%vin(j, plind(i,2:npl), index - 1) - pl%vin(j, i, index - 1)
                  end do
                  pl%plenc(i)%xh(:,:) = xbeg(:,:)
                  pl%plenc(i)%vh(:,:) = vbeg(:,:)
                  call pl%tpenc(i)%set_beg_end(xbeg = xbeg, xend = xend)
                  call pl%tpenc(i)%step(pl%cbenc(i), pl%plenc(i), config, time, dti)
                  do j = 1, NDIM
                     pl%tpenc(i)%xheliocen(j, :) = pl%tpenc(i)%xh(j, :) + pl%xin(j, i, index)
                  end do
                  time = config%t + j * dti
                  call pl%tpenc(i)%peri_pass(cb, pl, time, dti, .false., index, i, config) 
               end do
               where(pl%tpenc(i)%status(:) == ACTIVE) pl%tpenc(i)%status(:) = INACTIVE
            end associate
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
      integer(I4B)                                   :: i, j
      logical, dimension(:), allocatable             :: encmask

      associate(pl => self, npl => self%nbody, nenc => self%nenc, tpenc => self%tpenc, cbenc => self%cbenc, &
         plenc => self%plenc, GMpl => self%Gmass)

         do i = 1, npl
            if (nenc(i) == 0) cycle 
            ! There are inner encounters with this planet
            if (allocated(encmask)) deallocate(encmask)
            allocate(encmask(nenc(i)))
            encmask(:) = tp%plencP(:) == i

            ! Save the index value of the planet corresponding to this encounter 
            tpenc(i)%ipleP = i
            tpenc(i)%lplanetocentric = .true.
            tpenc(i)%nbody = nenc(i)
            !call tpenc(i)%setup(nenc(i))
            ! Create space for the heliocentric position values for those acceleration calculations that need them
            allocate(tpenc(i)%xheliocen(NDIM, nenc(i)))
            allocate(tpenc(i)%xh(NDIM, nenc(i)))
            allocate(tpenc(i)%vh(NDIM, nenc(i)))
            allocate(tpenc(i)%ah(NDIM, nenc(i)))
            allocate(tpenc(i)%status(nenc(i)))
            allocate(tpenc(i)%mu(nenc(i)))
            allocate(tpenc(i)%isperi(nenc(i)))
            allocate(tpenc(i)%name(nenc(i)))
            allocate(tpenc(i)%plperP(nenc(i)))
            allocate(tpenc(i)%lperi(nenc(i)))
            allocate(tpenc(i)%peri(nenc(i)))
            tpenc(i)%status(:) = ACTIVE
            tpenc(i)%lperi(:) = .false.
            tpenc(i)%plperP(:) = 0
            tpenc(i)%name(:) = pack(tp%name(:), encmask(:)) 
            !call tp%spill(tpenc(i), pl%encmask(:,i))
            ! Grab all the encountering test particles and convert them to a planetocentric frame
            do j = 1, NDIM 
               tpenc(i)%xheliocen(j, :) = pack(tp%xh(j,:), encmask(:)) 
               tpenc(i)%xh(j, :) = tpenc(i)%xheliocen(j, :) - pl%xin(j, i, 0)
               tpenc(i)%vh(j, :) = pack(tp%vh(j,:), encmask(:)) - pl%vin(j, i, 0)
            end do

            ! Make sure that the test particles get the planetocentric value of mu 
            call tpenc(i)%set_mu(pl%cbenc(i)) 

            ! Save the heliocentric position values of the encountering planet
            !allocate(tpenc(i)%xh_pl(NDIM,0:NTPHENC))
            ! Save the encountering planet's values of oblateness acceleration 
            if (config%loblatecb) then
               allocate(tpenc(i)%aoblin_pl(NDIM,0:NTPHENC))
               tpenc(i)%aoblin_pl(:,:) = pl%aoblin(:,i,0:NTPHENC)
            end if
         end do
      end associate
      return
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
      integer(I4B), dimension(:), allocatable :: tpind
      logical, dimension(:), allocatable :: encmask

      associate(pl => self, nenc => self%nenc, npl => self%nbody, ntp => tp%nbody, name => tp%name, &
         tpenc => self%tpenc, plenc => self%plenc, cbenc => self%cbenc)

         do i = 1, npl
            if (nenc(i) == 0) cycle
            allocate(tpind(nenc(i)))
            ! Index array of encountering test particles
            if (allocated(encmask)) deallocate(encmask)
            allocate(encmask(nenc(i)))
            encmask(:) = tp%plencP(:) == i
            tpind(:) = pack([(j,j=1,ntp)], encmask(:))

            ! Copy the results of the integration back over and shift back to heliocentric reference
            tp%status(tpind(1:nenc(i))) = tpenc(i)%status(1:nenc(i)) 
            do j = 1, NDIM
               tp%xh(j, tpind(1:nenc(i))) = tpenc(i)%xh(j,1:nenc(i)) + pl%xin(j, i, NTPHENC) 
               tp%vh(j, tpind(1:nenc(i))) = tpenc(i)%vh(j,1:nenc(i)) + pl%vin(j, i, NTPHENC)
            end do
            deallocate(tpenc(i)%xheliocen)
            deallocate(tpenc(i)%xh)
            deallocate(tpenc(i)%vh)
            deallocate(tpenc(i)%ah)
            deallocate(tpenc(i)%status)
            deallocate(tpenc(i)%mu)
            deallocate(tpenc(i)%isperi)
            deallocate(tpenc(i)%name)
            deallocate(tpenc(i)%plperP)
            deallocate(tpenc(i)%lperi)
            deallocate(tpenc(i)%peri)
            deallocate(tpind)
         end do
      end associate

      return
   end subroutine rmvs_step_end_planetocentric
   
end submodule s_rmvs_step