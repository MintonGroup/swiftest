submodule(rmvs_classes) s_rmvs_step
   use swiftest
contains
   module subroutine rmvs_step_system(self, config)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step.f90
      implicit none
      ! Arguments
      class(rmvs_nbody_system),          intent(inout)  :: self    !! RMVS nbody system object
      class(swiftest_configuration),     intent(in)     :: config  !! Input collection of  configuration parameters 
      ! Internals
      logical :: lencounter, lfirstpl, lfirsttp 
      real(DP) :: rts
      real(DP), dimension(:,:), allocatable :: xbeg, xend, vbeg
      integer(I4B) :: i
 
      select type(cb => self%cb)
      class is (rmvs_cb)
      select type(pl => self%pl)
      class is (rmvs_pl)
      select type(tp => self%tp)
      class is (rmvs_tp)
      associate(ntp => tp%nbody, npl => pl%nbody, t => config%t, dt => config%dt)
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
            pl%outer(0)%x(:,:) = xbeg(:,:)
            pl%outer(0)%v(:,:) = vbeg(:,:)
            call pl%step(cb, config, t, dt) 
            pl%outer(NTENC)%x(:,:) = pl%xh(:,:)
            pl%outer(NTENC)%v(:,:) = pl%vh(:,:)
            call tp%set_beg_end(xend = pl%xh)
            call rmvs_interp_out(pl,cb, dt, config)
            call rmvs_step_out(pl, cb, tp, dt, config)
            call tp%reverse_status()
            call tp%set_beg_end(xbeg = xbeg, xend = xend)
            tp%lfirst = .true.
            call tp%step(cb, pl, config, t, dt)
            where (tp%status(:) == INACTIVE) tp%status(:) = ACTIVE
            pl%lfirst = lfirstpl
            tp%lfirst = lfirsttp
         else
            call whm_step_system(self, config)
         end if
      end associate
      end select
      end select
      end select
      return

   end subroutine rmvs_step_system 

   subroutine rmvs_interp_in(pl, cb, dt, enc_index, config)
      !! author: David A. Minton
      !!
      !! Interpolate planet positions between two Keplerian orbits in inner encounter regio
      !!
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_interp_in.f90 
      !!
      !! Adapted from Hal Levison's Swift routine rmvs3_interp.f
      implicit none
      ! Arguments
      class(rmvs_pl),                intent(inout) :: pl     !! RMVS test particle object
      class(rmvs_cb),                intent(inout) :: cb     !! RMVS central body particle type
      real(DP),                      intent(in)    :: dt     !! Step size
      integer(I4B),                  intent(in)    :: enc_index !! Outer substep number within current se
      class(swiftest_configuration), intent(in)    :: config !! Swiftest configuration file
      ! Internals
      integer(I4B)                    :: i, j
      real(DP)                        :: frac, dntphenc
      real(DP), dimension(:,:), allocatable      :: xtmp, vtmp, xh_original
      real(DP), dimension(:), allocatable         :: msun, dti
      integer(I4B), dimension(:), allocatable      :: iflag

      associate (npl => pl%nbody)
         dntphenc = real(NTPHENC, kind=DP)

         ! Set the endpoints of the inner region from the outer region values in the current outer step index
         pl%inner(0)%x(:,:) = pl%outer(enc_index - 1)%x(:, :)
         pl%inner(0)%v(:,:) = pl%outer(enc_index - 1)%v(:, :)
         pl%inner(NTPHENC)%x(:,:) = pl%outer(enc_index)%x(:, :)
         pl%inner(NTPHENC)%v(:,:) = pl%outer(enc_index)%v(:, :)

         allocate(xtmp,mold=pl%xh)
         allocate(vtmp,mold=pl%vh)
         allocate(msun(npl))
         allocate(dti(npl))
         allocate(iflag(npl))
         dti(:) = dt / dntphenc
         msun(:) = cb%Gmass
         xtmp(:, :) = pl%inner(0)%x(:, :)
         vtmp(:, :) = pl%inner(0)%v(:, :)
         if (config%loblatecb) then
            allocate(xh_original,mold=pl%xh)
            pl%xh(:, :) = xtmp(:, :) ! Temporarily replace heliocentric position with inner substep values to calculate the oblateness terms
            call pl%obl_acc(cb)
            pl%inner(0)%aobl(:, :) = pl%aobl(:, :) ! Save the oblateness acceleration on the planet for this substep
         end if

         do j = 1, NTPHENC - 1
            !!$omp simd
            do i = 1, npl
               call drift_one(msun(i), xtmp(1,i), xtmp(2,i), xtmp(3,i), vtmp(1,i), vtmp(2,i), vtmp(3,i), dti(i), iflag(i))
            end do
            if (any(iflag(:) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /=0) then
                     write(*, *) " Planet ", pl%name(i), " is lost!!!!!!!!!!"
                     write(*, *) msun(i), dti(i)
                     write(*, *) xtmp(:,i)
                     write(*, *) vtmp(:,i)
                     write(*, *) " STOPPING "
                     call util_exit(failure)
                  end if
               end do
            end if
            frac = 1.0_DP - j / dntphenc
            pl%inner(j)%x(:, :) = frac * xtmp(:,:)
            pl%inner(j)%v(:, :) = frac * vtmp(:,:)
         end do

         xtmp(:,:) = pl%inner(NTPHENC)%x(:, :)
         vtmp(:,:) = pl%inner(NTPHENC)%v(:, :)
 
         do j = NTPHENC - 1, 1, -1
            !!$omp simd
            do i = 1, npl
               call drift_one(msun(i), xtmp(1,i), xtmp(2,i), xtmp(3,i), vtmp(1,i), vtmp(2,i), vtmp(3,i), -dti(i), iflag(i))
            end do
            if (any(iflag(:) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /=0) then
                     write(*, *) " Planet ", pl%name(i), " is lost!!!!!!!!!!"
                     write(*, *) msun(i), -dti(i)
                     write(*, *) xtmp(:,i)
                     write(*, *) vtmp(:,i)
                     write(*, *) " STOPPING "
                     call util_exit(failure)
                  end if
               end do
            end if
            frac = j / dntphenc
            pl%inner(j)%x(:, :) = pl%inner(j)%x(:, :) + frac * xtmp(:, :)
            pl%inner(j)%v(:, :) = pl%inner(j)%v(:, :) + frac * vtmp(:, :)

            if (config%loblatecb) then
               pl%xh(:,:) = pl%inner(j)%x(:, :)
               call pl%obl_acc(cb)
               pl%inner(j)%aobl(:, :) = pl%aobl(:, :) 
            end if
         end do
         ! Put the planet positions back into place
         if (config%loblatecb) call move_alloc(xh_original, pl%xh)
      end associate
      return

   end subroutine rmvs_interp_in

   subroutine rmvs_interp_out(pl, cb, dt, config)
      !! author: David A. Minton
      !!
      !! Interpolate planet positions between two Keplerian orbits in outer encounter region
      !!
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_interp_out.f90 
      !!
      !! Adapted from Hal Levison's Swift routine rmvs3_interp.f
      implicit none
      ! Arguments
      class(rmvs_pl), intent(inout)   :: pl !! RMVS test particle object
      class(rmvs_cb), intent(inout)   :: cb   !! RMVS central body particle type
      real(DP), intent(in)            :: dt   !! Step size
      class(swiftest_configuration), intent(in)    :: config !! Swiftest configuration file
      ! Internals
      integer(I4B)                    :: i, j
      real(DP)                        :: frac, dntenc
      real(DP), dimension(:,:), allocatable       :: xtmp, vtmp
      real(DP), dimension(:), allocatable         :: msun, dto
      integer(I4B), dimension(:), allocatable      :: iflag

   ! executable code
      dntenc = real(NTENC, DP)
      associate (npl => pl%nbody)
         allocate(xtmp, mold = pl%xh)
         allocate(vtmp, mold = pl%vh)
         allocate(msun(npl))
         allocate(dto(npl))
         allocate(iflag(npl))
         dto(:) = dt / dntenc
         msun(:) = cb%Gmass
         xtmp(:,:) = pl%outer(0)%x(:, :)
         vtmp(:,:) = pl%outer(0)%v(:, :)
         do j = 1, NTENC - 1
            !!$omp simd
            do i = 1, npl
               call drift_one(msun(i), xtmp(1,i), xtmp(2,i), xtmp(3,i), vtmp(1,i), vtmp(2,i), vtmp(3,i), dto(i), iflag(i))
            end do
            if (any(iflag(:) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /= 0) then
                     write(*, *) " Planet ", pl%name(i), " is lost!!!!!!!!!!"
                     write(*, *) msun(i), dto(i)
                     write(*, *) xtmp(:,i)
                     write(*, *) vtmp(:,i)
                     write(*, *) " STOPPING "
                     call util_exit(FAILURE)
                  end if
               end do
            end if
            frac = 1.0_DP - j / dntenc
            pl%outer(j)%x(:, :) = frac * xtmp(:,:)
            pl%outer(j)%v(:, :) = frac * vtmp(:,:)
         end do
         xtmp(:,:) = pl%outer(NTENC)%x(:, :)
         vtmp(:,:) = pl%outer(NTENC)%v(:, :)
         do j = NTENC - 1, 1, -1
            !!$omp simd
            do i = 1, npl
               call drift_one(msun(i), xtmp(1,i), xtmp(2,i), xtmp(3,i), vtmp(1,i), vtmp(2,i), vtmp(3,i), -dto(i), iflag(i))
            end do
            if (any(iflag(:) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /= 0) then
                     write(*, *) " Planet ", pl%name(i), " is lost!!!!!!!!!!"
                     write(*, *) msun(i), -dto(i)
                     write(*, *) xtmp(:,i)
                     write(*, *) vtmp(:,i)
                     write(*, *) " STOPPING "
                     call util_exit(FAILURE)
                  end if
               end do
            end if
            frac = j / dntenc
            pl%outer(j)%x(:, :) = pl%outer(j)%x(:, :) + frac * xtmp(:,:)
            pl%outer(j)%v(:, :) = pl%outer(j)%v(:, :) + frac * vtmp(:,:)
         end do
      end associate

      return

   end subroutine rmvs_interp_out   

   subroutine rmvs_peri_tp(tp, pl, t, dt, lfirst, enc_index, ipleP, config)
      !! author: David A. Minton
      !!
      !! Determine planetocentric pericenter passages for test particles in close encounters with a planet
      !! 
      !! Adapted from Hal Levison's Swift routine Adapted from Hal Levison's Swift routine util_peri.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_peri.f90
      implicit none
      ! Arguments
      class(rmvs_tp),                intent(inout) :: tp        !! RMVS test particle object (planetocentric) 
      class(rmvs_pl),                intent(inout) :: pl        !! RMVS massive body object (heliocentric)
      real(DP),                      intent(in)    :: t         !! current time
      real(DP),                      intent(in)    :: dt        !! step size
      logical,                       intent(in)    :: lfirst    !! Logical flag indicating whether current invocation is the first
      integer(I4B),                  intent(in)    :: enc_index !! Outer substep number within current set
      integer(I4B),                  intent(in)    :: ipleP     !!  index of RMVS planet being closely encountered
      class(swiftest_configuration), intent(in)    :: config    !! Input collection of  configuration parameters
      ! Internals
      integer(I4B)                                 :: i, id1, id2
      real(DP)                                     :: r2, mu, rhill2, vdotr, a, peri, capm, tperi, rpl
      real(DP), dimension(NDIM)                    :: xh1, xh2, vh1, vh2

      rhill2 = pl%rhill(ipleP)**2
      mu = pl%Gmass(ipleP)
      associate( nenc => tp%nbody, xpc => tp%xh, vpc => tp%vh)
         if (lfirst) then
            do i = 1, nenc
               if (tp%status(i) == ACTIVE) then
                  vdotr = dot_product(xpc(:, i), vpc(:, i))
                  if (vdotr > 0.0_DP) then
                     tp%isperi(i) = 1
                  else
                     tp%isperi(i) = -1
                  end if
               end if
            end do
         else
            do i = 1, nenc
               if (tp%status(i) == ACTIVE) then
                  vdotr = dot_product(xpc(:, i), vpc(:, i))
                  if (tp%isperi(i) == -1) then
                     if (vdotr >= 0.0_DP) then
                        tp%isperi(i) = 0
                        call orbel_xv2aqt(mu, xpc(:, i), vpc(:, i), a, peri, capm, tperi)
                        r2 = dot_product(xpc(:, i), xpc(:, i))
                        if ((abs(tperi) > FACQDT * dt) .or. (r2 > rhill2)) peri = sqrt(r2)
                        if (config%encounter_file /= "") then
                           id1 = pl%name(ipleP)
                           rpl = pl%radius(ipleP)
                           xh1(:) = pl%inner(enc_index)%x(:, ipleP)
                           vh1(:) = pl%inner(enc_index)%v(:, ipleP)
                           id2 = tp%name(i)
                           xh2(:) = xpc(:, i) + xh1(:)
                           vh2(:) = xpc(:, i) + vh1(:)
                           call io_write_encounter(t, id1, id2, mu, 0.0_DP, rpl, 0.0_DP, xh1(:), xh2(:), vh1(:), vh2(:),  &
                              config%encounter_file, config%out_type)
                        end if
                        if (tp%lperi(i)) then
                           if (peri < tp%peri(i)) then
                              tp%peri(i) = peri
                              tp%plperP(i) = ipleP
                           end if
                        else
                           tp%lperi(i) = .true.
                           tp%peri(i) = peri
                           tp%plperP(i) = ipleP
                        end if
                     end if
                  else
                     if (vdotr > 0.0_DP) then
                        tp%isperi(i) = 1
                     else
                        tp%isperi(i) = -1
                     end if
                  end if
               end if
            end do                   
         end if
         end associate
      return

   end subroutine rmvs_peri_tp

   subroutine rmvs_step_out(pl, cb, tp, dt, config)
      !! author: David A. Minton
      !!
      !! Step ACTIVE test particles ahead in the outer encounter region, setting up and calling the inner region
      !!    integration if necessar
      !! 
      !! Adapted from Hal Levison's Swift routines rmvs3_step_out.f and rmvs3_step_out2.f
      !! Adapted from David E. Kaufmann's Swifter routines rmvs_step_out.f90 and rmvs_step_out2.f90
      implicit none
      ! Arguments
      class(rmvs_pl),                    intent(inout)  :: pl !! RMVS massive body object
      class(rmvs_cb),                    intent(inout)  :: cb   !! RMVS central body object
      class(rmvs_tp),                    intent(inout)  :: tp   !! RMVS test particle object
      real(DP),                          intent(in)     :: dt   !! Step size
      class(swiftest_configuration),     intent(in)     :: config  !! Input collection of  configuration parameters
      ! Internals
      integer(I4B)                                      :: enc_index, j, k
      real(DP)                                          :: dto, outer_time, rts
      logical                                           :: lencounter, lfirsttp

      associate(npl => pl%nbody, ntp => tp%nbody, t => config%t)
         dto = dt / NTENC
         where(tp%plencP(:) == 0)
            tp%status(:) = INACTIVE
         elsewhere  
            tp%lperi(:) = .false.
         end where
         do enc_index = 1, NTENC
            outer_time = t + (enc_index - 1) * dto
            call tp%set_beg_end(xbeg = pl%outer(enc_index - 1)%x(:, :), &
                                vbeg = pl%outer(enc_index - 1)%v(:, :), &
                                xend = pl%outer(enc_index    )%x(:, :))
            rts = RHPSCALE
            lencounter = tp%encounter_check(cb, pl, dt, rts) 
            if (lencounter) then
               ! Interpolate planets in inner encounter region
               call rmvs_interp_in(pl, cb, dto, enc_index, config)
               ! Step through the inner region
               call rmvs_step_in(pl, cb, tp, config, outer_time, dto)
               lfirsttp = tp%lfirst
               tp%lfirst = .true.
               call tp%step(cb, pl, config, outer_time, dto)
               tp%lfirst = lfirsttp
            else
               call tp%step(cb, pl, config, outer_time, dto)
            end if
            do j = 1, npl
               if (pl%nenc(j) == 0) cycle
               where((tp%plencP(:) == j) .and. (tp%status(:) == INACTIVE)) tp%status(:) = ACTIVE
            end do
         end do
      end associate
      return

   end subroutine rmvs_step_out   

   subroutine rmvs_step_in(pl, cb, tp, config, outer_time, dto)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead in the inner encounter region
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step_in.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_in.f90
      implicit none
      ! Arguments
      class(rmvs_pl),                intent(inout)  :: pl   !! RMVS massive body object
      class(rmvs_cb),                intent(inout)  :: cb   !! RMVS central body object
      class(rmvs_tp),                intent(inout)  :: tp   !! RMVS test particle object
      class(swiftest_configuration), intent(in)     :: config  !! Input collection of configuration parameters 
      real(DP),                      intent(in)     :: outer_time    !! Current time
      real(DP),                      intent(in)     :: dto   !! Step size
      ! Internals
      logical                                        :: lfirsttp
      integer(I4B)                                   :: i, j, ipleP
      real(DP)                                       :: dti, inner_time
      real(DP), dimension(:, :), allocatable         :: xbeg, xend, vbeg


      associate(npl => pl%nbody, nenc => pl%nenc)
         dti = dto / NTPHENC
         allocate(xbeg, mold=pl%xh)
         allocate(xend, mold=pl%xh)
         allocate(vbeg, mold=pl%vh)
         if (config%loblatecb) call pl%obl_acc(cb)
         call rmvs_make_planetocentric(pl, cb, tp, config)
         do i = 1, npl
            if (nenc(i) == 0) cycle
            associate(cbenci => pl%planetocentric(i)%cb, plenci => pl%planetocentric(i)%pl, &
                        tpenci => pl%planetocentric(i)%tp)
               associate(enc_index => tpenci%index, plind => plenci%plind)
               ! There are inner encounters with this planet...switch to planetocentric coordinates to proceed
                  tpenci%lfirst = .true.
                  inner_time = outer_time
                  call rmvs_peri_tp(tpenci, pl, inner_time, dti, .true., 0, i, config) 
                  ! now step the encountering test particles fully through the inner encounter
                  lfirsttp = .true.
                  do enc_index = 1, NTPHENC ! Integrate over the encounter region, using the "substitute" planetocentric systems at each level
                     xbeg(:,1) = cb%xin(:) - pl%inner(enc_index - 1)%x(:, i)
                     xend(:,1) = cb%xin(:) - pl%inner(enc_index    )%x(:, i)
                     vbeg(:,1) = cb%vin(:) - pl%inner(enc_index - 1)%v(:, i)
                     do j = 2, npl
                        ipleP = plind(j)
                        xbeg(:,j) = pl%inner(enc_index - 1)%x(:,ipleP) - pl%inner(enc_index - 1)%x(:,i)  
                        xend(:,j) = pl%inner(enc_index    )%x(:,ipleP) - pl%inner(enc_index    )%x(:,i)  
                        vbeg(:,j) = pl%inner(enc_index - 1)%v(:,ipleP) - pl%inner(enc_index - 1)%v(:,i)
                     end do
                     plenci%xh(:,:) = xbeg(:,:)
                     plenci%vh(:,:) = vbeg(:,:)
                     call tpenci%set_beg_end(xbeg = xbeg, xend = xend)
                     call tpenci%step(cbenci, pl, config, inner_time, dti)
                     do j = 1, nenc(i)
                        tpenci%xheliocentric(:, j) = tpenci%xh(:, j) + pl%inner(enc_index    )%x(:,i)
                     end do
                     inner_time = outer_time + j * dti
                     call rmvs_peri_tp(tpenci, pl, inner_time, dti, .false., enc_index, i, config) 
                  end do
                  where(tpenci%status(:) == ACTIVE) tpenci%status(:) = INACTIVE
               end associate
            end associate
         end do
         call rmvs_end_planetocentric(pl, cb, tp)
      end associate
      return
   end subroutine rmvs_step_in

   subroutine rmvs_make_planetocentric(pl, cb, tp, config)
      !! author: David A. Minton
      !!
      !! When encounters are detected, this method will call the interpolation methods for the planets and 
      !! creates a Swiftest test particle structure for each planet's encountering test particles to simplify the 
      !! planetocentric calculations. This subroutine is not based on an existing one from Swift and Swifter
      !!
      implicit none
      ! Arguments
      class(rmvs_pl),                 intent(inout)  :: pl !! RMVS test particle object
      class(rmvs_cb),                 intent(inout)  :: cb   !! RMVS central body particle type
      class(rmvs_tp),                 intent(inout)  :: tp   !! RMVS test particle object
      class(swiftest_configuration),  intent(in)     :: config !! Input collection of configuration parameters 
      ! Internals
      integer(I4B)                                   :: i, j
      logical, dimension(:), allocatable             :: encmask

      associate(npl => pl%nbody, ntp => tp%nbody, GMpl => pl%Gmass, nenc => pl%nenc)
         do i = 1, npl
            if (nenc(i) == 0) cycle 
            ! There are inner encounters with this planet
            if (allocated(encmask)) deallocate(encmask)
            allocate(encmask(ntp))
            encmask(:) = tp%plencP(:) == i

            ! Create encountering test particle structure
            allocate(rmvs_tp :: pl%planetocentric(i)%tp)
            associate(tpenci => pl%planetocentric(i)%tp, cbenci => pl%planetocentric(i)%cb, plenci => pl%planetocentric(i)%pl)
               tpenci%lplanetocentric = .true.
               call tpenci%setup(nenc(i))
               tpenci%ipleP = i
               tpenci%status(:) = ACTIVE
               tpenci%name(:) = pack(tp%name(:), encmask(:)) 
               ! Grab all the encountering test particles and convert them to a planetocentric frame
               do j = 1, NDIM 
                  tpenci%xheliocentric(j, :) = pack(tp%xh(j,:), encmask(:)) 
                  tpenci%xh(j, :) = tpenci%xheliocentric(j, :) - pl%inner(0)%x(j, i)
                  tpenci%vh(j, :) = pack(tp%vh(j,:), encmask(:)) - pl%inner(0)%v(j, i)
               end do
               ! Make sure that the test particles get the planetocentric value of mu 
               call tpenci%set_mu(cbenci)
            end associate
         end do
      end associate
      return
   end subroutine rmvs_make_planetocentric

   subroutine rmvs_end_planetocentric(pl, cb, tp)
      !! author: David A. Minton
      !!
      !! Deallocates all of the encountering particle data structures for next time
      !!
      implicit none
      ! Arguments
      class(rmvs_pl),                 intent(inout)  :: pl !! RMVS test particle object
      class(rmvs_cb),                 intent(inout)  :: cb      !!  RMVS central body object
      class(rmvs_tp),                 intent(inout)  :: tp   !! RMVS test particle object
      ! Internals
      integer(I4B) :: i, j
      integer(I4B), dimension(:), allocatable :: tpind
      logical, dimension(:), allocatable :: encmask

      associate(nenc => pl%nenc, npl => pl%nbody, ntp => tp%nbody)
         do i = 1, npl
            if (nenc(i) == 0) cycle
            associate(tpenci => pl%planetocentric(i)%tp)
               allocate(tpind(nenc(i)))
               ! Index array of encountering test particles
               if (allocated(encmask)) deallocate(encmask)
               allocate(encmask(ntp))
               encmask(:) = tp%plencP(:) == i
               tpind(:) = pack([(j,j=1,ntp)], encmask(:))
      
               ! Copy the results of the integration back over and shift back to heliocentric reference
               tp%status(tpind(1:nenc(i))) = tpenci%status(1:nenc(i)) 
               do j = 1, NDIM
                  tp%xh(j, tpind(1:nenc(i))) = tpenci%xh(j,1:nenc(i)) + pl%inner(NTPHENC)%x(j, i)
                  tp%vh(j, tpind(1:nenc(i))) = tpenci%vh(j,1:nenc(i)) + pl%inner(NTPHENC)%v(j, i)
               end do
               deallocate(pl%planetocentric(i)%tp)
            end associate
         end do
      end associate
      return
   end subroutine rmvs_end_planetocentric
   
end submodule s_rmvs_step
