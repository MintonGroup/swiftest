!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(rmvs) s_rmvs_step
   use swiftest
contains

   module subroutine rmvs_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step.f90
      implicit none
      ! Arguments
      class(rmvs_nbody_system),   intent(inout)  :: self   !! RMVS nbody system object
      class(swiftest_parameters), intent(inout)  :: param  !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! Current simulation time
      real(DP),                   intent(in)    :: dt    !! Current stepsiz
      ! Internals
      logical :: lencounter, lfirstpl
      real(DP), dimension(:,:), allocatable :: rbeg, rend, vbeg

      if ((self%tp%nbody == 0) .or. (self%pl%nbody == 0)) then
         call whm_step_system(self, param, t, dt)
      else
         select type(cb => self%cb)
         class is (rmvs_cb)
            select type(pl => self%pl)
            class is (rmvs_pl)
               select type(tp => self%tp)
               class is (rmvs_tp)
                  associate(nbody_system => self, ntp => tp%nbody, npl => pl%nbody)
                     allocate(rbeg, source=pl%rh)
                     allocate(vbeg, source=pl%vh)
                     call pl%set_beg_end(rbeg = rbeg, vbeg = vbeg)
                     ! ****** Check for close encounters ***** !
                     call pl%set_renc(RHSCALE)
                     lencounter = tp%encounter_check(param, nbody_system, dt)
                     if (lencounter) then
                        lfirstpl = pl%lfirst
                        pl%outer(0)%x(:, 1:npl) = rbeg(:, 1:npl)
                        pl%outer(0)%v(:, 1:npl) = vbeg(:, 1:npl)
                        call pl%step(nbody_system, param, t, dt) 
                        pl%outer(NTENC)%x(:, 1:npl) = pl%rh(:, 1:npl)
                        pl%outer(NTENC)%v(:, 1:npl) = pl%vh(:, 1:npl)
                        call rmvs_interp_out(cb, pl, dt)
                        call rmvs_step_out(cb, pl, tp, nbody_system, param, t, dt) 
                        tp%lmask(1:ntp) = .not. tp%lmask(1:ntp)
                        call pl%set_beg_end(rbeg = rbeg, rend = rend)
                        tp%lfirst = .true.
                        call tp%step(nbody_system, param, t, dt)
                        tp%lmask(1:ntp) = .true.
                        pl%lfirst = lfirstpl
                        tp%lfirst = .true.
                        ! if (param%ltides) call nbody_system%step_spin(param, t, dt)
                     else
                        call whm_step_system(nbody_system, param, t, dt)
                     end if
                  end associate
               end select
            end select
         end select
      end if

      return
   end subroutine rmvs_step_system 


   subroutine rmvs_interp_out(cb, pl, dt)
      !! author: David A. Minton
      !!
      !! Interpolate planet positions between two Keplerian orbits in outer encounter region
      !!
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_interp_out.f90 
      !!
      !! Adapted from Hal Levison's Swift routine rmvs3_interp.f
      implicit none
      ! Arguments
      class(rmvs_cb),             intent(inout) :: cb      !! RMVS central body object
      class(rmvs_pl),             intent(inout) :: pl      !! RMVS massive body object
      real(DP),                   intent(in)    :: dt   !! Step size
      ! Internals
      integer(I4B)                              :: i, outer_index
      real(DP)                                  :: frac, dntenc
      real(DP),     dimension(:,:), allocatable :: xtmp, vtmp
      real(DP),     dimension(:),   allocatable :: GMcb, dto
      integer(I4B), dimension(:),   allocatable :: iflag
      character(len=STRMAX) :: message

      dntenc = real(NTENC, kind=DP)
      associate (npl => pl%nbody)
         allocate(xtmp, mold = pl%rh)
         allocate(vtmp, mold = pl%vh)
         allocate(GMcb(npl))
         allocate(dto(npl))
         allocate(iflag(npl))
         dto(1:npl) = dt / dntenc
         GMcb(1:npl) = cb%Gmass
         xtmp(:,1:npl) = pl%outer(0)%x(:, 1:npl)
         vtmp(:,1:npl) = pl%outer(0)%v(:, 1:npl)
         do outer_index = 1, NTENC - 1
            call swiftest_drift_one(GMcb(1:npl), xtmp(1,1:npl), xtmp(2,1:npl), xtmp(3,1:npl), &
                                        vtmp(1,1:npl), vtmp(2,1:npl), vtmp(3,1:npl), &
                                        dto(1:npl), iflag(1:npl))
            if (any(iflag(1:npl) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /= 0) then
                     write(message, *) " Planet ", pl%id(i), " is lost!!!!!!!!!!",new_line('a'), &
                                       GMcb(i), dto(i),new_line('a'), &
                                       xtmp(:,i),new_line('a'), &
                                       vtmp(:,i),new_line('a'), &
                                       " STOPPING "
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT,message)
                     call base_util_exit(FAILURE)
                  end if
               end do
            end if
            frac = 1.0_DP - outer_index / dntenc
            pl%outer(outer_index)%x(:, 1:npl) = frac * xtmp(:, 1:npl)
            pl%outer(outer_index)%v(:, 1:npl) = frac * vtmp(:, 1:npl)
         end do
         xtmp(:, 1:npl) = pl%outer(NTENC)%x(:, 1:npl)
         vtmp(:, 1:npl) = pl%outer(NTENC)%v(:, 1:npl)
         do outer_index = NTENC - 1, 1, -1
            call swiftest_drift_one(GMcb(1:npl), xtmp(1,1:npl), xtmp(2,1:npl), xtmp(3,1:npl), &
                                        vtmp(1,1:npl), vtmp(2,1:npl), vtmp(3,1:npl), &
                                       -dto(1:npl), iflag(1:npl))
            if (any(iflag(1:npl) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /= 0) then
                     write(message, *) " Planet ", pl%id(i), " is lost!!!!!!!!!!",new_line('a'), &
                                       GMcb(i), -dto(i), new_line('a'), &
                                       xtmp(:,i), new_line('a'), &
                                       vtmp(:,i), new_line('a'), &
                                       " STOPPING "
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT,message)
                     call base_util_exit(FAILURE)
                  end if
               end do
            end if
            frac = outer_index / dntenc
            pl%outer(outer_index)%x(:, 1:npl) = pl%outer(outer_index)%x(:, 1:npl) + frac * xtmp(:, 1:npl)
            pl%outer(outer_index)%v(:, 1:npl) = pl%outer(outer_index)%v(:, 1:npl) + frac * vtmp(:, 1:npl)
         end do
      end associate

      return
   end subroutine rmvs_interp_out   


   subroutine rmvs_step_out(cb, pl, tp, nbody_system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step ACTIVE test particles ahead in the outer encounter region, setting up and calling the inner region
      !!    integration if necessar
      !! 
      !! Adapted from Hal Levison's Swift routines rmvs3_step_out.f and rmvs3_step_out2.f
      !! Adapted from David E. Kaufmann's Swifter routines rmvs_step_out.f90 and rmvs_step_out2.f90
      implicit none
      ! Arguments
      class(rmvs_cb),             intent(inout) :: cb     !! RMVS central body object
      class(rmvs_pl),             intent(inout) :: pl     !! RMVS massive body object
      class(rmvs_tp),             intent(inout) :: tp     !! RMVS test particle object
      class(rmvs_nbody_system),   intent(inout) :: nbody_system !! RMVS nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t      !! Current simulation time
      real(DP),                   intent(in)    :: dt     !! Current stepsiz
      ! Internals
      integer(I4B)                              :: outer_index, j
      real(DP)                                  :: dto, outer_time
      logical                                   :: lencounter, lfirsttp

      associate(npl => pl%nbody, ntp => tp%nbody)
         dto = dt / NTENC
         where(tp%plencP(1:ntp) == 0)
            tp%lmask(1:ntp) = .false.
         elsewhere  
            tp%lperi(1:ntp) = .false.
         end where
         call pl%set_renc(RHPSCALE)
         do outer_index = 1, NTENC
            outer_time = t + (outer_index - 1) * dto
            call pl%set_beg_end(rbeg = pl%outer(outer_index - 1)%x(:, 1:npl), &
                                vbeg = pl%outer(outer_index - 1)%v(:, 1:npl), &
                                rend = pl%outer(outer_index    )%x(:, 1:npl))
            lencounter = tp%encounter_check(param, nbody_system, dto) 
            if (lencounter) then
               ! Interpolate planets in inner encounter region
               call rmvs_interp_in(cb, pl, nbody_system, param, dto, outer_index) 
               ! Step through the inner region
               call rmvs_step_in(cb, pl, tp, param, outer_time, dto)
               lfirsttp = tp%lfirst
               tp%lfirst = .true.
               call tp%step(nbody_system, param, outer_time, dto)
               tp%lfirst = lfirsttp
            else
               if (param%loblatecb) then
                  call swiftest_obl_acc(npl, cb%Gmass, cb%j2rp2, cb%j4rp4, pl%rbeg, pl%lmask, pl%outer(outer_index-1)%aobl, pl%Gmass, cb%aoblbeg)
                  call swiftest_obl_acc(npl, cb%Gmass, cb%j2rp2, cb%j4rp4, pl%rend, pl%lmask, pl%outer(outer_index)%aobl, pl%Gmass, cb%aoblend)
               end if
               call tp%step(nbody_system, param, outer_time, dto)
            end if
            do j = 1, npl
               if (pl%nenc(j) == 0) cycle
               tp%lfirst = .true.
               where((tp%plencP(1:ntp) == j) .and. (.not.tp%lmask(1:ntp)))
                  tp%lmask(1:ntp) = .true.
               end where 
            end do
         end do
      end associate

      return
   end subroutine rmvs_step_out


   subroutine rmvs_interp_in(cb, pl, nbody_system, param, dt, outer_index)
      !! author: David A. Minton
      !!
      !! Interpolate planet positions between two Keplerian orbits in inner encounter regio
      !!
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_interp_in.f90 
      !!
      !! Adapted from Hal Levison's Swift routine rmvs3_interp.f
      implicit none
      ! Arguments
      class(rmvs_cb),             intent(inout) :: cb          !! RMVS cenral body object
      class(rmvs_pl),             intent(inout) :: pl          !! RMVS massive body object
      class(rmvs_nbody_system),   intent(inout) :: nbody_system      !! RMVS nbody system object
      class(swiftest_parameters), intent(in)    :: param       !! Swiftest parameters file
      real(DP),                   intent(in)    :: dt          !! Step size
      integer(I4B),               intent(in)    :: outer_index !! Outer substep number within current set
      ! Internals
      integer(I4B)                              :: i, inner_index
      real(DP)                                  :: frac, dntphenc
      real(DP),     dimension(:,:), allocatable :: xtmp, vtmp, rh_original, ah_original
      real(DP),     dimension(:),   allocatable :: GMcb, dti
      integer(I4B), dimension(:),   allocatable :: iflag
      character(len=STRMAX) :: message

      associate (npl => nbody_system%pl%nbody)
         dntphenc = real(NTPHENC, kind=DP)

         ! Set the endpoints of the inner region from the outer region values in the current outer step index
         pl%inner(0)%x(:, 1:npl) = pl%outer(outer_index - 1)%x(:, 1:npl)
         pl%inner(0)%v(:, 1:npl) = pl%outer(outer_index - 1)%v(:, 1:npl)
         pl%inner(NTPHENC)%x(:, 1:npl) = pl%outer(outer_index)%x(:, 1:npl)
         pl%inner(NTPHENC)%v(:, 1:npl) = pl%outer(outer_index)%v(:, 1:npl)

         allocate(xtmp,mold=pl%rh)
         allocate(vtmp,mold=pl%vh)
         allocate(GMcb(npl))
         allocate(dti(npl))
         allocate(iflag(npl))
         dti(1:npl) = dt / dntphenc
         GMcb(1:npl) = cb%Gmass
         xtmp(:, 1:npl) = pl%inner(0)%x(:, 1:npl)
         vtmp(:, 1:npl) = pl%inner(0)%v(:, 1:npl)

         if ((param%loblatecb) .or. (param%ltides)) then
            allocate(rh_original, source=pl%rh)
            allocate(ah_original, source=pl%ah)
            pl%rh(:, 1:npl) = xtmp(:, 1:npl) ! Temporarily replace heliocentric position with inner substep values to calculate the oblateness terms
         end if
         if (param%loblatecb) then
            call pl%accel_obl(nbody_system)
            pl%inner(0)%aobl(:, 1:npl) = pl%aobl(:, 1:npl) ! Save the oblateness acceleration on the planet for this substep
         end if
         ! TODO: Implement tides
         ! if (param%ltides) then
         !    call pl%accel_tides(nbody_system)
         !    pl%inner(0)%atide(:, 1:npl) = pl%atide(:, 1:npl) ! Save the oblateness acceleration on the planet for this substep
         ! end if

         do inner_index = 1, NTPHENC - 1
            call swiftest_drift_one(GMcb(1:npl), xtmp(1,1:npl), xtmp(2,1:npl), xtmp(3,1:npl), &
                                        vtmp(1,1:npl), vtmp(2,1:npl), vtmp(3,1:npl), &
                                        dti(1:npl), iflag(1:npl))
            if (any(iflag(1:npl) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /=0) then
                     write(message, *) " Planet ", pl%id(i), " is lost!!!!!!!!!!", new_line('a'), &
                                       GMcb(i), dti(i), new_line('a'), &
                                       xtmp(:,i), new_line('a'), &
                                       vtmp(:,i), new_line('a'), &
                                       " STOPPING "
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                     call base_util_exit(FAILURE,param%display_unit)
                  end if
               end do
            end if
            frac = 1.0_DP - inner_index / dntphenc
            pl%inner(inner_index)%x(:, 1:npl) = frac * xtmp(:, 1:npl)
            pl%inner(inner_index)%v(:, 1:npl) = frac * vtmp(:, 1:npl)
         end do

         xtmp(:, 1:npl) = pl%inner(NTPHENC)%x(:, 1:npl)
         vtmp(:, 1:npl) = pl%inner(NTPHENC)%v(:, 1:npl)

         do inner_index = NTPHENC - 1, 1, -1
            call swiftest_drift_one(GMcb(1:npl), xtmp(1,1:npl), xtmp(2,1:npl), xtmp(3,1:npl), &
                                        vtmp(1,1:npl), vtmp(2,1:npl), vtmp(3,1:npl), &
                                        -dti(1:npl), iflag(1:npl))
            if (any(iflag(1:npl) /= 0)) then
               do i = 1, npl
                  if (iflag(i) /=0) then
                     write(*, *) " Planet ", pl%id(i), " is lost!!!!!!!!!!"
                     write(*, *) GMcb(i), -dti(i)
                     write(*, *) xtmp(:,i)
                     write(*, *) vtmp(:,i)
                     write(*, *) " STOPPING "
                     call base_util_exit(FAILURE,param%display_unit)
                  end if
               end do
            end if
            frac = inner_index / dntphenc
            pl%inner(inner_index)%x(:, 1:npl) = pl%inner(inner_index)%x(:, 1:npl) + frac * xtmp(:, 1:npl)
            pl%inner(inner_index)%v(:, 1:npl) = pl%inner(inner_index)%v(:, 1:npl) + frac * vtmp(:, 1:npl)

            if (param%loblatecb) then
               pl%rh(:,1:npl) = pl%inner(inner_index)%x(:, 1:npl)
               call pl%accel_obl(nbody_system)
               pl%inner(inner_index)%aobl(:, 1:npl) = pl%aobl(:, 1:npl) 
            end if
            ! TODO: Implement tides
            ! if (param%ltides) then 
            !    call pl%accel_tides(nbody_system)
            !    pl%inner(inner_index)%atide(:, 1:npl) = pl%atide(:, 1:npl)  
            ! end if
         end do
         if (param%loblatecb) then
            ! Calculate the final value of oblateness accelerations at the final inner substep
            pl%rh(:, 1:npl) = pl%inner(NTPHENC)%x(:, 1:npl)
            call pl%accel_obl(nbody_system)
            pl%inner(NTPHENC)%aobl(:, 1:npl) = pl%aobl(:, 1:npl) 
         end if
         ! TODO: Implement tides
         ! if (param%ltides) then
         !    call pl%accel_tides(nbody_system)
         !    pl%inner(NTPHENC)%atide(:, 1:npl) = pl%atide(:, 1:npl) 
         ! end if
         ! Put the planet positions and accelerations back into place 
         if (allocated(rh_original)) call move_alloc(rh_original, pl%rh)
         if (allocated(ah_original)) call move_alloc(ah_original, pl%ah)
      end associate
      return

   end subroutine rmvs_interp_in


   subroutine rmvs_step_in(cb, pl, tp, param, outer_time, dto)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead in the inner encounter region
      !! 
      !! Adapted from Hal Levison's Swift routine rmvs3_step_in.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_step_in.f90
      implicit none
      ! Arguments
      class(rmvs_cb),             intent(inout) :: cb         !! RMVS central body object
      class(rmvs_pl),             intent(inout) :: pl         !! RMVS massive body object
      class(rmvs_tp),             intent(inout) :: tp         !! RMVS test particle object
      class(swiftest_parameters), intent(inout) :: param      !! Current run configuration parameters 
      real(DP),                   intent(in)    :: outer_time !! Current time
      real(DP),                   intent(in)    :: dto        !! Outer step size
      ! Internals
      logical                                   :: lfirsttp
      integer(I4B)                              :: i, j
      real(DP)                                  :: dti, inner_time

      associate(npl => pl%nbody)
         dti = dto / NTPHENC
         call rmvs_make_planetocentric(param, cb, pl, tp)
         do i = 1, npl
            if (pl%nenc(i) == 0) cycle
            select type(planetocen_system => pl%planetocentric(i))
            class is (rmvs_nbody_system)
               select type(cbenci => planetocen_system%cb)
               class is (rmvs_cb)
                  select type(plenci => planetocen_system%pl)
                  class is (rmvs_pl)
                     select type(tpenci => planetocen_system%tp)
                     class is (rmvs_tp)
                        associate(inner_index => tpenci%index)
                           ! There are inner encounters with this planet...switch to planetocentric coordinates to proceed
                           tpenci%lfirst = .true.
                           inner_time = outer_time
                           call rmvs_peri_tp(tpenci, pl, inner_time, dti, .true., 0, i, param) 
                           ! now step the encountering test particles fully through the inner encounter
                           lfirsttp = .true.
                           do inner_index = 1, NTPHENC ! Integrate over the encounter region, using the "substitute" planetocentric systems at each level
                              plenci%rh(:, 1:npl) = plenci%inner(inner_index - 1)%x(:, 1:npl)
                              call plenci%set_beg_end(rbeg = plenci%inner(inner_index - 1)%x, &
                                                      rend = plenci%inner(inner_index)%x)

                              if (param%loblatecb) then
                                 cbenci%aoblbeg = cbenci%inner(inner_index - 1)%aobl(:, 1)
                                 cbenci%aoblend = cbenci%inner(inner_index    )%aobl(:, 1)
                              end if
                              if (param%ltides) then
                                 cbenci%atidebeg = cbenci%inner(inner_index - 1)%atide(:, 1)
                                 cbenci%atideend = cbenci%inner(inner_index    )%atide(:, 1)
                              end if

                              call tpenci%step(planetocen_system, param, inner_time, dti)
                              do j = 1, pl%nenc(i)
                                 tpenci%rheliocentric(:, j) = tpenci%rh(:, j) + pl%inner(inner_index)%x(:,i)
                              end do
                              inner_time = outer_time + j * dti
                              call rmvs_peri_tp(tpenci, pl, inner_time, dti, .false., inner_index, i, param) 
                           end do
                           tpenci%lmask(:) = .false.
                        end associate
                     end select
                  end select
               end select
            end select
         end do
         call rmvs_end_planetocentric(pl, tp)
      end associate
      return
   end subroutine rmvs_step_in


   subroutine rmvs_make_planetocentric(param, cb, pl, tp)
      !! author: David A. Minton
      !!
      !! When encounters are detected, this method will call the interpolation methods for the planets and 
      !! creates a Swiftest test particle structure for each planet's encountering test particles to simplify the 
      !! planetocentric calculations. This subroutine is not based on an existing one from Swift and Swifter
      !!
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration paramete
      class(rmvs_cb),             intent(inout) :: cb     !! RMVS central body object
      class(rmvs_pl),             intent(inout) :: pl     !! RMVS massive body object
      class(rmvs_tp),             intent(inout) :: tp     !! RMVS test particle object

      ! Internals
      integer(I4B)                        :: i, j, inner_index, ipc2hc
      logical, dimension(:), allocatable  :: encmask

      associate (npl => pl%nbody, ntp => tp%nbody)
         do i = 1, npl
            if (pl%nenc(i) == 0) cycle 
            ! There are inner encounters with this planet
            if (allocated(encmask)) deallocate(encmask)
            allocate(encmask(ntp))
            encmask(1:ntp) = tp%plencP(1:ntp) == i
            allocate(rmvs_tp :: pl%planetocentric(i)%tp)
            ! Create encountering test particle structure
            select type(cbenci => pl%planetocentric(i)%cb)
            class is (rmvs_cb)
               select type(plenci => pl%planetocentric(i)%pl)
               class is (rmvs_pl)
                  select type(tpenci => pl%planetocentric(i)%tp)
                  class is (rmvs_tp)
                     tpenci%lplanetocentric = .true.
                     associate(nenci => pl%nenc(i))
                        call tpenci%setup(nenci, param)
                        tpenci%cb_heliocentric = cb
                        tpenci%ipleP = i
                        tpenci%lmask(1:nenci) = .true.
                        tpenci%status(1:nenci) = ACTIVE
                        ! Grab all the encountering test particles and convert them to a planetocentric frame
                        tpenci%id(1:nenci) = pack(tp%id(1:ntp), encmask(1:ntp)) 
                        do j = 1, NDIM 
                           tpenci%rheliocentric(j, 1:nenci) = pack(tp%rh(j,1:ntp), encmask(:)) 
                           tpenci%rh(j, 1:nenci) = tpenci%rheliocentric(j, 1:nenci) - pl%inner(0)%x(j, i)
                           tpenci%vh(j, 1:nenci) = pack(tp%vh(j, 1:ntp), encmask(1:ntp)) - pl%inner(0)%v(j, i)
                        end do
                        tpenci%lperi(1:nenci) = pack(tp%lperi(1:ntp), encmask(1:ntp)) 
                        tpenci%plperP(1:nenci) = pack(tp%plperP(1:ntp), encmask(1:ntp)) 
                        ! Make sure that the test particles get the planetocentric value of mu 
                        allocate(cbenci%inner(0:NTPHENC))
                        do inner_index = 0, NTPHENC 
                           allocate(plenci%inner(inner_index)%x, mold=pl%inner(inner_index)%x)
                           allocate(plenci%inner(inner_index)%v, mold=pl%inner(inner_index)%x)
                           allocate(cbenci%inner(inner_index)%x(NDIM,1))
                           allocate(cbenci%inner(inner_index)%v(NDIM,1))
                           cbenci%inner(inner_index)%x(:,1)    =  pl%inner(inner_index)%x(:, i) 
                           cbenci%inner(inner_index)%v(:,1)    =  pl%inner(inner_index)%v(:, i) 
                           plenci%inner(inner_index)%x(:,1)    = -cbenci%inner(inner_index)%x(:,1)
                           plenci%inner(inner_index)%v(:,1)    = -cbenci%inner(inner_index)%v(:,1)
   
                           if (param%loblatecb) then
                              allocate(plenci%inner(inner_index)%aobl, mold=pl%inner(inner_index)%aobl)
                              allocate(cbenci%inner(inner_index)%aobl(NDIM,1))
                              cbenci%inner(inner_index)%aobl(:,1) =  pl%inner(inner_index)%aobl(:, i) 
                           end if
   
                           if (param%ltides) then  
                              allocate(plenci%inner(inner_index)%atide, mold=pl%inner(inner_index)%atide)
                              allocate(cbenci%inner(inner_index)%atide(NDIM,1))
                              cbenci%inner(inner_index)%atide(:,1) =  pl%inner(inner_index)%atide(:, i) 
                           end if
   
                           do j = 2, npl
                              ipc2hc = plenci%plind(j)
                              plenci%inner(inner_index)%x(:,j) = pl%inner(inner_index)%x(:, ipc2hc) &
                                                                 - cbenci%inner(inner_index)%x(:,1)
                              plenci%inner(inner_index)%v(:,j) = pl%inner(inner_index)%v(:, ipc2hc) &
                                                                 - cbenci%inner(inner_index)%v(:,1)
                           end do
                        end do
                        call tpenci%set_mu(cbenci)
                     end associate
                  end select
               end select
            end select
         end do
      end associate

      return
   end subroutine rmvs_make_planetocentric


   subroutine rmvs_peri_tp(tp, pl, t, dt, lfirst, inner_index, ipleP, param)
      !! author: David A. Minton
      !!
      !! Determine planetocentric pericenter passages for test particles in close encounters with a planet
      !! 
      !! Adapted from Hal Levison's Swift routine Adapted from Hal Levison's Swift routine util_peri.f
      !! Adapted from David E. Kaufmann's Swifter routine rmvs_peri.f90
      implicit none
      ! Arguments
      class(rmvs_tp),             intent(inout) :: tp        !! RMVS test particle object (planetocentric) 
      class(rmvs_pl),             intent(inout) :: pl        !! RMVS massive body object (heliocentric)
      real(DP),                   intent(in)    :: t         !! current time
      real(DP),                   intent(in)    :: dt        !! step size
      logical,                    intent(in)    :: lfirst    !! Logical flag indicating whether current invocation is the first
      integer(I4B),               intent(in)    :: inner_index !! Outer substep number within current set
      integer(I4B),               intent(in)    :: ipleP     !!  index of RMVS planet being closely encountered
      class(swiftest_parameters), intent(in)    :: param    !! Current run configuration parameters
      ! Internals
      integer(I4B)              :: i
      real(DP)                  :: r2, mu, rhill2, vdotr, a, peri, capm, tperi

      rhill2 = pl%rhill(ipleP)**2
      mu = pl%Gmass(ipleP)
      associate(nenc => tp%nbody, xpc => tp%rh, vpc => tp%vh)
         if (lfirst) then
            do i = 1, nenc
               if (tp%lmask(i)) then
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
               if (tp%lmask(i)) then
                  vdotr = dot_product(xpc(:, i), vpc(:, i))
                  if (tp%isperi(i) == -1) then
                     if (vdotr >= 0.0_DP) then
                        tp%isperi(i) = 0
                        call swiftest_orbel_xv2aqt(mu, xpc(1,i), xpc(2,i), xpc(3,i), vpc(1,i), vpc(2,i), vpc(3,i), &
                                          a, peri, capm, tperi)
                        r2 = dot_product(xpc(:, i), xpc(:, i))
                        if ((abs(tperi) > FACQDT * dt) .or. (r2 > rhill2)) peri = sqrt(r2)
                        ! TODO: write NetCDF encounter output writer
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


   subroutine rmvs_end_planetocentric(pl, tp)
      !! author: David A. Minton
      !!
      !! Deallocates all of the encountering particle data structures for next time
      !!
      implicit none
      ! Arguments
      class(rmvs_pl), intent(inout) :: pl     !! RMVS massive body object
      class(rmvs_tp), intent(inout) :: tp     !! RMVS test particle objec
      ! Internals
      integer(I4B) :: i, j, inner_index
      integer(I4B), dimension(:), allocatable :: tpind
      logical, dimension(:), allocatable :: encmask

      associate (npl => pl%nbody, ntp => tp%nbody)
         do i = 1, npl
            if (pl%nenc(i) == 0) cycle
            select type(cbenci => pl%planetocentric(i)%cb)
            class is (rmvs_cb)
               select type(plenci => pl%planetocentric(i)%pl)
               class is (rmvs_pl)
                  select type(tpenci => pl%planetocentric(i)%tp)
                  class is (rmvs_tp)
                     associate(nenci => pl%nenc(i))
                        if (allocated(tpind)) deallocate(tpind)
                        allocate(tpind(nenci))
                        ! Index array of encountering test particles
                        if (allocated(encmask)) deallocate(encmask)
                        allocate(encmask(ntp))
                        encmask(1:ntp) = tp%plencP(1:ntp) == i
                        tpind(1:nenci) = pack([(j,j=1,ntp)], encmask(1:ntp))
               
                        ! Copy the results of the integration back over and shift back to heliocentric reference
                        tp%status(tpind(1:nenci)) = tpenci%status(1:nenci) 
                        tp%lmask(tpind(1:nenci)) = tpenci%lmask(1:nenci) 
                        do j = 1, NDIM
                           tp%rh(j, tpind(1:nenci)) = tpenci%rh(j,1:nenci) + pl%inner(NTPHENC)%x(j, i)
                           tp%vh(j, tpind(1:nenci)) = tpenci%vh(j,1:nenci) + pl%inner(NTPHENC)%v(j, i)
                        end do
                        tp%lperi(tpind(1:nenci)) = tpenci%lperi(1:nenci)
                        tp%plperP(tpind(1:nenci)) = tpenci%plperP(1:nenci)
                        deallocate(pl%planetocentric(i)%tp)
                        deallocate(cbenci%inner)
                        do inner_index = 0, NTPHENC 
                           deallocate(plenci%inner(inner_index)%x) 
                           deallocate(plenci%inner(inner_index)%v) 
                           if (allocated(plenci%inner(inner_index)%aobl))  deallocate(plenci%inner(inner_index)%aobl)
                           if (allocated(plenci%inner(inner_index)%atide)) deallocate(plenci%inner(inner_index)%atide)
                        end do
                     end associate
                  end select
               end select
            end select
         end do
      end associate
      return
   end subroutine rmvs_end_planetocentric
   
end submodule s_rmvs_step
