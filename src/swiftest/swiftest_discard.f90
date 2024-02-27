! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_discard
contains

   module subroutine swiftest_discard_system(self, param)
      !! author: David A. Minton
      !!
      !! Calls the discard methods for each body class and then the write method if any discards were detected
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody_system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      logical :: lpl_discards, ltp_discards, lpl_check, ltp_check
      logical, dimension(:), allocatable :: ldiscard
      integer(I4B) :: i, nstart, nend, nsub
      character(len=STRMAX) :: idstr
      class(swiftest_pl), allocatable :: plsub
      class(swiftest_tp), allocatable :: tpsub

      lpl_check = allocated(self%pl_discards) .and. self%pl%nbody > 0
      ltp_check = allocated(self%tp_discards) .and. self%tp%nbody > 0

      associate(nbody_system => self,tp => self%tp,pl => self%pl,tp_discards => self%tp_discards,pl_discards => self%pl_discards, &
               npl => self%pl%nbody, ntp => self%tp%nbody, t => self%t, collision_history => self%collision_history, &
               collider => self%collider)
         lpl_discards = .false.
         ltp_discards = .false.
         if (lpl_check .and. pl%nbody > 0) then
            pl%ldiscard = pl%status(:) /= ACTIVE
            call pl%discard(nbody_system, param)
            if (npl > 0) lpl_discards = any(pl%ldiscard(1:npl))
         end if
            
         if (ltp_check .and. tp%nbody > 0) then
            tp%ldiscard = tp%status(:) /= ACTIVE
            call tp%discard(nbody_system, param)
            if (ntp > 0) ltp_discards = any(tp%ldiscard(1:ntp))
            if (npl > 0) lpl_discards = any(pl%ldiscard(1:npl))
         end if

         if (ltp_discards.or.lpl_discards) then
            ! Advance the collision id number and save it
            collider%maxid_collision = collider%maxid_collision + 1
            collider%collision_id = collider%maxid_collision
            collider%impactors%regime = COLLRESOLVE_REGIME_MERGE
            write(idstr,*) collider%collision_id
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, "collision_id " // trim(adjustl(idstr)))

            if (ltp_discards) then
               allocate(ldiscard, source=tp%ldiscard(:))
               do i = 1, ntp
                  if (ldiscard(i)) call tp%info(i)%set_value(collision_id=collider%collision_id)
               end do
               allocate(tpsub, mold=tp)
               call tp%spill(tpsub, ldiscard, ldestructive=.true.)
               nsub = tpsub%nbody
               nstart = tp_discards%nbody + 1
               nend = tp_discards%nbody + nsub
               call tp_discards%append(tpsub, lsource_mask=[(.true., i = 1, nsub)])
               deallocate(ldiscard)
               select type(before => collider%before)
               class is (swiftest_nbody_system)
                  if (allocated(before%tp)) deallocate(before%tp)
                  allocate(before%tp, source=tp_discards)
               end select
               call tp_discards%setup(0,param) 
            end if

            if (lpl_discards) then ! In the base integrators, massive bodies are not true discards. The discard is 
                                              ! simply used to trigger a snapshot.
               if (param%lenergy) call self%conservation_report(param, lterminal=.false.)
               allocate(ldiscard, source=pl%ldiscard(:))
               do i = 1, npl
                  if (ldiscard(i)) call pl%info(i)%set_value(collision_id=collider%collision_id)
               end do
               allocate(plsub, mold=pl)
               call pl%spill(plsub, ldiscard, ldestructive=.false.)
               nsub = plsub%nbody
               nstart = pl_discards%nbody + 1
               nend = pl_discards%nbody + nsub
               call pl_discards%append(plsub, lsource_mask=[(.true., i = 1, nsub)])
               deallocate(ldiscard)
               pl%ldiscard(1:npl) = .false.
               ! Save the before snapshots
               select type(before => collider%before)
               class is (swiftest_nbody_system)
                  if (allocated(before%pl)) deallocate(before%pl)
                  allocate(before%pl, source=pl_discards)
               end select
               call pl_discards%setup(0,param) 
            end if


            call collision_history%take_snapshot(param,nbody_system, t, "particle") 
         end if
         
      end associate

      return
   end subroutine swiftest_discard_system


   module subroutine swiftest_discard_pl(self, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Placeholder method for discarding massive bodies. This method does nothing except to ensure that the discard flag is set 
      !! to false. This method is intended to be overridden by more advanced integrators.
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameter

      if (self%nbody == 0) return
      self%ldiscard(1:self%nbody) = .false.

      return
   end subroutine swiftest_discard_pl


   module subroutine swiftest_discard_tp(self, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Check to see if particles should be discarded based on their positions relative to the massive bodies
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine: discard.f90
      !! Adapted from Hal Levison's Swift routine discard.
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameter

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, cb => nbody_system%cb, pl => nbody_system%pl, npl => nbody_system%pl%nbody, &
                tp_discards => nbody_system%tp_discards, pl_discards => nbody_system%pl_discards)

         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or. &
             (param%rmaxu >= 0.0_DP) .or. ((param%qmin >= 0.0_DP) .and. (param%qmin_coord == "BARY"))) then
            call pl%h2b(cb) 
            call tp%h2b(cb) 
         end if

         if ((param%rmin >= 0.0_DP) .or. &
             (param%rmax >= 0.0_DP) .or. & 
             (param%rmaxu >= 0.0_DP)) then
               call swiftest_discard_cb_tp(tp, nbody_system, param)
         end if
         if (param%qmin >= 0.0_DP) then
            call swiftest_discard_peri_tp(tp, nbody_system, param)
         end if
         if (param%lclose) then
            call swiftest_discard_pl_tp(tp, nbody_system, param)
         end if
      end associate

      return
   end subroutine swiftest_discard_tp


   subroutine swiftest_discard_cb_tp(tp, nbody_system, param)
      !! author: David A. Minton
      !!
      !!  Check to see if test particles should be discarded based on their positions relative to the Sun
      !!        or because they are unbound from the nbody_system
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_sun.f90
      !! Adapted from Hal Levison's Swift routine discard_sun.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp     !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)        :: i
      real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
      character(len=STRMAX) :: idstr, timestr, message

      associate(ntp => tp%nbody, cb => nbody_system%cb, Gmtot => nbody_system%Gmtot)
         rmin2 = max(param%rmin * param%rmin, cb%radius * cb%radius)
         rmax2 = param%rmax**2
         rmaxu2 = param%rmaxu**2
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               rh2 = dot_product(tp%rh(:, i), tp%rh(:, i))
               if ((param%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  tp%status(i) = DISCARDED_RMAX
                  write(idstr, *) tp%id(i)
                  write(timestr, *) nbody_system%t
                  write(message, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                              " too far from the central body at t = " // trim(adjustl(timestr))
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                  tp%ldiscard(i) = .true.
                  tp%lmask(i) = .false.
                  call tp%info(i)%set_value(status="DISCARDED_RMAX", discard_time=nbody_system%t, discard_rh=tp%rh(:,i), &
                                            discard_vh=tp%vh(:,i))
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  tp%status(i) = DISCARDED_RMIN
                  write(idstr, *) tp%id(i)
                  write(timestr, *) nbody_system%t
                  write(message, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                              " too close to the central body at t = " // trim(adjustl(timestr))
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                  tp%ldiscard(i) = .true.
                  tp%lmask(i) = .false.
                  call tp%info(i)%set_value(status="DISCARDED_RMIN", discard_time=nbody_system%t, discard_rh=tp%rh(:,i), &
                                            discard_vh=tp%vh(:,i), discard_body_id=cb%id)
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(tp%rb(:, i),  tp%rb(:, i))
                  vb2 = dot_product(tp%vb(:, i), tp%vb(:, i))
                  energy = 0.5_DP * vb2 - Gmtot / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     tp%status(i) = DISCARDED_RMAXU
                     write(idstr, *) tp%id(i)
                     write(timestr, *) nbody_system%t
                     write(message, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                                 " is unbound and too far from barycenter at t = " // trim(adjustl(timestr))
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                     tp%ldiscard(i) = .true.
                     tp%lmask(i) = .false.
                     call tp%info(i)%set_value(status="DISCARDED_RMAXU", discard_time=nbody_system%t, discard_rh=tp%rh(:,i), &
                                               discard_vh=tp%vh(:,i))
                  end if
               end if
            end if
         end do
      end associate

      return
   end subroutine swiftest_discard_cb_tp


   subroutine swiftest_discard_peri_tp(tp, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_peri.f90
      !! Adapted from Hal Levison's Swift routine discard_peri.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameterss
      ! Internals
      integer(I4B)              :: i, j, ih
      real(DP)                  :: r2
      real(DP), dimension(NDIM) :: dx
      character(len=STRMAX) :: idstr, timestr, message
   
      associate(cb => nbody_system%cb, ntp => tp%nbody, pl => nbody_system%pl, npl => nbody_system%pl%nbody, t => nbody_system%t)
         call tp%get_peri(nbody_system, param)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               if (tp%isperi(i) == 0) then
                  ih = 1
                  do j = 1, npl
                     dx(:) = tp%rh(:, i) - pl%rh(:, j)
                     r2 = dot_product(dx(:), dx(:))
                     if (r2 <= (pl%rhill(j))**2) ih = 0
                  end do
                  if (ih == 1) then
                     if ((tp%atp(i) >= param%qmin_alo) .and.    &
                        (tp%atp(i) <= param%qmin_ahi) .and.    &           
                        (tp%peri(i) <= param%qmin)) then
                        tp%status(i) = DISCARDED_PERI
                        write(idstr, *) tp%id(i)
                        write(timestr, *) nbody_system%t
                        write(message, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                                    " perihelion distance too small at t = " // trim(adjustl(timestr))
                        
                        call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                        tp%ldiscard(i) = .true.
                        call tp%info(i)%set_value(status="DISCARDED_PERI", discard_time=nbody_system%t, discard_rh=tp%rh(:,i), &
                                                  discard_vh=tp%vh(:,i), discard_body_id=cb%id)
                     end if
                  end if
               end if
            end if
         end do
      end associate

      return
   end subroutine swiftest_discard_peri_tp


   subroutine swiftest_discard_pl_tp(tp, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_pl.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp     !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals 
      integer(I4B)              :: i, j, isp
      real(DP)                  :: r2min, radius
      real(DP), dimension(NDIM) :: dx, dv
      character(len=STRMAX) :: idstri, idstrj, timestr, message
   
      associate(ntp => tp%nbody, pl => nbody_system%pl, npl => nbody_system%pl%nbody, t => nbody_system%t, dt => param%dt)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               do j = 1, npl
                  dx(:) = tp%rh(:, i) - pl%rh(:, j)
                  dv(:) = tp%vh(:, i) - pl%vh(:, j)
                  radius = pl%radius(j)
                  call swiftest_discard_pl_close(dx(:), dv(:), dt, radius**2, isp, r2min)
                  if (isp /= 0) then
                     tp%status(i) = DISCARDED_PLR
                     tp%lmask(i) = .false.
                     pl%ldiscard(j) = .true.
                     write(idstri, *) tp%id(i)
                     write(idstrj, *) pl%id(j)
                     write(timestr, *) nbody_system%t
                     write(message, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstri)) // ")" &
                                                   // "  too close to massive body " // trim(adjustl(pl%info(j)%name)) // " (" &
                                                   // trim(adjustl(idstrj)) // ")" &
                                                   // " at t = " // trim(adjustl(timestr))
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                     tp%ldiscard(i) = .true.
                     call tp%info(i)%set_value(status="DISCARDED_PLR", discard_time=nbody_system%t, discard_rh=tp%rh(:,i), &
                                               discard_vh=tp%vh(:,i), discard_body_id=pl%id(j))
                     exit
                  end if
               end do
            end if
         end do
      end associate

      return
   end subroutine swiftest_discard_pl_tp
   

   subroutine swiftest_discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
      !! author: David A. Minton
      !!
      !!  Check to see if a test particle and massive body are having, or will have within the next time step, an encounter such
      !!          that the separation distance r is less than some critical radius rcrit (or r**2 < rcrit**2 = r2crit)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl_close.f90
      !! Adapted from Hal Levison's Swift routine discard_pl_close.f
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(in)    :: dx, dv
      real(DP), intent(in)                  :: dt, r2crit
      integer(I4B), intent(out)             :: iflag
      real(DP), intent(out)                 :: r2min
      ! Internals
      real(DP) :: r2, v2, vdotr, tmin
      
      r2 = dot_product(dx(:), dx(:))
      if (r2 <= r2crit) then
         iflag = 1
      else
         vdotr = dot_product(dx(:), dv(:))
         if (vdotr > 0.0_DP) then
            iflag = 0
         else
            v2 = dot_product(dv(:), dv(:))
            tmin = -vdotr / v2
            if (tmin < dt) then
               r2min = r2 - vdotr * vdotr / v2
            else
               r2min = r2 + 2 * vdotr * dt + v2 * dt**2
            end if
            r2min = min(r2min, r2)
            if (r2min <= r2crit) then
               iflag = 1
            else
               iflag = 0
            end if
         end if
      end if
   
      return
   end subroutine swiftest_discard_pl_close

end submodule s_swiftest_discard
