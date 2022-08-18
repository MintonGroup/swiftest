submodule (swiftest_classes) s_discard
   use swiftest
contains

   module subroutine discard_system(self, param)
      !! author: David A. Minton
      !!
      !! Calls the discard methods for each body class and then the write method if any discards were detected
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      logical :: lany_discards, lpl_discards, ltp_discards, lpl_check, ltp_check

      lpl_check = allocated(self%pl_discards)
      ltp_check = allocated(self%tp_discards)

      associate(system => self, tp => self%tp, pl => self%pl, tp_discards => self%tp_discards, pl_discards => self%pl_discards)
         lpl_discards = .false.
         ltp_discards = .false.
         if (lpl_check) then
            call pl%discard(system, param)
            lpl_discards = (pl_discards%nbody > 0)
         end if
            
         if (ltp_check) then
            call tp%discard(system, param)
            ltp_discards = (tp_discards%nbody > 0)
         end if

         if (lpl_discards .or. ltp_discards) call system%write_discard(param)
         if (lpl_discards .and. param%lenergy) call self%conservation_report(param, lterminal=.false.)
         if (lpl_check) call pl_discards%setup(0,param) 
         if (ltp_check) call tp_discards%setup(0,param) 
         
      end associate

      return
   end subroutine discard_system


   module subroutine discard_pl(self, system, param)
      !! author: David A. Minton
      !!
      !!  Placeholder method for discarding massive bodies. This method does nothing except to ensure that the discard flag is set to false. 
      !!  This method is intended to be overridden by more advanced integrators.
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameter

      if (self%nbody == 0) return
      self%ldiscard(1:self%nbody) = .false.

      return
   end subroutine discard_pl


   module subroutine discard_tp(self, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if particles should be discarded based on their positions relative to the massive bodies
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine: discard.f90
      !! Adapted from Hal Levison's Swift routine discard.
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameter
      ! Internals
      logical, dimension(:), allocatable :: ldiscard
      integer(I4B) :: npl, ntp

      associate(tp => self, cb => system%cb, pl => system%pl)
         ntp = tp%nbody
         npl = pl%nbody

         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or. &
             (param%rmaxu >= 0.0_DP) .or. ((param%qmin >= 0.0_DP) .and. (param%qmin_coord == "BARY"))) then
            call pl%h2b(cb) 
            call tp%h2b(cb) 
         end if

         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or.  (param%rmaxu >= 0.0_DP)) call discard_cb_tp(tp, system, param)
         if (param%qmin >= 0.0_DP) call discard_peri_tp(tp, system, param)
         if (param%lclose) call discard_pl_tp(tp, system, param)
         if (any(tp%ldiscard(1:ntp))) then
            allocate(ldiscard, source=tp%ldiscard)
            call tp%spill(system%tp_discards, ldiscard(1:ntp), ldestructive=.true.)
         end if
      end associate

      return
   end subroutine discard_tp


   subroutine discard_cb_tp(tp, system, param)
      !! author: David A. Minton
      !!
      !!  Check to see if test particles should be discarded based on their positions relative to the Sun
      !!        or because they are unbound from the system
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_sun.f90
      !! Adapted from Hal Levison's Swift routine discard_sun.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp     !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)        :: i
      real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
      character(len=STRMAX) :: idstr, timestr

      associate(ntp => tp%nbody, cb => system%cb, Gmtot => system%Gmtot)
         rmin2 = max(param%rmin * param%rmin, cb%radius * cb%radius)
         rmax2 = param%rmax**2
         rmaxu2 = param%rmaxu**2
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               rh2 = dot_product(tp%xh(:, i), tp%xh(:, i))
               if ((param%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  tp%status(i) = DISCARDED_RMAX
                  write(idstr, *) tp%id(i)
                  write(timestr, *) param%t
                  write(*, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                              " too far from the central body at t = " // trim(adjustl(timestr))
                  tp%ldiscard(i) = .true.
                  tp%lmask(i) = .false.
                  call tp%info(i)%set_value(status="DISCARDED_RMAX", discard_time=param%t, discard_xh=tp%xh(:,i), &
                                            discard_vh=tp%vh(:,i))
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  tp%status(i) = DISCARDED_RMIN
                  write(idstr, *) tp%id(i)
                  write(timestr, *) param%t
                  write(*, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                              " too close to the central body at t = " // trim(adjustl(timestr))
                  tp%ldiscard(i) = .true.
                  tp%lmask(i) = .false.
                  call tp%info(i)%set_value(status="DISCARDED_RMIN", discard_time=param%t, discard_xh=tp%xh(:,i), &
                                            discard_vh=tp%vh(:,i), discard_body_id=cb%id)
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(tp%xb(:, i),  tp%xb(:, i))
                  vb2 = dot_product(tp%vb(:, i), tp%vb(:, i))
                  energy = 0.5_DP * vb2 - Gmtot / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     tp%status(i) = DISCARDED_RMAXU
                     write(idstr, *) tp%id(i)
                     write(timestr, *) param%t
                     write(*, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                                 " is unbound and too far from barycenter at t = " // trim(adjustl(timestr))
                     tp%ldiscard(i) = .true.
                     tp%lmask(i) = .false.
                     call tp%info(i)%set_value(status="DISCARDED_RMAXU", discard_time=param%t, discard_xh=tp%xh(:,i), &
                                               discard_vh=tp%vh(:,i))
                  end if
               end if
            end if
         end do
      end associate

      return
   end subroutine discard_cb_tp


   subroutine discard_peri_tp(tp, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_peri.f90
      !! Adapted from Hal Levison's Swift routine discard_peri.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameterss
      ! Internals
      logical, save             :: lfirst = .true.
      integer(I4B)              :: i, j, ih
      real(DP)                  :: r2
      real(DP), dimension(NDIM) :: dx
      character(len=STRMAX) :: idstr, timestr
   
      associate(cb => system%cb, ntp => tp%nbody, pl => system%pl, npl => system%pl%nbody, t => param%t)
         call tp%get_peri(system, param)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               if (tp%isperi(i) == 0) then
                  ih = 1
                  do j = 1, npl
                     dx(:) = tp%xh(:, i) - pl%xh(:, j)
                     r2 = dot_product(dx(:), dx(:))
                     if (r2 <= (pl%rhill(j))**2) ih = 0
                  end do
                  if (ih == 1) then
                     if ((tp%atp(i) >= param%qmin_alo) .and.    &
                        (tp%atp(i) <= param%qmin_ahi) .and.    &           
                        (tp%peri(i) <= param%qmin)) then
                        tp%status(i) = DISCARDED_PERI
                        write(idstr, *) tp%id(i)
                        write(timestr, *) param%t
                        write(*, *) "Particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstr)) // ")" // &
                                    " perihelion distance too small at t = " // trim(adjustl(timestr))
                        tp%ldiscard(i) = .true.
                        call tp%info(i)%set_value(status="DISCARDED_PERI", discard_time=param%t, discard_xh=tp%xh(:,i), &
                                                  discard_vh=tp%vh(:,i), discard_body_id=pl%id(j))
                     end if
                  end if
               end if
            end if
         end do
      end associate

      return
   end subroutine discard_peri_tp


   subroutine discard_pl_tp(tp, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_pl.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: tp     !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals 
      integer(I4B)              :: i, j, isp
      real(DP)                  :: r2min, radius
      real(DP), dimension(NDIM) :: dx, dv
      character(len=STRMAX) :: idstri, idstrj, timestr
   
      associate(ntp => tp%nbody, pl => system%pl, npl => system%pl%nbody, t => param%t, dt => param%dt)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               do j = 1, npl
                  dx(:) = tp%xh(:, i) - pl%xh(:, j)
                  dv(:) = tp%vh(:, i) - pl%vh(:, j)
                  radius = pl%radius(j)
                  call discard_pl_close(dx(:), dv(:), dt, radius**2, isp, r2min)
                  if (isp /= 0) then
                     tp%status(i) = DISCARDED_PLR
                     tp%lmask(i) = .false.
                     pl%ldiscard(j) = .true.
                     write(idstri, *) tp%id(i)
                     write(idstrj, *) pl%id(j)
                     write(timestr, *) param%t
                     write(*, *) "Test particle " // trim(adjustl(tp%info(i)%name)) // " ("  // trim(adjustl(idstri)) // ")" &
                                                  // "  too close to massive body " // trim(adjustl(pl%info(j)%name)) // " ("  // trim(adjustl(idstrj)) // ")" &
                                                  // " at t = " // trim(adjustl(timestr))
                     tp%ldiscard(i) = .true.
                     call tp%info(i)%set_value(status="DISCARDED_PLR", discard_time=param%t, discard_xh=tp%xh(:,i), &
                                               discard_vh=tp%vh(:,i), discard_body_id=pl%id(j))
                     exit
                  end if
               end do
            end if
         end do
      end associate

      return
   end subroutine discard_pl_tp
   

   subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
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
   end subroutine discard_pl_close

end submodule s_discard
