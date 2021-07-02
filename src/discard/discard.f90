submodule (swiftest_classes) s_discard
   use swiftest
contains
   module subroutine discard_system(self, param)
      !! author: David A. Minton
      !!
      !! Check to see if particles should be discarded based on their positions relative to the massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard.f90
      !! Adapted from Hal Levison's Swift routine discard.f
      implicit none
      class(swiftest_nbody_system),  intent(inout) :: self    !! Swiftest system object
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters of parameters

      if (self%tp%nbody == 0) return 
      select type(self)
      class is (whm_nbody_system)
         associate(cb => self%cb, pl => self%pl, tp => self%tp, t => param%t, dt => param%dt, &
                   msys => self%msys, discards => self%tp_discards, &
                   ntp => self%tp%nbody, ldiscard => self%tp%ldiscard)
            if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or. &
               (param%rmaxu >= 0.0_DP) .or. ((param%qmin >= 0.0_DP) .and. (param%qmin_coord == "BARY"))) then
                  call pl%h2b(cb) 
                  if (ntp > 0) call tp%h2b(cb) 
            end if
            if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or.  (param%rmaxu >= 0.0_DP)) then
               if (ntp > 0) call tp%discard_sun(cb, param, t, msys)
            end if
            if (param%qmin >= 0.0_DP .and. ntp > 0) call tp%discard_peri(cb, pl, param, t, msys)
            if (param%lclose .and. ntp > 0) call tp%discard_pl(pl, t, dt)
          
            if (any(tp%ldiscard(1:ntp))) then
               ! Spill the discards to the spill list
               call tp%spill(discards, ldiscard)
               call self%write_discard(param, discards)
            end if
         end associate  
         end select
      return
   end subroutine discard_system

   module subroutine discard_sun_tp(self, cb, param, t, msys)
      !! author: David A. Minton
      !!
      !!  Check to see if test particles should be discarded based on their positions relative to the Sun
      !!        or because they are unbound from the system
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_sun.f90
      !! Adapted from Hal Levison's Swift routine discard_sun.f
      implicit none
      ! Arguments
      class(swiftest_tp),            intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body object
      class(swiftest_parameters), intent(in)    :: param !! parameters
      real(DP),                      intent(in)    :: t      !! Current simulation tim
      real(DP),                      intent(in)    :: msys   !! Total system mass
      ! Internals
      integer(I4B)        :: i
      real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2

      associate(tp => self, ntp => self%nbody)
         rmin2 = max(param%rmin * param%rmin, cb%radius * cb%radius)
         rmax2 = param%rmax**2
         rmaxu2 = param%rmaxu**2
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               rh2 = dot_product(tp%xh(:, i), tp%xh(:, i))
               if ((param%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  tp%status(i) = DISCARDED_RMAX
                  write(*, *) "Particle ", tp%name(i), " too far from sun at t = ", t
                  tp%ldiscard(i) = .true.
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  tp%status(i) = DISCARDED_RMIN
                  write(*, *) "Particle ", tp%name(i), " too close to sun at t = ", t
                  tp%ldiscard(i) = .true.
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(tp%xb(:, i),  tp%xb(:, i))
                  vb2 = dot_product(tp%vb(:, i), tp%vb(:, i))
                  energy = 0.5_DP * vb2 - msys / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     tp%status(i) = DISCARDED_RMAXU
                     write(*, *) "Particle ", tp%name(i), " is unbound and too far from barycenter at t = ", t
                     tp%ldiscard(i) = .true.
                  end if
               end if
            end if
         end do
      end associate

      return
   end subroutine discard_sun_tp

   module subroutine discard_peri_tp(self, cb, pl, param, t, msys)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_peri.f90
      !! Adapted from Hal Levison's Swift routine discard_peri.f
      implicit none
      ! Arguments
      class(swiftest_tp),            intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body object
      class(swiftest_pl),            intent(inout) :: pl     !! Swiftest massive body object
      class(swiftest_parameters), intent(in)    :: param !! parameters
      real(DP),                      intent(in)    :: t      !! Current simulation tim
      real(DP),                      intent(in)    :: msys   !! Total system mass 
      ! Internals
      logical, save             :: lfirst = .true.
      integer(I4B)              :: i, j, ih
      real(DP)                  :: r2
      real(DP), dimension(NDIM) :: dx
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody, qmin_coord => param%qmin_coord)
         if (lfirst) then
            call util_hills(npl, pl)
            call util_peri(lfirst, ntp, tp, cb%Gmass, msys, param%qmin_coord)
            lfirst = .false.
         else
            call util_peri(lfirst, ntp, tp, cb%Gmass, msys, param%qmin_coord)
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
                           write(*, *) "Particle ", tp%name(i), " perihelion distance too small at t = ", t
                           tp%ldiscard(i) = .true.
                        end if
                     end if
                  end if
               end if
            end do
         end if
      end associate
      return
   
   end subroutine discard_peri_tp

   module subroutine discard_pl_tp(self, pl, t, dt)
      !! author: David A. Minton
      !!
      !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_pl.f
      implicit none
      ! Arguments
      class(swiftest_tp),            intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_pl),            intent(inout) :: pl     !! Swiftest massive body object
      real(DP),                      intent(in)    :: t      !! Current simulation tim
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals 
      integer(I4B)              :: i, j, isp
      real(DP)                  :: r2min, radius
      real(DP), dimension(NDIM) :: dx, dv
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               do j = 1, npl
                  dx(:) = tp%xh(:, i) - pl%xh(:, j)
                  dv(:) = tp%vh(:, i) - pl%vh(:, j)
                  radius = pl%radius(j)
                  call discard_pl_close(dx(:), dv(:), dt, radius**2, isp, r2min)
                  if (isp /= 0) then
                     tp%status(i) = DISCARDED_PLR
                     pl%ldiscard(j) = .true.
                     write(*, *) "Particle ", tp%name(i), " too close to massive body ", pl%name(j), " at t = ", t
                     tp%ldiscard(i) = .true.
                     exit
                  end if
               end do
            end if
         end do
   
      end associate
      return
   
   end subroutine discard_pl_tp

   module subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
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
