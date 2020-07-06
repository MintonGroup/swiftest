submodule (swiftest_classes) discard_implementations
contains
   module procedure discard_system
      !! author: David A. Minton
      !!
      !! Check to see if particles should be discarded based on their positions relative to the massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard.f90
      !! Adapted from Hal Levison's Swift routine discard.f
      use swiftest
      implicit none
      logical, dimension(:), allocatable :: lspill_list

      real(DP) :: msys
      if (self%tp%nbody == 0) return 
      select type(self)
      class is (whm_nbody_system)
         associate(cb => self%cb, pl => self%pl, tp => self%tp, &
                   t => config%t, msys => self%msys, discards => self%tp_discards, &
                   ntp => self%tp%nbody)
            if ((config%rmin >= 0.0_DP) .or. (config%rmax >= 0.0_DP) .or. &
               (config%rmaxu >= 0.0_DP) .or. ((config%qmin >= 0.0_DP) .and. (config%qmin_coord == "BARY"))) then
                  call pl%h2b(cb) 
                  if (ntp > 0) call tp%h2b(cb) 
            end if
            if ((config%rmin >= 0.0_DP) .or. (config%rmax >= 0.0_DP) .or.  (config%rmaxu >= 0.0_DP)) then
               if (ntp > 0) call tp%discard_sun(cb, config, t, msys)
            end if
            if (config%qmin >= 0.0_DP .and. ntp > 0) call tp%discard_peri(cb, pl, config, t, msys)
            if (config%lclose .and. ntp > 0) call tp%discard_pl(cb, pl, config, t, msys)
          
            if (any(tp%ldiscard(1:ntp))) then
               ! Spill the discards to the spill list
               allocate(lspill_list, source = tp%ldiscard)
               call whm_discard_spill(tp, discards, lspill_list) 
               call self%write_discard(config, discards)
               deallocate(lspill_list)
            end if
         end associate  
         end select
      return
   end procedure discard_system

   module procedure discard_sun_tp 
      !! author: David A. Minton
      !!
      !!  Check to see if test particles should be discarded based on their positions relative to the Sun
      !!        or because they are unbound from the system
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_sun.f90
      !! Adapted from Hal Levison's Swift routine discard_sun.f
      use swiftest
      integer(I4B)        :: i
      real(DP)            :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2

      associate(tp => self, ntp => self%nbody)
         rmin2 = max(config%rmin * config%rmin, cb%radius * cb%radius)
         rmax2 = config%rmax * config%rmax
         rmaxu2 = config%rmaxu * config%rmaxu
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               rh2 = tp%xh(i, :) .dot. tp%xh(i, :) 
               if ((config%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  tp%status(i) = DISCARDED_RMAX
                  write(*, *) "Particle ", tp%name(i), " too far from sun at t = ", t
                  tp%ldiscard(i) = .true.
               else if ((config%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  tp%status(i) = DISCARDED_RMIN
                  write(*, *) "Particle ", tp%name(i), " too close to sun at t = ", t
                  tp%ldiscard(i) = .true.
               else if (config%rmaxu >= 0.0_DP) then
                  rb2 = tp%xb(i, :) .dot. tp%xb(i,:) 
                  vb2 = tp%vb(i, :) .dot. tp%vb(i, :) 
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
   
   end procedure discard_sun_tp

   module procedure discard_peri_tp
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_peri.f90
      !! Adapted from Hal Levison's Swift routine discard_peri.f
      use swiftest
      logical, save             :: lfirst = .true.
      integer(I4B)              :: i, j, ih, ntp, npl
      real(DP)                  :: r2
      real(DP), dimension(NDIM) :: dx
   
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody)
         call util_hills(npl, pl)
         if (lfirst) then
            call util_peri(lfirst, ntp, tp, cb%Gmass, msys, config%qmin_coord)
            lfirst = .false.
         else
            call util_peri(lfirst, ntp, tp, pl%mass(1), msys, config%qmin_coord)
            do i = 1, ntp
               if (tp%status(i) == ACTIVE) then
                  if (tp%isperi(i) == 0) then
                     ih = 1
                     do j = 2, npl
                        dx(:) = tp%xh(:,i) - pl%xh(:,j)
                        r2 = dx(:) .dot. dx(:) 
                        if (r2 <= pl%rhill(j) * pl%rhill(j)) ih = 0
                     end do
                     if (ih == 1) then
                        if ((tp%atp(i) >= config%qmin_alo) .and.    &
                           (tp%atp(i) <= config%qmin_ahi) .and.    &           
                           (tp%peri(i) <= config%qmin)) then
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
   
   end procedure discard_peri_tp

   module procedure discard_pl_tp
      !! author: David A. Minton
      !!
      !! Check to see if test particles should be discarded based on their positions relative to the massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_pl.f
      use swiftest
      implicit none
   
      integer(I4B)              :: i, j, isp, ntp, npl
      real(DP)                  :: r2min, radius
      real(DP), dimension(NDIM) :: dx, dv
   
      associate(tp => self, ntp => self%nbody, npl => pl%nbody)
         do i = 1, ntp
            if (tp%status(i) == ACTIVE) then
               do j = 2, npl
                  dx(:) = tp%xh(i, :) - pl%xh(j, :)
                  dv(:) = tp%vh(i, :) - pl%vh(j, :)
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
   
   end procedure discard_pl_tp

   module procedure discard_pl_close
      !! author: David A. Minton
      !!
      !!  Check to see if a test particle and massive body are having, or will have within the next time step, an encounter such
      !!          that the separation distance r is less than some critical radius rcrit (or r**2 < rcrit**2 = r2crit)
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: discard_pl_close.f90
      !! Adapted from Hal Levison's Swift routine discard_pl_close.f
      use swiftest
      real(DP) :: r2, v2, vdotr, tmin
   
      r2 = dx(:) .dot. dx(:) 
      if (r2 <= r2crit) then
         iflag = 1
      else
         vdotr = dx(:) .dot. dv(:) 
         if (vdotr > 0.0_DP) then
            iflag = 0
         else
            v2 = dv(:) .dot. dv(:) 
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
   
      end procedure discard_pl_close

end submodule discard_implementations
