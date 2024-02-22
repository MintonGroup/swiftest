
! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (collision) s_collision_check
   use swiftest
   use symba, only : symba_pl, symba_tp
contains

   pure elemental subroutine collision_check_one(xr, yr, zr, vxr, vyr, vzr, Gmtot, rlim, dt, lvdotr, lcollision, lclosest)
      !! author: David A. Minton
      !! 
      !! Check for a merger between a single pair of particles
      !!
      !! Adapted from David E. Kaufmann's Swifter routines symba_merge_tp.f90 and symba_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      real(DP), intent(in) :: xr, yr, zr    !! Relative position vector components
      real(DP), intent(in) :: vxr, vyr, vzr !! Relative velocity vector components
      real(DP), intent(in) :: Gmtot         !! Sum of G*mass of colliding bodies
      real(DP), intent(in) :: rlim          !! Collision limit - Typically the sum of the radii of colliding bodies
      real(DP), intent(in) :: dt            !! Step size
      logical,  intent(in) :: lvdotr        !! Logical flag indicating that these two bodies are approaching in the current substep
      logical, intent(out) :: lcollision    !! Logical flag indicating whether these two bodies will collide or not
      logical, intent(out) :: lclosest      !! Logical flag indicating that, while not a collision, this is the closest approach for
                                            !!  this pair of bodies
      ! Internals
      real(DP) :: r2, rlim2, a, e, q, vdotr, tcr2, dt2

      r2 = xr**2 + yr**2 + zr**2
      rlim2 = rlim**2
      lclosest = .false. 
      if (r2 <= rlim2) then ! checks if bodies are actively colliding in this time step
         lcollision = .true.
      else ! if they are not actively colliding in  this time step, checks if they are going to collide next time step based on 
           ! velocities and q 
         lcollision = .false.
         vdotr = xr * vxr + yr * vyr + zr * vzr
         if (lvdotr .and. (vdotr > 0.0_DP)) then 
            tcr2 = r2 / (vxr**2 + vyr**2 + vzr**2)
            dt2 = dt**2
            if (tcr2 <= dt2) then
               call swiftest_orbel_xv2aeq(Gmtot, xr, yr, zr, vxr, vyr, vzr, a, e, q)
               lcollision = (q < rlim) 
            end if
            lclosest = .not. lcollision
         end if
      end if

      return
   end subroutine collision_check_one


   module subroutine collision_check_plpl(self, nbody_system, param, t, dt, irec, lany_collision)
      !! author: David A. Minton
      !!
      !! Check for merger between massive bodies and test particles in SyMBA
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge.f90 and symba_merge_tp.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      class(collision_list_plpl), intent(inout) :: self !! SyMBA pl-tp encounter list object
      class(base_nbody_system), intent(inout) :: nbody_system  !! SyMBA nbody system object
      class(base_parameters),intent(inout) :: param !! Current run configuration parameters 
      real(DP),  intent(in) :: t !! current time
      real(DP), intent(in) :: dt  !! step size
      integer(I4B), intent(in)  :: irec  !! Current recursion level
      logical, intent(out) :: lany_collision !! Returns true if cany pair of encounters resulted in a collision 
      ! Internals
      logical, dimension(:), allocatable        :: lcollision, lmask
      real(DP), dimension(NDIM)                 :: xr, vr
      integer(I4B)                              :: i, j
      integer(I8B)                              :: k, nenc
      real(DP)                                  :: rlim, Gmtot
      logical                                   :: lany_closest

      lany_collision = .false.
      if (self%nenc == 0) return

      select type(nbody_system)
      class is (swiftest_nbody_system)
         associate(pl => nbody_system%pl)

            nenc = self%nenc
            allocate(lmask(nenc))
            lmask(:) = (self%status(1:nenc) == ACTIVE)
            select type(pl)
            class is (symba_pl)
               lmask(:) = lmask(:).and. (pl%levelg(self%index1(1:nenc)) >= irec)
            end select
            if (.not.any(lmask(:))) return

            allocate(lcollision(nenc))
            lcollision(:) = .false.
            self%lclosest(:) = .false.

#ifdef DOCONLOC
            do concurrent(k = 1_I8B:nenc, lmask(k)) shared(self,pl,lmask, dt, lcollision) local(i,j,xr,vr,rlim,Gmtot)
#else
            do concurrent(k = 1_I8B:nenc, lmask(k))
#endif
               i = self%index1(k)
               j = self%index2(k)
               xr(:) = pl%rh(:, i) - pl%rh(:, j) 
               vr(:) = pl%vb(:, i) - pl%vb(:, j)
               rlim = pl%radius(i) + pl%radius(j)
               Gmtot = pl%Gmass(i) + pl%Gmass(j)
               call collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), Gmtot, rlim, dt, self%lvdotr(k), lcollision(k), &
                                        self%lclosest(k))
            end do

            lany_collision = any(lcollision(:))
            lany_closest = (param%lenc_save_closest .and. any(self%lclosest(:)))

            if (lany_collision .or. lany_closest) then
               call pl%rh2rb(nbody_system%cb) ! Update the central body barycenteric position vector to get us out of DH and into 
                                              ! barycentric coordinates
               do k = 1_I8B, nenc
                  if (.not.lcollision(k) .and. .not. self%lclosest(k)) cycle
                  i = self%index1(k)
                  j = self%index2(k)
                  self%r1(:,k) = pl%rh(:,i) + nbody_system%cb%rb(:)
                  self%v1(:,k) = pl%vb(:,i) 
                  if (lcollision(k)) then
                     self%status(k) = COLLIDED
                     self%tcollision(k) = t
                  end if
                  self%r2(:,k) = pl%rh(:,j) + nbody_system%cb%rb(:)
                  self%v2(:,k) = pl%vb(:,j) 
                  if (lcollision(k)) then
                     ! Check to see if either of these bodies has been involved with a collision before, and if so, make this a 
                     ! collider pair
                     if (pl%lcollision(i) .or. pl%lcollision(j)) call pl%make_impactors([i,j])

                     ! Set the collision flag for these to bodies to true in case they become involved in another collision later in
                     ! the step
                     pl%lcollision([i, j]) = .true.
                     pl%status([i, j]) = COLLIDED
                     call pl%info(i)%set_value(status="COLLIDED")
                     call pl%info(j)%set_value(status="COLLIDED")
                  end if

               end do

               ! Extract the pl-pl encounter list and return the pl-pl collision_list
               call self%extract_collisions(nbody_system, param)
            end if

            ! Take snapshots of pairs of bodies at close approach (but not collision) if requested
            if (lany_closest) call nbody_system%encounter_history%take_snapshot(param, nbody_system, t, "closest") 

         end associate
      end select
      return 
   end subroutine collision_check_plpl


   module subroutine collision_check_pltp(self, nbody_system, param, t, dt, irec, lany_collision)
      !! author: David A. Minton
      !!
      !! Check for merger between massive bodies and test particles in SyMBA
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge.f90 and symba_merge_tp.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      class(collision_list_pltp),intent(inout) :: self  !! SyMBA pl-tp encounter list object
      class(base_nbody_system), intent(inout) :: nbody_system !! SyMBA nbody system object
      class(base_parameters), intent(inout) :: param !! Current run configuration parameters 
      real(DP), intent(in) :: t !! current time
      real(DP), intent(in) :: dt  !! step size
      integer(I4B), intent(in)  :: irec  !! Current recursion level
      logical, intent(out)  :: lany_collision !! Returns true if cany pair of encounters resulted in a collision 
      ! Internals
      logical, dimension(:), allocatable        :: lcollision, lmask
      real(DP), dimension(NDIM)                 :: xr, vr
      integer(I4B)                              :: i, j
      integer(I8B)                              :: k, nenc
      logical                                   :: lany_closest
      character(len=STRMAX)                     :: timestr, idstri, idstrj, message
      class(collision_list_pltp), allocatable       :: tmp       

      lany_collision = .false.
      if (self%nenc == 0) return

      select type(nbody_system)
      class is (swiftest_nbody_system)
         associate(pl => nbody_system%pl, tp => nbody_system%tp)

            nenc = self%nenc
            allocate(lmask(nenc))
            lmask(:) = (self%status(1:nenc) == ACTIVE) 
            select type(pl)
            class is (symba_pl)
               select type(tp)
               class is (symba_tp)
                  lmask(:) = lmask(:) .and. (tp%levelg(self%index2(1:nenc)) >= irec)
               end select
            end select
            if (.not.any(lmask(:))) return

            allocate(lcollision(nenc))
            lcollision(:) = .false.
            self%lclosest(:) = .false.

#ifdef DOCONLOC
            do concurrent(k = 1_I8B:nenc, lmask(k)) shared(self,pl,tp,lmask, dt, lcollision) local(i,j,xr,vr)
#else
            do concurrent(k = 1_I8B:nenc, lmask(k))
#endif
               i = self%index1(k)
               j = self%index2(k)
               xr(:) = pl%rh(:, i) - tp%rh(:, j) 
               vr(:) = pl%vb(:, i) - tp%vb(:, j)
               call collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%Gmass(i), pl%radius(i), dt, self%lvdotr(k), &
                                        lcollision(k), self%lclosest(k))
            end do

            lany_collision = any(lcollision(:))
            lany_closest = (param%lenc_save_closest .and. any(self%lclosest(:)))


            if (lany_collision .or. lany_closest) then
               call pl%rh2rb(nbody_system%cb) ! Update the central body barycenteric position vector to get us out of DH and into 
                                              ! barycentric coordiantes
               do k = 1, nenc
                  if (.not.lcollision(k) .and. .not. self%lclosest(k)) cycle
                  i = self%index1(k)
                  j = self%index2(k)
                  self%r1(:,k) = pl%rh(:,i) + nbody_system%cb%rb(:)
                  self%v1(:,k) = pl%vb(:,i) 
                  if (lcollision(k)) then
                     self%status(k) = COLLIDED
                     self%tcollision(k) = t
                  end if

                  self%r2(:,k) = tp%rh(:,j) + nbody_system%cb%rb(:)
                  self%v2(:,k) = tp%vb(:,j) 
                  if (lcollision(k)) then
                     tp%status(j) = DISCARDED_PLR
                     tp%ldiscard(j) = .true.
                     write(idstri, *) pl%id(i)
                     write(idstrj, *) tp%id(j)
                     write(timestr, *) t
                     call tp%info(j)%set_value(status="DISCARDED_PLR", discard_time=t, discard_rh=tp%rh(:,j), discard_vh=tp%vh(:,j))
                     write(message, *) "Particle " // trim(adjustl(tp%info(j)%name)) // " ("  // trim(adjustl(idstrj)) // ")" &
                           //  " collided with massive body " // trim(adjustl(pl%info(i)%name)) // " (" // trim(adjustl(idstri)) &
                           // ")" //  " at t = " // trim(adjustl(timestr))
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
                  end if
               end do

               ! Extract the pl-tp encounter list and return the pl-tp collision_list
               call self%extract_collisions(nbody_system, param)
            end if

            ! Take snapshots of pairs of bodies at close approach (but not collision) if requested
            if (lany_closest) call nbody_system%encounter_history%take_snapshot(param, nbody_system, t, "closest") 
         end associate
      end select

      return 
   end subroutine collision_check_pltp

end submodule s_collision_check