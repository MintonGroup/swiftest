submodule (symba_classes) s_symba_collision
   use swiftest
contains
   module subroutine symba_collision_check_plplenc(self, system, param, t, dt, irec)
      !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
      !!
      !! Check for merger between massive bodies in SyMBA. If the user has turned on the FRAGMENTATION feature, it will call the 
      !! symba_regime subroutine to determine what kind of collision will occur.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      class(symba_plplenc),       intent(inout) :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t      !! current time
      real(DP),                   intent(in)    :: dt     !! step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
   end subroutine symba_collision_check_plplenc

   module subroutine symba_collision_check_pltpenc(self, system, param, t, dt, irec)
      !! author: David A. Minton
      !!
      !! Check for merger between massive bodies and test particles in SyMBA
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge_tp.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      class(symba_pltpenc),       intent(inout) :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t      !! current time
      real(DP),                   intent(in)    :: dt     !! step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
      ! Internals
      logical, dimension(:), allocatable        :: lcollision, mask
      real(DP), dimension(NDIM)                 :: xr, vr
      integer(I4B)                              :: k

      select type(pl => system%pl)
      class is (symba_pl)
         select type(tp => system%tp)
         class is (symba_tp)
            associate(pltpenc_list => self, npltpenc => self%nenc, plind => self%index1(1:self%nenc), tpind => self%index2(1:self%nenc))
               allocate(lcollision(npltpenc), mask(npltpenc))
               mask(:) = ((pltpenc_list%status(1:npltpenc) == ACTIVE) .and. (pl%levelg(plind) >= irec) .and. (tp%levelg(tpind) >= irec))
               lcollision(:) = .false.
               do concurrent(k = 1:npltpenc, mask(k))
                  associate(i => plind(k), j => tpind(k))
                     xr(:) = pl%xh(:, i) - tp%xh(:, j) 
                     vr(:) = pl%vb(:, i) - tp%vb(:, j)
                     lcollision(i) = symba_collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%Gmass(i), pl%radius(i), dt, pltpenc_list%lvdotr(k))
                  end associate
               end do

               if (any(lcollision(:))) then
                  where(lcollision(1:npltpenc))
                     pltpenc_list%status(1:npltpenc) = COLLISION
                     tp%status(tpind(1:npltpenc)) = DISCARDED_PLR
                     tp%ldiscard(tpind(1:npltpenc)) = .true.
                  end where
                  do k = 1, npltpenc
                     if (pltpenc_list%status(k) /= COLLISION) cycle
                     write(*,*) 'Test particle ',tp%id(tpind(k)), ' collided with massive body ',pl%id(plind(k)), ' at time ',t
                  end do
               end if
            end associate
         end select
      end select
      return 
   end subroutine symba_collision_check_pltpenc

   pure elemental function symba_collision_check_one(xr, yr, zr, vxr, vyr, vzr, Gmtot, rlim, dt, lvdotr) result(lcollision)
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
      ! Result
      logical              :: lcollision    !! Logical flag indicating whether these two bodies will collide or not
      ! Internals
      real(DP) :: r2, rlim2, a, e, q, vdotr, tcr2, dt2

      r2 = xr**2 + yr**2 + zr**2
      rlim2 = rlim**2
   
      if (r2 <= rlim2) then ! checks if bodies are actively colliding in this time step
         lcollision = .true.
      else ! if they are not actively colliding in  this time step, checks if they are going to collide next time step based on velocities and q 
         lcollision = .false.
         vdotr = xr * vxr + yr * vyr + zr * vzr
         if (lvdotr .and. (vdotr > 0.0_DP)) then 
            tcr2 = r2 / (vxr**2 + vyr**2 + vzr**2)
            dt2 = dt**2
            if (tcr2 <= dt2) then
               call orbel_xv2aeq(Gmtot, [xr, yr, zr], [vxr, vyr, vzr], a, e, q)
               lcollision = (q < rlim) 
            end if
         end if
      end if
      return
   end function symba_collision_check_one


end submodule s_symba_collision