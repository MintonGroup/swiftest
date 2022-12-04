!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (symba_classes) s_symba_collision
   use swiftest

contains

   module function symba_collision_casedisruption(system, param, colliders, frag)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic disruption collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system), intent(inout) :: system    !! SyMBA nbody system object
      class(symba_parameters),   intent(inout) :: param     !! Current run configuration parameters with SyMBA additions
      class(fraggle_colliders),  intent(inout) :: colliders !! Fraggle colliders object        
      class(fraggle_fragments),  intent(inout) :: frag      !! Fraggle fragmentation system object 
      ! Result
      integer(I4B)                             :: status    !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)          :: i, nfrag
      logical               :: lfailure
      character(len=STRMAX) :: message 

      select case(frag%regime)
      case(COLLRESOLVE_REGIME_DISRUPTION)
         message = "Disruption between"
      case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
         message = "Supercatastrophic disruption between"
      end select
      call symba_collision_collider_message(system%pl, colliders%idx, message)
      call io_log_one_message(FRAGGLE_LOG_OUT, message)

      ! Collisional fragments will be uniformly distributed around the pre-impact barycenter
      call frag%set_mass_dist(colliders, param)

      ! Generate the position and velocity distributions of the fragments
      call frag%generate_fragments(colliders, system, param, lfailure)

      if (lfailure) then
         call io_log_one_message(FRAGGLE_LOG_OUT, "No fragment solution found, so treat as a pure hit-and-run")
         status = ACTIVE 
         nfrag = 0
         select type(pl => system%pl)
         class is (symba_pl)
            pl%status(colliders%idx(:)) = status
            pl%ldiscard(colliders%idx(:)) = .false.
            pl%lcollision(colliders%idx(:)) = .false.
         end select
      else
         ! Populate the list of new bodies
         nfrag = frag%nbody
         write(message, *) nfrag
         call io_log_one_message(FRAGGLE_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")
         select case(frag%regime)
         case(COLLRESOLVE_REGIME_DISRUPTION)
            status = DISRUPTION
         case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
            status = SUPERCATASTROPHIC
         end select
         frag%id(1:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag)]
         param%maxid = frag%id(nfrag)
         call symba_collision_mergeaddsub(system, param, colliders, frag, status)
      end if

      return
   end function symba_collision_casedisruption


   module function symba_collision_casehitandrun(system, param, colliders, frag)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system), intent(inout) :: system    !! SyMBA nbody system object
      class(symba_parameters),   intent(inout) :: param     !! Current run configuration parameters with SyMBA additions
      class(fraggle_colliders),  intent(inout) :: colliders !! Fraggle colliders object        
      class(fraggle_fragments),  intent(inout) :: frag      !! Fraggle fragmentation system object 
      ! Result
      integer(I4B)                             :: status    !! Status flag assigned to this outcom
      ! Internals
      integer(I4B)                            :: i, ibiggest, nfrag, jtarg, jproj
      logical                                 :: lpure 
      character(len=STRMAX) :: message

      message = "Hit and run between"
      call symba_collision_collider_message(system%pl, colliders%idx, message)
      call io_log_one_message(FRAGGLE_LOG_OUT, trim(adjustl(message)))

      if (colliders%mass(1) > colliders%mass(2)) then
         jtarg = 1
         jproj = 2
      else
         jtarg = 2
         jproj = 1
      end if

      if (frag%mass_dist(2) > 0.9_DP * colliders%mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
         call io_log_one_message(FRAGGLE_LOG_OUT, "Pure hit and run. No new fragments generated.")
         nfrag = 0
         lpure = .true.
      else ! Imperfect hit and run, so we'll keep the largest body and destroy the other
         lpure = .false.
         call frag%set_mass_dist(colliders, param)

         ! Generate the position and velocity distributions of the fragments
         call frag%generate_fragments(colliders, system, param, lpure)

         if (lpure) then
            call io_log_one_message(FRAGGLE_LOG_OUT, "Should have been a pure hit and run instead")
            nfrag = 0
         else
            nfrag = frag%nbody
            write(message, *) nfrag
            call io_log_one_message(FRAGGLE_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")
         end if
      end if
      if (lpure) then ! Reset these bodies back to being active so that nothing further is done to them
         status = HIT_AND_RUN_PURE
         select type(pl => system%pl)
         class is (symba_pl)
            pl%status(colliders%idx(:)) = ACTIVE
            pl%ldiscard(colliders%idx(:)) = .false.
            pl%lcollision(colliders%idx(:)) = .false.
         end select
      else
         ibiggest = colliders%idx(maxloc(system%pl%Gmass(colliders%idx(:)), dim=1))
         frag%id(1) = system%pl%id(ibiggest)
         frag%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
         param%maxid = frag%id(nfrag)
         status = HIT_AND_RUN_DISRUPT
         call symba_collision_mergeaddsub(system, param, colliders, frag, status)
      end if

      return
   end function symba_collision_casehitandrun


   module function symba_collision_casemerge(system, param, colliders, frag)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Merge massive bodies.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routines symba_merge_pl.f90 and symba_discard_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routines symba5_merge.f and discard_mass_merge.f
      implicit none
      ! Arguments
      class(symba_nbody_system), intent(inout) :: system    !! SyMBA nbody system object
      class(symba_parameters),   intent(inout) :: param     !! Current run configuration parameters with SyMBA additions
      class(fraggle_colliders),  intent(inout) :: colliders !! Fraggle colliders object        
      class(fraggle_fragments),  intent(inout) :: frag      !! Fraggle fragmentation system object 
      ! Result
      integer(I4B)                             :: status    !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                              :: i, j, k, ibiggest
      real(DP)                                  :: pe
      real(DP), dimension(NDIM)                 :: L_spin_new
      character(len=STRMAX) :: message

      message = "Merging"
      call symba_collision_collider_message(system%pl, colliders%idx, message)
      call io_log_one_message(FRAGGLE_LOG_OUT, message)

      select type(pl => system%pl)
      class is (symba_pl)

         call frag%set_mass_dist(colliders, param)
         ibiggest = colliders%idx(maxloc(pl%Gmass(colliders%idx(:)), dim=1))
         frag%id(1) = pl%id(ibiggest)
         frag%xb(:,1) = frag%xbcom(:)
         frag%vb(:,1) = frag%vbcom(:)

         if (param%lrotation) then
            ! Conserve angular momentum by putting pre-impact orbital momentum into spin of the new body
            L_spin_new(:) = colliders%L_orbit(:,1) + colliders%L_orbit(:,2) + colliders%L_spin(:,1) + colliders%L_spin(:,2)  

            ! Assume prinicpal axis rotation on 3rd Ip axis
            frag%rot(:,1) = L_spin_new(:) / (frag%Ip(3,1) * frag%mass(1) * frag%radius(1)**2)
         else ! If spin is not enabled, we will consider the lost pre-collision angular momentum as "escaped" and add it to our bookkeeping variable
            param%Lescape(:) = param%Lescape(:) + colliders%L_orbit(:,1) + colliders%L_orbit(:,2) 
         end if

         ! Keep track of the component of potential energy due to the pre-impact colliders%idx for book-keeping
         pe = 0.0_DP
         do j = 1, colliders%ncoll
            do i = j + 1, colliders%ncoll
               pe = pe - pl%Gmass(i) * pl%mass(j) / norm2(pl%xb(:, i) - pl%xb(:, j))
            end do
         end do
         system%Ecollisions = system%Ecollisions + pe 
         system%Euntracked  = system%Euntracked - pe 

         ! Update any encounter lists that have the removed bodies in them so that they instead point to the new 
         do k = 1, system%plplenc_list%nenc
            do j = 1, colliders%ncoll
               i = colliders%idx(j)
               if (i == ibiggest) cycle
               if (system%plplenc_list%id1(k) == pl%id(i)) then
                  system%plplenc_list%id1(k) = pl%id(ibiggest)
                  system%plplenc_list%index1(k) = i
               end if
               if (system%plplenc_list%id2(k) == pl%id(i)) then
                  system%plplenc_list%id2(k) = pl%id(ibiggest)
                  system%plplenc_list%index2(k) = i
               end if
               if (system%plplenc_list%id1(k) == system%plplenc_list%id2(k)) system%plplenc_list%status(k) = INACTIVE
            end do
         end do

         status = MERGED
         
         call symba_collision_mergeaddsub(system, param, colliders, frag, status) 

      end select

      return 
   end function symba_collision_casemerge


   subroutine symba_collision_collider_message(pl, collidx, collider_message)
      !! author: David A. Minton
      !!
      !! Prints a nicely formatted message about which bodies collided, including their names and ids.
      !! This subroutine appends the body names and ids to an input message.
      implicit none
      ! Arguments
      class(swiftest_pl),            intent(in)    :: pl            !! Swiftest massive body object
      integer(I4B),    dimension(:), intent(in)    :: collidx           !! Index of collisional colliders%idx members
      character(*),                  intent(inout) :: collider_message !! The message to print to the screen.
      ! Internals
      integer(I4B) :: i, n
      character(len=STRMAX) :: idstr
      
      n = size(collidx)
      if (n == 0) return

      do i = 1, n
         if (i > 1) collider_message = trim(adjustl(collider_message)) // " and "
         collider_message = " " // trim(adjustl(collider_message)) // " " // trim(adjustl(pl%info(collidx(i))%name))
         write(idstr, '(I10)') pl%id(collidx(i))
         collider_message = trim(adjustl(collider_message)) // " (" // trim(adjustl(idstr)) // ") "
      end do

      return
   end subroutine symba_collision_collider_message


   module function symba_collision_check_encounter(self, system, param, t, dt, irec) result(lany_collision)
      !! author: David A. Minton
      !!
      !! Check for merger between massive bodies and test particles in SyMBA
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge.f90 and symba_merge_tp.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      class(symba_encounter),     intent(inout) :: self           !! SyMBA pl-tp encounter list object
      class(symba_nbody_system),  intent(inout) :: system         !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param          !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t              !! current time
      real(DP),                   intent(in)    :: dt             !! step size
      integer(I4B),               intent(in)    :: irec           !! Current recursion level
      ! Result
      logical                                   :: lany_collision !! Returns true if cany pair of encounters resulted in a collision 
      ! Internals
      logical, dimension(:), allocatable        :: lcollision, lmask
      real(DP), dimension(NDIM)                 :: xr, vr
      integer(I4B)                              :: i, j, k, nenc
      real(DP)                                  :: rlim, Gmtot
      logical                                   :: isplpl
      character(len=STRMAX)                     :: timestr, idstri, idstrj, message
      class(symba_encounter), allocatable       :: tmp       

      lany_collision = .false.
      if (self%nenc == 0) return

      select type(self)
      class is (symba_plplenc)
         isplpl = .true.
      class default
         isplpl = .false.
      end select

      select type(pl => system%pl)
      class is (symba_pl)
         select type(tp => system%tp)
         class is (symba_tp)
            nenc = self%nenc
            allocate(lmask(nenc))
            lmask(:) = ((self%status(1:nenc) == ACTIVE) .and. (pl%levelg(self%index1(1:nenc)) >= irec))
            if (isplpl) then
               lmask(:) = lmask(:) .and. (pl%levelg(self%index2(1:nenc)) >= irec)
            else
               lmask(:) = lmask(:) .and. (tp%levelg(self%index2(1:nenc)) >= irec)
            end if
            if (.not.any(lmask(:))) return

            allocate(lcollision(nenc))
            lcollision(:) = .false.

            if (isplpl) then
               do concurrent(k = 1:nenc, lmask(k))
                  i = self%index1(k)
                  j = self%index2(k)
                  xr(:) = pl%rh(:, i) - pl%rh(:, j) 
                  vr(:) = pl%vb(:, i) - pl%vb(:, j)
                  rlim = pl%radius(i) + pl%radius(j)
                  Gmtot = pl%Gmass(i) + pl%Gmass(j)
                  lcollision(k) = symba_collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), &
                                                            Gmtot, rlim, dt, self%lvdotr(k))
               end do
            else
               do concurrent(k = 1:nenc, lmask(k))
                  i = self%index1(k)
                  j = self%index2(k)
                  xr(:) = pl%rh(:, i) - tp%rh(:, j) 
                  vr(:) = pl%vb(:, i) - tp%vb(:, j)
                  lcollision(k) = symba_collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), &
                                                            pl%Gmass(i), pl%radius(i), dt, self%lvdotr(k))
               end do
            end if

            lany_collision = any(lcollision(:))
            if (lany_collision) then
               call pl%xh2xb(system%cb) ! Update the central body barycenteric position vector to get us out of DH and into bary
               do k = 1, nenc
                  i = self%index1(k)
                  j = self%index2(k)
                  if (lcollision(k)) self%status(k) = COLLISION
                  self%tcollision(k) = t
                  self%x1(:,k) = pl%rh(:,i) + system%cb%xb(:)
                  self%v1(:,k) = pl%vb(:,i) 
                  if (isplpl) then
                     self%x2(:,k) = pl%rh(:,j) + system%cb%xb(:)
                     self%v2(:,k) = pl%vb(:,j) 
                     if (lcollision(k)) then
                        ! Check to see if either of these bodies has been involved with a collision before, and if so, make this a collisional colliders%idx
                        if (pl%lcollision(i) .or. pl%lcollision(j)) call pl%make_colliders([i,j])
   
                        ! Set the collision flag for these to bodies to true in case they become involved in another collision later in the step
                        pl%lcollision([i, j]) = .true.
                        pl%status([i, j]) = COLLISION
                        call pl%info(i)%set_value(status="COLLISION", discard_time=t, discard_rh=pl%rh(:,i), discard_vh=pl%vh(:,i))
                        call pl%info(j)%set_value(status="COLLISION", discard_time=t, discard_rh=pl%rh(:,j), discard_vh=pl%vh(:,j))
                     end if
                  else
                     self%x2(:,k) = tp%rh(:,j) + system%cb%xb(:)
                     self%v2(:,k) = tp%vb(:,j) 
                     if (lcollision(k)) then
                        tp%status(j) = DISCARDED_PLR
                        tp%ldiscard(j) = .true.
                        write(idstri, *) pl%id(i)
                        write(idstrj, *) tp%id(j)
                        write(timestr, *) t
                        call tp%info(j)%set_value(status="DISCARDED_PLR", discard_time=t, discard_rh=tp%rh(:,j), discard_vh=tp%vh(:,j))
                        write(message, *) "Particle " // trim(adjustl(tp%info(j)%name)) // " ("  // trim(adjustl(idstrj)) // ")" &
                              //  " collided with massive body " // trim(adjustl(pl%info(i)%name)) // " (" // trim(adjustl(idstri)) // ")" &
                              //  " at t = " // trim(adjustl(timestr))
                        call io_log_one_message(FRAGGLE_LOG_OUT, message)
                     end if
                  end if
               end do
            end if
         end select
      end select

      ! Extract the pl-pl or pl-tp encounter list and return the pl-pl or pl-tp collision_list
      if (lany_collision) then
         select type(self)
         class is (symba_plplenc)
            call self%extract_collisions(system, param)
         class default
            allocate(tmp, mold=self)
            call self%spill(tmp, lcollision, ldestructive=.true.) ! Remove this encounter pair from the encounter list
         end select
      end if

      return 
   end function symba_collision_check_encounter


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
               call orbel_xv2aeq(Gmtot, xr, yr, zr, vxr, vyr, vzr, a, e, q)
               lcollision = (q < rlim) 
            end if
         end if
      end if

      return
   end function symba_collision_check_one


   function symba_collision_consolidate_colliders(pl, cb, param, idx_parent, colliders) result(lflag)
      !! author: David A. Minton
      !! 
      !! Loops through the pl-pl collision list and groups families together by index. Outputs the indices of all colliders%idx members, 
      !! and pairs of quantities (x and v vectors, mass, radius, L_spin, and Ip) that can be used to resolve the collisional outcome.
      implicit none
      ! Arguments
      class(symba_pl),                                 intent(inout) :: pl               !! SyMBA massive body object
      class(symba_cb),                                 intent(inout) :: cb               !! SyMBA central body object
      class(symba_parameters),                         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(2),                   intent(inout) :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      class(fraggle_colliders),                        intent(out)   :: colliders
      ! Result
      logical                                                        :: lflag            !! Logical flag indicating whether a colliders%idx was successfully created or not
      ! Internals
      type collidx_array
         integer(I4B), dimension(:), allocatable :: id
         integer(I4B), dimension(:), allocatable :: idx
      end type collidx_array
      type(collidx_array), dimension(2) :: parent_child_index_array
      integer(I4B), dimension(2)       :: nchild
      integer(I4B)                     :: i, j, ncolliders, idx_child
      real(DP), dimension(2)           :: volume, density
      real(DP)                         :: mchild, volchild
      real(DP), dimension(NDIM)        :: xc, vc, xcom, vcom, xchild, vchild, xcrossv
      real(DP), dimension(NDIM,2)      :: mxc, vcc

      nchild(:) = pl%kin(idx_parent(:))%nchild 
      ! If all of these bodies share a parent, but this is still a unique collision, move the last child
      ! out of the parent's position and make it the secondary body
      if (idx_parent(1) == idx_parent(2)) then
         if (nchild(1) == 0) then ! There is only one valid body recorded in this pair (this could happen due to restructuring of the kinship relationships, though it should be rare)
            lflag = .false. 
            call pl%reset_kinship([idx_parent(1)])
            return
         end if
         idx_parent(2) = pl%kin(idx_parent(1))%child(nchild(1))
         nchild(1) = nchild(1) - 1
         nchild(2) = 0
         pl%kin(idx_parent(:))%nchild = nchild(:)
         pl%kin(idx_parent(2))%parent = idx_parent(1)
      end if

      colliders%mass(:) = pl%mass(idx_parent(:)) ! Note: This is meant to mass, not G*mass, as the collisional regime determination uses mass values that will be converted to Si
      colliders%radius(:) = pl%radius(idx_parent(:))
      volume(:) =  (4.0_DP / 3.0_DP) * PI * colliders%radius(:)**3
 
      ! Group together the ids and indexes of each collisional parent and its children
      do j = 1, 2
         allocate(parent_child_index_array(j)%idx(nchild(j)+ 1))
         allocate(parent_child_index_array(j)%id(nchild(j)+ 1))
         associate(idx_arr => parent_child_index_array(j)%idx, &
                   id_arr => parent_child_index_array(j)%id, &
                   ncj => nchild(j), &
                   pl => pl, &
                   plkinj => pl%kin(idx_parent(j)))
            idx_arr(1) = idx_parent(j)
            if (ncj > 0) idx_arr(2:ncj + 1) = plkinj%child(1:ncj)
            id_arr(:) = pl%id(idx_arr(:))
         end associate
      end do

      ! Consolidate the groups of collsional parents with any children they may have into a single "colliders%idx" index array
      ncolliders = 2 + sum(nchild(:))
      allocate(colliders%idx(ncolliders))
      colliders%idx = [parent_child_index_array(1)%idx(:),parent_child_index_array(2)%idx(:)]

      colliders%ncoll = count(pl%lcollision(colliders%idx(:)))
      colliders%idx = pack(colliders%idx(:), pl%lcollision(colliders%idx(:)))
      colliders%L_spin(:,:) = 0.0_DP
      colliders%Ip(:,:) = 0.0_DP

      ! Find the barycenter of each body along with its children, if it has any
      do j = 1, 2
         colliders%xb(:, j)  = pl%rh(:, idx_parent(j)) + cb%xb(:)
         colliders%vb(:, j)  = pl%vb(:, idx_parent(j))
         ! Assume principal axis rotation about axis corresponding to highest moment of inertia (3rd Ip)
         if (param%lrotation) then
            colliders%Ip(:, j) = colliders%mass(j) * pl%Ip(:, idx_parent(j))
            colliders%L_spin(:, j) = colliders%Ip(3, j) * colliders%radius(j)**2 * pl%rot(:, idx_parent(j))
         end if

         if (nchild(j) > 0) then
            do i = 1, nchild(j) ! Loop over all children and take the mass weighted mean of the properties
               idx_child = parent_child_index_array(j)%idx(i + 1)
               if (.not. pl%lcollision(idx_child)) cycle
               mchild = pl%mass(idx_child)
               xchild(:) = pl%rh(:, idx_child) + cb%xb(:)
               vchild(:) = pl%vb(:, idx_child)
               volchild = (4.0_DP / 3.0_DP) * PI * pl%radius(idx_child)**3
               volume(j) = volume(j) + volchild
               ! Get angular momentum of the child-parent pair and add that to the spin
               ! Add the child's spin
               if (param%lrotation) then
                  xcom(:) = (colliders%mass(j) * colliders%xb(:,j) + mchild * xchild(:)) / (colliders%mass(j) + mchild)
                  vcom(:) = (colliders%mass(j) * colliders%vb(:,j) + mchild * vchild(:)) / (colliders%mass(j) + mchild)
                  xc(:) = colliders%xb(:, j) - xcom(:)
                  vc(:) = colliders%vb(:, j) - vcom(:)
                  xcrossv(:) = xc(:) .cross. vc(:) 
                  colliders%L_spin(:, j) = colliders%L_spin(:, j) + colliders%mass(j) * xcrossv(:)
   
                  xc(:) = xchild(:) - xcom(:)
                  vc(:) = vchild(:) - vcom(:)
                  xcrossv(:) = xc(:) .cross. vc(:) 
                  colliders%L_spin(:, j) = colliders%L_spin(:, j) + mchild * xcrossv(:)

                  colliders%L_spin(:, j) = colliders%L_spin(:, j) + mchild * pl%Ip(3, idx_child)  &
                                                                           * pl%radius(idx_child)**2 &
                                                                           * pl%rot(:, idx_child)
                  colliders%Ip(:, j) = colliders%Ip(:, j) + mchild * pl%Ip(:, idx_child)
               end if

               ! Merge the child and parent
               colliders%mass(j) = colliders%mass(j) + mchild
               colliders%xb(:, j) = xcom(:)
               colliders%vb(:, j) = vcom(:)
            end do
         end if
         density(j) = colliders%mass(j) / volume(j)
         colliders%radius(j) = (3 * volume(j) / (4 * PI))**(1.0_DP / 3.0_DP)
         if (param%lrotation) colliders%Ip(:, j) = colliders%Ip(:, j) / colliders%mass(j)
      end do
      lflag = .true.

      xcom(:) = (colliders%mass(1) * colliders%xb(:, 1) + colliders%mass(2) * colliders%xb(:, 2)) / sum(colliders%mass(:))
      vcom(:) = (colliders%mass(1) * colliders%vb(:, 1) + colliders%mass(2) * colliders%vb(:, 2)) / sum(colliders%mass(:))
      mxc(:, 1) = colliders%mass(1) * (colliders%xb(:, 1) - xcom(:))
      mxc(:, 2) = colliders%mass(2) * (colliders%xb(:, 2) - xcom(:))
      vcc(:, 1) = colliders%vb(:, 1) - vcom(:)
      vcc(:, 2) = colliders%vb(:, 2) - vcom(:)
      colliders%L_orbit(:,:) = mxc(:,:) .cross. vcc(:,:)

      ! Destroy the kinship relationships for all members of this colliders%idx
      call pl%reset_kinship(colliders%idx(:))

      return
   end function symba_collision_consolidate_colliders


   module subroutine symba_collision_encounter_extract_collisions(self, system, param)
      !! author: David A. Minton
      !! 
      !! Processes the pl-pl encounter list remove only those encounters that led to a collision
      !!
      implicit none
      ! Arguments
      class(symba_plplenc),       intent(inout) :: self   !! SyMBA pl-pl encounter list
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      logical,      dimension(:), allocatable :: lplpl_collision
      logical,      dimension(:), allocatable :: lplpl_unique_parent
      integer(I4B), dimension(:), pointer     :: plparent
      integer(I4B), dimension(:), allocatable :: collision_idx, unique_parent_idx
      integer(I4B)                            :: i, index_coll, ncollisions, nunique_parent, nplplenc

      select type (pl => system%pl)
      class is (symba_pl)
         associate(plplenc_list => self, idx1 => self%index1, idx2 => self%index2, plparent => pl%kin%parent)
            nplplenc = plplenc_list%nenc
            allocate(lplpl_collision(nplplenc))
            lplpl_collision(:) = plplenc_list%status(1:nplplenc) == COLLISION
            if (.not.any(lplpl_collision)) return 
            ! Collisions have been detected in this step. So we need to determine which of them are between unique bodies.

            ! Get the subset of pl-pl encounters that lead to a collision
            ncollisions = count(lplpl_collision(:))
            allocate(collision_idx(ncollisions))
            collision_idx = pack([(i, i=1, nplplenc)], lplpl_collision)

            ! Get the subset of collisions that involve a unique pair of parents
            allocate(lplpl_unique_parent(ncollisions))

            lplpl_unique_parent(:) = plparent(idx1(collision_idx(:))) /= plparent(idx2(collision_idx(:)))
            nunique_parent = count(lplpl_unique_parent(:))
            allocate(unique_parent_idx(nunique_parent))
            unique_parent_idx = pack(collision_idx(:), lplpl_unique_parent(:))

            ! Scrub all pl-pl collisions involving unique pairs of parents, which will remove all duplicates and leave behind
            ! all pairs that have themselves as parents but are not part of the unique parent list. This can hapepn in rare cases
            ! due to restructuring of parent/child relationships when there are large numbers of multi-body collisions in a single
            ! step
            lplpl_unique_parent(:) = .true.
            do index_coll = 1, ncollisions
               associate(ip1 => plparent(idx1(collision_idx(index_coll))), ip2 => plparent(idx2(collision_idx(index_coll))))
                  lplpl_unique_parent(:) = .not. ( any(plparent(idx1(unique_parent_idx(:))) == ip1) .or. &
                                                   any(plparent(idx2(unique_parent_idx(:))) == ip1) .or. &
                                                   any(plparent(idx1(unique_parent_idx(:))) == ip2) .or. &
                                                   any(plparent(idx2(unique_parent_idx(:))) == ip2) )
               end associate
            end do

            ! Reassemble collision index list to include only those containing the unique pairs of parents, plus all the non-unique pairs that don't
            ! contain a parent body on the unique parent list.
            ncollisions = nunique_parent + count(lplpl_unique_parent)
            collision_idx = [unique_parent_idx(:), pack(collision_idx(:), lplpl_unique_parent(:))]

            ! Create a mask that contains only the pl-pl encounters that did not result in a collision, and then discard them
            lplpl_collision(:) = .false.
            lplpl_collision(collision_idx(:)) = .true.
            call plplenc_list%spill(system%plplcollision_list, lplpl_collision, ldestructive=.true.) ! Extract any encounters that are not collisions from the list.
         end associate
      end select

      return
   end subroutine symba_collision_encounter_extract_collisions


   module subroutine symba_collision_make_colliders_pl(self, idx)
      !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
      !!
      !! When a single body is involved in more than one collision in a single step, it becomes part of a colliders%idx.
      !! The largest body involved in a multi-body collision is the "parent" and all bodies that collide with it are its "children,"
      !! including those that collide with the children.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self !! SyMBA massive body object
      integer(I4B), dimension(2), intent(in)    :: idx  !! Array holding the indices of the two bodies involved in the collision
      ! Internals
      integer(I4B)                              :: i, j, index_parent, index_child, p1, p2
      integer(I4B)                              :: nchild_inherit, nchild_orig, nchild_new
      integer(I4B), dimension(:), allocatable   :: temp

      associate(pl => self)
         p1 = pl%kin(idx(1))%parent
         p2 = pl%kin(idx(2))%parent
         if (p1 == p2) return ! This is a collision between to children of a shared parent. We will ignore it.

         if (pl%mass(p1) > pl%mass(p2)) then
            index_parent = p1
            index_child = p2
         else
            index_parent = p2
            index_child = p1
         end if

         ! Expand the child array (or create it if necessary) and copy over the previous lists of children
         nchild_orig = pl%kin(index_parent)%nchild
         nchild_inherit = pl%kin(index_child)%nchild
         nchild_new = nchild_orig + nchild_inherit + 1
         allocate(temp(nchild_new))

         if (nchild_orig > 0) temp(1:nchild_orig) = pl%kin(index_parent)%child(1:nchild_orig)
         ! Find out if the child body has any children of its own. The new parent wil inherit these children
         if (nchild_inherit > 0) then
            temp(nchild_orig+1:nchild_orig+nchild_inherit) = pl%kin(index_child)%child(1:nchild_inherit)
            do i = 1, nchild_inherit
               j = pl%kin(index_child)%child(i)
               ! Set the childrens' parent to the new parent
               pl%kin(j)%parent = index_parent
            end do
         end if
         call pl%reset_kinship([index_child])
         ! Add the new child to its parent
         pl%kin(index_child)%parent = index_parent
         temp(nchild_new) = index_child
         ! Save the new child array to the parent
         pl%kin(index_parent)%nchild = nchild_new
         call move_alloc(from=temp, to=pl%kin(index_parent)%child)
      end associate

      return
   end subroutine symba_collision_make_colliders_pl


   subroutine symba_collision_mergeaddsub(system, param, colliders, frag, status)
      !! author:  David A. Minton
      !!
      !! Fills the pl_discards and pl_adds with removed and added bodies
      !!  
      implicit none
      ! Arguments
      class(symba_nbody_system), intent(inout) :: system    !! SyMBA nbody system object
      class(symba_parameters),   intent(inout) :: param     !! Current run configuration parameters with SyMBA additions
      class(fraggle_colliders),  intent(inout) :: colliders !! Fraggle colliders object        
      class(fraggle_fragments),  intent(inout) :: frag      !! Fraggle fragmentation system object 
      integer(I4B),              intent(in)    :: status    !! Status flag to assign to adds
      ! Internals
      integer(I4B) :: i, ibiggest, ismallest, iother, nstart, nend, ncolliders, nfrag
      logical, dimension(system%pl%nbody)    :: lmask
      class(symba_pl), allocatable           :: plnew, plsub
      character(*), parameter :: FRAGFMT = '("Newbody",I0.7)'
      character(len=NAMELEN) :: newname
   
      select type(pl => system%pl)
      class is (symba_pl)
         select type(pl_discards => system%pl_discards)
         class is (symba_merger)
            associate(info => pl%info, pl_adds => system%pl_adds, cb => system%cb, npl => pl%nbody)
               ! Add the colliders%idx bodies to the subtraction list
               ncolliders = colliders%ncoll
               nfrag = frag%nbody

               param%maxid_collision = max(param%maxid_collision, maxval(system%pl%info(:)%collision_id))
               param%maxid_collision = param%maxid_collision + 1

               ! Setup new bodies
               allocate(plnew, mold=pl)
               call plnew%setup(nfrag, param)
               ibiggest  = colliders%idx(maxloc(pl%Gmass(colliders%idx(:)), dim=1))
               ismallest = colliders%idx(minloc(pl%Gmass(colliders%idx(:)), dim=1))

               ! Copy over identification, information, and physical properties of the new bodies from the fragment list
               plnew%id(1:nfrag) = frag%id(1:nfrag) 
               plnew%xb(:, 1:nfrag) = frag%xb(:, 1:nfrag) 
               plnew%vb(:, 1:nfrag) = frag%vb(:, 1:nfrag)
               call pl%vb2vh(cb)
               call pl%xh2xb(cb)
               do i = 1, nfrag
                  plnew%rh(:,i) = frag%xb(:, i) - cb%xb(:)
                  plnew%vh(:,i) = frag%vb(:, i) - cb%vb(:)
               end do
               plnew%mass(1:nfrag) = frag%mass(1:nfrag)
               plnew%Gmass(1:nfrag) = param%GU * frag%mass(1:nfrag)
               plnew%radius(1:nfrag) = frag%radius(1:nfrag)
               plnew%density(1:nfrag) = frag%mass(1:nfrag) / frag%radius(1:nfrag)
               call plnew%set_rhill(cb)

               select case(status)
               case(DISRUPTION)
                  plnew%status(1:nfrag) = NEW_PARTICLE
                  do i = 1, nfrag
                     write(newname, FRAGFMT) frag%id(i)
                     call plnew%info(i)%set_value(origin_type="Disruption", origin_time=system%t, name=newname, &
                                                  origin_rh=plnew%rh(:,i), &
                                                  origin_vh=plnew%vh(:,i), collision_id=param%maxid_collision)
                  end do
                  do i = 1, ncolliders
                     if (colliders%idx(i) == ibiggest) then
                        iother = ismallest
                     else
                        iother = ibiggest
                     end if
                     call pl%info(colliders%idx(i))%set_value(status="Disruption", discard_time=system%t, &
                                                              discard_rh=pl%rh(:,i), discard_vh=pl%vh(:,i), discard_body_id=iother)
                  end do
               case(SUPERCATASTROPHIC)
                  plnew%status(1:nfrag) = NEW_PARTICLE
                  do i = 1, nfrag
                     write(newname, FRAGFMT) frag%id(i)
                     call plnew%info(i)%set_value(origin_type="Supercatastrophic", origin_time=system%t, name=newname, &
                                                  origin_rh=plnew%rh(:,i), origin_vh=plnew%vh(:,i), &
                                                  collision_id=param%maxid_collision)
                  end do
                  do i = 1, ncolliders
                     if (colliders%idx(i) == ibiggest) then
                        iother = ismallest
                     else
                        iother = ibiggest
                     end if
                     call pl%info(colliders%idx(i))%set_value(status="Supercatastrophic", discard_time=system%t, &
                                                              discard_rh=pl%rh(:,i), discard_vh=pl%vh(:,i), &
                                                              discard_body_id=iother)
                  end do
               case(HIT_AND_RUN_DISRUPT)
                  call plnew%info(1)%copy(pl%info(ibiggest))
                  plnew%status(1) = OLD_PARTICLE
                  do i = 2, nfrag
                     write(newname, FRAGFMT) frag%id(i)
                     call plnew%info(i)%set_value(origin_type="Hit and run fragment", origin_time=system%t, name=newname, &
                                                  origin_rh=plnew%rh(:,i), origin_vh=plnew%vh(:,i), &
                                                  collision_id=param%maxid_collision)
                  end do
                  do i = 1, ncolliders
                     if (colliders%idx(i) == ibiggest) cycle
                     iother = ibiggest
                     call pl%info(colliders%idx(i))%set_value(status="Hit and run fragmention", discard_time=system%t, &
                                                              discard_rh=pl%rh(:,i), discard_vh=pl%vh(:,i), &
                                                              discard_body_id=iother)
                  end do 
               case(MERGED)
                  call plnew%info(1)%copy(pl%info(ibiggest))
                  plnew%status(1) = OLD_PARTICLE
                  do i = 1, ncolliders
                     if (colliders%idx(i) == ibiggest) cycle

                     iother = ibiggest
                     call pl%info(colliders%idx(i))%set_value(status="MERGED", discard_time=system%t, discard_rh=pl%rh(:,i), &
                                                              discard_vh=pl%vh(:,i), discard_body_id=iother)
                  end do 
               end select
   
               if (param%lrotation) then
                  plnew%Ip(:, 1:nfrag) = frag%Ip(:, 1:nfrag)
                  plnew%rot(:, 1:nfrag) = frag%rot(:, 1:nfrag)
               end if
   
               ! if (param%ltides) then
               !    plnew%Q = pl%Q(ibiggest)
               !    plnew%k2 = pl%k2(ibiggest)
               !    plnew%tlag = pl%tlag(ibiggest)
               ! end if

               !Copy over or set integration parameters for new bodies
               plnew%lcollision(1:nfrag) = .false.
               plnew%ldiscard(1:nfrag) = .false.
               plnew%levelg(1:nfrag) = pl%levelg(ibiggest)
               plnew%levelm(1:nfrag) = pl%levelm(ibiggest)

               ! Log the properties of the new bodies
               call fraggle_io_log_pl(plnew, param)
   
               ! Append the new merged body to the list 
               nstart = pl_adds%nbody + 1
               nend = pl_adds%nbody + nfrag
               call pl_adds%append(plnew, lsource_mask=[(.true., i=1, nfrag)])
               ! Record how many bodies were added in this event
               pl_adds%ncomp(nstart:nend) = plnew%nbody

               ! Add the discarded bodies to the discard list
               pl%status(colliders%idx(:)) = MERGED
               pl%ldiscard(colliders%idx(:)) = .true.
               pl%lcollision(colliders%idx(:)) = .true.
               lmask(:) = .false.
               lmask(colliders%idx(:)) = .true.
               
               call plnew%setup(0, param)
               deallocate(plnew)

               allocate(plsub, mold=pl)
               call pl%spill(plsub, lmask, ldestructive=.false.)
   
               nstart = pl_discards%nbody + 1
               nend = pl_discards%nbody + ncolliders
               call pl_discards%append(plsub, lsource_mask=[(.true., i = 1, ncolliders)])

               ! Record how many bodies were subtracted in this event
               pl_discards%ncomp(nstart:nend) = ncolliders 
   
               call plsub%setup(0, param)
               deallocate(plsub)
            end associate
         end select
      end select
   
      return
   end subroutine symba_collision_mergeaddsub


   module subroutine symba_collision_resolve_fragmentations(self, system, param)
      !! author: David A. Minton
      !! 
      !! Process list of collisions, determine the collisional regime, and then create fragments.
      !!
      implicit none
      ! Arguments
      class(symba_plplenc),      intent(inout) :: self   !! SyMBA pl-pl encounter list
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      class(symba_parameters),   intent(inout) :: param  !! Current run configuration parameters with SyMBA additions
      ! Internals
      ! Internals
      integer(I4B), dimension(2)                  :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      logical                                     :: lgoodcollision
      integer(I4B)                                :: i
      type(fraggle_colliders)                     :: colliders !! Fraggle colliders object
      type(fraggle_fragments)                     :: frag      !! Fraggle fragmentation system object

      associate(plplcollision_list => self, ncollisions => self%nenc, idx1 => self%index1, idx2 => self%index2)
         select type(pl => system%pl)
         class is (symba_pl)
            select type (cb => system%cb)
            class is (symba_cb)
               do i = 1, ncollisions
                  idx_parent(1) = pl%kin(idx1(i))%parent
                  idx_parent(2) = pl%kin(idx2(i))%parent
                  lgoodcollision = symba_collision_consolidate_colliders(pl, cb, param, idx_parent, colliders) 
                  if ((.not. lgoodcollision) .or. any(pl%status(idx_parent(:)) /= COLLISION)) cycle

                  call colliders%regime(frag, system, param)
   
                  select case (frag%regime)
                  case (COLLRESOLVE_REGIME_DISRUPTION, COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
                     plplcollision_list%status(i) = symba_collision_casedisruption(system, param, colliders, frag)
                  case (COLLRESOLVE_REGIME_HIT_AND_RUN)
                     plplcollision_list%status(i) = symba_collision_casehitandrun(system, param, colliders, frag)
                  case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
                     plplcollision_list%status(i) = symba_collision_casemerge(system, param, colliders, frag)
                  case default 
                     write(*,*) "Error in symba_collision, unrecognized collision regime"
                     call util_exit(FAILURE)
                  end select
               end do
            end select
         end select
      end associate

      return
   end subroutine symba_collision_resolve_fragmentations


   module subroutine symba_collision_resolve_mergers(self, system, param)
      !! author: David A. Minton
      !! 
      !! Process list of collisions and merge colliding bodies together.
      !!
      implicit none
      ! Arguments
      class(symba_plplenc),      intent(inout) :: self   !! SyMBA pl-pl encounter list
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      class(symba_parameters),   intent(inout) :: param  !! Current run configuration parameters with SyMBA additions
      ! Internals
      integer(I4B), dimension(2)                  :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      logical                                     :: lgoodcollision
      integer(I4B)                                :: i
      type(fraggle_colliders)                     :: colliders !! Fraggle colliders object
      type(fraggle_fragments)                     :: frag      !! Fraggle fragmentation system object

      associate(plplcollision_list => self, ncollisions => self%nenc, idx1 => self%index1, idx2 => self%index2)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(cb => system%cb)
            class is (symba_cb)
               do i = 1, ncollisions
                  idx_parent(1) = pl%kin(idx1(i))%parent
                  idx_parent(2) = pl%kin(idx2(i))%parent
                  lgoodcollision = symba_collision_consolidate_colliders(pl, cb, param, idx_parent, colliders) 
                  if (.not. lgoodcollision) cycle
                  if (any(pl%status(idx_parent(:)) /= COLLISION)) cycle ! One of these two bodies has already been resolved
 
                  frag%regime = COLLRESOLVE_REGIME_MERGE
                  frag%mtot = sum(colliders%mass(:))
                  frag%mass_dist(1) = frag%mtot
                  frag%mass_dist(2) = 0.0_DP
                  frag%mass_dist(3) = 0.0_DP
                  frag%xbcom(:) = (colliders%mass(1) * colliders%xb(:,1) + colliders%mass(2) * colliders%xb(:,2)) / frag%mtot 
                  frag%vbcom(:) = (colliders%mass(1) * colliders%vb(:,1) + colliders%mass(2) * colliders%vb(:,2)) / frag%mtot
                  plplcollision_list%status(i) = symba_collision_casemerge(system, param, colliders, frag)
               end do
            end select
         end select
      end associate

      return
   end subroutine symba_collision_resolve_mergers


   module subroutine symba_collision_resolve_plplenc(self, system, param, t, dt, irec)
      !! author: David A. Minton
      !! 
      !! Process the pl-pl collision list, then modifiy the massive bodies based on the outcome of the collision
      !! 
      implicit none
      ! Arguments
      class(symba_plplenc),       intent(inout) :: self   !! SyMBA pl-pl encounter list
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters with SyMBA additions
      real(DP),                   intent(in)    :: t      !! Current simulation time
      real(DP),                   intent(in)    :: dt     !! Current simulation step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
      ! Internals
      real(DP) :: Eorbit_before, Eorbit_after
      logical :: lplpl_collision
      character(len=STRMAX) :: timestr
      class(symba_parameters), allocatable :: tmp_param
   
      associate(plplenc_list => self, plplcollision_list => system%plplcollision_list)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(param)
            class is (symba_parameters)
               if (plplcollision_list%nenc == 0) return ! No collisions to resolve
               ! Make sure that the heliocentric and barycentric coordinates are consistent with each other
               call pl%vb2vh(system%cb) 
               call pl%xh2xb(system%cb)
   
               ! Get the energy before the collision is resolved
               if (param%lenergy) then
                  call system%get_energy_and_momentum(param)
                  Eorbit_before = system%te
               end if

               do
                  write(timestr,*) t
                  call io_log_one_message(FRAGGLE_LOG_OUT, "")
                  call io_log_one_message(FRAGGLE_LOG_OUT, "***********************************************************" // &
                                                           "***********************************************************")
                  call io_log_one_message(FRAGGLE_LOG_OUT, "Collision between massive bodies detected at time t = " // &
                                                           trim(adjustl(timestr)))
                  call io_log_one_message(FRAGGLE_LOG_OUT, "***********************************************************" // &
                                                           "***********************************************************")
                  allocate(tmp_param, source=param)
                  if (param%lfragmentation) then
                     call plplcollision_list%resolve_fragmentations(system, param)
                  else
                     call plplcollision_list%resolve_mergers(system, param)
                  end if

                  ! Destroy the collision list now that the collisions are resolved
                  call plplcollision_list%setup(0_I8B)

                  if ((system%pl_adds%nbody == 0) .and. (system%pl_discards%nbody == 0)) exit

                  ! Save the add/discard information to file
                  call system%write_discard(tmp_param)

                  ! Rearrange the arrays: Remove discarded bodies, add any new bodies, resort, and recompute all indices and encounter lists
                  call pl%rearray(system, tmp_param)

                  ! Destroy the add/discard list so that we don't append the same body multiple times if another collision is detected
                  call system%pl_discards%setup(0, param)
                  call system%pl_adds%setup(0, param)
                  deallocate(tmp_param)

                  ! Check whether or not any of the particles that were just added are themselves in a collision state. This will generate a new plplcollision_list 
                  lplpl_collision = plplenc_list%collision_check(system, param, t, dt, irec)

                  if (.not.lplpl_collision) exit
               end do

               if (param%lenergy) then
                  call system%get_energy_and_momentum(param)
                  Eorbit_after = system%te
                  system%Ecollisions = system%Ecollisions + (Eorbit_after - Eorbit_before)
               end if

            end select 
         end select
      end associate

      return
   end subroutine symba_collision_resolve_plplenc


   module subroutine symba_collision_resolve_pltpenc(self, system, param, t, dt, irec)
      !! author: David A. Minton
      !! 
      !! Process the pl-tp collision list, then modifiy the massive bodies based on the outcome of the collision
      !! 
      implicit none
      ! Arguments
      class(symba_pltpenc),       intent(inout) :: self   !! SyMBA pl-pl encounter list
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters with SyMBA additions
      real(DP),                   intent(in)    :: t      !! Current simulation tim
      real(DP),                   intent(in)    :: dt     !! Current simulation step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
     
      ! Make sure coordinate systems are all synced up due to being inside the recursion at this point
      call system%pl%vb2vh(system%cb)
      call system%tp%vb2vh(system%cb%vb)
      call system%pl%b2h(system%cb)
      call system%tp%b2h(system%cb)

      ! Discard the collider
      call system%tp%discard(system, param)

      return
   end subroutine symba_collision_resolve_pltpenc

end submodule s_symba_collision