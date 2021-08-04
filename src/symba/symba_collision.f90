submodule (symba_classes) s_symba_collision
   use swiftest
contains

   module subroutine symba_collision_check_pltpenc(self, system, param, t, dt, irec)
      !! author: David A. Minton
      !!
      !! Check for merger between massive bodies and test particles in SyMBA
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge.f90 and symba_merge_tp.f90
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
      logical, dimension(:), allocatable        :: lcollision, lmask
      real(DP), dimension(NDIM)                 :: xr, vr
      integer(I4B)                              :: k
      real(DP)                                  :: rlim, mtot
      logical                                   :: isplpl

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
            associate(nenc => self%nenc, ind1 => self%index1, ind2 => self%index2)
               allocate(lmask(nenc))
               lmask(:) = ((self%status(1:nenc) == ACTIVE) .and. (pl%levelg(ind1(1:nenc)) >= irec))
               if (isplpl) then
                  lmask(:) = lmask(:) .and. (pl%levelg(ind2(1:nenc)) >= irec)
               else
                  lmask(:) = lmask(:) .and. (tp%levelg(ind2(1:nenc)) >= irec)
               end if
               if (.not.any(lmask(:))) return

               allocate(lcollision(nenc))
               lcollision(:) = .false.

               if (isplpl) then
                  do concurrent(k = 1:nenc, lmask(k))
                     xr(:) = pl%xh(:, ind1(k)) - pl%xh(:, ind2(k)) 
                     vr(:) = pl%vb(:, ind1(k)) - pl%vb(:, ind2(k))
                     rlim = pl%radius(ind1(k)) + pl%radius(ind2(k))
                     mtot = pl%Gmass(ind1(k)) + pl%Gmass(ind2(k))
                     lcollision(k) = symba_collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), mtot, rlim, dt, self%lvdotr(k))
                  end do
               else
                  do concurrent(k = 1:nenc, lmask(k))
                     xr(:) = pl%xh(:, ind1(k)) - tp%xh(:, ind2(k)) 
                     vr(:) = pl%vb(:, ind1(k)) - tp%vb(:, ind2(k))
                     lcollision(k) = symba_collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%Gmass(ind1(k)), pl%radius(ind1(k)), dt, self%lvdotr(k))
                  end do
               end if

               if (any(lcollision(:))) then
                  do k = 1, nenc
                     if (.not.lcollision(k)) cycle 
                     self%status(k) = COLLISION
                     self%x1(:,k) = pl%xh(:,ind1(k)) 
                     self%v1(:,k) = pl%vb(:,ind1(k)) 
                     if (isplpl) then
                        self%x2(:,k) = pl%xh(:,ind2(k))
                        self%v2(:,k) = pl%vb(:,ind2(k))

                        ! Check to see if either of these bodies has been involved with a collision before, and if so, make this a collisional family
                        if (pl%lcollision(ind1(k)) .or. pl%lcollision(ind2(k))) call pl%make_family([ind1(k),ind2(k)])

                        ! Set the collision flag for these to bodies to true in case they become involved in another collision later in the step
                        pl%lcollision([ind1(k), ind2(k)]) = .true.
                        pl%ldiscard([ind1(k), ind2(k)]) = .true.
                        pl%status([ind1(k), ind2(k)]) = COLLISION
                     else
                        self%x2(:,k) = tp%xh(:,ind2(k))
                        self%v2(:,k) = tp%vb(:,ind2(k))
                        tp%status(ind2(k)) = DISCARDED_PLR
                        tp%ldiscard(ind2(k)) = .true.
                        write(*,*) 'Test particle ',tp%id(ind2(k)), ' collided with massive body ',pl%id(ind1(k)), ' at time ',t
                     end if
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


   module subroutine symba_collision_encounter_scrub(self, system, param)
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
      logical,      dimension(self%nenc)      :: lplpl_collision
      logical,      dimension(:), allocatable :: lplpl_unique_parent
      integer(I4B), dimension(:), pointer     :: plparent
      integer(I4B), dimension(:), allocatable :: collision_idx, unique_parent_idx
      integer(I4B)                            :: i, index_coll, ncollisions, nunique_parent
      type(symba_plplenc)                     :: plplenc_noncollision


      select type (pl => system%pl)
      class is (symba_pl)
         associate(plplenc_list => self, nplplenc => self%nenc, idx1 => self%index1, idx2 => self%index2, plparent => pl%kin%parent)
            lplpl_collision(:) = plplenc_list%status(1:nplplenc) == COLLISION
            if (.not.any(lplpl_collision)) return

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
            lplpl_collision(:) = .true.
            lplpl_collision(collision_idx(:)) = .false.
            call plplenc_list%spill(plplenc_noncollision, lplpl_collision, ldestructive = .true.)
         end associate
      end select

      return
   end subroutine symba_collision_encounter_scrub


   module subroutine symba_collision_make_family_pl(self, idx)
      !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
      !!
      !! When a single body is involved in more than one collision in a single step, it becomes part of a family.
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

         if (pl%Gmass(p1) > pl%Gmass(p2)) then
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
         if (allocated(pl%kin(index_child)%child)) deallocate(pl%kin(index_child)%child)
         pl%kin(index_child)%nchild = 0
         ! Add the new child to its parent
         pl%kin(index_child)%parent = index_parent
         temp(nchild_new) = index_child
         ! Save the new child array to the parent
         pl%kin(index_parent)%nchild = nchild_new
         call move_alloc(from=temp, to=pl%kin(index_parent)%child)
      end associate

      return
   end subroutine symba_collision_make_family_pl

end submodule s_symba_collision