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


   function symba_collision_consolidate_familes(pl, param, idx_parent, family, x, v, mass, radius, L_spin, Ip) result(lflag)
      !! author: David A. Minton
      !! 
      !! Loops through the pl-pl collision list and groups families together by index. Outputs the indices of all family members, 
      !! and pairs of quantities (x and v vectors, mass, radius, L_spin, and Ip) that can be used to resolve the collisional outcome.
      implicit none
      ! Arguments
      class(symba_pl),                                 intent(inout) :: pl               !! SyMBA massive body object
      class(symba_parameters),                         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(2),                   intent(inout) :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      integer(I4B),    dimension(:),      allocatable, intent(out)   :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(NDIM,2),              intent(out)   :: x, v, L_spin, Ip !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(2),                   intent(out)   :: mass, radius     !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      ! Result
      logical                                                        :: lflag            !! Logical flag indicating whether a family was successfully created or not
      ! Internals
      type family_array
         integer(I4B), dimension(:), allocatable :: id
         integer(I4B), dimension(:), allocatable :: idx
      end type family_array
      type(family_array), dimension(2) :: parent_child_index_array
      integer(I4B), dimension(2)       :: nchild
      integer(I4B)                     :: i, j, fam_size, idx_child
      real(DP), dimension(2)           :: volume, density
      real(DP)                         :: mchild, mtot, volchild
      real(DP), dimension(NDIM)        :: xc, vc, xcom, vcom, xchild, vchild, xcrossv

      nchild(:) = pl%kin(idx_parent(:))%nchild 
      ! If all of these bodies share a parent, but this is still a unique collision, move the last child
      ! out of the parent's position and make it the secondary body
      if (idx_parent(1) == idx_parent(2)) then
         if (nchild(1) == 0) then ! There is only one valid body recorded in this pair (this could happen due to restructuring of the kinship relationships, though it should be rare)
            lflag = .false. 
            return
         end if
         idx_parent(2) = pl%kin(idx_parent(1))%child(nchild(1))
         nchild(1) = nchild(1) - 1
         nchild(2) = 0
         pl%kin(idx_parent(:))%nchild = nchild(:)
         pl%kin(idx_parent(2))%parent = idx_parent(1)
      end if

      mass(:) = pl%mass(idx_parent(:)) ! Note: This is meant to mass, not G*mass, as the collisional regime determination uses mass values that will be converted to Si
      radius(:) = pl%radius(idx_parent(:))
      volume(:) =  (4.0_DP / 3.0_DP) * PI * radius(:)**3
 
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

      ! Consolidate the groups of collsional parents with any children they may have into a single "family" index array
      fam_size = 2 + sum(nchild(:))
      allocate(family(fam_size))
      family = [parent_child_index_array(1)%idx(:),parent_child_index_array(2)%idx(:)]
      fam_size = count(pl%lcollision(family(:)))
      family = pack(family(:), pl%lcollision(family(:)))
      L_spin(:,:) = 0.0_DP
      Ip(:,:) = 0.0_DP

      ! Find the barycenter of each body along with its children, if it has any
      do j = 1, 2
         x(:, j)  = pl%xb(:, idx_parent(j))
         v(:, j)  = pl%vb(:, idx_parent(j))
         ! Assume principal axis rotation about axis corresponding to highest moment of inertia (3rd Ip)
         if (param%lrotation) then
            Ip(:, j) = mass(j) * pl%Ip(:, idx_parent(j))
            L_spin(:, j) = Ip(3, j) * radius(j)**2 * pl%rot(:, idx_parent(j))
         end if

         if (nchild(j) > 0) then
            do i = 1, nchild(j) ! Loop over all children and take the mass weighted mean of the properties
               idx_child = parent_child_index_array(j)%idx(i + 1)
               if (.not. pl%lcollision(idx_child)) cycle
               mchild = pl%mass(idx_child)
               xchild(:) = pl%xb(:, idx_child)
               vchild(:) = pl%vb(:, idx_child)
               volchild = (4.0_DP / 3.0_DP) * PI * pl%radius(idx_child)**3
               volume(j) = volume(j) + volchild
               ! Get angular momentum of the child-parent pair and add that to the spin
               ! Add the child's spin
               if (param%lrotation) then
                  xcom(:) = (mass(j) * x(:,j) + mchild * xchild(:)) / (mass(j) + mchild)
                  vcom(:) = (mass(j) * v(:,j) + mchild * vchild(:)) / (mass(j) + mchild)
                  xc(:) = x(:, j) - xcom(:)
                  vc(:) = v(:, j) - vcom(:)
                  xcrossv(:) = xc(:) .cross. vc(:) 
                  L_spin(:, j) = L_spin(:, j) + mass(j) * xcrossv(:)
   
                  xc(:) = xchild(:) - xcom(:)
                  vc(:) = vchild(:) - vcom(:)
                  xcrossv(:) = xc(:) .cross. vc(:) 
                  L_spin(:, j) = L_spin(:, j) + mchild * xcrossv(:)

                  L_spin(:, j) = L_spin(:, j) + mchild * pl%Ip(3, idx_child) * pl%radius(idx_child)**2 * pl%rot(:, idx_child)
                  Ip(:, j) = Ip(:, j) + mchild * pl%Ip(:, idx_child)
               end if

               ! Merge the child and parent
               mass(j) = mass(j) + mchild
               x(:, j) = xcom(:)
               v(:, j) = vcom(:)
            end do
         end if
         density(j) =  mass(j) / volume(j)
         radius(j) = ((3 * mass(j)) / (density(j) * 4 * pi))**(1.0_DP / 3.0_DP)
         if (param%lrotation) Ip(:, j) = Ip(:, j) / mass(j)
      end do
      lflag = .true.

      return
   end function symba_collision_consolidate_familes


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
            if (any(lplpl_collision)) then ! Collisions have been detected in this step. So we need to determine which of them are between unique bodies.

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
            end if
            call plplenc_list%spill(plplenc_noncollision, .not.lplpl_collision, ldestructive=.true.) ! Remove any encounters that are not collisions from the list.
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


   module subroutine symba_collision_resolve_fragmentations(self, system, param)
      !! author: David A. Minton
      !! 
      !! Process list of collisions, determine the collisional regime, and then create fragments.
      !!
      implicit none
      ! Arguments
      class(symba_plplenc),      intent(inout) :: self   !! SyMBA pl-pl encounter list
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      class(symba_parameters),   intent(in)    :: param  !! Current run configuration parameters with SyMBA additions
      ! Internals
      ! Internals
      integer(I4B), dimension(:),     allocatable :: family           !! List of indices of all bodies inovlved in the collision
      integer(I4B), dimension(2)                  :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      real(DP),     dimension(NDIM,2)             :: x, v, L_spin, Ip !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),     dimension(2)                  :: mass, radius     !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      logical                                     :: lgoodcollision
      integer(I4B)                                :: i, status, jtarg, jproj, regime
      real(DP), dimension(2)                      :: radius_si, mass_si, density_si
      real(DP)                                    :: mtiny_si, Mcb_si
      real(DP), dimension(NDIM)                   :: x1_si, v1_si, x2_si, v2_si
      real(DP)                                    :: mlr, mslr, mtot, dentot, msys, msys_new, Qloss, impact_parameter
      integer(I4B), parameter                     :: NRES = 3   !! Number of collisional product results
      real(DP), dimension(NRES)                   :: mass_res   

      associate(plpl_collisions => self, ncollisions => self%nenc, idx1 => self%index1, idx2 => self%index2, cb => system%cb)
         select type(pl => system%pl)
         class is (symba_pl)
            do i = 1, ncollisions
               idx_parent(1) = pl%kin(idx1(i))%parent
               idx_parent(2) = pl%kin(idx2(i))%parent
               lgoodcollision = symba_collision_consolidate_familes(pl, param, idx_parent, family, x, v, mass, radius, L_spin, Ip)
               if (.not. lgoodcollision) cycle
               if (any(pl%status(idx_parent(:)) /= COLLISION)) cycle ! One of these two bodies has already been resolved

               ! Convert all quantities to SI units and determine which of the pair is the projectile vs. target before sending them 
               ! to symba_regime
               if (mass(1) > mass(2)) then
                  jtarg = 1
                  jproj = 2
               else
                  jtarg = 2
                  jproj = 1
               end if
               mass_si(:)    = (mass(:)) * param%MU2KG                              !! The collective mass of the parent and its children
               radius_si(:)  = radius(:) * param%DU2M                               !! The collective radius of the parent and its children
               x1_si(:)      = plpl_collisions%x1(:,i) * param%DU2M                 !! The position of the parent from inside the step (at collision)
               v1_si(:)      = plpl_collisions%v1(:,i) * param%DU2M / param%TU2S    !! The velocity of the parent from inside the step (at collision)
               x2_si(:)      = plpl_collisions%x2(:,i) * param%DU2M                 !! The position of the parent from inside the step (at collision)
               v2_si(:)      = plpl_collisions%v2(:,i) * param%DU2M / param%TU2S    !! The velocity of the parent from inside the step (at collision)
               density_si(:) = mass_si(:) / (4.0_DP / 3._DP * PI * radius_si(:)**3) !! The collective density of the parent and its children
               Mcb_si        = cb%mass * param%MU2KG 
               mtiny_si      = (param%MTINY / param%GU) * param%MU2KG
            
               mass_res(:) = 0.0_DP
         
               mtot = sum(mass_si(:)) 
               dentot = sum(mass_si(:) * density_si(:)) / mtot 

               !! Use the positions and velocities of the parents from indside the step (at collision) to calculate the collisional regime
               call fragmentation_regime(Mcb_si, mass_si(jtarg), mass_si(jproj), radius_si(jtarg), radius_si(jproj), x1_si(:), x2_si(:),& 
                     v1_si(:), v2_si(:), density_si(jtarg), density_si(jproj), regime, mlr, mslr, mtiny_si, Qloss)

               mass_res(1) = min(max(mlr, 0.0_DP), mtot)
               mass_res(2) = min(max(mslr, 0.0_DP), mtot)
               mass_res(3) = min(max(mtot - mlr - mslr, 0.0_DP), mtot)
               mass_res(:) = (mass_res(:) / param%MU2KG) * param%GU
               Qloss = Qloss * (param%GU / param%MU2KG) * (param%TU2S / param%DU2M)**2

               select case (regime)
               case (COLLRESOLVE_REGIME_DISRUPTION)
                  status = symba_fragmentation_casedisruption(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)
               case (COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
                  status = symba_fragmentation_casesupercatastrophic(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)
               case (COLLRESOLVE_REGIME_HIT_AND_RUN)
                  status = symba_fragmentation_casehitandrun(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)
               case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
                  status = symba_fragmentation_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip) 
               case default 
                  write(*,*) "Error in symba_collision, unrecognized collision regime"
                  call util_exit(FAILURE)
               end select
            end do
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
      class(symba_parameters),   intent(in)    :: param  !! Current run configuration parameters with SyMBA additions
      ! Internals
      integer(I4B), dimension(:),     allocatable :: family           !! List of indices of all bodies inovlved in the collision
      integer(I4B), dimension(2)                  :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      real(DP),     dimension(NDIM,2)             :: x, v, L_spin, Ip !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),     dimension(2)                  :: mass, radius     !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      logical                                     :: lgoodcollision
      integer(I4B)                                :: i, status

      associate(plpl_collisions => self, ncollisions => self%nenc, idx1 => self%index1, idx2 => self%index2)
         select type(pl => system%pl)
         class is (symba_pl)
            do i = 1, ncollisions
               idx_parent(1) = pl%kin(idx1(i))%parent
               idx_parent(2) = pl%kin(idx2(i))%parent
               lgoodcollision = symba_collision_consolidate_familes(pl, param, idx_parent, family, x, v, mass, radius, L_spin, Ip)
               if (.not. lgoodcollision) cycle
               if (any(pl%status(idx_parent(:)) /= COLLISION)) cycle ! One of these two bodies has already been resolved
               status = symba_fragmentation_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip) 
            end do
         end select
      end associate

      return
   end subroutine symba_collision_resolve_mergers

end submodule s_symba_collision