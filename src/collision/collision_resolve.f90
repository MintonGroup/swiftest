! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (collision) s_collision_resolve
   use swiftest
contains

   module subroutine collision_resolve_consolidate_impactors(self, nbody_system, param, idx_parent, lflag)
      !! author: David A. Minton
      !! 
      !! Loops through the pl-pl collision list and groups families together by index. Outputs the indices of all impactors%id 
      !! members, and pairs of quantities (x and v vectors, mass, radius, L_spin, and Ip) that can be used to resolve the 
      !! collisional outcome.
      implicit none
      ! Arguments
      class(collision_impactors),intent(out) :: self !! Collision impactors object
      class(base_nbody_system),intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters), intent(in) :: param !! Current run configuration parameters with Swiftest additions
      integer(I4B), dimension(:), intent(inout) :: idx_parent !! Index of the two bodies considered the "parents" of the collision
      logical, intent(out) :: lflag !! Logical flag indicating whether a impactors%id was successfully created or not
      ! Internals
      type collidx_array
         integer(I4B), dimension(:), allocatable :: id
         integer(I4B), dimension(:), allocatable :: idx
      end type collidx_array
      type(collidx_array), dimension(2) :: parent_child_index_array
      integer(I4B), dimension(2)       :: nchild
      integer(I4B)                     :: i, j, nimpactors, idx_child
      real(DP), dimension(2)           :: volume, density
      real(DP)                         :: mchild, volchild, rrel_mag, rlim, mtot, vdotr
      real(DP), dimension(NDIM)        :: xc, vc, xcom, vcom, xchild, vchild, xcrossv, rrel, vrel, rrel_unit, vrel_unit, dr
      real(DP), dimension(NDIM,2)      :: mxc, vcc

      select type(nbody_system)
      class is (swiftest_nbody_system)
         associate(impactors => self, pl => nbody_system%pl, cb => nbody_system%cb)

            nchild(:) = pl%kin(idx_parent(:))%nchild 
            ! If all of these bodies share a parent, but this is still a unique collision, move the last child
            ! out of the parent's position and make it the secondary body
            if (idx_parent(1) == idx_parent(2)) then
               if (nchild(1) == 0) then ! There is only one valid body recorded in this pair (this could happen due to restructuring
                                        ! of the kinship relationships, though it should be rare)
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

            impactors%Gmass(:) = pl%Gmass(idx_parent(:)) 
            impactors%mass(:)  = pl%mass(idx_parent(:)) 
            impactors%radius(:) = pl%radius(idx_parent(:))
            volume(:) =  (4.0_DP / 3.0_DP) * PI * impactors%radius(:)**3
      
            ! Group together the ids and indexes of each collisional parent and its children
            do j = 1, 2
               allocate(parent_child_index_array(j)%idx(nchild(j)+ 1))
               allocate(parent_child_index_array(j)%id(nchild(j)+ 1))
               associate(idx_arr => parent_child_index_array(j)%idx, &
                        id_arr => parent_child_index_array(j)%id, &
                        ncj => nchild(j), &
                        plkinj => pl%kin(idx_parent(j)))
                  idx_arr(1) = idx_parent(j)
                  if (ncj > 0) idx_arr(2:ncj + 1) = plkinj%child(1:ncj)
                  id_arr(:) = pl%id(idx_arr(:))
               end associate
            end do

            ! Consolidate the groups of collsional parents with any children they may have into a single "impactors%id" index array
            nimpactors = 2 + sum(nchild(:))
            allocate(impactors%id(nimpactors))
            impactors%id = [parent_child_index_array(1)%idx(:),parent_child_index_array(2)%idx(:)]

            impactors%ncoll = count(pl%lcollision(impactors%id(:)))
            impactors%id = pack(impactors%id(:), pl%lcollision(impactors%id(:)))
            impactors%L_spin(:,:) = 0.0_DP
            impactors%Ip(:,:) = 0.0_DP

            ! Find the barycenter of each body along with its children, if it has any
            do j = 1, 2
               impactors%rb(:, j)  = pl%rh(:, idx_parent(j)) + cb%rb(:)
               impactors%vb(:, j)  = pl%vb(:, idx_parent(j))
               ! Assume principal axis rotation about axis corresponding to highest moment of inertia (3rd Ip)
               if (param%lrotation) then
                  impactors%Ip(:, j) = impactors%mass(j) * pl%Ip(:, idx_parent(j))
                  impactors%L_spin(:, j) = impactors%Ip(3, j) * impactors%radius(j)**2 * pl%rot(:, idx_parent(j))
               end if

               if (nchild(j) > 0) then
                  do i = 1, nchild(j) ! Loop over all children and take the mass weighted mean of the properties
                     idx_child = parent_child_index_array(j)%idx(i + 1)
                     if (.not. pl%lcollision(idx_child)) cycle
                     mchild = pl%mass(idx_child)
                     xchild(:) = pl%rh(:, idx_child) + cb%rb(:)
                     vchild(:) = pl%vb(:, idx_child)
                     volchild = (4.0_DP / 3.0_DP) * PI * pl%radius(idx_child)**3
                     volume(j) = volume(j) + volchild
                     ! Get angular momentum of the child-parent pair and add that to the spin
                     ! Add the child's spin
                     if (param%lrotation) then
                        xcom(:) = (impactors%mass(j) * impactors%rb(:,j) + mchild * xchild(:)) / (impactors%mass(j) + mchild)
                        vcom(:) = (impactors%mass(j) * impactors%vb(:,j) + mchild * vchild(:)) / (impactors%mass(j) + mchild)
                        xc(:) = impactors%rb(:, j) - xcom(:)
                        vc(:) = impactors%vb(:, j) - vcom(:)
                        xcrossv(:) = xc(:) .cross. vc(:) 
                        impactors%L_spin(:, j) = impactors%L_spin(:, j) + impactors%mass(j) * xcrossv(:)
         
                        xc(:) = xchild(:) - xcom(:)
                        vc(:) = vchild(:) - vcom(:)
                        xcrossv(:) = xc(:) .cross. vc(:) 
                        impactors%L_spin(:, j) = impactors%L_spin(:, j) + mchild * xcrossv(:)

                        impactors%L_spin(:, j) = impactors%L_spin(:, j) + mchild * pl%Ip(3, idx_child)  &
                                                                                 * pl%radius(idx_child)**2 &
                                                                                 * pl%rot(:, idx_child)
                        impactors%Ip(:, j) = impactors%Ip(:, j) + mchild * pl%Ip(:, idx_child)
                     end if

                     ! Merge the child and parent
                     impactors%mass(j) = impactors%mass(j) + mchild
                     impactors%rb(:, j) = xcom(:)
                     impactors%vb(:, j) = vcom(:)
                  end do
               end if
               density(j) = impactors%mass(j) / volume(j)
               impactors%radius(j) = (3 * volume(j) / (4 * PI))**(1.0_DP / 3.0_DP)
               if (param%lrotation) then
                  impactors%Ip(:, j) = impactors%Ip(:, j) / impactors%mass(j)
                  impactors%rot(:,j) = impactors%L_spin(:, j) / (impactors%Ip(3,j) * impactors%mass(j) * impactors%radius(j)**2)
               end if
            end do
            lflag = .true.

            mtot = sum(impactors%mass(:))
            xcom(:) = (impactors%mass(1) * impactors%rb(:, 1) + impactors%mass(2) * impactors%rb(:, 2)) / mtot
            vcom(:) = (impactors%mass(1) * impactors%vb(:, 1) + impactors%mass(2) * impactors%vb(:, 2)) / mtot

            ! Shift the impactors so that they are not overlapping
            rlim = sum(impactors%radius(1:2))
            rrel = impactors%rb(:,2) - impactors%rb(:,1)
            rrel_mag = .mag. rrel 
            if (rrel_mag < rlim) then
               rrel_unit = .unit.rrel
               vrel = impactors%vb(:,2) - impactors%vb(:,1)
               vrel_unit = .unit.vrel
               vdotr = dot_product(vrel_unit, rrel)
               dr(:) = -(vdotr - sign(1.0_DP, vdotr) * sqrt(rlim**2 - rrel_mag**2 + vdotr**2)) * vrel_unit(:)
               dr(:) = (1.0_DP + 2*epsilon(1.0_DP)) * dr(:)
               impactors%rb(:,1) = impactors%rb(:,1) - dr(:) *  impactors%mass(2) / mtot
               impactors%rb(:,2) = impactors%rb(:,2) + dr(:) *  impactors%mass(1) / mtot
               rrel = impactors%rb(:,2) - impactors%rb(:,1)
               rrel_mag = .mag. rrel 
            end if

            mxc(:, 1) = impactors%mass(1) * (impactors%rb(:, 1) - xcom(:))
            mxc(:, 2) = impactors%mass(2) * (impactors%rb(:, 2) - xcom(:))
            vcc(:, 1) = impactors%vb(:, 1) - vcom(:)
            vcc(:, 2) = impactors%vb(:, 2) - vcom(:)
            impactors%L_orbit(:,:) = mxc(:,:) .cross. vcc(:,:)

            ! Destroy the kinship relationships for all members of this impactors%id
            call pl%reset_kinship(impactors%id(:))

         end associate
      end select
      return
   end subroutine collision_resolve_consolidate_impactors


   module subroutine collision_resolve_extract_plpl(self, nbody_system, param)
      !! author: David A. Minton
      !! 
      !! Processes the pl-pl encounter list remove only those encounters that led to a collision
      !!
      implicit none
      ! Arguments
      class(collision_list_plpl), intent(inout) :: self         !! pl-pl encounter list
      class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),     intent(in)    :: param        !! Current run configuration parameters
      ! Internals
      logical,      dimension(:), allocatable :: lplpl_collision
      logical,      dimension(:), allocatable :: lplpl_unique_parent
      integer(I4B), dimension(:), pointer     :: plparent
      integer(I4B)                            :: nunique_parent
      integer(I8B)                            :: ncollisions, index_coll, k, nplplenc
      integer(I8B), dimension(:), allocatable :: unique_parent_idx, collision_idx

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type (pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(idx1 => self%index1, idx2 => self%index2, plparent => pl%kin%parent)
            nplplenc = self%nenc
            allocate(lplpl_collision(nplplenc))
            lplpl_collision(:) = self%status(1_I8B:nplplenc) == COLLIDED
            if (.not.any(lplpl_collision)) return 
            ! Collisions have been detected in this step. So we need to determine which of them are between unique bodies.

            ! Get the subset of pl-pl encounters that lead to a collision
            ncollisions = count(lplpl_collision(:))
            allocate(collision_idx(ncollisions))
            collision_idx = pack([(k, k=1_I8B, nplplenc)], lplpl_collision)

            ! Get the subset of collisions that involve a unique pair of parents
            allocate(lplpl_unique_parent(ncollisions))

            lplpl_unique_parent(:) = plparent(idx1(collision_idx(:))) /= plparent(idx2(collision_idx(:)))
            nunique_parent = count(lplpl_unique_parent(:))
            allocate(unique_parent_idx(nunique_parent))
            unique_parent_idx = pack(collision_idx(:), lplpl_unique_parent(:))

            ! Scrub all pl-pl collisions involving unique pairs of parents, which will remove all duplicates and leave behind
            ! all pairs that have themselves as parents but are not part of the unique parent list. This can happen in rare cases
            ! due to restructuring of parent/child relationships when there are large numbers of multi-body collisions in a single
            ! step
            lplpl_unique_parent(:) = .true.
            do index_coll = 1_I8B, ncollisions
               associate(ip1 => plparent(idx1(collision_idx(index_coll))), ip2 => plparent(idx2(collision_idx(index_coll))))
                  lplpl_unique_parent(:) = .not. ( any(plparent(idx1(unique_parent_idx(:))) == ip1) .or. &
                                                   any(plparent(idx2(unique_parent_idx(:))) == ip1) .or. &
                                                   any(plparent(idx1(unique_parent_idx(:))) == ip2) .or. &
                                                   any(plparent(idx2(unique_parent_idx(:))) == ip2) )
               end associate
            end do

            ! Reassemble collision index list to include only those containing the unique pairs of parents, plus all the non-unique
            ! pairs that don't contain a parent body on the unique parent list.
            ncollisions = nunique_parent + count(lplpl_unique_parent)
            collision_idx = [unique_parent_idx(:), pack(collision_idx(:), lplpl_unique_parent(:))]

            ! Create a mask that contains only the pl-pl encounters that did not result in a collision, and then discard them
            lplpl_collision(:) = .false.
            lplpl_collision(collision_idx(:)) = .true.
            ! Extract any encounters that are not collisions from the list.
            call self%spill(nbody_system%plpl_collision, lplpl_collision, ldestructive=.true.) 
         end associate
      end select
      end select

      return
   end subroutine collision_resolve_extract_plpl


   module subroutine collision_resolve_extract_pltp(self, nbody_system, param)
      implicit none
      class(collision_list_pltp), intent(inout) :: self   !! pl-tp encounter list
      class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),     intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      logical,      dimension(:), allocatable :: lpltp_collision
      integer(I8B)                            :: ncollisions, index_coll, k, npltpenc
      integer(I8B), dimension(:), allocatable :: collision_idx

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type (pl => nbody_system%pl)
      class is (swiftest_pl)
      select type (tp => nbody_system%tp)
      class is (swiftest_tp)
         associate(idx1 => self%index1, idx2 => self%index2)
            npltpenc = self%nenc
            allocate(lpltp_collision(npltpenc))
            lpltp_collision(:) = self%status(1_I8B:npltpenc) == COLLIDED
            if (.not.any(lpltp_collision)) return 
            ! Collisions have been detected in this step. So we need to determine which of them are between unique bodies.

            ! Get the subset of pl-tp encounters that lead to a collision
            ncollisions = count(lpltp_collision(:))
            allocate(collision_idx(ncollisions))
            collision_idx = pack([(k, k=1_I8B, npltpenc)], lpltp_collision)

            ! Create a mask that contains only the pl-tp encounters that did not result in a collision, and then discard them
            lpltp_collision(:) = .false.
            lpltp_collision(collision_idx(:)) = .true.
            ! Extract any encounters that are not collisions from the list.
            call self%spill(nbody_system%pltp_collision, lpltp_collision, ldestructive=.true.) 
         end associate
      end select
      end select
      end select


      return
   end subroutine collision_resolve_extract_pltp


   module subroutine collision_resolve_make_impactors_pl(pl, idx)
      !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
      !!
      !! When a single body is involved in more than one collision in a single step, it becomes part of a collision family
      !! The largest body involved in a multi-body collision is the "parent" and all bodies that collide with it are its "children,"
      !! including those that collide with the children.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      ! Arguments
      class(base_object),           intent(inout) :: pl !! Swiftest massive body object
      integer(I4B), dimension(:), intent(in)      :: idx  !! Array holding the indices of the two bodies involved in the collision
      ! Internals
      integer(I4B)                              :: i, j, index_parent, index_child, p1, p2
      integer(I4B)                              :: nchild_inherit, nchild_orig, nchild_new
      integer(I4B), dimension(:), allocatable   :: temp

      select type(pl)
      class is (swiftest_pl)

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
      end select

      return
   end subroutine collision_resolve_make_impactors_pl


   module subroutine collision_resolve_mergeaddsub(nbody_system, param, t, status)
      !! author:  David A. Minton
      !!
      !! Fills the pl_discards and pl_adds with removed and added bodies
      !!  
      use symba, only : symba_pl
      implicit none
      ! Arguments
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters with Swiftest additions
      real(DP),                 intent(in)    :: t            !! Time of collision
      integer(I4B),             intent(in)    :: status       !! Status flag to assign to adds
      ! Internals
      integer(I4B) :: i, ibiggest, ismallest, iother, nimpactors, nfrag, nameidx
      logical, dimension(:), allocatable  :: lmask
      class(swiftest_pl), allocatable :: plnew, plsub
      character(*), parameter :: FRAGFMT = '("Newbody",I0.7)'
      character(*), parameter :: MERGEFMT = '(A,I0.7)'
      character(*), parameter :: MERGE_PREPEND_TEXT = "_MERGE"
      integer(I4B) :: merge_text_length 
      character(len=NAMELEN) :: merge_text
      character(len=NAMELEN) :: newname, origin_type
      real(DP) :: volume
  
      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         associate(pl => nbody_system%pl, pl_discards => nbody_system%pl_discards, info => nbody_system%pl%info, &
                   pl_adds => nbody_system%pl_adds, cb => nbody_system%cb, collider => nbody_system%collider,  &
                   impactors => nbody_system%collider%impactors,fragments => nbody_system%collider%fragments)

            ! Add the impactors%id bodies to the subtraction list
            nimpactors = impactors%ncoll
            nfrag = fragments%nbody

            ! Setup new bodies
            allocate(plnew, mold=pl)
            call plnew%setup(nfrag, param)
            ibiggest  = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
            ismallest = impactors%id(minloc(pl%Gmass(impactors%id(:)), dim=1))

            ! Copy over identification, information, and physical properties of the new bodies from the fragment list
            plnew%id(1:nfrag) = fragments%id(1:nfrag) 
            plnew%rb(:, 1:nfrag) = fragments%rb(:, 1:nfrag) 
            plnew%vb(:, 1:nfrag) = fragments%vb(:, 1:nfrag)
            plnew%status(1:nfrag) = ACTIVE
            call plnew%b2h(cb)
            plnew%mass(1:nfrag) = fragments%mass(1:nfrag)
            plnew%Gmass(1:nfrag) = param%GU * fragments%mass(1:nfrag)
            plnew%radius(1:nfrag) = fragments%radius(1:nfrag)
            if (allocated(pl%a)) call plnew%xv2el(cb)

            if (param%lrotation) then
               plnew%Ip(:, 1:nfrag) = fragments%Ip(:, 1:nfrag)
               plnew%rot(:, 1:nfrag) = fragments%rot(:, 1:nfrag)
            end if

            call plnew%set_rhill(cb)

            ! if (param%ltides) then
            !    plnew%Q = pl%Q(ibiggest)
            !    plnew%k2 = pl%k2(ibiggest)
            !    plnew%tlag = pl%tlag(ibiggest)
            ! end if

#ifdef DOCONLOC
            do concurrent(i = 1:nfrag) shared(plnew,fragments) local(volume)
#else
            do concurrent(i = 1:nfrag)
#endif
               volume = 4.0_DP/3.0_DP * PI * plnew%radius(i)**3
               plnew%density(i) = fragments%mass(i) / volume
            end do

            select case(status)
            case(SUPERCATASTROPHIC)
               plnew%status(1:nfrag) = NEW_PARTICLE
               do i = 1, nfrag
                  write(newname, FRAGFMT) fragments%id(i)
                  call plnew%info(i)%set_value(origin_type="Supercatastrophic", origin_time=t, name=newname, &
                                                origin_rh=plnew%rh(:,i), origin_vh=plnew%vh(:,i), &
                                                collision_id=collider%maxid_collision)
               end do
               do i = 1, nimpactors
                  if (impactors%id(i) == ibiggest) then
                     iother = ismallest
                  else
                     iother = ibiggest
                  end if
                  call pl%info(impactors%id(i))%set_value(status="Supercatastrophic", discard_time=t, &
                                                            discard_rh=pl%rh(:,i), discard_vh=pl%vh(:,i), &
                                                            discard_body_id=iother)
               end do
            case(DISRUPTED,HIT_AND_RUN_DISRUPT)
               if (status == DISRUPTED) then
                  write(origin_type,*) "Disruption"
               else if (status == HIT_AND_RUN_DISRUPT) then
                  write(origin_type,*) "Hit and run fragmentation"
               end if
               call plnew%info(1)%copy(pl%info(ibiggest))
               plnew%status(1) = OLD_PARTICLE
               do i = 2, nfrag
                  write(newname, FRAGFMT) fragments%id(i)
                  call plnew%info(i)%set_value(origin_type=origin_type, origin_time=t, name=newname, &
                                                origin_rh=plnew%rh(:,i), origin_vh=plnew%vh(:,i), &
                                                collision_id=collider%maxid_collision)
               end do
               do i = 1, nimpactors
                  if (impactors%id(i) == ibiggest) cycle
                  iother = ibiggest
                  call pl%info(impactors%id(i))%set_value(status=origin_type, discard_time=t, &
                                                            discard_rh=pl%rh(:,i), discard_vh=pl%vh(:,i), &
                                                            discard_body_id=iother)
               end do 
            case(MERGED)
               write(origin_type,*) "Merger"
               call plnew%info(1)%copy(pl%info(ibiggest))
               nbody_system%maxid = nbody_system%maxid + 1
               plnew%id(1) = nbody_system%maxid

               ! Appends an index number to the end of the original name to make it unique, but still identifiable as the original.
               ! If there is already an index number appended, replace it
               write(merge_text,MERGEFMT) MERGE_PREPEND_TEXT,plnew%id(1)
               merge_text_length = len(trim(adjustl(merge_text)))
               nameidx = index(plnew%info(1)%name, MERGE_PREPEND_TEXT) - 1
               if (nameidx < 0) nameidx = min(len(trim(adjustl(plnew%info(1)%name))), NAMELEN - merge_text_length)
               write(newname,*) trim(adjustl(plnew%info(1)%name(1:nameidx))) // trim(adjustl(merge_text))
               plnew%status(1) = NEW_PARTICLE
               call plnew%info(1)%set_value(origin_type=origin_type, origin_time=t, name=newname, &
                                            origin_rh=plnew%rh(:,1), origin_vh=plnew%vh(:,1), &
                                            collision_id=collider%maxid_collision)
               do i = 1, nimpactors
                  if (impactors%id(i) == ibiggest) cycle

                  iother = ibiggest
                  call pl%info(impactors%id(i))%set_value(status="MERGED", discard_time=t, discard_rh=pl%rh(:,i), &
                                                            discard_vh=pl%vh(:,i), discard_body_id=iother)
               end do 
            end select

            !Copy over or set integration parameters for new bodies
            plnew%lcollision(1:nfrag) = .false.
            plnew%ldiscard(1:nfrag) = .false.
            select type(pl)
            class is (symba_pl)
            select type(plnew)
            class is (symba_pl)
               plnew%levelg(1:nfrag) = pl%levelg(ibiggest)
               plnew%levelm(1:nfrag) = pl%levelm(ibiggest)
            end select
            end select

            plnew%lmtiny(1:nfrag) = plnew%Gmass(1:nfrag) < param%GMTINY
            where(plnew%lmtiny(1:nfrag))
               plnew%info(1:nfrag)%particle_type = PL_TINY_TYPE_NAME 
            elsewhere
               plnew%info(1:nfrag)%particle_type = PL_TYPE_NAME 
            end where

            ! Append the new merged body to the list 
            call pl_adds%append(plnew, lsource_mask=[(.true., i=1, nfrag)])

            ! Add the discarded bodies to the discard list
            pl%status(impactors%id(:)) = MERGED
            pl%ldiscard(impactors%id(:)) = .true.
            pl%lcollision(impactors%id(:)) = .true.
            allocate(lmask, mold=pl%lmask)
            lmask(:) = .false.
            lmask(impactors%id(:)) = .true.
            
            allocate(plsub, mold=pl)
            call pl%spill(plsub, lmask, ldestructive=.false.)

            ! call pl_discards%append(plsub, lsource_mask=[(.true., i = 1, nimpactors)])

            ! Save the before/after snapshots
            select type(before => collider%before)
            class is (swiftest_nbody_system)
               call move_alloc(plsub, before%pl)
            end select

            select type(after => collider%after)
            class is (swiftest_nbody_system)
               call move_alloc(plnew, after%pl)
            end select

         end associate

      end select
      end select 
   
      return
   end subroutine collision_resolve_mergeaddsub


   module subroutine collision_resolve_plpl(self, nbody_system, param, t, dt, irec)
      !! author: David A. Minton
      !! 
      !! Process the pl-pl collision list, then modifiy the massive bodies based on the outcome of the collision
      !! 
      implicit none
      ! Arguments
      class(collision_list_plpl), intent(inout) :: self   !! Swiftest pl-pl encounter list
      class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),     intent(inout) :: param  !! Current run configuration parameters with Swiftest additions
      real(DP),                   intent(in)    :: t      !! Current simulation time
      real(DP),                   intent(in)    :: dt     !! Current simulation step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
      ! Internals
      real(DP) :: E_before, E_after, mnew
      real(DP), dimension(NDIM) ::L_before, L_after, dL
      logical :: lplpl_collision
      character(len=STRMAX) :: timestr, idstr
      integer(I4B), dimension(2) :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      logical  :: lgoodcollision
      integer(I4B) :: i, j, nnew, loop
      integer(I8B) :: k, ncollisions
      integer(I4B), dimension(:), allocatable :: idnew
      integer(I4B), parameter :: MAXCASCADE = 1000
  
      select type (nbody_system)
      class is (swiftest_nbody_system)
      select type(pl => nbody_system%pl)
      class is (swiftest_pl)
      select type(param)
      class is (swiftest_parameters)
         associate(plpl_collision => nbody_system%plpl_collision, &
                   collision_history => nbody_system%collision_history, pl => nbody_system%pl, cb => nbody_system%cb, &
                   collider => nbody_system%collider, fragments => nbody_system%collider%fragments, &
                   impactors => nbody_system%collider%impactors)
            if (plpl_collision%nenc == 0) return ! No collisions to resolve


            ! Make sure that the heliocentric and barycentric coordinates are consistent with each other
            call pl%vb2vh(nbody_system%cb) 
            call pl%rh2rb(nbody_system%cb)

            ! Get the energy before the collision is resolved
            if (param%lenergy) then
               call nbody_system%get_energy_and_momentum(param)
               E_before = nbody_system%te
               L_before(:) = nbody_system%L_total(:)
            end if

            do loop = 1, MAXCASCADE
               associate( idx1 => plpl_collision%index1, idx2 => plpl_collision%index2)
                  ncollisions = plpl_collision%nenc
                  write(timestr,*) t
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT,&
                                                            "***********************************************************" // &
                                                            "***********************************************************")
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Collision between massive bodies detected at time t = " // &
                                                            trim(adjustl(timestr)))
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, &
                                                            "***********************************************************" // &
                                                            "***********************************************************")

                  do k = 1_I8B, ncollisions
                     idx_parent(1) = pl%kin(idx1(k))%parent
                     idx_parent(2) = pl%kin(idx2(k))%parent
                     call impactors%consolidate(nbody_system, param, idx_parent, lgoodcollision)
                     if ((.not. lgoodcollision) .or. any(pl%status(idx_parent(:)) /= COLLIDED)) cycle

                     ! Advance the collision id number and save it
                     collider%maxid_collision = max(collider%maxid_collision, maxval(nbody_system%pl%info(:)%collision_id))
                     collider%maxid_collision = collider%maxid_collision + 1
                     collider%collision_id = collider%maxid_collision
                     write(idstr,*) collider%collision_id
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT, "collision_id " // trim(adjustl(idstr)))

                     ! Get the collision regime
                     call collider%get_regime(nbody_system, param)

                     call collision_history%take_snapshot(param,nbody_system, t, "before") 

                     ! Generate the new bodies resulting from the collision
                     call collider%generate(nbody_system, param, t)

                     call collision_history%take_snapshot(param,nbody_system, t, "after") 

                     plpl_collision%status(k) = collider%status
                     call impactors%dealloc()
                  end do

                  ! Destroy the collision list now that the collisions are resolved
                  call plpl_collision%setup(0_I8B)

                  if ((nbody_system%pl_adds%nbody == 0) .and. (.not.any(pl%ldiscard(:)))) exit
                  if (allocated(idnew)) deallocate(idnew)
                  nnew = nbody_system%pl_adds%nbody
                  allocate(idnew, source=nbody_system%pl_adds%id)
                  mnew = sum(nbody_system%pl_adds%mass(:))

                  ! Rearrange the arrays: Remove discarded bodies, add any new bodies, re-sort, and recompute all indices and 
                  ! encounter lists
                  call pl%rearray(nbody_system, param)

                  ! Destroy the add/discard list so that we don't append the same body multiple times if another collision 
                  ! is detected
                  call nbody_system%pl_discards%setup(0, param)
                  call nbody_system%pl_adds%setup(0, param)

                  if (param%lenergy) then
                     call nbody_system%get_energy_and_momentum(param)
                     L_after(:) = nbody_system%L_total(:)
                     dL = L_after(:) - L_before(:)

                     ! Add some velocity torque to the new bodies to remove residual angular momentum difference
                     do j = 1, nnew
                        i = findloc(pl%id,idnew(j),dim=1)
                        if (i == 0) cycle
                        call collision_util_velocity_torque(-dL * pl%mass(i)/mnew, pl%mass(i), pl%rb(:,i), pl%vb(:,i)) 
                     end do

                     call nbody_system%get_energy_and_momentum(param)
                     E_after = nbody_system%te
                     nbody_system%E_collisions = nbody_system%E_collisions + (E_after - E_before)
                     L_after(:) = nbody_system%L_total(:)
                     dL = L_after(:) - L_before(:)
                  end if


                  ! Check whether or not any of the particles that were just added are themselves in a collision state. This will 
                  ! generate a new plpl_collision 
                  call self%collision_check(nbody_system, param, t, dt, irec, lplpl_collision)

                  if (.not.lplpl_collision) exit
                  if (loop == MAXCASCADE) then
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT,"A runaway collisional cascade has been detected in " // &
                                                                        "collision_resolve_plpl.")
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT,"Consider reducing the step size or changing the " // &
                                                                        "parameters in the collisional model to reduce the " // &
                                                                        "number of fragments.")
                     call base_util_exit(FAILURE,unit=param%display_unit)
                  end if
               end associate
            end do

         end associate
      end select
      end select
      end select

      return
   end subroutine collision_resolve_plpl


   module subroutine collision_resolve_pltp(self, nbody_system, param, t, dt, irec)
      !! author: David A. Minton
      !! 
      !! Process the pl-tp collision list, then modifiy the massive bodies based on the outcome of the collision
      !! 
      implicit none
      ! Arguments
      class(collision_list_pltp), intent(inout) :: self   !! Swiftest pl-pl encounter list
      class(base_nbody_system),   intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),     intent(inout) :: param  !! Current run configuration parameters with Swiftest additions
      real(DP),                   intent(in)    :: t      !! Current simulation tim
      real(DP),                   intent(in)    :: dt     !! Current simulation step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
      ! Internals
      class(swiftest_pl), allocatable :: plsub
      class(swiftest_tp), allocatable :: tpsub
      logical :: lpltp_collision
      character(len=STRMAX) :: timestr, idstr
      integer(I4B) :: i, j, nnew, loop
      integer(I8B) :: k, ncollisions
      integer(I4B), dimension(:), allocatable :: idnew      
      logical, dimension(:), allocatable :: lmask
     
      ! Make sure coordinate systems are all synced up due to being inside the recursion at this point
      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
         associate(pltp_collision => nbody_system%pltp_collision, &
            collision_history => nbody_system%collision_history, pl => nbody_system%pl, cb => nbody_system%cb, &
            tp => nbody_system%tp, collider => nbody_system%collider, impactors => nbody_system%collider%impactors)
            call pl%vb2vh(nbody_system%cb)
            call tp%vb2vh(nbody_system%cb%vb)
            call pl%b2h(nbody_system%cb)
            call tp%b2h(nbody_system%cb)

            ! Restructure the massive bodies based on the outcome of the collision
            call tp%rearray(nbody_system, param)

            ! Check for discards
            call nbody_system%tp%discard(nbody_system, param)

            associate(idx1 => pltp_collision%index1, idx2 => pltp_collision%index2)
               ncollisions = pltp_collision%nenc
               write(timestr,*) t
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "")
               call swiftest_io_log_one_message(COLLISION_LOG_OUT,"***********************************************************" // &
                                                                  "***********************************************************")
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Collision between test particle and massive body detected " // &
                                                                   "at time t = " // trim(adjustl(timestr)))
               call swiftest_io_log_one_message(COLLISION_LOG_OUT,"***********************************************************" // &
                                                                  "***********************************************************")

               do k = 1_I8B, ncollisions
                  ! Advance the collision id number and save it
                  collider%maxid_collision = max(collider%maxid_collision, maxval(nbody_system%pl%info(:)%collision_id))
                  collider%maxid_collision = collider%maxid_collision + 1
                  collider%collision_id = collider%maxid_collision
                  write(idstr,*) collider%collision_id
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "collision_id " // trim(adjustl(idstr)))
                  collider%impactors%regime = COLLRESOLVE_REGIME_MERGE
                  allocate(lmask, mold=pl%lmask)
                  lmask(:) = .false.
                  lmask(idx1(k)) = .true.
                  
                  allocate(plsub, mold=pl)
                  call pl%spill(plsub, lmask, ldestructive=.false.)
      
                  ! Save the before snapshots
                  select type(before => collider%before)
                  class is (swiftest_nbody_system)
                     call move_alloc(plsub, before%pl)
                  end select

                  deallocate(lmask)
                  allocate(lmask, mold=tp%lmask)
                  lmask(:) = .false.
                  lmask(idx2(k)) = .true.
                  
                  allocate(tpsub, mold=tp)
                  call tp%spill(tpsub, lmask, ldestructive=.false.)
      
                  ! Save the before snapshots
                  select type(before => collider%before)
                  class is (swiftest_nbody_system)
                     call move_alloc(tpsub, before%tp)
                  end select

                  call collision_history%take_snapshot(param,nbody_system, t, "particle") 

                  call impactors%dealloc()
               end do

               ! Destroy the collision list now that the collisions are resolved
               call pltp_collision%setup(0_I8B)

            end associate

         end associate
      end select
      end select

      return
   end subroutine collision_resolve_pltp

end submodule s_collision_resolve
