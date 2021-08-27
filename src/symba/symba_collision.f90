submodule (symba_classes) s_symba_collision
   use swiftest

   integer(I4B), parameter :: NFRAG_DISRUPT = 12
   integer(I4B), parameter :: NFRAG_SUPERCAT = 20
contains

   module function symba_collision_casedisruption(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic disruption collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
      real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, istart, nfrag, nfamily, nstart, nend
      real(DP)                                :: mtot, avg_dens
      real(DP), dimension(NDIM)               :: xcom, vcom, Ip_new
      real(DP), dimension(2)                  :: vol
      real(DP), dimension(:, :), allocatable  :: vb_frag, xb_frag, rot_frag, Ip_frag
      real(DP), dimension(:), allocatable     :: m_frag, rad_frag
      integer(I4B), dimension(:), allocatable :: id_frag
      logical                                 :: lfailure

      write(*, '("Disruption between bodies ",I8,99(:,",",I8))') system%pl%id(family(:))

      ! Collisional fragments will be uniformly distributed around the pre-impact barycenter
      nfrag = NFRAG_DISRUPT 
      allocate(m_frag(nfrag))
      allocate(rad_frag(nfrag))
      allocate(xb_frag(NDIM, nfrag))
      allocate(vb_frag(NDIM, nfrag))
      allocate(rot_frag(NDIM, nfrag))
      allocate(Ip_frag(NDIM, nfrag))
      allocate(id_frag(nfrag))

      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot

      ! Get mass weighted mean of Ip and average density
      Ip_new(:) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / mtot
      vol(:) = 4._DP / 3._DP * PI * radius(:)**3
      avg_dens = mtot / sum(vol(:))

      ! Distribute the mass among fragments, with a branch to check for the size of the second largest fragment
      m_frag(1) = mass_res(1)
      if (mass_res(2) > mass_res(1) / 3._DP) then
         m_frag(2) = mass_res(2)
         istart = 3
      else
         istart = 2
      end if
      ! Distribute remaining mass among the remaining bodies
      do i = istart, nfrag
         m_frag(i) = (mtot - sum(m_frag(1:istart - 1))) / (nfrag - istart + 1) 
      end do

      ! Distribute any residual mass if there is any and set the radius
      m_frag(nfrag) = m_frag(nfrag) + (mtot - sum(m_frag(:)))
      rad_frag(1:nfrag) = (3 * m_frag(:) / (4 * PI * avg_dens))**(1.0_DP / 3.0_DP)
      id_frag(1:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag)]
      param%maxid = id_frag(nfrag)

      do i = 1, nfrag
         Ip_frag(:, i) = Ip_new(:)
      end do

      call fragmentation_initialize(system, param, family, x, v, L_spin, Ip, mass, radius, &
                                    nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lfailure)

      if (lfailure) then
         write(*,*) 'No fragment solution found, so treat as a pure hit-and-run'
         status = ACTIVE 
         nfrag = 0
         select type(pl => system%pl)
         class is (symba_pl)
            pl%status(family(:)) = status
            pl%ldiscard(family(:)) = .false.
            pl%lcollision(family(:)) = .false.
         end select
      else
         ! Populate the list of new bodies
         write(*,'("Generating ",I2.0," fragments")') nfrag
         status = DISRUPTION
         call symba_collision_mergeaddsub(system, param, family, id_frag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, status)
      end if

      return
   end function symba_collision_casedisruption


   module function symba_collision_casehitandrun(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
      real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, nfrag, jproj, jtarg, idstart, ibiggest, nfamily
      real(DP)                                :: mtot, avg_dens
      real(DP), dimension(NDIM)               :: xcom, vcom
      real(DP), dimension(2)                  :: vol
      real(DP), dimension(:, :), allocatable  :: vb_frag, xb_frag, rot_frag, Ip_frag
      real(DP), dimension(:), allocatable     :: m_frag, rad_frag
      integer(I4B), dimension(:), allocatable :: id_frag
      logical                                 :: lpure
      logical,  dimension(system%pl%nbody)    :: lmask

      write(*, '("Hit and run between bodies ",I8,99(:,",",I8))')  system%pl%id(family(:))

      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot
      lpure = .false.

      ! The largest body will stay untouched
      if (mass(1) > mass(2)) then
         jtarg = 1
         jproj = 2
      else
         jtarg = 2
         jproj = 1
      end if

      if (mass_res(2) > 0.9_DP * mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
         write(*,*) 'Pure hit and run. No new fragments generated.'
         nfrag = 0
         lpure = .true.
      else ! Imperfect hit and run, so we'll keep the largest body and destroy the other
         nfrag = NFRAG_DISRUPT - 1
         lpure = .false.
         allocate(m_frag(nfrag))
         allocate(id_frag(nfrag))
         allocate(rad_frag(nfrag))
         allocate(xb_frag(NDIM, nfrag))
         allocate(vb_frag(NDIM, nfrag))
         allocate(rot_frag(NDIM, nfrag))
         allocate(Ip_frag(NDIM, nfrag))
         m_frag(1) = mass(jtarg)
         ibiggest = family(maxloc(system%pl%Gmass(family(:)), dim=1))
         id_frag(1) = system%pl%id(ibiggest)
         rad_frag(1) = radius(jtarg)
         xb_frag(:, 1) = x(:, jtarg) 
         vb_frag(:, 1) = v(:, jtarg)
         Ip_frag(:,1) = Ip(:, jtarg)

         ! Get mass weighted mean of Ip and average density
         vol(:) = 4._DP / 3._DP * pi * radius(:)**3
         avg_dens = mass(jproj) / vol(jproj)
         m_frag(2:nfrag) = (mtot - m_frag(1)) / (nfrag - 1) 
         rad_frag(2:nfrag) = (3 * m_frag(2:nfrag) / (4 * PI * avg_dens))**(1.0_DP / 3.0_DP)
         m_frag(nfrag) = m_frag(nfrag) + (mtot - sum(m_frag(:)))
         id_frag(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
         param%maxid = id_frag(nfrag)

         do i = 1, nfrag
            Ip_frag(:, i) = Ip(:, jproj)
         end do

         ! Put the fragments on the circle surrounding the center of mass of the system
         call fragmentation_initialize(system, param, family, x, v, L_spin, Ip, mass, radius, &
                           nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lpure)
         if (lpure) then
            write(*,*) 'Should have been a pure hit and run instead'
            nfrag = 0
         else
            write(*,'("Generating ",I2.0," fragments")') nfrag
         end if
      end if
      if (lpure) then ! Reset these bodies back to being active so that nothing further is done to them
         status = HIT_AND_RUN_PURE
         select type(pl => system%pl)
         class is (symba_pl)
            pl%status(family(:)) = ACTIVE
            pl%ldiscard(family(:)) = .false.
            pl%lcollision(family(:)) = .false.
         end select
      else
         status = HIT_AND_RUN_DISRUPT
         call symba_collision_mergeaddsub(system, param, family, id_frag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, status)
      end if

      return
   end function symba_collision_casehitandrun


   module function symba_collision_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Merge planets.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routines symba_merge_pl.f90 and symba_discard_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routines symba5_merge.f and discard_mass_merge.f
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(in)    :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(in)    :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                              :: i, j, k, ibiggest, nfamily
      real(DP)                                  :: volume_new, pe
      real(DP), dimension(NDIM)                 :: xc, vc, xcrossv
      real(DP), dimension(2)                    :: vol
      real(DP), dimension(NDIM)                 :: L_orb_old, L_spin_old
      real(DP), dimension(NDIM)                 :: L_spin_new
      logical,  dimension(system%pl%nbody)      :: lmask
      real(DP), dimension(NDIM, 1)              :: vb_frag, xb_frag, rot_frag, Ip_frag
      real(DP), dimension(1)                    :: m_frag, rad_frag
      integer(I4B), dimension(1)                :: id_frag

      select type(pl => system%pl)
      class is (symba_pl)
         write(*, '("Merging bodies ",I8,99(:,",",I8))') pl%id(family(:))

         ibiggest = family(maxloc(pl%Gmass(family(:)), dim=1))
         id_frag(1) = pl%id(ibiggest)

         m_frag(1) = sum(mass(:))

         ! Merged body is created at the barycenter of the original bodies
         xb_frag(:,1) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / m_frag(1)
         vb_frag(:,1) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / m_frag(1)

         ! Get mass weighted mean of Ip and 
         vol(:) = 4._DP / 3._DP * PI * radius(:)**3
         volume_new = sum(vol(:))
         rad_frag(1) = (3 * volume_new / (4 * PI))**(1._DP / 3._DP)

         L_orb_old(:) = 0.0_DP

         ! Compute orbital angular momentum of pre-impact system
         do i = 1, 2
            xc(:) = x(:, i) - xb_frag(:,1)
            vc(:) = v(:, i) - vb_frag(:,1)
            xcrossv(:) = xc(:) .cross. vc(:)
            L_orb_old(:) = L_orb_old(:) + mass(i) * xcrossv(:)
         end do
      
         if (param%lrotation) then
            Ip_frag(:,1) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / m_frag(1)
            L_spin_old(:) = L_spin(:,1) + L_spin(:,2)

            ! Conserve angular momentum by putting pre-impact orbital momentum into spin of the new body
            L_spin_new(:) = L_orb_old(:) + L_spin_old(:) 

            ! Assume prinicpal axis rotation on 3rd Ip axis
            rot_frag(:,1) = L_spin_new(:) / (Ip_frag(3,1) * m_frag(1) * rad_frag(1)**2)
         else ! If spin is not enabled, we will consider the lost pre-collision angular momentum as "escaped" and add it to our bookkeeping variable
            param%Lescape(:) = param%Lescape(:) + L_orb_old(:) 
         end if

         ! Keep track of the component of potential energy due to the pre-impact family for book-keeping
         nfamily = size(family(:))
         pe = 0.0_DP
         do j = 1, nfamily
            do i = j + 1, nfamily
               pe = pe - pl%Gmass(i) * pl%mass(j) / norm2(pl%xb(:, i) - pl%xb(:, j))
            end do
         end do
         param%Ecollisions = param%Ecollisions + pe 
         param%Euntracked  = param%Euntracked - pe 

         ! Update any encounter lists that have the removed bodies in them so that they instead point to the new 
         do k = 1, system%plplenc_list%nenc
            do j = 1, nfamily
               i = family(j)
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
         call symba_collision_mergeaddsub(system, param, family, id_frag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, status)

      end select

      return 
   end function symba_collision_casemerge


   module function symba_collision_casesupercatastrophic(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a supercatastrophic collision
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
      real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, nfrag, ibiggest, nfamily, nstart, nend
      real(DP)                                :: mtot, avg_dens, min_frag_mass
      real(DP), dimension(NDIM)               :: xcom, vcom
      real(DP), dimension(2)                  :: vol
      real(DP), dimension(NDIM)               :: Ip_new
      real(DP), dimension(:, :), allocatable  :: vb_frag, xb_frag, rot_frag, Ip_frag
      real(DP), dimension(:), allocatable     :: m_frag, rad_frag
      integer(I4B), dimension(:), allocatable :: id_frag
      logical                                 :: lfailure
      logical,  dimension(system%pl%nbody)    :: lmask

      write(*, '("Supercatastrophic disruption between bodies ",I8,99(:,",",I8))')  system%pl%id(family(:))

      ! Collisional fragments will be uniformly distributed around the pre-impact barycenter
      nfrag = NFRAG_SUPERCAT
      allocate(m_frag(nfrag))
      allocate(rad_frag(nfrag))
      allocate(id_frag(nfrag))
      allocate(xb_frag(NDIM, nfrag))
      allocate(vb_frag(NDIM, nfrag))
      allocate(rot_frag(NDIM, nfrag))
      allocate(Ip_frag(NDIM, nfrag))

      mtot = sum(mass(:))
      xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot 
      vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot

      ! Get mass weighted mean of Ip and average density
      Ip_new(:) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / mtot
      vol(:) = 4._DP / 3._DP * pi * radius(:)**3
      avg_dens = mtot / sum(vol(:))

      ! If we are adding the first and largest fragment (lr), check to see if its mass is SMALLER than an equal distribution of 
      ! mass between all fragments. If so, we will just distribute the mass equally between the fragments
      min_frag_mass = mtot / nfrag
      if (mass_res(1) < min_frag_mass) then
         m_frag(:) = min_frag_mass
      else
         m_frag(1) = mass_res(1)
         m_frag(2:nfrag) = (mtot - mass_res(1)) / (nfrag - 1)
      end if
      ! Distribute any residual mass if there is any and set the radius
      m_frag(nfrag) = m_frag(nfrag) + (mtot - sum(m_frag(:)))
      rad_frag(1:nfrag) = (3 * m_frag(:) / (4 * PI * avg_dens))**(1.0_DP / 3.0_DP)
      id_frag(1:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag)]
      param%maxid = id_frag(nfrag)

      do i = 1, nfrag
         Ip_frag(:, i) = Ip_new(:)
      end do

      call fragmentation_initialize(system, param, family, x, v, L_spin, Ip, mass, radius, &
                                    nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lfailure)

      if (lfailure) then
         write(*,*) 'No fragment solution found, so treat as a pure hit-and-run'
         status = ACTIVE 
         nfrag = 0
         select type(pl => system%pl)
         class is (symba_pl)
            pl%status(family(:)) = status
            pl%ldiscard(family(:)) = .false.
            pl%lcollision(family(:)) = .false.
         end select
      else
         ! Populate the list of new bodies
         write(*,'("Generating ",I2.0," fragments")') nfrag
         status = SUPERCATASTROPHIC
         call symba_collision_mergeaddsub(system, param, family, id_frag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, status)
      end if

      return
   end function symba_collision_casesupercatastrophic


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
      integer(I4B)                              :: k
      real(DP)                                  :: rlim, Gmtot
      logical                                   :: isplpl

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
                     Gmtot = pl%Gmass(ind1(k)) + pl%Gmass(ind2(k))
                     lcollision(k) = symba_collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), Gmtot, rlim, dt, self%lvdotr(k))
                  end do
               else
                  do concurrent(k = 1:nenc, lmask(k))
                     xr(:) = pl%xh(:, ind1(k)) - tp%xh(:, ind2(k)) 
                     vr(:) = pl%vb(:, ind1(k)) - tp%vb(:, ind2(k))
                     lcollision(k) = symba_collision_check_one(xr(1), xr(2), xr(3), vr(1), vr(2), vr(3), pl%Gmass(ind1(k)), pl%radius(ind1(k)), dt, self%lvdotr(k))
                  end do
               end if

               if (any(lcollision(1:nenc))) call pl%xh2xb(system%cb) ! Update the central body barycenteric position vector to get us out of DH and into bary
               do k = 1, nenc
                  if (lcollision(k)) self%status(k) = COLLISION
                  self%t(k) = t
                  self%x1(:,k) = pl%xh(:,ind1(k)) + system%cb%xb(:)
                  self%v1(:,k) = pl%vb(:,ind1(k)) 
                  if (isplpl) then
                     self%x2(:,k) = pl%xh(:,ind2(k)) + system%cb%xb(:)
                     self%v2(:,k) = pl%vb(:,ind2(k)) 
                     if (lcollision(k)) then
                        ! Check to see if either of these bodies has been involved with a collision before, and if so, make this a collisional family
                        if (pl%lcollision(ind1(k)) .or. pl%lcollision(ind2(k))) call pl%make_family([ind1(k),ind2(k)])
   
                        ! Set the collision flag for these to bodies to true in case they become involved in another collision later in the step
                        pl%lcollision([ind1(k), ind2(k)]) = .true.
                        pl%ldiscard([ind1(k), ind2(k)]) = .true.
                        pl%status([ind1(k), ind2(k)]) = COLLISION
                     end if
                  else
                     self%x2(:,k) = tp%xh(:,ind2(k)) + system%cb%xb(:)
                     self%v2(:,k) = tp%vb(:,ind2(k)) 
                     if (lcollision(k)) then
                        tp%status(ind2(k)) = DISCARDED_PLR
                        tp%ldiscard(ind2(k)) = .true.
                        write(*,*) 'Test particle ',tp%id(ind2(k)), ' collided with massive body ',pl%id(ind1(k)), ' at time ',t
                     end if
                  end if
               end do
            end associate
         end select
      end select

      lany_collision = any(lcollision(:))

      ! Extract the pl-pl encounter list and return the plplcollision_list
      if (lany_collision) then
         select type(plplenc_list => self)
         class is (symba_plplenc)
            call plplenc_list%extract_collisions(system, param)
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
               call orbel_xv2aeq(Gmtot, [xr, yr, zr], [vxr, vyr, vzr], a, e, q)
               lcollision = (q < rlim) 
            end if
         end if
      end if

      return
   end function symba_collision_check_one


   function symba_collision_consolidate_familes(pl, cb, param, idx_parent, family, x, v, mass, radius, L_spin, Ip) result(lflag)
      !! author: David A. Minton
      !! 
      !! Loops through the pl-pl collision list and groups families together by index. Outputs the indices of all family members, 
      !! and pairs of quantities (x and v vectors, mass, radius, L_spin, and Ip) that can be used to resolve the collisional outcome.
      implicit none
      ! Arguments
      class(symba_pl),                                 intent(inout) :: pl               !! SyMBA massive body object
      class(symba_cb),                                 intent(inout) :: cb               !! SyMBA central body object
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
         x(:, j)  = pl%xh(:, idx_parent(j)) + cb%xb(:)
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
               xchild(:) = pl%xh(:, idx_child) + cb%xb(:)
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
      logical,      dimension(self%nenc)      :: lplpl_collision
      logical,      dimension(:), allocatable :: lplpl_unique_parent
      integer(I4B), dimension(:), pointer     :: plparent
      integer(I4B), dimension(:), allocatable :: collision_idx, unique_parent_idx
      integer(I4B)                            :: i, index_coll, ncollisions, nunique_parent

      select type (pl => system%pl)
      class is (symba_pl)
         associate(plplenc_list => self, nplplenc => self%nenc, idx1 => self%index1, idx2 => self%index2, plparent => pl%kin%parent)
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


   subroutine symba_collision_mergeaddsub(system, param, family, id_frag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, status)
      !! author:  David A. Minton
      !!
      !! Fills the pl_discards and pl_adds with removed and added bodies
      !!  
      implicit none
      ! Arguments
      class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
      class(symba_parameters),         intent(inout) :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      integer(I4B),    dimension(:),   intent(in)    :: id_frag          !! List of fragment ids
      real(DP),        dimension(:),   intent(in)    :: m_frag, rad_frag !! Distribution of fragment mass and radii
      real(DP),        dimension(:,:), intent(in)    :: Ip_frag          !! Fragment rotational inertia vectors
      real(DP),        dimension(:,:), intent(in)    :: xb_frag, vb_frag, rot_frag !! Fragment barycentric position, barycentric velocity, and rotation vectors
      integer(I4B),                    intent(in)    :: status           !! Status flag to assign to adds
      ! Internals
      integer(I4B) :: i, ibiggest, nstart, nend, nfamily, nfrag
      logical, dimension(system%pl%nbody)    :: lmask
      class(symba_pl), allocatable            :: plnew
   
      select type(pl => system%pl)
      class is (symba_pl)
         select type(pl_discards => system%pl_discards)
         class is (symba_merger)
            select type(info => pl%info)
            class is (symba_particle_info)
               associate(pl_adds => system%pl_adds, cb => system%cb)
      
                  ! Add the family bodies to the subtraction list
                  nfamily = size(family(:))
                  nfrag   = size(m_frag(:))
                  lmask(:) = .false.
                  lmask(family(:)) = .true.
                  pl%status(family(:)) = MERGED
                  nstart = pl_discards%nbody + 1
                  nend = pl_discards%nbody + nfamily
                  call pl_discards%append(pl, lmask)
                  pl%ldiscard(family(:)) = .true.
                  pl%lcollision(family(:)) = .true.
      
                  ! Record how many bodies were subtracted in this event
                  pl_discards%ncomp(nstart:nend) = nfamily 
      
                  ! Setup new bodies
                  allocate(plnew, mold=pl)
                  call plnew%setup(nfrag, param)
                  ibiggest = family(maxloc(pl%Gmass(family(:)), dim=1))
   
                  ! Copy over identification, information, and physical properties of the new bodies from the fragment list
                  plnew%id(1:nfrag) = id_frag(1:nfrag) 
                  param%maxid = param%maxid + nfrag
                  plnew%xb(:, 1:nfrag) = xb_frag(:, 1:nfrag) 
                  plnew%vb(:, 1:nfrag) = vb_frag(:, 1:nfrag)
                  do i = 1, nfrag
                     plnew%xh(:,i) = xb_frag(:, i) - cb%xb(:)
                     plnew%vh(:,i) = vb_frag(:, i) - cb%vb(:)
                  end do
                  plnew%mass(1:nfrag) = m_frag(1:nfrag)
                  plnew%Gmass(1:nfrag) = param%GU * m_frag(1:nfrag)
                  plnew%radius(1:nfrag) = rad_frag(1:nfrag)
                  plnew%density(1:nfrag) = m_frag(1:nfrag) / rad_frag(1:nfrag)
   
                  select type(info => plnew%info)
                  class is (symba_particle_info)
                     select case(status)
                     case(DISRUPTION)
                        info(1:nfrag)%origin_type = "Disruption"
                        plnew%status(1:nfrag) = NEW_PARTICLE
                        info(1:nfrag)%origin_time = param%t
                        do i = 1, nfrag
                           info(i)%origin_xh(:) = plnew%xh(:,i)
                           info(i)%origin_vh(:) = plnew%vh(:,i)
                        end do
                     case(SUPERCATASTROPHIC)
                        info(1:nfrag)%origin_type = "Supercatastrophic"
                        plnew%status(1:nfrag) = NEW_PARTICLE
                        info(1:nfrag)%origin_time = param%t
                        do i = 1, nfrag
                           info(i)%origin_xh(:) = plnew%xh(:,i)
                           info(i)%origin_vh(:) = plnew%vh(:,i)
                        end do
                     case(HIT_AND_RUN_DISRUPT)
                        select type(plinfo => pl%info)
                        class is (symba_particle_info)
                           info(1)%name = plinfo(ibiggest)%name
                           info(1)%origin_xh(:) = plinfo(ibiggest)%origin_xh(:)
                           info(1)%origin_vh(:) = plinfo(ibiggest)%origin_vh(:)
                        end select
                        plnew%status(1) = OLD_PARTICLE
                        plnew%status(2:nfrag) = NEW_PARTICLE
                        info(2:nfrag)%origin_type = "Hit and run fragment"
                        info(2:nfrag)%origin_time = param%t
                        do i = 2, nfrag
                           info(i)%origin_xh(:) = plnew%xh(:,i)
                           info(i)%origin_vh(:) = plnew%vh(:,i)
                        end do
                     case(MERGED)
                        select type(plinfo => pl%info)
                        class is (symba_particle_info)
                           info(1)%name = plinfo(ibiggest)%name
                           info(1)%origin_xh(:) = plinfo(ibiggest)%origin_xh(:)
                           info(1)%origin_vh(:) = plinfo(ibiggest)%origin_vh(:)
                        end select
                        plnew%status(1) = OLD_PARTICLE
                     end select
                  end select
      
                  if (param%lrotation) then
                     plnew%Ip(:, 1:nfrag) = Ip_frag(:, 1:nfrag)
                     plnew%rot(:, 1:nfrag) = rot_frag(:, 1:nfrag)
                  end if
      
                  if (param%ltides) then
                     plnew%Q = pl%Q(ibiggest)
                     plnew%k2 = pl%k2(ibiggest)
                     plnew%tlag = pl%tlag(ibiggest)
                  end if

                  call plnew%set_mu(cb)
                  !Copy over or set integration parameters for new bodies
                  plnew%lcollision(1:nfrag) = .false.
                  plnew%ldiscard(1:nfrag) = .false.
                  plnew%levelg(1:nfrag) = pl%levelg(ibiggest)
                  plnew%levelm(1:nfrag) = pl%levelm(ibiggest)
      
                  ! Append the new merged body to the list and record how many we made
                  nstart = pl_adds%nbody + 1
                  nend = pl_adds%nbody + plnew%nbody
                  call pl_adds%append(plnew, lsource_mask=[(.true., i=1, plnew%nbody)])
                  pl_adds%ncomp(nstart:nend) = plnew%nbody
      
                  call plnew%setup(0, param)
                  deallocate(plnew)
               end associate
            end select
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
      integer(I4B), dimension(:),     allocatable :: family           !! List of indices of all bodies inovlved in the collision
      integer(I4B), dimension(2)                  :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      real(DP),     dimension(NDIM,2)             :: x, v, L_spin, Ip !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),     dimension(2)                  :: mass, radius     !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      logical                                     :: lgoodcollision
      integer(I4B)                                :: i, jtarg, jproj, regime
      real(DP), dimension(2)                      :: radius_si, mass_si, density_si
      real(DP)                                    :: min_mfrag_si, Mcb_si
      real(DP), dimension(NDIM)                   :: x1_si, v1_si, x2_si, v2_si
      real(DP)                                    :: mlr, mslr, mtot, dentot, msys, msys_new, Qloss, impact_parameter
      integer(I4B), parameter                     :: NRES = 3   !! Number of collisional product results
      real(DP), dimension(NRES)                   :: mass_res   

      associate(plplcollision_list => self, ncollisions => self%nenc, idx1 => self%index1, idx2 => self%index2)
         select type(pl => system%pl)
         class is (symba_pl)
            select type (cb => system%cb)
            class is (symba_cb)
               do i = 1, ncollisions
                  idx_parent(1) = pl%kin(idx1(i))%parent
                  idx_parent(2) = pl%kin(idx2(i))%parent
                  lgoodcollision = symba_collision_consolidate_familes(pl, cb, param, idx_parent, family, x, v, mass, radius, L_spin, Ip)
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
                  x1_si(:)      = plplcollision_list%x1(:,i) * param%DU2M                 !! The position of the parent from inside the step (at collision)
                  v1_si(:)      = plplcollision_list%v1(:,i) * param%DU2M / param%TU2S    !! The velocity of the parent from inside the step (at collision)
                  x2_si(:)      = plplcollision_list%x2(:,i) * param%DU2M                 !! The position of the parent from inside the step (at collision)
                  v2_si(:)      = plplcollision_list%v2(:,i) * param%DU2M / param%TU2S    !! The velocity of the parent from inside the step (at collision)
                  density_si(:) = mass_si(:) / (4.0_DP / 3._DP * PI * radius_si(:)**3) !! The collective density of the parent and its children
                  Mcb_si        = cb%mass * param%MU2KG 
                  min_mfrag_si  = (param%min_GMfrag / param%GU) * param%MU2KG
               
                  mass_res(:) = 0.0_DP
            
                  mtot = sum(mass_si(:)) 
                  dentot = sum(mass_si(:) * density_si(:)) / mtot 
   
                  !! Use the positions and velocities of the parents from indside the step (at collision) to calculate the collisional regime
                  call fragmentation_regime(Mcb_si, mass_si(jtarg), mass_si(jproj), radius_si(jtarg), radius_si(jproj), x1_si(:), x2_si(:),& 
                        v1_si(:), v2_si(:), density_si(jtarg), density_si(jproj), regime, mlr, mslr, min_mfrag_si, Qloss)
   
                  mass_res(1) = min(max(mlr, 0.0_DP), mtot)
                  mass_res(2) = min(max(mslr, 0.0_DP), mtot)
                  mass_res(3) = min(max(mtot - mlr - mslr, 0.0_DP), mtot)
                  mass_res(:) = (mass_res(:) / param%MU2KG) 
                  Qloss = Qloss * (param%TU2S / param%DU2M)**2 / param%MU2KG
   
                  select case (regime)
                  case (COLLRESOLVE_REGIME_DISRUPTION)
                     plplcollision_list%status(i) = symba_collision_casedisruption(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)
                  case (COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
                     plplcollision_list%status(i) = symba_collision_casesupercatastrophic(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)
                  case (COLLRESOLVE_REGIME_HIT_AND_RUN)
                     plplcollision_list%status(i) = symba_collision_casehitandrun(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)
                  case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
                     plplcollision_list%status(i) = symba_collision_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip) 
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
      integer(I4B), dimension(:),     allocatable :: family           !! List of indices of all bodies inovlved in the collision
      integer(I4B), dimension(2)                  :: idx_parent       !! Index of the two bodies considered the "parents" of the collision
      real(DP),     dimension(NDIM,2)             :: x, v, L_spin, Ip !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),     dimension(2)                  :: mass, radius     !! Output values that represent a 2-body equivalent of a possibly 2+ body collision
      logical                                     :: lgoodcollision
      integer(I4B)                                :: i

      associate(plplcollision_list => self, ncollisions => self%nenc, idx1 => self%index1, idx2 => self%index2)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(cb => system%cb)
            class is (symba_cb)
               do i = 1, ncollisions
                  idx_parent(1) = pl%kin(idx1(i))%parent
                  idx_parent(2) = pl%kin(idx2(i))%parent
                  lgoodcollision = symba_collision_consolidate_familes(pl, cb, param, idx_parent, family, x, v, mass, radius, L_spin, Ip)
                  if (.not. lgoodcollision) cycle
                  if (any(pl%status(idx_parent(:)) /= COLLISION)) cycle ! One of these two bodies has already been resolved
   
                  plplcollision_list%status(i) = symba_collision_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip) 
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
                  write(*, *) "Collision between massive bodies detected at time t = ", t
                  if (param%lfragmentation) then
                     call plplcollision_list%resolve_fragmentations(system, param)
                  else
                     call plplcollision_list%resolve_mergers(system, param)
                  end if

                  ! Destroy the collision list now that the collisions are resolved
                  call plplcollision_list%setup(0)

                  if ((system%pl_adds%nbody == 0) .and. (system%pl_discards%nbody == 0)) exit

                  ! Save the add/discard information to file
                  call system%write_discard(param)

                  ! Rearrange the arrays: Remove discarded bodies, add any new bodies, resort, and recompute all indices and encounter lists
                  call pl%rearray(system, param)

                  ! Destroy the add/discard list so that we don't append the same body multiple times if another collision is detected
                  call system%pl_discards%setup(0, param)
                  call system%pl_adds%setup(0, param)

                  ! Check whether or not any of the particles that were just added are themselves in a collision state. This will generate a new plplcollision_list 
                  lplpl_collision = plplenc_list%collision_check(system, param, t, dt, irec)

                  if (.not.lplpl_collision) exit
               end do

               if (param%lenergy) then
                  call system%get_energy_and_momentum(param)
                  Eorbit_after = system%te
                  param%Ecollisions = param%Ecollisions + (Eorbit_after - Eorbit_before)
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
      
      call system%tp%xh2xb(system%cb)
      call system%tp%discard(system, param)

      return
   end subroutine symba_collision_resolve_pltpenc

end submodule s_symba_collision