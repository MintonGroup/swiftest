submodule (symba_classes) s_symba_fragmentation
   use swiftest
contains

   module function symba_fragmentation_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip)  result(status)
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
      class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
      integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
      real(DP),        dimension(:,:), intent(in)    :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      real(DP),        dimension(:),   intent(in)    :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
      ! Result
      integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                              :: i, j, ibiggest, nfamily, nstart, nend
      real(DP)                                  :: mass_new, radius_new, volume_new, pe
      real(DP), dimension(NDIM)                 :: xcom, vcom, xc, vc, xcrossv
      real(DP), dimension(2)                    :: vol
      real(DP), dimension(NDIM)                 :: L_orb_old, L_spin_old
      real(DP), dimension(NDIM)                 :: L_spin_new, rot_new, Ip_new
      logical,  dimension(system%pl%nbody)      :: lmask
      class(symba_pl), allocatable              :: plnew
   
      select type(pl => system%pl)
      class is (symba_pl)
         associate(mergeadd_list => system%mergeadd_list, mergesub_list => system%mergesub_list, cb => system%cb)
            status = MERGED
            write(*, '("Merging bodies ",99(I8,",",:))') pl%id(family(:))
            mass_new = sum(mass(:))
      
            ! Merged body is created at the barycenter of the original bodies
            xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mass_new
            vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mass_new
      
            ! Get mass weighted mean of Ip and 
            vol(:) = 4._DP / 3._DP * PI * radius(:)**3
            volume_new = sum(vol(:))
            radius_new = (3 * volume_new / (4 * PI))**(1._DP / 3._DP)

            L_orb_old(:) = 0.0_DP

            ! Compute orbital angular momentum of pre-impact system
            do i = 1, 2
               xc(:) = x(:, i) - xcom(:)
               vc(:) = v(:, i) - vcom(:)
               xcrossv(:) = xc(:) .cross. vc(:)
               L_orb_old(:) = L_orb_old(:) + mass(i) * xcrossv(:)
            end do
         
            if (param%lrotation) then
               Ip_new(:) = (mass(1) * Ip(:,1) + mass(2) * Ip(:,2)) / mass_new
               L_spin_old(:) = L_spin(:,1) + L_spin(:,2)

               ! Conserve angular momentum by putting pre-impact orbital momentum into spin of the new body
               L_spin_new(:) = L_orb_old(:) + L_spin_old(:) 
      
               ! Assume prinicpal axis rotation on 3rd Ip axis
               rot_new(:) = L_spin_new(:) / (Ip_new(3) * mass_new * radius_new**2)
            else ! If spin is not enabled, we will consider the lost pre-collision angular momentum as "escaped" and add it to our bookkeeping variable
               system%Lescape(:) = system%Lescape(:) + L_orb_old(:) 
            end if
      
            ! Keep track of the component of potential energy due to the pre-impact family for book-keeping
            nfamily = size(family(:))
            pe = 0.0_DP
            do j = 1, nfamily
               do i = j + 1, nfamily
                  pe = pe - pl%mass(i) * pl%mass(j) / norm2(pl%xb(:, i) - pl%xb(:, j))
               end do
            end do
            system%Ecollisions  = system%Ecollisions + pe 
            system%Euntracked = system%Euntracked - pe 
     
            ! Add the family bodies to the subtraction list
            lmask(:) = .false.
            lmask(family(:)) = .true.
            pl%status(family(:)) = MERGED
            nstart = mergesub_list%nbody + 1
            nend = mergesub_list%nbody + nfamily
            call mergesub_list%append(pl, lmask)
            ! Record how many bodies were subtracted in this event
            mergesub_list%ncomp(nstart:nend) = nfamily 

            ! Create the new merged body 
            allocate(plnew, mold=pl)
            call plnew%setup(1, param)

            ! The merged body's name will be that of the largest of the two parents 
            ibiggest = maxloc(pl%Gmass(family(:)), dim=1)
            plnew%id(1) = pl%id(family(ibiggest))
            plnew%status(1) = ACTIVE
            plnew%xb(:,1) = xcom(:)
            plnew%vb(:,1) = vcom(:)
            plnew%xh(:,1) = xcom(:) - cb%xb(:)
            plnew%vh(:,1) = vcom(:) - cb%vb(:)
            plnew%mass(1) = mass_new
            plnew%Gmass(1) = param%GU * mass_new
            plnew%density(1) = mass_new / volume_new
            plnew%radius(1) = radius_new
            plnew%info(1) = pl%info(family(ibiggest)) 
            if (param%lrotation) then
               pl%Ip(:,1) = Ip_new(:)
               pl%rot(:,1) = rot_new(:)
            end if
            if (param%ltides) then
               plnew%Q = pl%Q(ibiggest)
               plnew%k2 = pl%k2(ibiggest)
               plnew%tlag = pl%tlag(ibiggest)
            end if

            ! Append the new merged body to the list and record how many we made
            nstart = mergeadd_list%nbody + 1
            nend = mergeadd_list%nbody + plnew%nbody
            call mergeadd_list%append(plnew)
            mergeadd_list%ncomp(nstart:nend) = plnew%nbody

            call plnew%setup(0, param)
            deallocate(plnew)

         end associate
      end select
   
      return 

   end function symba_fragmentation_casemerge

end submodule s_symba_fragmentation
