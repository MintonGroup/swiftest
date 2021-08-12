submodule (symba_classes) s_symba_discard
   use swiftest
contains

   subroutine symba_discard_cb_pl(pl, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if planets should be discarded based on their positions relative to the central body.
      !! If a body gets flagged here when it has also been previously flagged for a collision with another massive body,
      !! its collisional status will be revoked. Discards due to colliding with or escaping the central body take precedence 
      !! over pl-pl collisions
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_sun.f90
      !! Adapted from Hal Levison's Swift routine discard_massive5.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i, j
      real(DP)     :: energy, vb2, rb2, rh2, rmin2, rmax2, rmaxu2
   
      associate(npl => pl%nbody, cb => system%cb)
         call system%set_msys()
         rmin2 = param%rmin**2
         rmax2 = param%rmax**2
         rmaxu2 = param%rmaxu**2
         do i = 1, npl
            if (pl%status(i) == ACTIVE) then
               rh2 = dot_product(pl%xh(:,i), pl%xh(:,i))
               if ((param%rmax >= 0.0_DP) .and. (rh2 > rmax2)) then
                  pl%ldiscard(i) = .true.
                  pl%lcollision(i) = .false. 
                  pl%status(i) = DISCARDED_RMAX
                  write(*, *) "Massive body ",  pl%id(i), " too far from the central body at t = ", param%t
               else if ((param%rmin >= 0.0_DP) .and. (rh2 < rmin2)) then
                  pl%ldiscard(i) = .true.
                  pl%lcollision(i) = .false. 
                  pl%status(i) = DISCARDED_RMIN
                  write(*, *) "Massive body ", pl%id(i), " too close to the central body at t = ", param%t
               else if (param%rmaxu >= 0.0_DP) then
                  rb2 = dot_product(pl%xb(:,i), pl%xb(:,i))
                  vb2 = dot_product(pl%vb(:,i), pl%vb(:,i))
                  energy = 0.5_DP * vb2 - system%Gmtot / sqrt(rb2)
                  if ((energy > 0.0_DP) .and. (rb2 > rmaxu2)) then
                     pl%ldiscard(i) = .true.
                     pl%lcollision(i) = .false. 
                     pl%status(i) = DISCARDED_RMAXU
                     write(*, *) "Massive body ", pl%id(i), " is unbound and too far from barycenter at t = ", param%t
                  end if
               end if
            end if
         end do
      end associate
   
      return
   end subroutine symba_discard_cb_pl


   subroutine symba_discard_conserve_mtm(pl, system, param, ipl, lescape_body)
      !! author: David A. Minton
      !! 
      !! Conserves system momentum when a body is lost from the system or collides with central body
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: pl
      class(symba_nbody_system), intent(inout) :: system
      class(symba_parameters),   intent(inout) :: param
      integer(I4B),              intent(in)    :: ipl
      logical,                   intent(in)         :: lescape_body
      ! Internals
      real(DP), dimension(NDIM) :: Lpl, Ltot, Lcb, xcom, vcom
      real(DP)                  :: pe, ke_orbit, ke_spin
      integer(I4B)              :: i, oldstat
   
      select type(cb => system%cb)
      class is (symba_cb)
   
         ! Add the potential and kinetic energy of the lost body to the records
         pe = -cb%mass * pl%mass(ipl) / norm2(pl%xb(:, ipl) - cb%xb(:))
         ke_orbit = 0.5_DP * pl%mass(ipl) * dot_product(pl%vb(:, ipl), pl%vb(:, ipl)) 
         if (param%lrotation) then
            ke_spin  = 0.5_DP * pl%mass(ipl) * pl%radius(ipl)**2 * pl%Ip(3, ipl) * dot_product(pl%rot(:, ipl), pl%rot(:, ipl))
         else
            ke_spin = 0.0_DP
         end if
   
         ! Add the pre-collision ke of the central body to the records
         ! Add planet mass to central body accumulator
         if (lescape_body) then
            system%Mescape = system%Mescape + pl%mass(ipl)
            do i = 1, pl%nbody
               if (i == ipl) cycle
               pe = pe - pl%mass(i) * pl%mass(ipl) / norm2(pl%xb(:, ipl) - pl%xb(:, i))
            end do
   
            Ltot(:) = 0.0_DP
            do i = 1, pl%nbody
               Lpl(:) = pL%mass(i) * pl%xb(:,i) .cross. pl%vb(:, i)
               Ltot(:) = Ltot(:) + Lpl(:)
            end do
            Ltot(:) = Ltot(:) + cb%mass * cb%xb(:) .cross. cb%vb(:)
            call pl%b2h(cb)
            oldstat = pl%status(ipl)
            pl%status(ipl) = INACTIVE
            call pl%h2b(cb)
            pl%status(ipl) = oldstat
            do i = 1, pl%nbody
               if (i == ipl) cycle
               Lpl(:) = pl%mass(i) * pl%xb(:,i) .cross. pl%vb(:, i)
               Ltot(:) = Ltot(:) - Lpl(:) 
            end do 
            Ltot(:) = Ltot(:) - cb%mass * cb%xb(:) .cross. cb%vb(:)
            system%Lescape(:) = system%Lescape(:) + Ltot(:)
            if (param%lrotation) system%Lescape(:) = system%Lescape + pl%mass(ipl) * pl%radius(ipl)**2 * pl%Ip(3, ipl) * pl%rot(:, ipl)
   
         else
            xcom(:) = (pl%mass(ipl) * pl%xb(:, ipl) + cb%mass * cb%xb(:)) / (cb%mass + pl%mass(ipl))
            vcom(:) = (pl%mass(ipl) * pl%vb(:, ipl) + cb%mass * cb%vb(:)) / (cb%mass + pl%mass(ipl))
            Lpl(:) = (pl%xb(:,ipl) - xcom(:)) .cross. pL%vb(:,ipl) - vcom(:)
            if (param%lrotation) Lpl(:) = pl%mass(ipl) * (Lpl(:) + pl%radius(ipl)**2 * pl%Ip(3,ipl) * pl%rot(:, ipl))
     
            Lcb(:) = cb%mass * (cb%xb(:) - xcom(:)) .cross. (cb%vb(:) - vcom(:))
   
            ke_orbit = ke_orbit + 0.5_DP * cb%mass * dot_product(cb%vb(:), cb%vb(:)) 
            if (param%lrotation) ke_spin = ke_spin + 0.5_DP * cb%mass * cb%radius**2 * cb%Ip(3) * dot_product(cb%rot(:), cb%rot(:))
            ! Update mass of central body to be consistent with its total mass
            cb%dM = cb%dM + pl%mass(ipl)
            cb%dR = cb%dR + 1.0_DP / 3.0_DP * (pl%radius(ipl) / cb%radius)**3 - 2.0_DP / 9.0_DP * (pl%radius(ipl) / cb%radius)**6
            cb%mass = cb%M0 + cb%dM
            cb%Gmass = param%GU * cb%mass 
            cb%radius = cb%R0 + cb%dR
            param%rmin = cb%radius
            ! Add planet angular momentum to central body accumulator
            cb%dL(:) = Lpl(:) + Lcb(:) + cb%dL(:)
            ! Update rotation of central body to by consistent with its angular momentum 
            if (param%lrotation) then
               cb%rot(:) = (cb%L0(:) + cb%dL(:)) / (cb%Ip(3) * cb%mass * cb%radius**2)        
               ke_spin  = ke_spin - 0.5_DP * cb%mass * cb%radius**2 * cb%Ip(3) * dot_product(cb%rot(:), cb%rot(:)) 
            end if
            cb%xb(:) = xcom(:)
            cb%vb(:) = vcom(:)
            ke_orbit = ke_orbit - 0.5_DP * cb%mass * dot_product(cb%vb(:), cb%vb(:)) 
         end if
         call pl%b2h(cb)
   
         ! We must do this for proper book-keeping, since we can no longer track this body's contribution to energy directly
         if (lescape_body) then
            system%Ecollisions  = system%Ecollisions + ke_orbit + ke_spin + pe
            system%Euntracked  = system%Euntracked - (ke_orbit + ke_spin + pe)
         else
            system%Ecollisions  = system%Ecollisions + pe 
            system%Euntracked = system%Euntracked - pe
         end if
   
      end select
      return
   end subroutine symba_discard_conserve_mtm


   subroutine symba_discard_nonplpl(pl, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if planets should be discarded based on their positions or because they are unbound
      !s
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_massive5.f 
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
    
      ! First check for collisions with the central body
      associate(npl => pl%nbody, cb => system%cb)
         if (npl == 0) return 
         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or. &
             (param%rmaxu >= 0.0_DP) .or. ((param%qmin >= 0.0_DP) .and. (param%qmin_coord == "BARY"))) then
            call pl%h2b(cb) 
         end if
         if ((param%rmin >= 0.0_DP) .or. (param%rmax >= 0.0_DP) .or.  (param%rmaxu >= 0.0_DP)) then
            call symba_discard_cb_pl(pl, system, param)
         end if
         if (param%qmin >= 0.0_DP .and. npl > 0) call symba_discard_peri_pl(pl, system, param)
      end associate

      return
   end subroutine symba_discard_nonplpl


   subroutine symba_discard_nonplpl_conservation(pl, system, param)
      !! author: David A. Minton
      !!
      !! If there are any bodies that are removed due to either colliding with the central body or escaping the systme,
      !! we need to track the conserved quantities with the system bookkeeping terms.
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: pl     !! SyMBA test particle object
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      class(symba_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B)                            :: i, ndiscard, dstat
      logical                                 :: lescape
      logical, dimension(pl%nbody)            :: discard_l_pl
      integer(I4B), dimension(:), allocatable :: discard_index_list

      associate(npl => pl%nbody)
         discard_l_pl(1:npl) = pl%ldiscard(1:npl) .and. .not. pl%lcollision(1:npl) ! These are bodies that are discarded but not flagged as pl-pl collision
         ndiscard = count(discard_l_pl(:)) 
         allocate(discard_index_list(ndiscard))
         discard_index_list(:) = pack([(i, i = 1, npl)], discard_l_pl(1:npl))
         do i = 1, ndiscard
            dstat = pl%status(discard_index_list(i)) 
            if ((dstat == DISCARDED_RMIN) .or. (dstat == DISCARDED_PERI)) then
               lescape = .false.
            else if ((dstat == DISCARDED_RMAX) .or. (dstat == DISCARDED_RMAXU)) then
               lescape = .true.
            else 
               cycle
            end if
            ! Conserve all the quantities
            call symba_discard_conserve_mtm(pl, system, param, discard_index_list(i), lescape)
         end do
      end associate

      return
   end subroutine symba_discard_nonplpl_conservation


   subroutine symba_discard_peri_pl(pl, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if a test particle should be discarded because its perihelion distance becomes too small
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_peri_pl.f90
      !! Adapted from Hal Levison's Swift routine discard_mass_peri.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: pl     !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      logical, save      :: lfirst = .true.
      logical            :: lfirst_orig
      integer(I4B)       :: i


      lfirst_orig = pl%lfirst
      pl%lfirst = lfirst
      if (lfirst) then
         call pl%get_peri(system, param)          
         lfirst = .false.
      else
         call pl%get_peri(system, param)          
         do i = 1, pl%nbody
            if (pl%status(i) == ACTIVE) then
               if ((pl%isperi(i) == 0) .and. (pl%nplenc(i)== 0)) then
                  if ((pl%atp(i) >= param%qmin_alo) .and. (pl%atp(i) <= param%qmin_ahi) .and. (pl%peri(i) <= param%qmin)) then
                     pl%ldiscard(i) = .true.
                     pl%lcollision(i) = .false.
                     pl%status(i) = DISCARDED_PERI
                     write(*, *) "Particle ", pl%id(i), " perihelion distance too small at t = ", param%t
                  end if
               end if
            end if
         end do
      end if
      pl%lfirst = lfirst_orig
   
      return
   end subroutine symba_discard_peri_pl


   module subroutine symba_discard_pl(self, system, param)
      !! author: David A. Minton
      !!
      !! Call the various flavors of discards for massive bodies in SyMBA runs, including discards due to colling with the central body, 
      !! escaping the system, or colliding with each other.
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      real(DP) :: Eorbit_before, Eorbit_after
   
      select type(system)
      class is (symba_nbody_system)
         select type(param)
         class is (symba_parameters)
            associate(pl => self, plplenc_list => system%plplenc_list, plplcollision_list => system%plplcollision_list)
               call pl%h2b(system%cb) 

               ! First deal with the non pl-pl collisions
               call symba_discard_nonplpl(self, system, param)

               ! Extract the pl-pl encounter list and return the plplcollision_list
               call plplenc_list%extract_collisions(system, param)
               call plplenc_list%write(pl, pl, param)

               if ((plplcollision_list%nenc > 0) .and. any(pl%lcollision(:))) then 
                  write(*, *) "Collision between massive bodies detected at time t = ",param%t
                  if (param%lfragmentation) then
                     call plplcollision_list%resolve_fragmentations(system, param)
                  else
                     call plplcollision_list%resolve_mergers(system, param)
                  end if
                  ! Destroy the collision list now that the collisions are resolved
                  call plplcollision_list%setup(0)
               end if

               if (any(pl%ldiscard(:))) then
                  if (param%lenergy) then
                     call system%get_energy_and_momentum(param)
                     Eorbit_before = system%te
                  end if
                  call symba_discard_nonplpl_conservation(self, system, param)
                  call pl%rearray(system, param)
                  if (param%lenergy) then
                     call system%get_energy_and_momentum(param)
                     Eorbit_after = system%te
                     system%Ecollisions = system%Ecollisions + (Eorbit_after - Eorbit_before)
                  end if
               end if

            end associate
         end select 
      end select

      return
   end subroutine symba_discard_pl

end submodule s_symba_discard