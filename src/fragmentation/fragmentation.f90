submodule(swiftest_classes) s_fragmentation
   use swiftest
contains

   module subroutine fragmentation_initialize(system, param, family, x, v, L_spin, Ip, mass, radius, &
                           nfrag, Ip_frag, m_frag, rad_frag, xb_frag, vb_frag, rot_frag, Qloss, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initialize the position and velocity of fragments to conserve energy and momentum.
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(swiftest_nbody_system),                              intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),                                intent(in)    :: param  !! Current run configuration parameters 
      integer(I4B),                 dimension(:),                intent(in)    :: family !! Index of bodies involved in the collision
      real(DP),                     dimension(:,:),              intent(inout) :: x, v, L_spin, Ip !! Two-body equivalent position, vector, spin momentum, and rotational inertia values for the collision
      real(DP),                     dimension(:),                intent(inout) :: mass, radius     !! Two-body equivalent mass and radii for the bodies in the collision
      integer(I4B),                                              intent(inout) :: nfrag            !! Number of fragments to generate
      real(DP),                     dimension(:),   allocatable, intent(inout) :: m_frag, rad_frag !! Distribution of fragment mass and radii
      real(DP),                     dimension(:,:), allocatable, intent(inout) :: Ip_frag          !! Fragment rotational inertia vectors
      real(DP),                     dimension(:,:), allocatable, intent(inout) :: xb_frag, vb_frag, rot_frag !! Fragment barycentric position, barycentric velocity, and rotation vectors
      real(DP),                                                  intent(inout) :: Qloss !! Energy lost during the collision
      logical,                                                   intent(out)   :: lfailure !! Answers the question: Should this have been a merger instead?
      ! Internals
      real(DP)                                :: mscale, dscale, vscale, tscale, Lscale, Escale ! Scale factors that reduce quantities to O(~1) in the collisional system
      real(DP)                                :: mtot 
      real(DP), dimension(NDIM)               :: xcom, vcom
      integer(I4B)                            :: ii, npl_new
      logical, dimension(:), allocatable      :: lexclude
      real(DP), dimension(NDIM, 2)            :: rot, L_orb, mxc, vc
      real(DP), dimension(:,:), allocatable   :: x_frag, v_frag, v_r_unit, v_t_unit, v_h_unit
      real(DP), dimension(:), allocatable     :: rmag, rotmag, v_r_mag, v_t_mag
      real(DP), dimension(NDIM)               :: Ltot_before
      real(DP), dimension(NDIM)               :: Ltot_after
      real(DP)                                :: Etot_before, ke_orbit_before, ke_spin_before, pe_before, Lmag_before
      real(DP)                                :: Etot_after,  ke_orbit_after,  ke_spin_after,  pe_after,  Lmag_after, dEtot, dLmag
      real(DP), dimension(NDIM)               :: L_frag_tot, L_frag_orb, L_frag_spin, L_frag_budget, Lorbit_before, Lorbit_after, Lspin_before, Lspin_after, dL
      real(DP)                                :: ke_frag_budget, ke_frag_orbit, ke_frag_spin, ke_tot_deficit, ke_avg_deficit, ke_avg_deficit_old
      real(DP), dimension(NDIM)               :: x_col_unit, y_col_unit, z_col_unit
      character(len=*), parameter             :: fmtlabel = "(A14,10(ES11.4,1X,:))"
      integer(I4B)                            :: try
      integer(I4B), parameter                 :: NFRAG_MIN = 7 !! The minimum allowable number of fragments (set to 6 because that's how many unknowns are needed in the tangential velocity calculation)
      real(DP)                                :: r_max_start, r_max_start_old, r_max, f_spin 
      real(DP), parameter                     :: Ltol = 10 * epsilon(1.0_DP)
      real(DP), parameter                     :: Etol = 1e-8_DP
      integer(I4B), parameter                 :: MAXTRY = 3000
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes, fpe_quiet_modes

      if (nfrag < NFRAG_MIN) then
         write(*,*) "symba_frag_pos needs at least ",NFRAG_MIN," fragments, but only ",nfrag," were given."
         lfailure = .true.
         return
      end if

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      fpe_quiet_modes(:) = .false.
      call ieee_set_halting_mode(IEEE_ALL,fpe_quiet_modes)

      allocate(x_frag, source=xb_frag)
      allocate(v_frag, source=vb_frag)
      allocate(rmag(nfrag))
      allocate(rotmag(nfrag))
      allocate(v_r_mag(nfrag))
      allocate(v_t_mag(nfrag))
      allocate(v_r_unit(NDIM,nfrag))
      allocate(v_t_unit(NDIM,nfrag))
      allocate(v_h_unit(NDIM,nfrag))

      rmag(:) = 0.0_DP
      rotmag(:) = 0.0_DP
      v_r_mag(:) = 0.0_DP
      v_t_mag(:) = 0.0_DP
      v_r_unit(:,:) = 0.0_DP
      v_t_unit(:,:) = 0.0_DP
      v_h_unit(:,:) = 0.0_DP

      associate(pl => system%pl, npl => system%pl%nbody)
         npl_new = npl + nfrag
         allocate(lexclude(npl_new))
         lexclude(1:npl) = pl%status(1:npl) == INACTIVE
         lexclude(npl+1:npl_new) = .true.
      end associate

      call set_scale_factors()

      ! Compute orbital angular momentum of pre-impact system
      mxc(:, 1) = mass(1) * (x(:, 1) - xcom(:))
      mxc(:, 2) = mass(2) * (x(:, 2) - xcom(:))
      vc(:, 1) = v(:, 1) - vcom(:)
      vc(:, 2) = v(:, 2) - vcom(:)
      L_orb(:,:) = mxc(:,:) .cross. vc(:,:)

      ! Compute orbital angular momentum of pre-impact system. We'll use this to start the coordinate system, but it will get updated as we divide up the angular momentum
      L_frag_orb(:) = L_orb(:, 1) + L_orb(:, 2)
      L_frag_spin(:) = L_spin(:, 1) + L_spin(:, 2)
      L_frag_budget(:) = L_frag_orb(:) + L_frag_spin(:)
      f_spin = 0.05_DP

      call reset_fragments()
      call define_coordinate_system()

      ! Calculate the initial energy of the system without the collisional family
      call calculate_system_energy(linclude_fragments=.false.)
      
      r_max_start = 1 * norm2(x(:,2) - x(:,1))
      try = 1
      lfailure = .false.
      ke_tot_deficit = 0.0_DP
      do while (try < MAXTRY)
         lfailure = .false.

         call set_fragment_position_vectors()

         do concurrent (ii = 1:nfrag)
            vb_frag(:, ii) = vcom(:)
         end do

         call calculate_system_energy(linclude_fragments=.true.)
         L_frag_budget(:) = -dL(:)
         ! The ke constraints are calcualted in the collision frame, so undo the barycentric velocity component
         ke_frag_budget = -(dEtot - 0.5_DP * mtot * dot_product(vcom(:), vcom(:))) - Qloss 

         call set_fragment_spin(lfailure)
         if (.not.lfailure) call set_fragment_tan_vel(lfailure)
            
         if (lfailure) then
            ! write(*,*) 'Failed to find tangential velocities'
         else 
            call set_fragment_radial_velocities(lfailure)
            ! if (lfailure) write(*,*) 'Failed to find radial velocities'
            if (.not.lfailure) then
               call calculate_system_energy(linclude_fragments=.true.)
               if ((abs(dEtot + Qloss) > Etol) .or. (dEtot > 0.0_DP)) then
                  ! write(*,*) 'Failed due to high energy error: ',dEtot, abs(dEtot + Qloss) / Etol
                  lfailure = .true.
               else if (abs(dLmag) / Lmag_before > Ltol) then
                  ! write(*,*) 'Failed due to high angular momentum error: ', dLmag / Lmag_before
                  lfailure = .true.
               end if
            end if
         end if
   
         if (.not.lfailure) exit
         call restructure_failed_fragments()
         call reset_fragments()
         try = try + 1
      end do
      call restore_scale_factors()

      ! write(*,        "(' -------------------------------------------------------------------------------------')")
      ! write(*,        "('  Final diagnostic')")
      ! write(*,        "(' -------------------------------------------------------------------------------------')")
      ! call calculate_system_energy(linclude_fragments=.true.)
      if (lfailure) then
         ! write(*,*) "symba_frag_pos failed after: ",try," tries"
         do ii = 1, nfrag
            vb_frag(:, ii) = vcom(:)
         end do
      ! else
      !    write(*,*) "symba_frag_pos succeeded after: ",try," tries"
      !    write(*,        "(' dL_tot should be very small' )")
      !    write(*,fmtlabel) ' dL_tot      |', dLmag / Lmag_before
      !    write(*,        "(' dE_tot should be negative and equal to Qloss' )")
      !    write(*,fmtlabel) ' dE_tot      |', dEtot / abs(Etot_before)
      !    write(*,fmtlabel) ' Qloss       |', -Qloss / abs(Etot_before)
      !    write(*,fmtlabel) ' dE - Qloss  |', (Etot_after - Etot_before + Qloss) / abs(Etot_before)
      end if
      ! write(*,        "(' -------------------------------------------------------------------------------------')")

      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily

      return 

      contains

         ! Because of the complexity of this procedure, we have chosen to break it up into a series of nested subroutines.
         subroutine set_scale_factors()
            !! author: David A. Minton
            !!
            !! Scales dimenional quantities to ~O(1) with respect to the collisional system. This scaling makes it easier for the non-linear minimization 
            !! to converge on a solution
            implicit none
            integer(I4B) :: i

            ! Find the center of mass of the collisional system	
            mtot = sum(mass(:)) 
            xcom(:) = (mass(1) * x(:,1) + mass(2) * x(:,2)) / mtot
            vcom(:) = (mass(1) * v(:,1) + mass(2) * v(:,2)) / mtot

            ! Set scale factors
            Escale = 0.5_DP * (mass(1) * dot_product(v(:,1), v(:,1)) + mass(2) * dot_product(v(:,2), v(:,2)))
            dscale = sum(radius(:))
            mscale = mtot 
            vscale = sqrt(Escale / mscale) 
            tscale = dscale / vscale 
            Lscale = mscale * dscale * vscale

            xcom(:) = xcom(:) / dscale
            vcom(:) = vcom(:) / vscale

            mtot = mtot / mscale
            mass = mass / mscale
            radius = radius / dscale
            x = x / dscale
            v = v / vscale
            L_spin = L_spin / Lscale
            do i = 1, 2
               rot(:,i) = L_spin(:,i) / (mass(i) * radius(i)**2 * Ip(3, i))
            end do

            m_frag = m_frag / mscale
            rad_frag = rad_frag / dscale
            Qloss = Qloss / Escale

            return
         end subroutine set_scale_factors


         subroutine restore_scale_factors()
            !! author: David A. Minton
            !!
            !! Restores dimenional quantities back to the system units
            implicit none
            integer(I4B) :: i

            call ieee_set_halting_mode(IEEE_ALL,.false.)
            ! Restore scale factors
            xcom(:) = xcom(:) * dscale
            vcom(:) = vcom(:) * vscale

            mtot = mtot * mscale
            mass = mass * mscale
            radius = radius * dscale
            x = x * dscale
            v = v * vscale
            L_spin = L_spin * Lscale
            do i = 1, 2
               rot(:,i) = L_spin(:,i) * (mass(i) * radius(i)**2 * Ip(3, i))
            end do

            m_frag = m_frag * mscale
            rad_frag = rad_frag * dscale
            rot_frag = rot_frag / tscale
            x_frag = x_frag * dscale
            v_frag = v_frag * vscale
            Qloss = Qloss * Escale

            do i = 1, nfrag
               xb_frag(:, i) = x_frag(:, i) + xcom(:)
               vb_frag(:, i) = v_frag(:, i) + vcom(:)
            end do

            Etot_before = Etot_before * Escale
            pe_before = pe_before * Escale
            ke_spin_before = ke_spin_before * Escale
            ke_orbit_before = ke_orbit_before * Escale
            Ltot_before = Ltot_before * Lscale
            Lmag_before = Lmag_before * Lscale 
            Etot_after = Etot_after * Escale
            pe_after = pe_after * Escale
            ke_spin_after = ke_spin_after * Escale
            ke_orbit_after = ke_orbit_after * Escale
            Ltot_after = Ltot_after * Lscale
            Lmag_after = Lmag_after * Lscale 

            dL(:) = Ltot_after(:) - Ltot_before(:)
            dLmag = .mag.dL(:)
            dEtot = Etot_after - Etot_before

            mscale = 1.0_DP
            dscale = 1.0_DP
            vscale = 1.0_DP
            tscale = 1.0_DP
            Lscale = 1.0_DP
            Escale = 1.0_DP

            return
         end subroutine restore_scale_factors

         subroutine reset_fragments()
            !! author: David A. Minton
            !!
            !! Resets all tracked fragment quantities in order to do a fresh calculation
            !! Initialize the fragments with 0 velocity and spin so we can divide up the balance between the tangential, radial, and spin components while conserving momentum
            implicit none

            xb_frag(:,:) = 0.0_DP
            vb_frag(:,:) = 0.0_DP
            x_frag(:,:) = 0.0_DP
            v_frag(:,:) = 0.0_DP
            rot_frag(:,:) = 0.0_DP
            v_t_mag(:) = 0.0_DP
            v_r_mag(:) = 0.0_DP
            ke_frag_orbit = 0.0_DP
            ke_frag_spin = 0.0_DP
            L_frag_orb(:) = 0.0_DP
            L_frag_spin(:) = 0.0_DP

            return
         end subroutine reset_fragments


         subroutine define_coordinate_system()
            !! author: David A. Minton
            !!
            !! Defines the collisional coordinate system, including the unit vectors of both the system and individual fragments.
            implicit none
            integer(I4B) :: i
            real(DP), dimension(NDIM) ::  x_cross_v, delta_r, delta_v
            real(DP)   :: r_col_norm, v_col_norm
            real(DP), dimension(NDIM, nfrag) :: L_sigma

            delta_v(:) = v(:, 2) - v(:, 1)
            v_col_norm = .mag. delta_v(:)
            delta_r(:) = x(:, 2) - x(:, 1)
            r_col_norm = .mag. delta_r(:)

            ! We will initialize fragments on a plane defined by the pre-impact system, with the z-axis aligned with the angular momentum vector
            ! and the y-axis aligned with the pre-impact distance vector.
            y_col_unit(:) = delta_r(:) / r_col_norm 
            z_col_unit(:) = (L_frag_budget(:) - L_frag_spin(:)) / (.mag. (L_frag_budget(:) - L_frag_spin(:)))
            ! The cross product of the y- by z-axis will give us the x-axis
            x_col_unit(:) = y_col_unit(:) .cross. z_col_unit(:)

            if (.not.any(x_frag(:,:) > 0.0_DP)) return
            rmag(:) = .mag. x_frag(:,:)

            call random_number(L_sigma(:,:)) ! Randomize the tangential velocity direction. This helps to ensure that the tangential velocity doesn't completely line up with the angular momentum vector,
                                             ! otherwise we can get an ill-conditioned system
            do concurrent(i = 1:nfrag, rmag(i) > 0.0_DP)
               v_r_unit(:, i) = x_frag(:, i) / rmag(i)
               v_h_unit(:, i) = z_col_unit(:) + 2e-1_DP * (L_sigma(:,i) - 0.5_DP)
               v_h_unit(:, i) = v_h_unit(:, i) / (.mag. v_h_unit(:, i))
               v_t_unit(:, i) = v_h_unit(:, i) .cross. v_r_unit(:, i)
               v_t_unit(:, i) = v_t_unit(:, i) / (.mag. v_t_unit(:, i))
            end do

            return
         end subroutine define_coordinate_system


         subroutine construct_temporary_system(tmpsys, tmpparam)
            !! Author: David A. Minton
            !!
            !! Constructs a temporary internal system consisting of active bodies and additional fragments. This internal temporary system is used to calculate system energy with and without fragments
            !! and optionally including fragments.
            implicit none
            ! Arguments
            class(swiftest_nbody_system), allocatable, intent(inout) :: tmpsys
            class(swiftest_parameters), allocatable, intent(inout) :: tmpparam
            ! Internals
            logical, dimension(:), allocatable :: lexclude_tmp

            associate(pl => system%pl, npl => system%pl%nbody, cb => system%cb)
               if (size(lexclude) /= npl + nfrag) then 
                  allocate(lexclude_tmp(npl + nfrag))
                  lexclude_tmp(1:npl) = lexclude(1:npl)
                  call move_alloc(lexclude_tmp, lexclude)
               end if
               where (pl%status(1:npl) == INACTIVE) ! Safety check in case one of the included bodies has been previously deactivated 
                  lexclude(1:npl) = .true.  
               elsewhere
                  lexclude(1:npl) = .false. 
               end where
               lexclude(npl+1:(npl + nfrag)) = .true.
               allocate(tmpparam, source=param)
               call setup_construct_system(tmpsys, param)
               call tmpsys%tp%setup(0, param)
               deallocate(tmpsys%cb)
               allocate(tmpsys%cb, source=cb)
               call tmpsys%pl%setup(npl + nfrag, tmpparam)
               call tmpsys%pl%fill(pl, .not.lexclude)
               call tmpsys%rescale(tmpparam, mscale, dscale, tscale)

            end associate

         return
         end subroutine construct_temporary_system


         subroutine add_fragments_to_tmpsys(tmpsys, tmpparam)
            !! Author: David A. Minton
            !!
            !! Adds fragments to the temporary system pl object
            implicit none
            ! Arguments
            class(swiftest_nbody_system), intent(inout) :: tmpsys
            class(swiftest_parameters), intent(inout) :: tmpparam
            ! Internals
            integer(I4B) :: i
            class(swiftest_pl), allocatable :: pl_discards
            logical, dimension(:), allocatable :: lexclude_tmp

            associate(pl => system%pl, npl => system%pl%nbody)
               npl_new = npl + nfrag

               tmpsys%pl%mass(npl+1:npl_new) = m_frag(1:nfrag)
               tmpsys%pl%Gmass(npl+1:npl_new) = m_frag(1:nfrag) * tmpparam%GU
               tmpsys%pl%radius(npl+1:npl_new) = rad_frag(1:nfrag)
               do concurrent (i = 1:nfrag)
                  tmpsys%pl%xb(:,npl+i) =  xb_frag(:,i) 
                  tmpsys%pl%vb(:,npl+i) =  vb_frag(:,i) 
                  tmpsys%pl%xh(:,npl+i) =  xb_frag(:,i) - tmpsys%cb%xb(:)
                  tmpsys%pl%vh(:,npl+i) =  vb_frag(:,i) - tmpsys%cb%vb(:)
               end do
               if (tmpparam%lrotation) then
                  tmpsys%pl%Ip(:,npl+1:npl_new) = Ip_frag(:,1:nfrag)
                  tmpsys%pl%rot(:,npl+1:npl_new) = rot_frag(:,1:nfrag)
               end if
               ! Disable the collisional family for subsequent energy calculations and coordinate shifts
               lexclude(family(:)) = .true.
               lexclude(npl+1:npl_new) = .false.
               where(lexclude(1:npl_new)) 
                  tmpsys%pl%status(1:npl_new) = INACTIVE
               elsewhere
                  tmpsys%pl%status(1:npl_new) = ACTIVE
               end where
               allocate(pl_discards, mold=tmpsys%pl)
               call tmpsys%pl%spill(pl_discards, lspill_list=lexclude(1:npl_new), ldestructive=.true.)
               npl_new = count(.not.lexclude(:))

               if (size(lexclude) /= npl_new) then 
                  allocate(lexclude_tmp(npl_new))
                  call move_alloc(lexclude_tmp, lexclude)
               end if
               lexclude(1:npl_new) = .false.

            end associate

            return
         end subroutine add_fragments_to_tmpsys


         subroutine calculate_system_energy(linclude_fragments)
            !! Author: David A. Minton
            !!
            !! Calculates total system energy, including all bodies in the pl list that do not have a corresponding value of the lexclude array that is true
            !! and optionally including fragments.
            implicit none
            ! Arguments
            logical,                intent(in) :: linclude_fragments
            ! Internals
            integer(I4B) :: i, nplm
            logical, dimension(:), allocatable :: lexclude_tmp
            logical :: lk_plpl
            class(swiftest_nbody_system), allocatable :: tmpsys
            class(swiftest_parameters), allocatable   :: tmpparam

            ! Build the internal planet list out of the non-excluded bodies and optionally with fragments appended. This
            ! will get passed to the energy calculation subroutine so that energy is computed exactly the same way is it
            ! is in the main program. This will temporarily expand the planet list in a temporary system object called tmpsys to feed it into symba_energy
            associate(pl => system%pl, npl => system%pl%nbody, cb => system%cb)

               ! Because we're making a copy of symba_pl with the excludes/fragments appended, we need to deallocate the
               ! big k_plpl array and recreate it when we're done, otherwise we run the risk of blowing up the memory by
               ! allocating two of these ginormous arrays simulteouously. This is not particularly efficient, but as this
               ! subroutine should be called relatively infrequently, it shouldn't matter too much.
               lk_plpl = allocated(pl%k_plpl)
               if (lk_plpl) deallocate(pl%k_plpl)

               call construct_temporary_system(tmpsys, tmpparam)
               if (linclude_fragments) call add_fragments_to_tmpsys(tmpsys, tmpparam)

               call tmpsys%pl%index(param)

               call tmpsys%get_energy_and_momentum(param) 

               ! Restore the big array
               deallocate(tmpsys%pl%k_plpl) 

               if (lk_plpl) call pl%index(param)

               ! Calculate the current fragment energy and momentum balances
               if (linclude_fragments) then
                  Lorbit_after(:) = tmpsys%Lorbit
                  Lspin_after(:) = tmpsys%Lspin
                  Ltot_after(:) = tmpsys%Lorbit(:) + tmpsys%Lspin(:)
                  Lmag_after = norm2(Ltot_after(:))
                  ke_orbit_after = tmpsys%ke_orbit
                  ke_spin_after = tmpsys%ke_spin
                  pe_after = tmpsys%pe
                  Etot_after = tmpsys%te
                  dEtot = Etot_after - Etot_before 
                  dL(:) = Ltot_after(:) - Ltot_before(:)
                  dLmag = .mag.dL(:)
               else
                  Lorbit_before(:) = tmpsys%Lorbit
                  Lspin_before(:) = tmpsys%Lspin
                  Ltot_before(:) = tmpsys%Lorbit(:) + tmpsys%Lspin(:)
                  Lmag_before = norm2(Ltot_before(:))
                  ke_orbit_before = tmpsys%ke_orbit
                  ke_spin_before = tmpsys%ke_spin
                  pe_before = tmpsys%pe
                  Etot_before = tmpsys%te
               end if
            end associate

            return
         end subroutine calculate_system_energy


         subroutine calculate_fragment_ang_mtm() 
            !! Author: David A. Minton
            !!
            !! Calcualtes the current angular momentum of the fragments
            implicit none
            integer(I4B) :: i

            L_frag_orb(:) = 0.0_DP
            L_frag_spin(:) = 0.0_DP

            do i = 1, nfrag
               L_frag_orb(:) = L_frag_orb(:) + m_frag(i) * (x_frag(:, i) .cross. v_frag(:, i))
               L_frag_spin(:) = L_frag_spin(:) + m_frag(i) * rad_frag(i)**2 * Ip_frag(:, i) * rot_frag(:, i)
            end do

            L_frag_tot(:) = L_frag_orb(:) + L_frag_spin(:)

            return
         end subroutine calculate_fragment_ang_mtm


         subroutine shift_vector_to_origin(m_frag, vec_frag)
            !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! Adjusts the position or velocity of the fragments as needed to align them with the origin
            implicit none
            ! Arguments
            real(DP), dimension(:),   intent(in)    :: m_frag    !! Fragment masses
            real(DP), dimension(:,:), intent(inout) :: vec_frag  !! Fragment positions or velocities in the center of mass frame

            ! Internals
            real(DP), dimension(NDIM)               :: mvec_frag, COM_offset
            integer(I4B)                            :: i

            mvec_frag(:) = 0.0_DP

            do i = 1, nfrag
               mvec_frag = mvec_frag(:) + vec_frag(:,i) * m_frag(i)
            end do
            COM_offset(:) = -mvec_frag(:) / mtot
            do i = 1, nfrag 
               vec_frag(:, i) = vec_frag(:, i) + COM_offset(:)
            end do

            return
         end subroutine shift_vector_to_origin


         subroutine set_fragment_position_vectors()
            !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! Initializes the orbits of the fragments around the center of mass. The fragments are initially placed on a plane defined by the 
            !! pre-impact angular momentum. They are distributed on an ellipse surrounding the center of mass.
            !! The initial positions do not conserve energy or momentum, so these need to be adjusted later.

            implicit none
            real(DP)  :: dis, rad
            logical, dimension(:), allocatable :: loverlap
            integer(I4B) :: i, j

            allocate(loverlap(nfrag))

            ! Place the fragments into a region that is big enough that we should usually not have overlapping bodies
            ! An overlapping bodies will collide in the next time step, so it's not a major problem if they do (it just slows the run down)
            r_max = r_max_start
            rad = sum(radius(:))

            ! We will treat the first two fragments of the list as special cases. They get initialized the maximum distances apart along the original impactor distance vector.
            ! This is done because in a regular disruption, the first body is the largest, the second the second largest, and the rest are smaller equal-mass fragments.

            call random_number(x_frag(:,3:nfrag))
            loverlap(:) = .true.
            do while (any(loverlap(3:nfrag)))
               x_frag(:, 1) = x(:, 1) - xcom(:) 
               x_frag(:, 2) = x(:, 2) - xcom(:)
               r_max = r_max + 0.1_DP * rad
               do i = 3, nfrag
                  if (loverlap(i)) then
                     call random_number(x_frag(:,i))
                     x_frag(:, i) = 2 * (x_frag(:, i) - 0.5_DP) * r_max 
                  end if
               end do
               loverlap(:) = .false.
               do j = 1, nfrag
                  do i = j + 1, nfrag
                     dis = norm2(x_frag(:,j) - x_frag(:,i))
                     loverlap(i) = loverlap(i) .or. (dis <= (rad_frag(i) + rad_frag(j))) 
                  end do
               end do
            end do
            call shift_vector_to_origin(m_frag, x_frag)
            call define_coordinate_system()

            do i = 1, nfrag
               xb_frag(:,i) = x_frag(:,i) + xcom(:)
            end do

            xcom(:) = 0.0_DP
            do i = 1, nfrag
               xcom(:) = xcom(:) + m_frag(i) * xb_frag(:,i) 
            end do
            xcom(:) = xcom(:) / mtot

            return
         end subroutine set_fragment_position_vectors


         subroutine set_fragment_spin(lerr)
            !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! Sets the spins of a collection of fragments such that they conserve angular momentum without blowing the fragment kinetic energy budget.
            !!
            !! A failure will trigger a restructuring of the fragments so we will try new values of the radial position distribution.
            implicit none
            ! Arguments
            logical, intent(out) :: lerr
            ! Internals
            real(DP), dimension(NDIM) :: L_remainder, rot_L, rot_ke
            integer(I4B) :: i

            lerr = .false.

            ! Start the first two bodies with the same rotation as the original two impactors, then distribute the remaining angular momentum among the rest
            rot_frag(:,1:2) = rot(:, :)
            rot_frag(:,3:nfrag) = 0.0_DP
            call calculate_fragment_ang_mtm()
            L_remainder(:) = L_frag_budget(:) - L_frag_spin(:)

            ke_frag_spin = 0.0_DP
            do i = 1, nfrag
               ! Convert a fraction (f_spin) of either the remaining angular momentum or kinetic energy budget into spin, whichever gives the smaller rotation so as not to blow any budgets
               rot_ke(:) = sqrt(2 * f_spin * ke_frag_budget / (nfrag * m_frag(i) * rad_frag(i)**2 * Ip_frag(3, i))) * L_remainder(:) / norm2(L_remainder(:))
               rot_L(:) = f_spin * L_remainder(:) / (nfrag * m_frag(i) * rad_frag(i)**2 * Ip_frag(3, i))
               if (norm2(rot_ke) < norm2(rot_L)) then
                  rot_frag(:,i) = rot_frag(:, i) + rot_ke(:)
               else
                  rot_frag(:, i) = rot_frag(:, i) + rot_L(:)
               end if
               ke_frag_spin = ke_frag_spin + m_frag(i) * Ip_frag(3, i) * rad_frag(i)**2 * dot_product(rot_frag(:, i), rot_frag(:, i))
            end do
            ke_frag_spin = 0.5_DP * ke_frag_spin

            lerr = ((ke_frag_budget - ke_frag_spin - ke_frag_orbit) < 0.0_DP)

            return
         end subroutine set_fragment_spin


         subroutine set_fragment_tan_vel(lerr)
            !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! Adjusts the tangential velocities and spins of a collection of fragments such that they conserve angular momentum without blowing the fragment kinetic energy budget.
            !! This procedure works in several stages, with a goal to solve the angular and linear momentum constraints on the fragments, while still leaving a positive balance of
            !! our fragment kinetic energy (ke_frag_budget) that we can put into the radial velocity distribution.
            !!
            !! The first thing we'll try to do is solve for the tangential velocities of the first 6 fragments, using angular and linear momentum as constraints and an initial
            !! tangential velocity distribution for the remaining bodies (if there are any) that distributes their angular momentum equally between them.
            !! If that doesn't work and we blow our kinetic energy budget, we will attempt to find a tangential velocity distribution that minimizes the kinetic energy while
            !! conserving momentum. 
            !!
            !! A failure will trigger a restructuring of the fragments so we will try new values of the radial position distribution.
            implicit none
            ! Arguments
            logical, intent(out)  :: lerr
            ! Internals
            integer(I4B) :: i
            real(DP), parameter                  :: TOL_MIN = 1e-1_DP ! This doesn't have to be very accurate, as we really just want a tangential velocity distribution with less kinetic energy than our initial guess.
            real(DP), parameter                  :: TOL_INIT = 1e-14_DP
            integer(I4B), parameter              :: MAXLOOP = 10
            real(DP)                             :: tol
            real(DP), dimension(:), allocatable  :: v_t_initial
            real(DP), dimension(nfrag)           :: kefrag
            type(lambda_obj)                     :: spinfunc
            type(lambda_obj_err)                 :: objective_function
            real(DP), dimension(NDIM)            :: Li, L_remainder

            lerr = .false.

            ! write(*,*) '***************************************************'
            ! write(*,*) 'Original dis   : ',norm2(x(:,2) - x(:,1))
            ! write(*,*) 'r_max          : ',r_max
            ! write(*,*) 'f_spin         : ',f_spin
            ! write(*,*) '***************************************************'
            ! write(*,*) 'Energy balance so far: '
            ! write(*,*) 'ke_frag_budget : ',ke_frag_budget
            ! write(*,*) 'ke_orbit_before: ',ke_orbit_before 
            ! write(*,*) 'ke_orbit_after : ',ke_orbit_after  
            ! write(*,*) 'ke_spin_before : ',ke_spin_before 
            ! write(*,*) 'ke_spin_after  : ',ke_spin_after  
            ! write(*,*) 'pe_before      : ',pe_before 
            ! write(*,*) 'pe_after       : ',pe_after  
            ! write(*,*) 'Qloss          : ',Qloss
            ! write(*,*) '***************************************************'

            allocate(v_t_initial, mold=v_t_mag)
            v_t_initial(:) = 0.0_DP
            v_frag(:,:) = 0.0_DP

            ! Next we will solve for the tangential component of the velocities that both conserves linear momentum and uses the remaining angular momentum not used in spin.
            ! This will be done using a linear solver that solves for the tangential velocities of the first 6 fragments, constrained by the linear and angular momentum vectors, 
            ! which is embedded in a non-linear minimizer that will adjust the tangential velocities of the remaining i>6 fragments to minimize kinetic energy for a given momentum solution
            ! The initial conditions fed to the minimizer for the fragments will be the remaining angular momentum distributed between the fragments.
            call calculate_fragment_ang_mtm()
            call define_coordinate_system() ! Make sure that the tangential velocity components are defined properly
            L_remainder(:) = L_frag_budget(:) - L_frag_spin(:)
            do i = 1, nfrag
               v_t_initial(i) = norm2(L_remainder(:)) / ((nfrag - i + 1) * m_frag(i) * norm2(x_frag(:,i)))
               Li(:) = m_frag(i) * (x_frag(:,i) .cross. (v_t_initial(i) * v_t_unit(:, i)))
               L_remainder(:) = L_remainder(:) - Li(:)
            end do

            ! Find the local kinetic energy minimum for the system that conserves linear and angular momentum
            objective_function = lambda_obj(tangential_objective_function, lerr)

            tol = TOL_INIT
            do while(tol < TOL_MIN)
               v_t_mag(7:nfrag) = util_minimize_bfgs(objective_function, nfrag-6, v_t_initial(7:nfrag), tol, MAXLOOP, lerr)
               ! Now that the KE-minimized values of the i>6 fragments are found, calculate the momentum-conserving solution for tangential velociteis
               v_t_initial(7:nfrag) = v_t_mag(7:nfrag)
               if (.not.lerr) exit
               tol = tol * 2_DP ! Keep increasing the tolerance until we converge on a solution
            end do
            v_t_mag(1:nfrag) = solve_fragment_tan_vel(v_t_mag_input=v_t_initial(7:nfrag), lerr=lerr)

            ! Perform one final shift of the radial velocity vectors to align with the center of mass of the collisional system (the origin)
            vb_frag(:,1:nfrag) = vmag_to_vb(v_r_mag(1:nfrag), v_r_unit(:,1:nfrag), v_t_mag(1:nfrag), v_t_unit(:,1:nfrag), m_frag(1:nfrag), vcom(:)) 
            do concurrent (i = 1:nfrag)
               v_frag(:,i) = vb_frag(:,i) - vcom(:)
            end do

            ! Now do a kinetic energy budget check to make sure we are still within the budget.
            kefrag = 0.0_DP
            do concurrent(i = 1:nfrag)
               kefrag(i) = m_frag(i) * dot_product(vb_frag(:, i), vb_frag(:, i))
            end do
            ke_frag_orbit = 0.5_DP * sum(kefrag(:))

            ! If we are over the energy budget, flag this as a failure so we can try again
            lerr = ((ke_frag_budget - ke_frag_spin - ke_frag_orbit) < 0.0_DP)
            ! write(*,*) 'Tangential'
            ! write(*,*) 'Failure? ',lerr
            ! call calculate_fragment_ang_mtm()
            ! write(*,*) '|L_remainder| : ',.mag.(L_frag_budget(:) - L_frag_tot(:)) / Lmag_before
            ! write(*,*) 'ke_frag_budget: ',ke_frag_budget
            ! write(*,*) 'ke_frag_spin  : ',ke_frag_spin
            ! write(*,*) 'ke_tangential : ',ke_frag_orbit
            ! write(*,*) 'ke_radial     : ',ke_frag_budget - ke_frag_spin - ke_frag_orbit

            return
         end subroutine set_fragment_tan_vel


         function tangential_objective_function(v_t_mag_input, lerr) result(fval)
            !! Author: David A. Minton
            !!
            !! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
            implicit none
            ! Arguments
            real(DP), dimension(:),   intent(in)  :: v_t_mag_input   !! Unknown tangential component of velocity vector set previously by angular momentum constraint
            logical,                  intent(out) :: lerr            !! Error flag
            ! Result
            real(DP)                              :: fval
            ! Internals
            integer(I4B) :: i
            real(DP), dimension(NDIM,nfrag) :: v_shift
            real(DP), dimension(nfrag) :: v_t_new, kearr
            real(DP) :: keo

            lerr = .false.

            v_t_new(:) = solve_fragment_tan_vel(v_t_mag_input=v_t_mag_input(:), lerr=lerr)
            v_shift(:,:) = vmag_to_vb(v_r_mag, v_r_unit, v_t_new, v_t_unit, m_frag, vcom) 

            kearr = 0.0_DP
            do concurrent(i = 1:nfrag)
               kearr(i) = m_frag(i) * dot_product(v_shift(:, i), v_shift(:, i))
            end do
            keo = 0.5_DP * sum(kearr(:))
            fval = keo 
            lerr = .false.

            return
         end function tangential_objective_function


         function solve_fragment_tan_vel(lerr, v_t_mag_input) result(v_t_mag_output)
            !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! Adjusts the positions, velocities, and spins of a collection of fragments such that they conserve angular momentum
            implicit none
            ! Arguments
            logical,                          intent(out) :: lerr            !! Error flag 
            real(DP), dimension(:), optional, intent(in)  :: v_t_mag_input   !! Unknown tangential velocities for fragments 7:nfrag
            ! Internals
            integer(I4B)                            :: i
            ! Result
            real(DP), dimension(:), allocatable     :: v_t_mag_output

            real(DP), dimension(2 * NDIM, 2 * NDIM) :: A ! LHS of linear equation used to solve for momentum constraint in Gauss elimination code
            real(DP), dimension(2 * NDIM)           :: b  ! RHS of linear equation used to solve for momentum constraint in Gauss elimination code
            real(DP), dimension(NDIM)               :: L_lin_others, L_orb_others, L, vtmp

            lerr = .false.
            ! We have 6 constraint equations (2 vector constraints in 3 dimensions each)
            ! The first 3 are that the linear momentum of the fragments is zero with respect to the collisional barycenter
            ! The second 3 are that the sum of the angular momentum of the fragments is conserved from the pre-impact state
            L_lin_others(:) = 0.0_DP
            L_orb_others(:) = 0.0_DP
            do i = 1, nfrag
               if (i <= 2 * NDIM) then ! The tangential velocities of the first set of bodies will be the unknowns we will solve for to satisfy the constraints
                  A(1:3, i) = m_frag(i) * v_t_unit(:, i) 
                  A(4:6, i) = m_frag(i) * rmag(i) * (v_r_unit(:, i) .cross. v_t_unit(:, i))
               else if (present(v_t_mag_input)) then
                  vtmp(:) = v_t_mag_input(i - 6) * v_t_unit(:, i)
                  L_lin_others(:) = L_lin_others(:) + m_frag(i) * vtmp(:)
                  L(:) = m_frag(i) * (x_frag(:, i) .cross. vtmp(:)) 
                  L_orb_others(:) = L_orb_others(:) + L(:)
               end if
            end do
            b(1:3) = -L_lin_others(:)
            b(4:6) = L_frag_budget(:) - L_frag_spin(:) - L_orb_others(:)
            allocate(v_t_mag_output(nfrag))
            v_t_mag_output(1:6) = util_solve_linear_system(A, b, 6, lerr)
            if (present(v_t_mag_input)) v_t_mag_output(7:nfrag) = v_t_mag_input(:)

            return 
         end function solve_fragment_tan_vel


         subroutine set_fragment_radial_velocities(lerr)
            !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! 
            !! Adjust the fragment velocities to set the fragment orbital kinetic energy. This will minimize the difference between the fragment kinetic energy and the energy budget
            implicit none
            ! Arguments
            logical,                intent(out)   :: lerr
            ! Internals
            real(DP), parameter                   :: TOL_MIN = Etol   ! This needs to be more accurate than the tangential step, as we are trying to minimize the total residual energy
            real(DP), parameter                   :: TOL_INIT = 1e-14_DP
            integer(I4B), parameter               :: MAXLOOP = 100
            real(DP)                              :: ke_radial, tol 
            integer(I4B)                          :: i, j
            real(DP), dimension(:), allocatable   :: v_r_initial, v_r_sigma
            real(DP), dimension(:,:), allocatable :: v_r
            type(lambda_obj)                      :: objective_function

            ! Set the "target" ke for the radial component
            ke_radial = ke_frag_budget - ke_frag_spin - ke_frag_orbit
            
            allocate(v_r_initial, source=v_r_mag)
            ! Initialize radial velocity magnitudes with a random value that is approximately 10% of that found by distributing the kinetic energy equally
            allocate(v_r_sigma, source=v_r_mag)
            call random_number(v_r_sigma(1:nfrag))
            v_r_sigma(1:nfrag) = sqrt(1.0_DP + 2 * (v_r_sigma(1:nfrag) - 0.5_DP) * 1e-4_DP) 
            v_r_initial(1:nfrag) = v_r_sigma(1:nfrag) * sqrt(abs(2 * ke_radial) / (m_frag(1:nfrag) * nfrag)) 

            ! Initialize the lambda function using a structure constructor that calls the init method
            ! Minimize the ke objective function using the BFGS optimizer
            objective_function = lambda_obj(radial_objective_function)
            tol = TOL_INIT
            do while(tol < TOL_MIN)
               v_r_mag = util_minimize_bfgs(objective_function, nfrag, v_r_initial, tol, MAXLOOP, lerr)
               if (.not.lerr) exit
               tol = tol * 2 ! Keep increasing the tolerance until we converge on a solution
               v_r_initial(:) = v_r_mag(:)
            end do
            
            ! Shift the radial velocity vectors to align with the center of mass of the collisional system (the origin)
            ke_frag_orbit = 0.0_DP
            vb_frag(:,1:nfrag) = vmag_to_vb(v_r_mag(1:nfrag), v_r_unit(:,1:nfrag), v_t_mag(1:nfrag), v_t_unit(:,1:nfrag), m_frag(1:nfrag), vcom(:)) 
            do i = 1, nfrag
               v_frag(:, i) = vb_frag(:, i) - vcom(:)
               ke_frag_orbit = ke_frag_orbit + m_frag(i) * dot_product(vb_frag(:, i), vb_frag(:, i))
            end do
            ke_frag_orbit = 0.5_DP * ke_frag_orbit

            ! write(*,*) 'Radial'
            ! write(*,*) 'Failure? ',lerr 
            ! write(*,*) 'ke_frag_budget: ',ke_frag_budget
            ! write(*,*) 'ke_frag_spin  : ',ke_frag_spin
            ! write(*,*) 'ke_frag_orbit : ',ke_frag_orbit
            ! write(*,*) 'ke_remainder  : ',ke_frag_budget - (ke_frag_orbit + ke_frag_spin)
            lerr = .false.

            return
         end subroutine set_fragment_radial_velocities


         function radial_objective_function(v_r_mag_input) result(fval) 
            !! Author: David A. Minton
            !!
            !! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
            implicit none
            ! Arguments
            real(DP), dimension(:),   intent(in)  :: v_r_mag_input   !! Unknown radial component of fragment velocity vector
            ! Result
            real(DP)                              :: fval      !! The objective function result, which is the square of the difference between the calculated fragment kinetic energy and our target
                                                               !! Minimizing this brings us closer to our objective
            ! Internals
            integer(I4B)                         :: i
            real(DP), dimension(:,:), allocatable :: v_shift
            real(DP), dimension(nfrag)             :: kearr
            real(DP)                               :: keo, ke_radial

            allocate(v_shift, mold=vb_frag)
            v_shift(:,:) = vmag_to_vb(v_r_mag_input, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) 
            do concurrent(i = 1:nfrag)
               kearr(i) = m_frag(i) * (Ip_frag(3, i) * rad_frag(i)**2 * dot_product(rot_frag(:, i), rot_frag(:, i)) + dot_product(v_shift(:, i), v_shift(:, i)))
            end do
            keo = 2 * ke_frag_budget - sum(kearr(:))
            ke_radial = ke_frag_budget - ke_frag_orbit - ke_frag_spin
            ! The following ensures that fval = 0 is a local minimum, which is what the BFGS method is searching for
            fval = (keo / (2 * ke_radial))**2

            return
         end function radial_objective_function


         function vmag_to_vb(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) result(vb) 
            !! Author: David A. Minton
            !!
            !! Converts radial and tangential velocity magnitudes into barycentric velocity
            implicit none
            ! Arguments
            real(DP), dimension(:),   intent(in)  :: v_r_mag   !! Unknown radial component of fragment velocity vector
            real(DP), dimension(:),   intent(in)  :: v_t_mag   !! Tangential component of velocity vector set previously by angular momentum constraint
            real(DP), dimension(:,:), intent(in)  :: v_r_unit, v_t_unit !! Radial and tangential unit vectors for each fragment
            real(DP), dimension(:),   intent(in)  :: m_frag    !! Fragment masses
            real(DP), dimension(:),   intent(in)  :: vcom      !! Barycentric velocity of collisional system center of mass
            ! Result
            real(DP), dimension(:,:), allocatable   :: vb
            ! Internals
            integer(I4B) :: i

            allocate(vb, mold=v_r_unit)
            ! Make sure the velocity magnitude stays positive
            do i = 1, nfrag
               vb(:,i) = abs(v_r_mag(i)) * v_r_unit(:, i)
            end do
            ! In order to keep satisfying the kinetic energy constraint, we must shift the origin of the radial component of the velocities to the center of mass
            call shift_vector_to_origin(m_frag, vb)
            
            do i = 1, nfrag
               vb(:, i) = vb(:, i) + v_t_mag(i) * v_t_unit(:, i) + vcom(:)
            end do

            return
         end function vmag_to_vb


         subroutine restructure_failed_fragments()
            !! Author: David A. Minton
            !!
            !! We failed to find a set of positions and velocities that satisfy all the constraints, and so we will alter the fragments and try again.
            implicit none
            integer(I4B) :: i
            real(DP), dimension(:), allocatable  :: m_frag_new, rad_frag_new
            real(DP), dimension(:,:), allocatable  :: xb_frag_new, vb_frag_new, Ip_frag_new, rot_frag_new
            real(DP) :: delta_r, delta_r_max
            real(DP), parameter :: ke_avg_deficit_target = 0.0_DP 

            ke_tot_deficit = ke_tot_deficit - (ke_frag_budget - ke_frag_orbit - ke_frag_spin)
            ke_avg_deficit = ke_tot_deficit / try
            ! Introduce a bit of noise in the radius determination so we don't just flip flop between similar failed positions
            call random_number(delta_r_max)
            delta_r_max = sum(radius(:)) * (1.0_DP + 2e-1_DP * (delta_r_max - 0.5_DP))
            if (try > 1) then
               ! Linearly interpolate the last two failed solution ke deficits to find a new distance value to try
               delta_r = (r_max_start - r_max_start_old) * (ke_avg_deficit_target - ke_avg_deficit_old) / (ke_avg_deficit - ke_avg_deficit_old)
               if (abs(delta_r) > delta_r_max) delta_r = sign(delta_r_max, delta_r)
            else
               delta_r = delta_r_max
            end if
            r_max_start_old = r_max_start
            r_max_start = r_max_start + delta_r ! The larger lever arm can help if the problem is in the angular momentum step
            ke_avg_deficit_old = ke_avg_deficit

            if (f_spin > epsilon(1.0_DP)) then ! Try reducing the fraction in spin
               f_spin = f_spin / 2
            else
               f_spin = 0.0_DP
            end if

            return
         end subroutine restructure_failed_fragments
   end subroutine fragmentation_initialize


   module subroutine fragmentation_regime(Mcb, m1, m2, rad1, rad2, xh1, xh2, vb1, vb2, den1, den2, regime, Mlr, Mslr, min_mfrag, Qloss)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Determine the collisional regime of two colliding bodies. 
      !! Current version requires all values to be converted to SI units prior to calling the function
      !!       References:
      !!       Kokubo, E., Genda, H., 2010. Formation of Terrestrial Planets from Protoplanets Under a Realistic Accretion 
      !!          Condition. ApJL 714, L21. https://doi.org/10.1088/2041-8205/714/1/L21
      !!       Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
      !!          Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
      !!       Mustill, A.J., Davies, M.B., Johansen, A., 2018. The dynamical evolution of transiting planetary systems including 
      !!          a realistic collision prescription. Mon Not R Astron Soc 478, 2896–2908. https://doi.org/10.1093/mnras/sty1273
      !!       Rufu, R., Aharonson, O., 2019. Impact Dynamics of Moons Within a Planetary Potential. J. Geophys. Res. Planets 124, 
      !!          1008–1019. https://doi.org/10.1029/2018JE005798
      !!       Stewart, S.T., Leinhardt, Z.M., 2012. Collisions between Gravity-dominated Bodies. II. The Diversity of Impact 
      !!          Outcomes during the End Stage of Planet Formation. ApJ 751, 32. https://doi.org/10.1088/0004-637X/751/1/32
      !!
      implicit none
      ! Arguments
      integer(I4B), intent(out)         :: regime
      real(DP), intent(out)          :: Mlr, Mslr
      real(DP), intent(in)           :: Mcb, m1, m2, rad1, rad2, den1, den2, min_mfrag 
      real(DP), dimension(:), intent(in)   :: xh1, xh2, vb1, vb2
      real(DP), intent(out)          :: Qloss !! The residual energy after the collision 
      ! Constants
      integer(I4B), parameter :: N1 = 1  !number of objects with mass equal to the largest remnant from LS12
      integer(I4B), parameter :: N2 = 2  !number of objects with mass larger than second largest remnant from LS12
      real(DP), parameter   :: DENSITY1 = 1000.0_DP !standard density parameter from LS12 [kg/m3]
      real(DP), parameter   :: MU_BAR = 0.37_DP !0.385#0.37#0.3333# 3.978 # 1/3 material parameter for hydrodynamic planet-size bodies (LS12)
      real(DP), parameter   :: BETA = 2.85_DP !slope of sfd for remnants from LS12 2.85
      real(DP), parameter   :: C1 = 2.43_DP  !! Kokubo & Genda (2010) eq. (3)
      real(DP), parameter   :: C2 = -0.0408_DP !! Kokubo & Genda (2010) eq. (3)
      real(DP), parameter   :: C3 = 1.86_DP !! Kokubo & Genda (2010) eq. (3)
      real(DP), parameter   :: C4 = 1.08_DP !! Kokubo & Genda (2010) eq. (3)
      real(DP), parameter   :: CRUFU = 2.0_DP - 3 * MU_BAR ! central potential variable from Rufu and Aharonson (2019)
      real(DP), parameter   :: SUPERCAT_QRATIO = 1.8_DP ! See Section 4.1 of LS12
      ! Internals
      real(DP)           :: a1, alpha, aint, b, bcrit, c_star, egy, zeta, l, lint, mu, phi, theta
      real(DP)           :: Qr, Qrd_pstar, Qr_erosion, Qr_supercat
      real(DP)           :: Vhr, Verosion, Vescp, Vhill, Vimp, Vsupercat
      real(DP)           :: Mint, Mtot
      real(DP)           :: Rp, rhill 
      real(DP)           :: Mresidual
      real(DP)           :: U_binding

      Vimp = norm2(vb2(:) - vb1(:))
      b = calc_b(xh2, vb2, xh1, vb1)
      l = (rad1 + rad2) * (1 - b)
      egy = 0.5_DP * dot_product(vb1, vb1) - GC * Mcb / norm2(xh1)
      a1 = - GC * Mcb / 2.0_DP / egy
      Mtot = m1 + m2 
      mu = (m1 * m2) / Mtot
      if (l < 2 * rad2) then
         !calculate Mint
         phi = 2 * acos((l - rad2) / rad2)
         aint = rad2**2 * (PI - (phi - sin(phi)) / 2.0_DP)
         lint = 2 * sqrt(rad2**2 - (rad2 - l / 2.0_DP) ** 2) 
         Mint = aint * lint  ![kg]
         alpha = (l**2) * (3 * rad2 - l) / (4 * (rad2**3))
      else
         alpha = 1.0_DP
         Mint = m2
      end if 
      Rp = (3 * (m1 / den1 + alpha * m2 / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! (Mustill et al. 2018)
      c_star = calc_c_star(Rp)
      !calculate Vescp
      Vescp = sqrt(2 * GC * Mtot / Rp) !Mustill et al. 2018 eq 6 
      !calculate rhill
      rhill = a1 * (m1 / 3.0_DP / (Mcb + m1))**(1.0_DP/3.0_DP)
      !calculate Vhill
      if ((rad2 + rad1) < rhill) then 
         Vhill = sqrt(2 * GC * m1 * ((rhill**2 - rhill * (rad1 + rad2)) / &
         (rhill**2 - 0.5_DP * (rad1 + rad2)**2)) / (rad1 + rad2))
      else
         Vhill = Vescp
      end if 
      !calculate Qr_pstar
      Qrd_pstar = calc_Qrd_pstar(m1, m2, alpha, c_star) * (Vhill / Vescp)**CRUFU !Rufu and Aharaonson eq (3)
      !calculate Verosion
      Qr_erosion = 2 * (1.0_DP - m1 / Mtot) * Qrd_pstar
      Verosion = (2 * Qr_erosion * Mtot / mu)** (1.0_DP / 2.0_DP)
      Qr = mu*(Vimp**2) / Mtot / 2.0_DP
      !calculate mass largest remnant Mlr 
      Mlr = (1.0_DP - Qr / Qrd_pstar / 2.0_DP) * Mtot  ! [kg] # LS12 eq (5)
      !calculate Vsupercat
      Qr_supercat = SUPERCAT_QRATIO * Qrd_pstar ! See LS12 Section 4.1 
      Vsupercat = sqrt(2 * Qr_supercat * Mtot / mu)
      !calculate Vhr
      zeta = (m1 - m2) / Mtot
      theta = 1.0_DP - b
      Vhr = Vescp * (C1 * zeta**2 * theta**(2.5_DP) + C2 * zeta**2 + C3 * theta**(2.5_DP) + C4) ! Kokubo & Genda (2010) eq. (3)
      bcrit = rad1 / (rad1 + rad2)
      Qloss = 0.0_DP
      U_binding = (3.0_DP * Mtot) / (5.0_DP * Rp) ! LS12 eq. 27

      if ((m1 < min_mfrag).or.(m2 < min_mfrag)) then 
         regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
         Mlr = Mtot
         Mslr = 0.0_DP
         Qloss = 0.0_DP
         write(*,*) "FORCE MERGE"
      else 
         if( Vimp < Vescp) then
            regime = COLLRESOLVE_REGIME_MERGE !perfect merging regime
            Mlr = Mtot
            Mslr = 0.0_DP
            Qloss = 0.0_DP
         else if (Vimp < Verosion) then 
            if (b < bcrit) then
               regime = COLLRESOLVE_REGIME_MERGE !partial accretion regime"
               Mlr = Mtot
               Mslr = 0.0_DP
               Qloss = 0.0_DP
            else if ((b > bcrit) .and. (Vimp < Vhr)) then
               regime = COLLRESOLVE_REGIME_MERGE ! graze and merge
               Mlr = Mtot
               Mslr = 0.0_DP
               Qloss = 0.0_DP
            else
               Mlr = m1
               Mslr = calc_Qrd_rev(m2, m1, Mint, den1, den2, Vimp, c_star)
               regime = COLLRESOLVE_REGIME_HIT_AND_RUN !hit and run
               Qloss = (c_star + 1.0_DP) * U_binding ! Qr 
            end if 
         else if (Vimp > Verosion .and. Vimp < Vsupercat) then
            if (m2 < 0.001_DP * m1) then 
               regime = COLLRESOLVE_REGIME_MERGE !cratering regime"
               Mlr = Mtot
               Mslr = 0.0_DP
               Qloss = 0.0_DP
            else 
               Mslr = Mtot * (3.0_DP - BETA) * (1.0_DP - N1 * Mlr / Mtot) / (N2 * BETA)  ! LS12 eq (37)
               regime = COLLRESOLVE_REGIME_DISRUPTION !disruption
               Qloss = (c_star + 1.0_DP) * U_binding ! Qr - Qr_erosion
            end if 
         else if (Vimp > Vsupercat) then 
            Mlr = Mtot * 0.1_DP * (Qr / (Qrd_pstar * SUPERCAT_QRATIO))**(-1.5_DP)   !LS12 eq (44)
            Mslr = Mtot * (3.0_DP - BETA) * (1.0_DP - N1 * Mlr / Mtot) / (N2 * BETA)  !LS12 eq (37)
            regime = COLLRESOLVE_REGIME_SUPERCATASTROPHIC ! supercatastrophic
            Qloss = (c_star + 1.0_DP) * U_binding ! Qr - Qr_supercat
         else 
            write(*,*) "Error no regime found in symba_regime"
         end if 
      end if 
      Mresidual = Mtot - Mlr - Mslr
      if (Mresidual < 0.0_DP) then ! prevents final masses from going negative
         Mlr = Mlr + Mresidual
      end if
         
      return 

      ! Internal functions
      contains
         function calc_Qrd_pstar(Mtarg, Mp, alpha, c_star) result(Qrd_pstar)
            !! author: Jennifer L.L. Pouplin and Carlisle A. Wishard
            !!
            !! Calculates the corrected Q* for oblique impacts. See Eq. (15) of LS12.
            !!       Reference:
            !!       Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
            !!          Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
            !! 
            implicit none
            ! Arguments
            real(DP),intent(in) :: Mtarg, Mp, alpha, c_star
            ! Result
            real(DP)      :: Qrd_pstar
            ! Internals
            real(DP)      :: Qrd_star1, mu_alpha, mu, Qrd_star

            ! calc mu, mu_alpha
            mu = (Mtarg * Mp) / (Mtarg + Mp)  ! [kg]
            mu_alpha = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)  ! [kg]
            ! calc Qrd_star1
            Qrd_star1 = (c_star * 4 * PI * DENSITY1 * GC * Rp**2) / 5.0_DP
            ! calc Qrd_star
            Qrd_star = Qrd_star1 * (((Mp / Mtarg + 1.0_DP)**2) / (4 * Mp / Mtarg))**(2.0_DP / (3.0_DP * MU_BAR) - 1.0_DP)  !(eq 23)
            ! calc Qrd_pstar, v_pstar
            Qrd_pstar = ((mu / mu_alpha)**(2.0_DP - 3.0_DP * MU_BAR / 2.0_DP)) * Qrd_star  ! (eq 15)

            return
         end function calc_Qrd_pstar

         function calc_Qrd_rev(Mp, Mtarg, Mint, den1, den2, Vimp, c_star) result(Mslr)
            !! author: Jennifer L.L. Pouplin and Carlisle A. Wishard
            !!
            !! Calculates mass of second largest fragment.
            !! 
            implicit none
            ! Arguments
            real(DP),intent(in) :: Mp, Mtarg, Mint, den1, den2, Vimp, c_star
            ! Result
            real(DP) :: Mslr
            ! Internals
            real(DP) :: mtot_rev, mu_rev, gamma_rev, Qrd_star1, Qrd_star, mu_alpha_rev
            real(DP) :: Qrd_pstar, Rc1, Qr_rev, Qrd_pstar_rev, Qr_supercat_rev

            ! calc Mslr, Rc1, mu, gammalr
            mtot_rev =  Mint + Mp
            Rc1 = (3 * (Mint / den1 + Mp / den2) / (4 * PI))**(1.0_DP/3.0_DP) ! [m] Mustill et al 2018
            mu_rev = (Mint * Mp) / mtot_rev ! [kg] eq 49 LS12
            mu_alpha_rev = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)
            gamma_rev = Mint / Mp ! eq 50 LS12
            !calc Qr_rev
            Qr_rev = mu_rev * (Vimp**2) / (2 * mtot_rev)
            ! calc Qrd_star1, v_star1
            Qrd_star1 = (c_star * 4 * PI * mtot_rev * GC ) / Rc1 / 5.0_DP
            ! calc Qrd_pstar_rev
            Qrd_star = Qrd_star1 * (((gamma_rev + 1.0_DP)**2) / (4 * gamma_rev)) ** (2.0_DP / (3.0_DP * MU_BAR) - 1.0_DP) !(eq 52)
            Qrd_pstar = Qrd_star * ((mu_rev / mu_alpha_rev)**(2.0_DP - 3.0_DP * MU_BAR / 2.0_DP))
            Qrd_pstar_rev = Qrd_pstar * (Vhill / Vescp)**CRUFU !Rufu and Aharaonson eq (3)
            !calc Qr_supercat_rev
            Qr_supercat_rev = 1.8_DP * Qrd_pstar_rev 
            if (Qr_rev > Qr_supercat_rev ) then 
               Mslr = mtot_rev * (0.1_DP * ((Qr_rev / (Qrd_pstar_rev * 1.8_DP))**(-1.5_DP)))   !eq (44)
            else if ( Qr_rev < Qrd_pstar_rev ) then 
               Mslr = Mp 
            else 
               Mslr = (1.0_DP - Qr_rev / Qrd_pstar_rev / 2.0_DP) * (mtot_rev)  ! [kg] #(eq 5)
            end if 

            if ( Mslr > Mp ) Mslr = Mp !check conservation of mass

            return
         end function calc_Qrd_rev

         function calc_b(proj_pos, proj_vel, targ_pos, targ_vel) result(sintheta)
            !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! Calculates the impact factor b = sin(theta), where theta is the angle between the relative velocity
            !! and distance vectors of the target and projectile bodies. See Fig. 2 of Leinhardt and Stewart (2012)
            !! 
            implicit none
            !! Arguments
            real(DP), dimension(:), intent(in) :: proj_pos, proj_vel, targ_pos, targ_vel
            !! Result
            real(DP)             :: sintheta
            !! Internals
            real(DP), dimension(NDIM)     :: imp_vel, distance, x_cross_v      

            imp_vel(:) = proj_vel(:) - targ_vel(:)
            distance(:) = proj_pos(:) - targ_pos(:)
            x_cross_v(:) = distance(:) .cross. imp_vel(:) 
            sintheta = norm2(x_cross_v(:)) / norm2(distance(:)) / norm2(imp_vel(:))
            return 
         end function calc_b

         function calc_c_star(Rc1) result(c_star)
            !! author: David A. Minton
            !!
            !! Calculates c_star as a function of impact equivalent radius. It inteRpolates between 5 for ~1 km sized bodies to
            !! 1.8 for ~10000 km sized bodies. See LS12 Fig. 4 for details.
            !! 
            implicit none
            !! Arguments
            real(DP), intent(in) :: Rc1
            !! Result
            real(DP)             :: c_star
            !! Internals
            real(DP), parameter  :: loR   = 1.0e3_DP ! Lower bound of inteRpolation size (m)
            real(DP), parameter  :: hiR   = 1.0e7_DP ! Upper bound of inteRpolation size (m)
            real(DP), parameter  :: loval = 5.0_DP   ! Value of C* at lower bound
            real(DP), parameter  :: hival = 1.9_DP   ! Value of C* at upper bound

            if (Rc1 < loR) then
               c_star = loval
            else if (Rc1 < hiR) then
               c_star = loval + (hival - loval) * log(Rc1 / loR) / log(hiR /loR)
            else
               c_star = hival
            end if
            return
         end function calc_c_star 

   end subroutine fragmentation_regime

end submodule s_fragmentation