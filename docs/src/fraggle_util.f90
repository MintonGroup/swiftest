submodule(fraggle_classes) s_fraggle_util
   use swiftest
contains

   module subroutine fraggle_util_add_fragments_to_system(frag, colliders, system, param)
      !! Author: David A. Minton
      !!
      !! Adds fragments to the temporary system pl object
      implicit none
      ! Arguments
      class(fraggle_fragments),     intent(in)    :: frag      !! Fraggle fragment system object
      class(fraggle_colliders),     intent(in)    :: colliders !! Fraggle collider system object
      class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param     !! Current swiftest run configuration parameters
      ! Internals
      integer(I4B) :: i, npl_before, npl_after
      logical, dimension(:), allocatable :: lexclude

      associate(nfrag => frag%nbody, pl => system%pl, cb => system%cb)
         npl_after = pl%nbody
         npl_before = npl_after - nfrag
         allocate(lexclude(npl_after))

         pl%status(npl_before+1:npl_after) = ACTIVE
         pl%mass(npl_before+1:npl_after) = frag%mass(1:nfrag)
         pl%Gmass(npl_before+1:npl_after) = frag%mass(1:nfrag) * param%GU
         pl%radius(npl_before+1:npl_after) = frag%radius(1:nfrag)
         do concurrent (i = 1:nfrag)
            pl%xb(:,npl_before+i) =  frag%xb(:,i) 
            pl%vb(:,npl_before+i) =  frag%vb(:,i) 
            pl%xh(:,npl_before+i) =  frag%xb(:,i) - cb%xb(:)
            pl%vh(:,npl_before+i) =  frag%vb(:,i) - cb%vb(:)
         end do
         if (param%lrotation) then
            pl%Ip(:,npl_before+1:npl_after) = frag%Ip(:,1:nfrag)
            pl%rot(:,npl_before+1:npl_after) = frag%rot(:,1:nfrag)
         end if
         ! This will remove the colliders from the system since we've replaced them with fragments
         lexclude(1:npl_after) = .false.
         lexclude(colliders%idx(1:colliders%ncoll)) = .true.
         where(lexclude(1:npl_after)) 
            pl%status(1:npl_after) = INACTIVE
         elsewhere
            pl%status(1:npl_after) = ACTIVE
         endwhere

      end associate

      return
   end subroutine fraggle_util_add_fragments_to_system
   

   module subroutine fraggle_util_ang_mtm(self) 
      !! Author: David A. Minton
      !!
      !! Calcualtes the current angular momentum of the fragments
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout)  :: self !! Fraggle fragment system object
      ! Internals
      integer(I4B) :: i

      associate(frag => self, nfrag => self%nbody)
         frag%L_orbit(:) = 0.0_DP
         frag%L_spin(:) = 0.0_DP
   
         do i = 1, nfrag
            frag%L_orbit(:) = frag%L_orbit(:) + frag%mass(i) * (frag%x_coll(:, i) .cross. frag%v_coll(:, i))
            frag%L_spin(:) = frag%L_spin(:) + frag%mass(i) * frag%radius(i)**2 * frag%Ip(:, i) * frag%rot(:, i)
         end do
      end associate

      return
   end subroutine fraggle_util_ang_mtm


   module subroutine fraggle_util_construct_temporary_system(frag, system, param, tmpsys, tmpparam)
      !! Author: David A. Minton
      !!
      !! Constructs a temporary internal system consisting of active bodies and additional fragments. This internal temporary system is used to calculate system energy with and without fragments
      !! and optionally including fragments.
      implicit none
      ! Arguments
      class(fraggle_fragments),                   intent(in)  :: frag     !! Fraggle fragment system object
      class(swiftest_nbody_system),               intent(in)  :: system   !! Original swiftest nbody system object
      class(swiftest_parameters),                 intent(in)  :: param    !! Current swiftest run configuration parameters
      class(swiftest_nbody_system), allocatable,  intent(out) :: tmpsys   !! Output temporary swiftest nbody system object
      class(swiftest_parameters),   allocatable,  intent(out) :: tmpparam !! Output temporary configuration run parameters
      ! Internals
      logical, dimension(:), allocatable :: linclude
      integer(I4B) :: npl_tot

      associate(nfrag => frag%nbody, pl => system%pl, npl => system%pl%nbody, cb => system%cb)
         ! Set up a new system based on the original
         if (allocated(tmpparam)) deallocate(tmpparam)
         if (allocated(tmpsys)) deallocate(tmpsys)
         allocate(tmpparam, source=param)
         call setup_construct_system(tmpsys, tmpparam)

         ! No test particles necessary for energy/momentum calcs
         call tmpsys%tp%setup(0, param)

         ! Replace the empty central body object with a copy of the original
         deallocate(tmpsys%cb)
         allocate(tmpsys%cb, source=cb)

         ! Make space for the fragments
         npl_tot = npl + nfrag
         call tmpsys%pl%setup(npl_tot, tmpparam)
         allocate(linclude(npl_tot))

         ! Fill up the temporary system with all of the original bodies, leaving the spaces for fragments empty until we add them in later
         linclude(1:npl) = .true.
         linclude(npl+1:npl_tot) = .false.
         call tmpsys%pl%fill(pl, linclude)

         ! Scale the temporary system to the natural units of the current Fraggle calculation
         call tmpsys%rescale(tmpparam, frag%mscale, frag%dscale, frag%tscale)

      end associate

      return
   end subroutine fraggle_util_construct_temporary_system


   module subroutine fraggle_util_get_energy_momentum(self, colliders, system, param, lbefore)
      !! Author: David A. Minton
      !!
      !! Calculates total system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision outcome state (lbefore = .false.)
      !! This subrourtine works by building a temporary internal massive body object out of the non-excluded bodies and optionally with fragments appended. 
      !! This will get passed to the energy calculation subroutine so that energy is computed exactly the same way is it is in the main program. 
      !! This will temporarily expand the massive body object in a temporary system object called tmpsys to feed it into symba_energy
      implicit none
      ! Arguments
      class(fraggle_fragments),     intent(inout) :: self      !! Fraggle fragment system object
      class(fraggle_colliders),     intent(inout) :: colliders !! Fraggle collider system object
      class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param     !! Current swiftest run configuration parameters
      logical,                      intent(in)    :: lbefore   !! Flag indicating that this the "before" state of the system, with colliders included and fragments excluded or vice versa
      ! Internals
      integer(I4B) :: i, nplm
      logical, dimension(:), allocatable :: lexclude
      logical :: lk_plpl
      logical, save :: ladd_frags
      class(swiftest_nbody_system), allocatable, save :: tmpsys
      class(swiftest_parameters), allocatable, save   :: tmpparam
      integer(I4B)  :: npl_before, npl_after

      associate(frag => self, nfrag => self%nbody, pl => system%pl, cb => system%cb)

         ! Because we're making a copy of the massive body object with the excludes/fragments appended, we need to deallocate the
         ! big k_plpl array and recreate it when we're done, otherwise we run the risk of blowing up the memory by
         ! allocating two of these ginormous arrays simulteouously. This is not particularly efficient, but as this
         ! subroutine should be called relatively infrequently, it shouldn't matter too much.

         npl_before = pl%nbody
         npl_after = npl_before + nfrag

         ! Build the exluded body logical mask
         allocate(lexclude(npl_after))
         if (lbefore) then
            lexclude(1:npl_before) = .false.
            lexclude(npl_before+1:npl_after) = .true.
            call fraggle_util_construct_temporary_system(frag, system, param, tmpsys, tmpparam)
         else
            lexclude(1:npl_after) = .false.
            lexclude(colliders%idx(1:colliders%ncoll)) = .true.
            if (.not.allocated(tmpsys)) then
               write(*,*) "Error in fraggle_util_get_energy_momentum. This must be called with lbefore=.true. at least once before calling it with lbefore=.false."
               call util_exit(FAILURE)
            end if
            call fraggle_util_add_fragments_to_system(frag, colliders, tmpsys, tmpparam)
         end if 

         call tmpsys%pl%index(param)

         call tmpsys%get_energy_and_momentum(param) 


         ! Calculate the current fragment energy and momentum balances
         if (lbefore) then
            frag%Lorbit_before(:) = tmpsys%Lorbit(:)
            frag%Lspin_before(:) = tmpsys%Lspin(:)
            frag%Ltot_before(:) = tmpsys%Ltot(:)
            frag%ke_orbit_before = tmpsys%ke_orbit
            frag%ke_spin_before = tmpsys%ke_spin
            frag%pe_before = tmpsys%pe
            frag%Etot_before = tmpsys%te
         else
            frag%Lorbit_after(:) = tmpsys%Lorbit(:)
            frag%Lspin_after(:) = tmpsys%Lspin(:)
            frag%Ltot_after(:) = tmpsys%Ltot(:)
            frag%ke_orbit_after = tmpsys%ke_orbit
            frag%ke_spin_after = tmpsys%ke_spin
            frag%pe_after = tmpsys%pe
            frag%Etot_after = tmpsys%te
         end if
      end associate

      return
   end subroutine fraggle_util_get_energy_momentum


   module subroutine fraggle_util_restructure(self, colliders, try, f_spin, r_max_start)
      !! Author: David A. Minton
      !!
      !! Restructure the inputs after a failed attempt failed to find a set of positions and velocities that satisfy the energy and momentum constraints
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: self        !! Fraggle fragment system object
      class(fraggle_colliders), intent(in)    :: colliders   !! Fraggle collider system object
      integer(I4B),             intent(in)    :: try         !! The current number of times Fraggle has tried to find a solution
      real(DP),                 intent(inout) :: f_spin      !! Fraction of energy/momentum that goes into spin. This decreases ater a failed attempt
      real(DP),                 intent(inout) :: r_max_start !! The maximum radial distance that the position calculation starts with. This increases after a failed attempt
      ! Internals
      integer(I4B) :: i
      real(DP), save :: ke_tot_deficit, r_max_start_old, ke_avg_deficit_old
      real(DP), dimension(:), allocatable  :: m_frag_new, rad_frag_new
      real(DP), dimension(:,:), allocatable  :: xb_frag_new, vb_frag_new, Ip_frag_new, rot_frag_new
      real(DP) :: delta_r, delta_r_max, ke_avg_deficit
      real(DP), parameter :: ke_avg_deficit_target = 0.0_DP 

      ! Introduce a bit of noise in the radius determination so we don't just flip flop between similar failed positions
      associate(frag => self)
         call random_number(delta_r_max)
         delta_r_max = sum(colliders%radius(:)) * (1.0_DP + 2e-1_DP * (delta_r_max - 0.5_DP))
         if (try == 1) then
            ke_tot_deficit = - (frag%ke_budget - frag%ke_orbit - frag%ke_spin)
            ke_avg_deficit = ke_tot_deficit
            delta_r = delta_r_max
         else
            ! Linearly interpolate the last two failed solution ke deficits to find a new distance value to try
            ke_tot_deficit = ke_tot_deficit - (frag%ke_budget - frag%ke_orbit - frag%ke_spin)
            ke_avg_deficit = ke_tot_deficit / try
            delta_r = (r_max_start - r_max_start_old) * (ke_avg_deficit_target - ke_avg_deficit_old) / (ke_avg_deficit - ke_avg_deficit_old)
            if (abs(delta_r) > delta_r_max) delta_r = sign(delta_r_max, delta_r)
         end if
         r_max_start_old = r_max_start
         r_max_start = r_max_start + delta_r ! The larger lever arm can help if the problem is in the angular momentum step
         ke_avg_deficit_old = ke_avg_deficit
   
         if (f_spin > epsilon(1.0_DP)) then ! Try reducing the fraction in spin
            f_spin = f_spin / 2
         else
            f_spin = 0.0_DP
         end if
      end associate 

      return
   end subroutine fraggle_util_restructure


   module subroutine fraggle_util_shift_vector_to_origin(m_frag, vec_frag)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the position or velocity of the fragments as needed to align them with the origin
      implicit none
      ! Arguments
      real(DP), dimension(:),   intent(in)    :: m_frag    !! Fragment masses
      real(DP), dimension(:,:), intent(inout) :: vec_frag  !! Fragment positions or velocities in the center of mass frame

      ! Internals
      real(DP), dimension(NDIM) :: mvec_frag, COM_offset
      integer(I4B) :: i, nfrag
      real(DP) :: mtot

      mvec_frag(:) = 0.0_DP
      mtot = sum(m_frag)
      nfrag = size(m_frag)

      do i = 1, nfrag
         mvec_frag = mvec_frag(:) + vec_frag(:,i) * m_frag(i)
      end do
      COM_offset(:) = -mvec_frag(:) / mtot
      do i = 1, nfrag 
         vec_frag(:, i) = vec_frag(:, i) + COM_offset(:)
      end do

      return
   end subroutine fraggle_util_shift_vector_to_origin


   module function fraggle_util_vmag_to_vb(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) result(vb) 
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
      integer(I4B) :: i, nfrag

      allocate(vb, mold=v_r_unit)
      ! Make sure the velocity magnitude stays positive
      nfrag = size(m_frag)
      do i = 1, nfrag
         vb(:,i) = abs(v_r_mag(i)) * v_r_unit(:, i)
      end do
      ! In order to keep satisfying the kinetic energy constraint, we must shift the origin of the radial component of the velocities to the center of mass
      call fraggle_util_shift_vector_to_origin(m_frag, vb)
      
      do i = 1, nfrag
         vb(:, i) = vb(:, i) + v_t_mag(i) * v_t_unit(:, i) + vcom(:)
      end do

      return
   end function fraggle_util_vmag_to_vb


end submodule s_fraggle_util
