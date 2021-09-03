submodule(fraggle_classes) s_fraggle_generate
   use swiftest

   integer(I4B), parameter :: NFRAG_MIN = 7 !! The minimum allowable number of fragments (set to 6 because that's how many unknowns are needed in the tangential velocity calculation)
   real(DP),     parameter :: F_SPIN_FIRST = 0.05_DP !! The initial try value of the fraction of energy or momenum in spin (whichever has the lowest kinetic energy)
   real(DP), parameter     :: FRAGGLE_LTOL = 10 * epsilon(1.0_DP)
   real(DP), parameter     :: FRAGGLE_ETOL = 1e-8_DP

contains

   module subroutine fraggle_generate_fragments(self, colliders, system, param, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Generates a system of fragments in barycentric coordinates that conserves energy and momentum.
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(fraggle_fragments),     intent(inout) :: self      !! Fraggle system object the outputs will be the fragmentation 
      class(fraggle_colliders),     intent(inout) :: colliders !! Fraggle colliders object containing the two-body equivalent values of the colliding bodies 
      class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param     !! Current run configuration parameters 
      logical,                      intent(out)   :: lfailure  !! Answers the question: Should this have been a merger instead?
      ! Internals
      integer(I4B)                         :: i
      integer(I4B)                         :: try
      real(DP)                             :: r_max_start, f_spin, dEtot, dLmag
      integer(I4B), parameter              :: MAXTRY = 100
      logical                              :: lk_plpl
      logical, dimension(size(IEEE_ALL))   :: fpe_halting_modes, fpe_quiet_modes
      logical, dimension(size(IEEE_USUAL)) :: fpe_flag 
      character(len=STRMAX)                :: message

      ! The minimization and linear solvers can sometimes lead to floating point exceptions. Rather than halting the code entirely if this occurs, we
      ! can simply fail the attempt and try again. So we need to turn off any floating point exception halting modes temporarily 
      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      fpe_quiet_modes(:) = .false.
      call ieee_set_halting_mode(IEEE_ALL,fpe_quiet_modes)

      associate(frag => self, nfrag => self%nbody, pl => system%pl)

         if (nfrag < NFRAG_MIN) then
            write(message,*) "Fraggle needs at least ",NFRAG_MIN," fragments, but only ",nfrag," were given."
            call fraggle_io_log_one_message(message)
            lfailure = .true.
            return
         end if
         f_spin = F_SPIN_FIRST

         lk_plpl = allocated(pl%k_plpl)
         if (lk_plpl) deallocate(pl%k_plpl)

         call frag%set_natural_scale(colliders)

         call frag%reset()

         ! Calculate the initial energy of the system without the collisional family
         call frag%get_energy_and_momentum(colliders, system, param, lbefore=.true.)
        
         ! Start out the fragments close to the initial separation distance. This will be increased if there is any overlap or we fail to find a solution
         r_max_start = 1 * norm2(colliders%xb(:,2) - colliders%xb(:,1))
         lfailure = .false.
         try = 1
         do while (try < MAXTRY)
            write(message,*) try
            call fraggle_io_log_one_message("Fraggle try " // trim(adjustl(message)))
            if (lfailure) then
               call frag%restructure(colliders, try, f_spin, r_max_start)
               call frag%reset()
               try = try + 1
            end if

            lfailure = .false.
            call ieee_set_flag(ieee_all, .false.) ! Set all fpe flags to quiet

            call fraggle_generate_pos_vec(frag, colliders, r_max_start)
            call frag%set_coordinate_system(colliders)

            ! Initial velocity guess will be the barycentric velocity of the colliding system so that the budgets are based on the much smaller collisional-frame velocities
            do concurrent (i = 1:nfrag)
               frag%vb(:, i) = frag%vbcom(:)
            end do

            call frag%get_energy_and_momentum(colliders, system, param, lbefore=.false.)
            call frag%set_budgets(colliders)

            call fraggle_generate_spins(frag, colliders, f_spin, lfailure)
            if (lfailure) then
               call fraggle_io_log_one_message("Fraggle failed to find spins")
               cycle
            end if

            call fraggle_generate_tan_vel(frag, colliders, lfailure)
            if (lfailure) then
               call fraggle_io_log_one_message("Fraggle failed to find tangential velocities")
               cycle
            end if

            call fraggle_generate_rad_vel(frag, colliders, lfailure)
            if (lfailure) then
               call fraggle_io_log_one_message("Fraggle failed to find radial velocities")
               cycle
            end if

            call frag%get_energy_and_momentum(colliders, system, param, lbefore=.false.)
            dEtot = frag%Etot_after - frag%Etot_before 
            dLmag = .mag. (frag%Ltot_after(:) - frag%Ltot_before(:))

            lfailure = ((abs(dEtot + frag%Qloss) > FRAGGLE_ETOL) .or. (dEtot > 0.0_DP)) 
            if (lfailure) then
               write(message, *) dEtot, abs(dEtot + frag%Qloss) / FRAGGLE_ETOL
               call fraggle_io_log_one_message("Fraggle failed due to high energy error: " // trim(adjustl(message)))
               cycle
            end if

            lfailure = ((abs(dLmag) / (.mag.frag%Ltot_before)) > FRAGGLE_LTOL) 
            if (lfailure) then
               write(message,*) dLmag / (.mag.frag%Ltot_before(:))
               call fraggle_io_log_one_message("Fraggle failed due to high angular momentum error: " // trim(adjustl(message)))
               cycle
            end if

            ! Check if any of the usual floating point exceptions happened, and fail the try if so
            call ieee_get_flag(ieee_usual, fpe_flag)
            lfailure = any(fpe_flag) 
            if (.not.lfailure) exit
            write(message,*) "Fraggle failed due to a floating point exception: ", fpe_flag
            call fraggle_io_log_one_message(message)
         end do

         write(message,*) try
         if (lfailure) then
            call fraggle_io_log_one_message("Fraggle fragment generation failed after " // trim(adjustl(message)) // " tries")
         else
            call fraggle_io_log_one_message("Fraggle fragment generation succeeded after " // trim(adjustl(message)) // " tries")
            call fraggle_io_log_generate(frag)
         end if

         call frag%set_original_scale(colliders)

         ! Restore the big array
         if (lk_plpl) call pl%index(param)
      end associate
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily

      return 
   end subroutine fraggle_generate_fragments


   subroutine fraggle_generate_pos_vec(frag, colliders, r_max_start)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the orbits of the fragments around the center of mass. The fragments are initially placed on a plane defined by the 
      !! pre-impact angular momentum. They are distributed on an ellipse surrounding the center of mass.
      !! The initial positions do not conserve energy or momentum, so these need to be adjusted later.
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: frag        !! Fraggle fragment system object
      class(fraggle_colliders), intent(inout) :: colliders   !! Fraggle collider system object
      real(DP),                 intent(in)    :: r_max_start !! Initial guess for the starting maximum radial distance of fragments
      ! Internals
      real(DP)  :: dis, rad, r_max
      logical, dimension(:), allocatable :: loverlap
      integer(I4B) :: i, j

      associate(nfrag => frag%nbody)
         allocate(loverlap(nfrag))

         ! Place the fragments into a region that is big enough that we should usually not have overlapping bodies
         ! An overlapping bodies will collide in the next time step, so it's not a major problem if they do (it just slows the run down)
         r_max = r_max_start
         rad = sum(colliders%radius(:))

         ! We will treat the first two fragments of the list as special cases. They get initialized the maximum distances apart along the original impactor distance vector.
         ! This is done because in a regular disruption, the first body is the largest, the second the second largest, and the rest are smaller equal-mass fragments.

         call random_number(frag%x_coll(:,3:nfrag))
         loverlap(:) = .true.
         do while (any(loverlap(3:nfrag)))
            frag%x_coll(:, 1) = colliders%xb(:, 1) - frag%xbcom(:) 
            frag%x_coll(:, 2) = colliders%xb(:, 2) - frag%xbcom(:)
            r_max = r_max + 0.1_DP * rad
            do i = 3, nfrag
               if (loverlap(i)) then
                  call random_number(frag%x_coll(:,i))
                  frag%x_coll(:, i) = 2 * (frag%x_coll(:, i) - 0.5_DP) * r_max 
               end if
            end do
            loverlap(:) = .false.
            do j = 1, nfrag
               do i = j + 1, nfrag
                  dis = norm2(frag%x_coll(:,j) - frag%x_coll(:,i))
                  loverlap(i) = loverlap(i) .or. (dis <= (frag%radius(i) + frag%radius(j))) 
               end do
            end do
         end do
         call fraggle_util_shift_vector_to_origin(frag%mass, frag%x_coll)
         call frag%set_coordinate_system(colliders)

         do i = 1, nfrag
            frag%xb(:,i) = frag%x_coll(:,i) + frag%xbcom(:)
         end do

         frag%xbcom(:) = 0.0_DP
         do i = 1, nfrag
            frag%xbcom(:) = frag%xbcom(:) + frag%mass(i) * frag%xb(:,i) 
         end do
         frag%xbcom(:) = frag%xbcom(:) / frag%mtot
      end associate

      return
   end subroutine fraggle_generate_pos_vec


   subroutine fraggle_generate_spins(frag, colliders, f_spin, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Calculates the spins of a collection of fragments such that they conserve angular momentum without blowing the fragment kinetic energy budget.
      !!
      !! A failure will trigger a restructuring of the fragments so we will try new values of the radial position distribution.
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: frag      !! Fraggle fragment system object
      class(fraggle_colliders), intent(in)    :: colliders !! Fraggle collider system object
      real(DP),                 intent(in)    :: f_spin    !! Fraction of energy or momentum that goes into spin (whichever gives the lowest kinetic energy)
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds! 
      ! Internals
      real(DP), dimension(NDIM) :: L_remainder, rot_L, rot_ke
      integer(I4B)              :: i
      character(len=STRMAX)     :: message

      associate(nfrag => frag%nbody)
         lfailure = .false.

         ! Start the first two bodies with the same rotation as the original two impactors, then distribute the remaining angular momentum among the rest
         frag%rot(:,1:2) = colliders%rot(:, :)
         frag%rot(:,3:nfrag) = 0.0_DP
         call frag%get_ang_mtm()
         L_remainder(:) = frag%L_budget(:) - frag%L_spin(:)

         frag%ke_spin = 0.0_DP
         do i = 1, nfrag
            ! Convert a fraction (f_spin) of either the remaining angular momentum or kinetic energy budget into spin, whichever gives the smaller rotation so as not to blow any budgets
            rot_ke(:) = sqrt(2 * f_spin * frag%ke_budget / (nfrag * frag%mass(i) * frag%radius(i)**2 * frag%Ip(3, i))) * L_remainder(:) / norm2(L_remainder(:))
            rot_L(:) = f_spin * L_remainder(:) / (nfrag * frag%mass(i) * frag%radius(i)**2 * frag%Ip(3, i))
            if (norm2(rot_ke) < norm2(rot_L)) then
               frag%rot(:,i) = frag%rot(:, i) + rot_ke(:)
            else
               frag%rot(:, i) = frag%rot(:, i) + rot_L(:)
            end if
            frag%ke_spin = frag%ke_spin + frag%mass(i) * frag%Ip(3, i) * frag%radius(i)**2 * dot_product(frag%rot(:, i), frag%rot(:, i))
         end do
         frag%ke_spin = 0.5_DP * frag%ke_spin

         lfailure = ((frag%ke_budget - frag%ke_spin - frag%ke_orbit) < 0.0_DP)

         if (lfailure) then
            call fraggle_io_log_one_message(" ")
            call fraggle_io_log_one_message("Spin failure diagnostics")
            write(message, *) frag%ke_budget
            call fraggle_io_log_one_message("ke_budget     : " // trim(adjustl(message)))
            write(message, *) frag%ke_spin
            call fraggle_io_log_one_message("ke_spin       : " // trim(adjustl(message)))
            write(message, *) frag%ke_orbit
            call fraggle_io_log_one_message("ke_orbit      : " // trim(adjustl(message)))
            write(message, *) frag%ke_budget - frag%ke_spin - frag%ke_orbit
            call fraggle_io_log_one_message("ke_remainder  : " // trim(adjustl(message)))
         end if

      end associate

      return
   end subroutine fraggle_generate_spins


   subroutine fraggle_generate_tan_vel(frag, colliders, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the tangential velocities and spins of a collection of fragments such that they conserve angular momentum without blowing the fragment kinetic energy budget.
      !! This procedure works in several stages, with a goal to solve the angular and linear momentum constraints on the fragments, while still leaving a positive balance of
      !! our fragment kinetic energy (frag%ke_budget) that we can put into the radial velocity distribution.
      !!
      !! The first thing we'll try to do is solve for the tangential velocities of the first 6 fragments, using angular and linear momentum as constraints and an initial
      !! tangential velocity distribution for the remaining bodies (if there are any) that distributes their angular momentum equally between them.
      !! If that doesn't work and we blow our kinetic energy budget, we will attempt to find a tangential velocity distribution that minimizes the kinetic energy while
      !! conserving momentum. 
      !!
      !! A failure will trigger a restructuring of the fragments so we will try new values of the radial position distribution.
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: frag      !! Fraggle fragment system object
      class(fraggle_colliders), intent(in)    :: colliders !! Fraggle collider system object
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds
      ! Internals
      integer(I4B) :: i
      real(DP), parameter                  :: TOL_MIN = 1e-1_DP ! This doesn't have to be very accurate, as we really just want a tangential velocity distribution with less kinetic energy than our initial guess.
      real(DP), parameter                  :: TOL_INIT = 1e-14_DP
      integer(I4B), parameter              :: MAXLOOP = 10
      real(DP)                             :: tol
      real(DP), dimension(:), allocatable  :: v_t_initial
      real(DP), dimension(frag%nbody)      :: kefrag
      type(lambda_obj)                     :: spinfunc
      type(lambda_obj_err)                 :: objective_function
      real(DP), dimension(NDIM)            :: Li, L_remainder, L_frag_tot
      character(len=STRMAX)                :: message

      associate(nfrag => frag%nbody)
         lfailure = .false.

         allocate(v_t_initial, mold=frag%v_t_mag)
         v_t_initial(:) = 0.0_DP
         frag%v_coll(:,:) = 0.0_DP

         ! Next we will solve for the tangential component of the velocities that both conserves linear momentum and uses the remaining angular momentum not used in spin.
         ! This will be done using a linear solver that solves for the tangential velocities of the first 6 fragments, constrained by the linear and angular momentum vectors, 
         ! which is embedded in a non-linear minimizer that will adjust the tangential velocities of the remaining i>6 fragments to minimize kinetic energy for a given momentum solution
         ! The initial conditions fed to the minimizer for the fragments will be the remaining angular momentum distributed between the fragments.
         call frag%get_ang_mtm()
         L_remainder(:) = frag%L_budget(:) - frag%L_spin(:)
         do i = 1, nfrag
            v_t_initial(i) = norm2(L_remainder(:)) / ((nfrag - i + 1) * frag%mass(i) * norm2(frag%x_coll(:,i)))
            Li(:) = frag%mass(i) * (frag%x_coll(:,i) .cross. (v_t_initial(i) * frag%v_t_unit(:, i)))
            L_remainder(:) = L_remainder(:) - Li(:)
         end do

         ! Find the local kinetic energy minimum for the system that conserves linear and angular momentum
         objective_function = lambda_obj(tangential_objective_function, lfailure)

         tol = TOL_INIT
         do while(tol < TOL_MIN)
            frag%v_t_mag(7:nfrag) = util_minimize_bfgs(objective_function, nfrag-6, v_t_initial(7:nfrag), tol, MAXLOOP, lfailure)
            ! Now that the KE-minimized values of the i>6 fragments are found, calculate the momentum-conserving solution for tangential velociteis
            v_t_initial(7:nfrag) = frag%v_t_mag(7:nfrag)
            if (.not.lfailure) exit
            tol = tol * 2_DP ! Keep increasing the tolerance until we converge on a solution
         end do
         frag%v_t_mag(1:nfrag) = solve_fragment_tan_vel(v_t_mag_input=v_t_initial(7:nfrag), lfailure=lfailure)

         ! Perform one final shift of the radial velocity vectors to align with the center of mass of the collisional system (the origin)
         frag%vb(:,1:nfrag) = fraggle_util_vmag_to_vb(frag%v_r_mag(1:nfrag), frag%v_r_unit(:,1:nfrag), frag%v_t_mag(1:nfrag), frag%v_t_unit(:,1:nfrag), frag%mass(1:nfrag), frag%vbcom(:)) 
         do concurrent (i = 1:nfrag)
            frag%v_coll(:,i) = frag%vb(:,i) - frag%vbcom(:)
         end do

         ! Now do a kinetic energy budget check to make sure we are still within the budget.
         kefrag = 0.0_DP
         do concurrent(i = 1:nfrag)
            kefrag(i) = frag%mass(i) * dot_product(frag%vb(:, i), frag%vb(:, i))
         end do
         frag%ke_orbit = 0.5_DP * sum(kefrag(:))

         ! If we are over the energy budget, flag this as a failure so we can try again
         lfailure = ((frag%ke_budget - frag%ke_spin - frag%ke_orbit) < 0.0_DP)
         if (lfailure) then
            call fraggle_io_log_one_message(" ")
            call fraggle_io_log_one_message("Tangential velocity failure diagnostics")
            call frag%get_ang_mtm()
            L_frag_tot = frag%L_spin(:) + frag%L_orbit(:)
            write(message, *) .mag.(frag%L_budget(:) - L_frag_tot(:)) / (.mag.frag%Ltot_before(:))
            call fraggle_io_log_one_message("|L_remainder| : " // trim(adjustl(message)))
            write(message, *) frag%ke_budget
            call fraggle_io_log_one_message("ke_budget     : " // trim(adjustl(message)))
            write(message, *) frag%ke_spin
            call fraggle_io_log_one_message("ke_spin       : " // trim(adjustl(message)))
            write(message, *) frag%ke_orbit
            call fraggle_io_log_one_message("ke_tangential : " // trim(adjustl(message)))
            write(message, *) frag%ke_budget - frag%ke_spin - frag%ke_orbit
            call fraggle_io_log_one_message("ke_radial     : " // trim(adjustl(message)))
         end if
      end associate

      return

      contains
         function solve_fragment_tan_vel(lfailure, v_t_mag_input) result(v_t_mag_output)
            !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
            !!
            !! Adjusts the positions, velocities, and spins of a collection of fragments such that they conserve angular momentum
            implicit none
            ! Arguments
            logical,                          intent(out) :: lfailure            !! Error flag 
            real(DP), dimension(:), optional, intent(in)  :: v_t_mag_input   !! Unknown tangential velocities for fragments 7:nfrag
            ! Internals
            integer(I4B)                            :: i
            ! Result
            real(DP), dimension(:), allocatable     :: v_t_mag_output
      
            real(DP), dimension(2 * NDIM, 2 * NDIM) :: A ! LHS of linear equation used to solve for momentum constraint in Gauss elimination code
            real(DP), dimension(2 * NDIM)           :: b  ! RHS of linear equation used to solve for momentum constraint in Gauss elimination code
            real(DP), dimension(NDIM)               :: L_lin_others, L_orb_others, L, vtmp
      
            associate(nfrag => frag%nbody)
               lfailure = .false.
               ! We have 6 constraint equations (2 vector constraints in 3 dimensions each)
               ! The first 3 are that the linear momentum of the fragments is zero with respect to the collisional barycenter
               ! The second 3 are that the sum of the angular momentum of the fragments is conserved from the pre-impact state
               L_lin_others(:) = 0.0_DP
               L_orb_others(:) = 0.0_DP
               do i = 1, nfrag
                  if (i <= 2 * NDIM) then ! The tangential velocities of the first set of bodies will be the unknowns we will solve for to satisfy the constraints
                     A(1:3, i) = frag%mass(i) * frag%v_t_unit(:, i) 
                     A(4:6, i) = frag%mass(i) * frag%rmag(i) * (frag%v_r_unit(:, i) .cross. frag%v_t_unit(:, i))
                  else if (present(v_t_mag_input)) then
                     vtmp(:) = v_t_mag_input(i - 6) * frag%v_t_unit(:, i)
                     L_lin_others(:) = L_lin_others(:) + frag%mass(i) * vtmp(:)
                     L(:) = frag%mass(i) * (frag%x_coll(:, i) .cross. vtmp(:)) 
                     L_orb_others(:) = L_orb_others(:) + L(:)
                  end if
               end do
               b(1:3) = -L_lin_others(:)
               b(4:6) = frag%L_budget(:) - frag%L_spin(:) - L_orb_others(:)
               allocate(v_t_mag_output(nfrag))
               v_t_mag_output(1:6) = util_solve_linear_system(A, b, 6, lfailure)
               if (present(v_t_mag_input)) v_t_mag_output(7:nfrag) = v_t_mag_input(:)
            end associate
            return 
         end function solve_fragment_tan_vel


         function tangential_objective_function(v_t_mag_input, lfailure) result(fval)
            !! Author: David A. Minton
            !!
            !! Objective function for evaluating how close our fragment velocities get to minimizing KE error from our required value
            implicit none
            ! Arguments
            real(DP), dimension(:),   intent(in)  :: v_t_mag_input   !! Unknown tangential component of velocity vector set previously by angular momentum constraint
            logical,                  intent(out) :: lfailure            !! Error flag
            ! Result
            real(DP)                              :: fval
            ! Internals
            integer(I4B) :: i
            real(DP), dimension(NDIM,frag%nbody) :: v_shift
            real(DP), dimension(frag%nbody) :: v_t_new, kearr
            real(DP) :: keo
      
            associate(nfrag => frag%nbody)
               lfailure = .false.
         
               v_t_new(:) = solve_fragment_tan_vel(v_t_mag_input=v_t_mag_input(:), lfailure=lfailure)
               v_shift(:,:) = fraggle_util_vmag_to_vb(frag%v_r_mag, frag%v_r_unit, v_t_new, frag%v_t_unit, frag%mass, frag%vbcom) 
         
               kearr = 0.0_DP
               do concurrent(i = 1:nfrag)
                  kearr(i) = frag%mass(i) * dot_product(v_shift(:, i), v_shift(:, i))
               end do
               keo = 0.5_DP * sum(kearr(:))
               fval = keo 
               lfailure = .false.
            end associate

            return
         end function tangential_objective_function

   end subroutine fraggle_generate_tan_vel


   subroutine fraggle_generate_rad_vel(frag, colliders, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! 
      !! Adjust the fragment velocities to set the fragment orbital kinetic energy. This will minimize the difference between the fragment kinetic energy and the energy budget
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: frag      !! Fraggle fragment system object
      class(fraggle_colliders), intent(in)    :: colliders !! Fraggle collider system object
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds! 
      ! Internals
      real(DP), parameter                   :: TOL_MIN = FRAGGLE_ETOL   ! This needs to be more accurate than the tangential step, as we are trying to minimize the total residual energy
      real(DP), parameter                   :: TOL_INIT = 1e-14_DP
      integer(I4B), parameter               :: MAXLOOP = 100
      real(DP)                              :: ke_radial, tol 
      integer(I4B)                          :: i, j
      real(DP), dimension(:), allocatable   :: v_r_initial, v_r_sigma
      real(DP), dimension(:,:), allocatable :: v_r
      type(lambda_obj)                      :: objective_function
      character(len=STRMAX)                 :: message

      associate(nfrag => frag%nbody)
         ! Set the "target" ke for the radial component
         ke_radial = frag%ke_budget - frag%ke_spin - frag%ke_orbit
         
         allocate(v_r_initial, source=frag%v_r_mag)
         ! Initialize radial velocity magnitudes with a random value that is approximately 10% of that found by distributing the kinetic energy equally
         allocate(v_r_sigma, source=frag%v_r_mag)
         call random_number(v_r_sigma(1:nfrag))
         v_r_sigma(1:nfrag) = sqrt(1.0_DP + 2 * (v_r_sigma(1:nfrag) - 0.5_DP) * 1e-4_DP) 
         v_r_initial(1:nfrag) = v_r_sigma(1:nfrag) * sqrt(abs(2 * ke_radial) / (frag%mass(1:nfrag) * nfrag)) 

         ! Initialize the lambda function using a structure constructor that calls the init method
         ! Minimize the ke objective function using the BFGS optimizer
         objective_function = lambda_obj(radial_objective_function)
         tol = TOL_INIT
         do while(tol < TOL_MIN)
            frag%v_r_mag = util_minimize_bfgs(objective_function, nfrag, v_r_initial, tol, MAXLOOP, lfailure)
            if (.not.lfailure) exit
            tol = tol * 2 ! Keep increasing the tolerance until we converge on a solution
            v_r_initial(:) = frag%v_r_mag(:)
         end do
         
         ! Shift the radial velocity vectors to align with the center of mass of the collisional system (the origin)
         frag%ke_orbit = 0.0_DP
         frag%vb(:,1:nfrag) = fraggle_util_vmag_to_vb(frag%v_r_mag(1:nfrag), frag%v_r_unit(:,1:nfrag), frag%v_t_mag(1:nfrag), frag%v_t_unit(:,1:nfrag), frag%mass(1:nfrag), frag%vbcom(:)) 
      do i = 1, nfrag
         frag%v_coll(:, i) = frag%vb(:, i) - frag%vbcom(:)
         frag%ke_orbit = frag%ke_orbit + frag%mass(i) * dot_product(frag%vb(:, i), frag%vb(:, i))
      end do
      frag%ke_orbit = 0.5_DP * frag%ke_orbit

      lfailure = abs((frag%ke_budget - (frag%ke_orbit + frag%ke_spin)) / frag%ke_budget) > FRAGGLE_ETOL
      if (lfailure) then
         call fraggle_io_log_one_message(" ")
         call fraggle_io_log_one_message("Radial velocity failure diagnostics")
         write(message, *) frag%ke_budget
         call fraggle_io_log_one_message("ke_budget     : " // trim(adjustl(message)))
         write(message, *) frag%ke_spin
         call fraggle_io_log_one_message("ke_spin       : " // trim(adjustl(message)))
         write(message, *) frag%ke_orbit
         call fraggle_io_log_one_message("ke_orbit : " // trim(adjustl(message)))
         write(message, *) frag%ke_budget - (frag%ke_orbit + frag%ke_spin)
         call fraggle_io_log_one_message("ke_remainder  : " // trim(adjustl(message)))
      end if

      end associate
      return

      contains
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
            integer(I4B)                          :: i
            real(DP), dimension(:,:), allocatable :: v_shift
            real(DP), dimension(frag%nbody)       :: kearr
            real(DP)                              :: keo, ke_radial
      
            associate(nfrag => frag%nbody)
               allocate(v_shift, mold=frag%vb)
               v_shift(:,:) = fraggle_util_vmag_to_vb(v_r_mag_input, frag%v_r_unit, frag%v_t_mag, frag%v_t_unit, frag%mass, frag%vbcom) 
               do concurrent(i = 1:nfrag)
                  kearr(i) = frag%mass(i) * (frag%Ip(3, i) * frag%radius(i)**2 * dot_product(frag%rot(:, i), frag%rot(:, i)) + dot_product(v_shift(:, i), v_shift(:, i)))
               end do
               keo = 2 * frag%ke_budget - sum(kearr(:))
               ke_radial = frag%ke_budget - frag%ke_orbit - frag%ke_spin
               ! The following ensures that fval = 0 is a local minimum, which is what the BFGS method is searching for
               fval = (keo / (2 * ke_radial))**2
            end associate
      
            return
         end function radial_objective_function

      end subroutine fraggle_generate_rad_vel

end submodule s_fraggle_generate