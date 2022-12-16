!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle_classes) s_fraggle_generate
   use swiftest

   integer(I4B), parameter :: NFRAG_MIN = 7 !! The minimum allowable number of fragments (set to 6 because that's how many unknowns are needed in the tangential velocity calculation)
   real(DP),     parameter :: F_SPIN_FIRST = 0.05_DP !! The initial try value of the fraction of energy or momenum in spin (whichever has the lowest kinetic energy)
   real(DP), parameter     :: FRAGGLE_LTOL = 10 * epsilon(1.0_DP)
   real(DP), parameter     :: FRAGGLE_ETOL = 1e-8_DP

contains

   module subroutine fraggle_generate_fragments(self, system, param, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Generates a system of fragments in barycentric coordinates that conserves energy and momentum.
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(fraggle_system),        intent(inout) :: self      !! Fraggle system object the outputs will be the fragmentation 
      class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param     !! Current run configuration parameters 
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
      real(DP), dimension(NDIM)            :: runit, vunit

      ! The minimization and linear solvers can sometimes lead to floating point exceptions. Rather than halting the code entirely if this occurs, we
      ! can simply fail the attempt and try again. So we need to turn off any floating point exception halting modes temporarily 
      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      fpe_quiet_modes(:) = .false.
      call ieee_set_halting_mode(IEEE_ALL,fpe_quiet_modes)

      select type(fragments => self%fragments)
      class is (fraggle_fragments)
      associate(collision => self, impactors => self%impactors, nfrag => fragments%nbody, pl => system%pl)

         write(message,*) nfrag
         call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle generating " // trim(adjustl(message)) // " fragments.")
         if (nfrag < NFRAG_MIN) then
            write(message,*) "Fraggle needs at least ",NFRAG_MIN," fragments, but only ",nfrag," were given."
            call io_log_one_message(FRAGGLE_LOG_OUT, message)
            lfailure = .true.
            return
         end if


         ! This is a factor that will "distort" the shape of the frgment cloud in the direction of the impact velocity 
         f_spin= .mag. (runit(:) .cross. vunit(:)) 

         if (param%lflatten_interactions) then
            lk_plpl = allocated(pl%k_plpl)
            if (lk_plpl) deallocate(pl%k_plpl)
         else 
            lk_plpl = .false.
         end if

         call fragments%set_natural_scale(impactors)

         call fragments%reset()

         ! Calculate the initial energy of the system without the collisional family
         call collision%get_energy_and_momentum(system, param, lbefore=.true.)
        
         ! Start out the fragments close to the initial separation distance. This will be increased if there is any overlap or we fail to find a solution
         r_max_start = 1.2_DP * .mag.(impactors%rb(:,2) - impactors%rb(:,1))
         lfailure = .false.
         try = 1
         do while (try < MAXTRY)
            write(message,*) try
            call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle try " // trim(adjustl(message)))
            if (lfailure) then
               call fragments%restructure(impactors, try, f_spin, r_max_start)
               call fragments%reset()
               try = try + 1
            end if

            lfailure = .false.
            call ieee_set_flag(ieee_all, .false.) ! Set all fpe flags to quiet

            call fraggle_generate_pos_vec(fragments, impactors, r_max_start)
            call fragments%set_coordinate_system(impactors)

            ! Initial velocity guess will be the barycentric velocity of the colliding system so that the budgets are based on the much smaller collisional-frame velocities
            do concurrent (i = 1:nfrag)
               fragments%vb(:, i) = fragments%vbcom(:)
            end do

            call fragments%get_energy_and_momentum(impactors, system, param, lbefore=.false.)
            call fragments%set_budgets()

            call fraggle_generate_spins(fragments, impactors, f_spin, lfailure)
            if (lfailure) then
               call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle failed to find spins")
               cycle
            end if

            call fraggle_generate_tan_vel(fragments, impactors, lfailure)
            if (lfailure) then
               call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle failed to find tangential velocities")
               cycle
            end if

            call fraggle_generate_rad_vel(fragments, impactors, lfailure)
            if (lfailure) then
               call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle failed to find radial velocities")
               cycle
            end if

            call fragments%get_energy_and_momentum(impactors, system, param, lbefore=.false.)
            dEtot = fragments%Etot_after - fragments%Etot_before 
            dLmag = .mag. (fragments%Ltot_after(:) - fragments%Ltot_before(:))
            exit

            lfailure = ((abs(dEtot + impactors%Qloss) > FRAGGLE_ETOL) .or. (dEtot > 0.0_DP)) 
            if (lfailure) then
               write(message, *) dEtot, abs(dEtot + impactors%Qloss) / FRAGGLE_ETOL
               call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle failed due to high energy error: " // &
                                                        trim(adjustl(message)))
               cycle
            end if

            lfailure = ((abs(dLmag) / (.mag.fragments%Ltot_before)) > FRAGGLE_LTOL) 
            if (lfailure) then
               write(message,*) dLmag / (.mag.fragments%Ltot_before(:))
               call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle failed due to high angular momentum error: " // &
                                                        trim(adjustl(message)))
               cycle
            end if

            ! Check if any of the usual floating point exceptions happened, and fail the try if so
            call ieee_get_flag(ieee_usual, fpe_flag)
            lfailure = any(fpe_flag) 
            if (.not.lfailure) exit
            write(message,*) "Fraggle failed due to a floating point exception: ", fpe_flag
            call io_log_one_message(FRAGGLE_LOG_OUT, message)

         end do
         write(message,*) try
         if (lfailure) then
            call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle fragment generation failed after " // &
                                                      trim(adjustl(message)) // " tries")
         else
            call io_log_one_message(FRAGGLE_LOG_OUT, "Fraggle fragment generation succeeded after " // &
                                                       trim(adjustl(message)) // " tries")
         end if

         call fragments%set_original_scale(impactors)

         ! Restore the big array
         if (lk_plpl) call pl%flatten(param)
      end associate
      end select
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily

      return 
   end subroutine fraggle_generate_fragments


   subroutine fraggle_generate_pos_vec(fragments, impactors, r_max_start)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the orbits of the fragments around the center of mass. The fragments are initially placed on a plane defined by the 
      !! pre-impact angular momentum. They are distributed on an ellipse surrounding the center of mass.
      !! The initial positions do not conserve energy or momentum, so these need to be adjusted later.
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: fragments        !! Fraggle fragment system object
      class(collision_impactors), intent(inout) :: impactors   !! Fraggle collider system object
      real(DP),                 intent(in)    :: r_max_start !! Initial guess for the starting maximum radial distance of fragments
      ! Internals
      real(DP)  :: dis, rad, r_max, fdistort
      real(DP), dimension(NDIM) :: runit, vunit
      logical, dimension(:), allocatable :: loverlap
      integer(I4B) :: i, j
      logical :: lnoncat, lhitandrun

      associate(nfrag => fragments%nbody)
         allocate(loverlap(nfrag))

         lnoncat = (impactors%regime /= COLLRESOLVE_REGIME_SUPERCATASTROPHIC) ! For non-catastrophic impacts, make the fragments act like ejecta and point away from the impact point
         lhitandrun = (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) ! Disruptive hit and runs have their own fragment distribution

         ! Place the fragments into a region that is big enough that we should usually not have overlapping bodies
         ! An overlapping bodies will collide in the next time step, so it's not a major problem if they do (it just slows the run down)
         r_max = r_max_start
         rad = sum(impactors%radius(:))


         ! Get the unit vectors for the relative position and velocity vectors. These are used to shift the fragment cloud depending on the 
         runit(:) = impactors%rb(:,2) - impactors%rb(:,1) 
         runit(:) = runit(:) / (.mag. runit(:))

         vunit(:) = impactors%vb(:,2) - impactors%vb(:,1) 
         vunit(:) = vunit(:) / (.mag. vunit(:))

         ! This is a factor that will "distort" the shape of the frgment cloud in the direction of the impact velocity 
         fdistort = .mag. (runit(:) .cross. vunit(:)) 

         ! We will treat the first two fragments of the list as special cases. They get initialized the maximum distances apart along the original impactor distance vector.
         ! This is done because in a regular disruption, the first body is the largest, the second the second largest, and the rest are smaller equal-mass fragments.
         call random_number(fragments%rc(:,3:nfrag))
         loverlap(:) = .true.
         do while (any(loverlap(3:nfrag)))
            fragments%rc(:, 1) = impactors%rb(:, 1) - fragments%rbcom(:) 
            fragments%rc(:, 2) = impactors%rb(:, 2) - fragments%rbcom(:)
            r_max = r_max + 0.1_DP * rad
            do i = 3, nfrag
               if (loverlap(i)) then
                  call random_number(fragments%rc(:,i))
                  fragments%rc(:,i) = 2 * (fragments%rc(:, i) - 0.5_DP)  
                  fragments%rc(:, i) = fragments%rc(:,i) + fdistort * vunit(:) 
                  fragments%rc(:, i) = r_max * fragments%rc(:, i) 
                  fragments%rc(:, i) = fragments%rc(:, i) + (fragments%rbimp(:) - fragments%rbcom(:)) ! Shift the center of the fragment cloud to the impact point rather than the CoM
                  !if (lnoncat .and. dot_product(fragments%rc(:,i), runit(:)) < 0.0_DP) fragments%rc(:, i) = -fragments%rc(:, i) ! Make sure the fragment cloud points away from the impact point
               end if
            end do
            loverlap(:) = .false.
            do j = 1, nfrag
               do i = j + 1, nfrag
                  dis = .mag.(fragments%rc(:,j) - fragments%rc(:,i))
                  loverlap(i) = loverlap(i) .or. (dis <= (fragments%radius(i) + fragments%radius(j))) 
               end do
            end do
         end do
         call fraggle_util_shift_vector_to_origin(fragments%mass, fragments%rc)
         call fragments%set_coordinate_system(impactors)

         do concurrent(i = 1:nfrag)
            fragments%rb(:,i) = fragments%rc(:,i) + fragments%rbcom(:)
         end do

         fragments%rbcom(:) = 0.0_DP
         do concurrent(i = 1:nfrag)
            fragments%rbcom(:) = fragments%rbcom(:) + fragments%mass(i) * fragments%rb(:,i) 
         end do
         fragments%rbcom(:) = fragments%rbcom(:) / fragments%mtot
      end associate

      return
   end subroutine fraggle_generate_pos_vec


   subroutine fraggle_generate_spins(fragments, impactors, f_spin, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Calculates the spins of a collection of fragments such that they conserve angular momentum without blowing the fragment kinetic energy budget.
      !!
      !! A failure will trigger a restructuring of the fragments so we will try new values of the radial position distribution.
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: fragments      !! Fraggle fragment system object
      class(collision_impactors), intent(inout) :: impactors   !! Fraggle collider system object
      real(DP),                 intent(in)    :: f_spin    !! Fraction of energy or momentum that goes into spin (whichever gives the lowest kinetic energy)
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds! 
      ! Internals
      real(DP), dimension(NDIM) :: L_remainder, rot_L, rot_ke, L
      real(DP), dimension(NDIM,fragments%nbody) :: frot_rand ! The random rotation factor applied to fragments
      real(DP), parameter :: frot_rand_mag = 1.50_DP ! The magnitude of the rotation variation to apply to the fragments
      integer(I4B)              :: i
      character(len=STRMAX)     :: message
      real(DP) :: ke_remainder, ke

      associate(nfrag => fragments%nbody)
         lfailure = .false.
         L_remainder(:) = fragments%L_budget(:)
         ke_remainder = fragments%ke_budget

         ! Add a fraction of the orbit angular momentum of the second body to the spin of the first body to kick things off
         L(:) = impactors%Lspin(:,1) + f_spin * (impactors%Lorbit(:,2) + impactors%Lspin(:,2))
         fragments%rot(:,1) = L(:) / (fragments%mass(1) * fragments%radius(1)**2 * fragments%Ip(3,1))
         L_remainder(:) = L_remainder(:) - L(:)

         ! Partition the spin momentum of the second body into the spin of the fragments, with some random variation
         L(:) = impactors%Lspin(:,2) / (nfrag - 1) 

         call random_number(frot_rand(:,2:nfrag))
         frot_rand(:,2:nfrag) = 2 * (frot_rand(:,2:nfrag) - 0.5_DP) * frot_rand_mag

         do i = 2, nfrag
            rot_L(:) = L(:) + frot_rand(:,i) * .mag.L(:)
            fragments%rot(:,i) = rot_L(:) / (fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,1))
            L_remainder(:) = L_remainder(:) - rot_L(:)
         end do

         ! Make sure we didn't blow our kinetic energy budget
         do i = 1, nfrag
            ke_remainder = ke_remainder - 0.5_DP * fragments%mass(i) * fragments%Ip(3, i) * fragments%radius(i)**2 * .mag.(fragments%rot(:, i))
         end do

         ! Distributed most of the remaining angular momentum amongst all the particles
         fragments%ke_spin = 0.0_DP
         if (.mag.(L_remainder(:)) > FRAGGLE_LTOL) then 
            do i = nfrag, 1, -1
               ! Convert a fraction (f_spin) of either the remaining angular momentum or kinetic energy budget into spin, whichever gives the smaller rotation so as not to blow any budgets
               rot_ke(:) = sqrt(2 * f_spin * ke_remainder / (i * fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3, i))) * L_remainder(:) / .mag.(L_remainder(:))
               rot_L(:) = f_spin * L_remainder(:) / (i * fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3, i))
               if (.mag.(rot_ke) < .mag.(rot_L)) then
                  fragments%rot(:,i) = fragments%rot(:,i) + rot_ke(:)
               else
                  fragments%rot(:, i) = fragments%rot(:,i) + rot_L(:)
               end if
               ke = 0.5_DP * fragments%mass(i) * fragments%Ip(3, i) * fragments%radius(i)**2 * norm2(fragments%rot(:, i))
               L(:) = fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3, i) * fragments%rot(:, i)
               ke_remainder = ke_remainder - ke
               L_remainder(:) = L_remainder(:) - L(:)
               fragments%ke_spin = fragments%ke_spin + ke
            end do
         end if

         lfailure = ((fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit) < 0.0_DP)

         if (lfailure) then
            call io_log_one_message(FRAGGLE_LOG_OUT, " ")
            call io_log_one_message(FRAGGLE_LOG_OUT, "Spin failure diagnostics")
            write(message, *) fragments%ke_budget
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_budget     : " // trim(adjustl(message)))
            write(message, *) fragments%ke_spin
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_spin       : " // trim(adjustl(message)))
            write(message, *) fragments%ke_orbit
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_orbit      : " // trim(adjustl(message)))
            write(message, *) fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_remainder  : " // trim(adjustl(message)))
         end if

      end associate

      return
   end subroutine fraggle_generate_spins


   subroutine fraggle_generate_tan_vel(fragments, impactors, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Adjusts the tangential velocities and spins of a collection of fragments such that they conserve angular momentum without blowing the fragment kinetic energy budget.
      !! This procedure works in several stages, with a goal to solve the angular and linear momentum constraints on the fragments, while still leaving a positive balance of
      !! our fragment kinetic energy (fragments%ke_budget) that we can put into the radial velocity distribution.
      !!
      !! The first thing we'll try to do is solve for the tangential velocities of the first 6 fragments, using angular and linear momentum as constraints and an initial
      !! tangential velocity distribution for the remaining bodies (if there are any) that distributes their angular momentum equally between them.
      !! If that doesn't work and we blow our kinetic energy budget, we will attempt to find a tangential velocity distribution that minimizes the kinetic energy while
      !! conserving momentum. 
      !!
      !! A failure will trigger a restructuring of the fragments so we will try new values of the radial position distribution.
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: fragments      !! Fraggle fragment system object
      class(collision_impactors), intent(inout) :: impactors   !! Fraggle collider system object
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds
      ! Internals
      integer(I4B) :: i, try
      real(DP), parameter                  :: TOL_MIN = 1e0_DP ! This doesn't have to be very accurate, as we really just want a tangential velocity distribution with less kinetic energy than our initial guess.
      real(DP), parameter                  :: TOL_INIT = 1e-12_DP
      real(DP), parameter                  :: VNOISE_MAG = 1e-1_DP !! Magnitude of the noise to apply to initial conditions to help minimizer find a solution in case of failure
      integer(I4B), parameter              :: MAXLOOP = 10
      integer(I4B), parameter              :: MAXTRY = 100
      real(DP)                             :: tol
      real(DP), dimension(:), allocatable  :: v_t_initial, v_t_output
      real(DP), dimension(fragments%nbody) :: kefrag, vnoise
      type(lambda_obj_err)                 :: objective_function
      real(DP), dimension(NDIM)            :: L_frag_tot
      character(len=STRMAX)                :: message
      real(DP)                             :: ke_diff

      associate(nfrag => fragments%nbody)
         lfailure = .false.

         allocate(v_t_initial, mold=fragments%v_t_mag)
         do try = 1, MAXTRY
            v_t_initial(1) = dot_product(impactors%vb(:,1),fragments%v_t_unit(:,1))
            do i = 2, nfrag
               v_t_initial(i) = dot_product(impactors%vb(:,2), fragments%v_t_unit(:,i)) 
            end do
            fragments%v_t_mag(:) = v_t_initial

            ! Find the local kinetic energy minimum for the system that conserves linear and angular momentum
            objective_function = lambda_obj(tangential_objective_function, lfailure)

            tol = TOL_INIT
            do while(tol < TOL_MIN)
               call util_minimize_bfgs(objective_function, nfrag-6, v_t_initial(7:nfrag), tol, MAXLOOP, lfailure, v_t_output)
               fragments%v_t_mag(7:nfrag) = v_t_output(:)
               ! Now that the KE-minimized values of the i>6 fragments are found, calculate the momentum-conserving solution for tangential velociteis
               v_t_initial(7:nfrag) = fragments%v_t_mag(7:nfrag)
               if (.not.lfailure) exit
               tol = tol * 2_DP ! Keep increasing the tolerance until we converge on a solution
               call random_number(vnoise(1:nfrag)) ! Adding a bit of noise to the initial conditions helps it find a solution more often
               vnoise(:) = 1.0_DP + VNOISE_MAG * (2 * vnoise(:) - 1._DP)
               v_t_initial(:) = v_t_initial(:) * vnoise(:)
            end do
            fragments%v_t_mag(1:nfrag) = solve_fragment_tan_vel(v_t_mag_input=v_t_initial(7:nfrag), lfailure=lfailure)

            ! Perform one final shift of the radial velocity vectors to align with the center of mass of the collisional system (the origin)
            fragments%vb(:,1:nfrag) = fraggle_util_vmag_to_vb(fragments%v_r_mag(1:nfrag), fragments%v_r_unit(:,1:nfrag), fragments%v_t_mag(1:nfrag), &
                                                         fragments%v_t_unit(:,1:nfrag), fragments%mass(1:nfrag), fragments%vbcom(:)) 
            do concurrent (i = 1:nfrag)
               fragments%vc(:,i) = fragments%vb(:,i) - fragments%vbcom(:)
            end do

            ! Now do a kinetic energy budget check to make sure we are still within the budget.
            kefrag = 0.0_DP
            do concurrent(i = 1:nfrag)
               kefrag(i) = fragments%mass(i) * dot_product(fragments%vb(:, i), fragments%vb(:, i))
            end do
            fragments%ke_orbit = 0.5_DP * sum(kefrag(:))

            ! If we are over the energy budget, flag this as a failure so we can try again
            ke_diff = fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit
            lfailure = ke_diff < 0.0_DP
            if (.not.lfailure) exit
            fragments%rc(:,:) = fragments%rc(:,:) * 1.1_DP
         end do
         if (lfailure) then
            call io_log_one_message(FRAGGLE_LOG_OUT, " ")
            call io_log_one_message(FRAGGLE_LOG_OUT, "Tangential velocity failure diagnostics")
            call fragments%get_angular_momentum()
            L_frag_tot = fragments%Lspin(:) + fragments%Lorbit(:)
            write(message, *) .mag.(fragments%L_budget(:) - L_frag_tot(:)) / (.mag.fragments%Ltot_before(:))
            call io_log_one_message(FRAGGLE_LOG_OUT, "|L_remainder| : " // trim(adjustl(message)))
            write(message, *) fragments%ke_budget
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_budget     : " // trim(adjustl(message)))
            write(message, *) fragments%ke_spin
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_spin       : " // trim(adjustl(message)))
            write(message, *) fragments%ke_orbit
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_tangential : " // trim(adjustl(message)))
            write(message, *) fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_radial     : " // trim(adjustl(message)))
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
      
            associate(nfrag => fragments%nbody)
               lfailure = .false.
               ! We have 6 constraint equations (2 vector constraints in 3 dimensions each)
               ! The first 3 are that the linear momentum of the fragments is zero with respect to the collisional barycenter
               ! The second 3 are that the sum of the angular momentum of the fragments is conserved from the pre-impact state
               L_lin_others(:) = 0.0_DP
               L_orb_others(:) = 0.0_DP
               do i = 1, nfrag
                  if (i <= 2 * NDIM) then ! The tangential velocities of the first set of bodies will be the unknowns we will solve for to satisfy the constraints
                     A(1:3, i) = fragments%mass(i) * fragments%v_t_unit(:, i) 
                     A(4:6, i) = fragments%mass(i) * fragments%rmag(i) * (fragments%v_r_unit(:, i) .cross. fragments%v_t_unit(:, i))
                  else if (present(v_t_mag_input)) then
                     vtmp(:) = v_t_mag_input(i - 6) * fragments%v_t_unit(:, i)
                     L_lin_others(:) = L_lin_others(:) + fragments%mass(i) * vtmp(:)
                     L(:) = fragments%mass(i) * (fragments%rc(:, i) .cross. vtmp(:)) 
                     L_orb_others(:) = L_orb_others(:) + L(:)
                  end if
               end do
               b(1:3) = -L_lin_others(:)
               b(4:6) = fragments%L_budget(:) - fragments%Lspin(:) - L_orb_others(:)
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
            real(DP), dimension(NDIM,fragments%nbody) :: v_shift
            real(DP), dimension(fragments%nbody) :: v_t_new, kearr
            real(DP) :: keo
      
            associate(nfrag => fragments%nbody)
               lfailure = .false.
         
               v_t_new(:) = solve_fragment_tan_vel(v_t_mag_input=v_t_mag_input(:), lfailure=lfailure)
               v_shift(:,:) = fraggle_util_vmag_to_vb(fragments%v_r_mag, fragments%v_r_unit, v_t_new, fragments%v_t_unit, fragments%mass, fragments%vbcom) 
         
               kearr = 0.0_DP
               do concurrent(i = 1:nfrag)
                  kearr(i) = fragments%mass(i) * dot_product(v_shift(:, i), v_shift(:, i))
               end do
               keo = 0.5_DP * sum(kearr(:))
               fval = keo 
               lfailure = .false.
            end associate

            return
         end function tangential_objective_function

   end subroutine fraggle_generate_tan_vel


   subroutine fraggle_generate_rad_vel(fragments, impactors, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! 
      !! Adjust the fragment velocities to set the fragment orbital kinetic energy. This will minimize the difference between the fragment kinetic energy and the energy budget
      implicit none
      ! Arguments
      class(fraggle_fragments), intent(inout) :: fragments      !! Fraggle fragment system object
      class(collision_impactors), intent(inout) :: impactors   !! Fraggle collider system object
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds! 
      ! Internals
      real(DP), parameter                   :: TOL_MIN = FRAGGLE_ETOL   ! This needs to be more accurate than the tangential step, as we are trying to minimize the total residual energy
      real(DP), parameter                   :: TOL_INIT = 1e-14_DP
      real(DP), parameter                   :: VNOISE_MAG = 1e-10_DP !! Magnitude of the noise to apply to initial conditions to help minimizer find a solution in case of failure
      integer(I4B), parameter               :: MAXLOOP = 100
      real(DP)                              :: ke_radial, tol 
      integer(I4B)                          :: i
      real(DP), dimension(:), allocatable   :: v_r_initial
      real(DP), dimension(fragments%nbody)       :: vnoise
      type(lambda_obj)                      :: objective_function
      character(len=STRMAX)                 :: message

      associate(nfrag => fragments%nbody)
         ! Set the "target" ke for the radial component
         
         allocate(v_r_initial, source=fragments%v_r_mag)
         ! Initialize radial velocity magnitudes with a random value that related to equipartition of kinetic energy with some noise and scaled with respect to the initial distance
         v_r_initial(1) = dot_product(impactors%vb(:,1),fragments%v_r_unit(:,1))
         fragments%ke_orbit = 0.5_DP * fragments%mass(1) * v_r_initial(1)**2

         ke_radial = fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit
         call random_number(vnoise(1:nfrag))
         vnoise(:) = 1.0_DP + VNOISE_MAG * (2 * vnoise(:) - 1.0_DP) 
         v_r_initial(2:nfrag) = -sqrt(abs(2 * ke_radial) / (fragments%mass(1:nfrag) * nfrag)) * vnoise(1:nfrag)

         ! Initialize the lambda function using a structure constructor that calls the init method
         ! Minimize the ke objective function using the BFGS optimizer
         objective_function = lambda_obj(radial_objective_function)
         tol = TOL_INIT
         do while(tol < TOL_MIN)
            call util_minimize_bfgs(objective_function, nfrag, v_r_initial, tol, MAXLOOP, lfailure, fragments%v_r_mag)
            if (.not.lfailure) exit
            tol = tol * 2 ! Keep increasing the tolerance until we converge on a solution
            v_r_initial(:) = fragments%v_r_mag(:)
            call random_number(vnoise(1:nfrag)) ! Adding a bit of noise to the initial conditions helps it find a solution more often
            vnoise(:) = 1.0_DP + VNOISE_MAG * (2 * vnoise(:) - 1._DP)
            v_r_initial(:) = v_r_initial(:) * vnoise(:)
         end do
         
         ! Shift the radial velocity vectors to align with the center of mass of the collisional system (the origin)
         fragments%ke_orbit = 0.0_DP
         fragments%ke_spin = 0.0_DP
         fragments%vb(:,1:nfrag) = fraggle_util_vmag_to_vb(fragments%v_r_mag(1:nfrag), fragments%v_r_unit(:,1:nfrag), &
                              fragments%v_t_mag(1:nfrag), fragments%v_t_unit(:,1:nfrag), fragments%mass(1:nfrag), fragments%vbcom(:)) 
         do i = 1, nfrag
            fragments%vc(:, i) = fragments%vb(:, i) - fragments%vbcom(:)
            fragments%ke_orbit = fragments%ke_orbit + fragments%mass(i) * norm2(fragments%vb(:, i))
            fragments%ke_spin = fragments%ke_spin + fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i) * norm2(fragments%rot(:,i))
         end do
         fragments%ke_orbit = 0.5_DP * fragments%ke_orbit
         fragments%ke_spin = 0.5_DP * fragments%ke_spin

         lfailure = abs((fragments%ke_budget - (fragments%ke_orbit + fragments%ke_spin)) / fragments%ke_budget) > FRAGGLE_ETOL
         if (lfailure) then
            call io_log_one_message(FRAGGLE_LOG_OUT, " ")
            call io_log_one_message(FRAGGLE_LOG_OUT, "Radial velocity failure diagnostics")
            write(message, *) fragments%ke_budget
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_budget     : " // trim(adjustl(message)))
            write(message, *) fragments%ke_spin
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_spin       : " // trim(adjustl(message)))
            write(message, *) fragments%ke_orbit
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_orbit : " // trim(adjustl(message)))
            write(message, *) fragments%ke_budget - (fragments%ke_orbit + fragments%ke_spin)
            call io_log_one_message(FRAGGLE_LOG_OUT, "ke_remainder  : " // trim(adjustl(message)))
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
            real(DP), dimension(fragments%nbody)       :: kearr
            real(DP)                              :: keo, ke_radial, rotmag2, vmag2
      
            allocate(v_shift, mold=fragments%vb)
            v_shift(:,:) = fraggle_util_vmag_to_vb(v_r_mag_input, fragments%v_r_unit, fragments%v_t_mag, fragments%v_t_unit, fragments%mass, fragments%vbcom) 
            !$omp do simd firstprivate(fragments)
            do i = 1,fragments%nbody
               rotmag2 = fragments%rot(1,i)**2 + fragments%rot(2,i)**2 + fragments%rot(3,i)**2
               vmag2 = v_shift(1,i)**2 + v_shift(2,i)**2 + v_shift(3,i)**2
               kearr(i) = fragments%mass(i) * (fragments%Ip(3, i) * fragments%radius(i)**2 * rotmag2 + vmag2) 
            end do
            !$omp end do simd
            keo = 2 * fragments%ke_budget - sum(kearr(:))
            ke_radial = fragments%ke_budget - fragments%ke_orbit - fragments%ke_spin
            ! The following ensures that fval = 0 is a local minimum, which is what the BFGS method is searching for
            fval = (keo / (2 * ke_radial))**2
      
            return
         end function radial_objective_function

      end subroutine fraggle_generate_rad_vel

end submodule s_fraggle_generate
