!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle) s_fraggle_generate
   use swiftest
   use symba

   integer(I4B), parameter :: NFRAG_MIN = 7 !! The minimum allowable number of fragments (set to 6 because that's how many unknowns are needed in the tangential velocity calculation)
   real(DP), parameter     :: FRAGGLE_LTOL = 1e-4_DP !10 * epsilon(1.0_DP)
   real(DP), parameter     :: FRAGGLE_ETOL = 1e-1_DP

contains

   module subroutine fraggle_generate_disrupt(self, nbody_system, param, t, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Generates a nbody_system of fragments in barycentric coordinates that conserves energy and momentum.
      use, intrinsic :: ieee_exceptions
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self         !! Fraggle system object the outputs will be the fragmentation 
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters 
      real(DP),                 intent(in)    :: t            !! Time of collision 
      logical, optional,        intent(out)   :: lfailure     !! Answers the question: Should this have been a merger instead?
       ! Internals
      integer(I4B)                         :: try
      real(DP)                             :: dEtot, dLmag
      integer(I4B), parameter              :: MAXTRY = 10
      logical                              :: lk_plpl, lfailure_local
      logical, dimension(size(IEEE_ALL))   :: fpe_halting_modes, fpe_quiet_modes
      logical, dimension(size(IEEE_USUAL)) :: fpe_flag 
      character(len=STRMAX)                :: message

      ! The minimization and linear solvers can sometimes lead to floating point exceptions. Rather than halting the code entirely if this occurs, we
      ! can simply fail the attempt and try again. So we need to turn off any floating point exception halting modes temporarily 
      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      fpe_quiet_modes(:) = .false.
      call ieee_set_halting_mode(IEEE_ALL,fpe_quiet_modes)

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(param)
      class is (swiftest_parameters)
      select type(fragments => self%fragments)
      class is (fraggle_fragments(*))
      associate(impactors => self%impactors, nfrag => fragments%nbody, pl => nbody_system%pl)

         write(message,*) nfrag
         call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle generating " // trim(adjustl(message)) // " fragments.")
         if (nfrag < NFRAG_MIN) then
            write(message,*) "Fraggle needs at least ",NFRAG_MIN," fragments, but only ",nfrag," were given."
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
            lfailure_local = .true.
            return
         end if

         if (param%lflatten_interactions) then
            lk_plpl = allocated(pl%k_plpl)
            if (lk_plpl) deallocate(pl%k_plpl)
         else 
            lk_plpl = .false.
         end if
         call ieee_set_flag(ieee_all, .false.) ! Set all fpe flags to quiet

         call self%set_natural_scale()

         call fragments%reset()

         lfailure_local = .false.
         try = 1

         do while (try < MAXTRY)
            write(message,*) try
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle try " // trim(adjustl(message)))
            if (lfailure_local) then
               call fragments%reset()
               try = try + 1
            end if
            lfailure_local = .false.

            ! Use the disruption collision model to generate initial conditions
            ! Compute the "before" energy/momentum and the budgets
            call self%get_energy_and_momentum(nbody_system, param, lbefore=.true.)
            call self%collision_disruption%disrupt(nbody_system, param, t)
            call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)
            call self%set_budgets()

            ! Minimize difference between energy/momentum and budgets
            ! call fraggle_generate_minimize(self, lfailure_local)
            ! call fraggle_generate_tan_vel(self, lfailure_local)
            call fraggle_generate_rad_vel(self, lfailure_local)

            call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)
            dEtot = self%Etot(2) - self%Etot(1)
            dLmag = .mag. (self%Ltot(:,2) - self%Ltot(:,1))

            lfailure_local = (dEtot > 0.0_DP) 
            if (lfailure_local) then
               write(message, *) "dEtot: ",dEtot, "dEtot/Qloss", dEtot / impactors%Qloss
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle failed due to energy gain: " // &
                                                        trim(adjustl(message)))
               !cycle
               lfailure_local = .false.
               exit
            end if

            lfailure_local = ((abs(dLmag) / (.mag.self%Ltot(:,1))) > FRAGGLE_LTOL) 
            if (lfailure_local) then
               write(message,*) dLmag / (.mag.self%Ltot(:,1))
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle failed due to high angular momentum error: " // &
                                                        trim(adjustl(message)))
               cycle
            end if

            ! Check if any of the usual floating point exceptions happened, and fail the try if so
            call ieee_get_flag(ieee_usual, fpe_flag)
            lfailure_local = any(fpe_flag) 
            if (.not.lfailure_local) exit
            write(message,*) "Fraggle failed due to a floating point exception: ", fpe_flag
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)

         end do
         write(message,*) try
         if (lfailure_local) then
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle fragment generation failed after " // &
                                                      trim(adjustl(message)) // " tries")
         else
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle fragment generation succeeded after " // &
                                                       trim(adjustl(message)) // " tries")
         end if

         call self%set_original_scale()

         ! Restore the big array
         if (lk_plpl) call pl%flatten(param)
         if (present(lfailure)) lfailure = lfailure_local
      end associate
      end select
      end select
      end select
      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily

      return 
   end subroutine fraggle_generate_disrupt


   subroutine fraggle_generate_minimize(collider, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Minimizes the differences between the energy and momentum and the budget
      !!
      !! A failure will trigger a restructuring of the fragments so we will try a new random distribution
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: collider !! Fraggle collision system object
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds
      ! Internals
      real(DP), parameter :: TOL_INIT = 1e-5_DP
      integer(I4B), parameter :: MAXLOOP = 100
      real(DP), dimension(collider%fragments%nbody) :: input_v
      real(DP), dimension(:), allocatable :: output_v
      type(lambda_obj) :: Efunc
      real(DP) :: tol, fval
      integer(I4B) :: loop,i, nelem
      logical :: lfirst_Efunc

      associate(impactors => collider%impactors,  nfrag => collider%fragments%nbody)
         select type(fragments => collider%fragments)
         class is (fraggle_fragments(*))

            nelem = nfrag 
            lfailure = .false.
            ! Find the local kinetic energy minimum for the nbody_system that conserves linear and angular momentum
            Efunc = lambda_obj(E_objective_function)
            tol = TOL_INIT

            fragments%r_unit(:,:) = .unit. fragments%vc(:,:)
            fragments%vmag(:) = .mag. fragments%vc(:,1:nfrag) 
            input_v(:) = fragments%vmag(1:nfrag)
            lfirst_Efunc = .true.
            fval = E_objective_function(input_v)

            call minimize_bfgs(Efunc, nelem, input_v, tol, MAXLOOP, lfailure, output_v)
            fval = E_objective_function(output_v)
            fragments%vmag(1:nfrag) = output_v(1:nfrag)
            do concurrent(i=1:nfrag)
               fragments%vc(:,i) = abs(fragments%vmag(i)) * fragments%r_unit(:,i)
            end do

            ! Set spins in order to conserve angular momentum
            call fragments%set_spins() 
            call collision_util_shift_vector_to_origin(fragments%mass, fragments%vc)
         end select
      end associate

      return

      contains

         function E_objective_function(val_input) result(fval) 
            !! Author: David A. Minton
            !!
            !! Objective function for evaluating how close our fragment trajectories are to the energy budget
            implicit none
            ! Arguments
            real(DP), dimension(:), intent(in)  :: val_input !! Flattened velocity and rotation arrays
            ! Result
            real(DP)                            :: fval      !! The objective function result, which is the sum of the squares of the difference between the calculated fragment kinetic energy and the components of angular and linear momentum
                                                             !! Minimizing this brings us closer to our objective
            ! Internals
            integer(I4B) :: i
            type(fraggle_fragments(:)), allocatable :: tmp_frag
            real(DP) :: deltaE
            real(DP), save :: deltaE_scale = 1.0_DP

            associate(impactors => collider%impactors, nfrag => collider%fragments%nbody)
               select type(fragments => collider%fragments)
               class is (fraggle_fragments(*))
                  allocate(tmp_frag, source=fragments)
                  tmp_frag%vmag(1:nfrag) = val_input(1:nfrag)
                  do concurrent(i=1:nfrag)
                     tmp_frag%vc(:,i) = abs(tmp_frag%vmag(i)) * tmp_frag%r_unit(:,i)
                  end do

                  call collision_util_shift_vector_to_origin(tmp_frag%mass, tmp_frag%vc)
                  ! Set spins in order to conserve angular momentum
                  call collision_util_set_spins(tmp_frag)

                  ! Get the current kinetic energy of the system
                  call tmp_frag%get_kinetic_energy()
                  deltaE = (tmp_frag%ke_budget - (tmp_frag%ke_orbit + tmp_frag%ke_spin)) / (tmp_frag%ke_budget)

                  ! The result will be scaled so such that the initial deltaE is 1.0. This helps improve the success of the
                  ! minimizer, by keeping values from getting to large
                  if (lfirst_Efunc) then 
                     deltaE_scale = deltaE
                     deltaE = 1.0_DP
                     lfirst_Efunc = .false.
                  else
                     deltaE = deltaE/deltaE_scale
                  end if
                  fval = deltaE**2
                  deallocate(tmp_frag)
                  
               end select
            end associate
         
            return
         end function E_objective_function


   end subroutine fraggle_generate_minimize


   subroutine fraggle_generate_tan_vel(collider, lfailure)
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
      class(collision_fraggle), intent(inout) :: collider      !! Fraggle fragment system object
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds
      ! Internals
      integer(I4B) :: i
      real(DP), parameter                  :: TOL_MIN = 1e-1_DP ! This doesn't have to be very accurate, as we really just want a tangential velocity distribution with less kinetic energy than our initial guess.
      real(DP), parameter                  :: TOL_INIT = 1e-6_DP
      real(DP), parameter                  :: VNOISE_MAG = 1e-3_DP !! Magnitude of the noise to apply to initial conditions to help minimizer find a solution in case of failure
      integer(I4B), parameter              :: MAXLOOP = 10
      real(DP)                             :: tol, fval
      real(DP), dimension(:), allocatable  :: v_t_initial
      real(DP), dimension(collider%fragments%nbody)      :: kefrag, vnoise
      type(lambda_obj_err)                 :: objective_function
      real(DP), dimension(NDIM)            :: Li, L_remainder, L_frag_tot
      character(len=STRMAX)                :: message
      real(DP), dimension(:), allocatable :: v_t_output
      logical                             :: lfirst_func

      associate(impactors => collider%impactors,  nfrag => collider%fragments%nbody)
         select type(fragments => collider%fragments)
         class is (fraggle_fragments(*))
            lfailure = .false.

            ! Solve for the tangential component of the velocities that both conserves linear momentum and uses the remaining angular momentum not used in spin.
            ! This will be done using a linear solver that solves for the tangential velocities of the first 6 fragments, constrained by the linear and angular momentum vectors, 
            ! which is embedded in a non-linear minimizer that will adjust the tangential velocities of the remaining i>6 fragments to minimize kinetic energy for a given momentum solution
            ! The initial conditions fed to the minimizer for the fragments will be the remaining angular momentum distributed between the fragments.
            tol = TOL_INIT
            lfirst_func = .true.
            do i = 1, nfrag
               fragments%v_t_mag(i) = dot_product(fragments%vc(:,i), fragments%t_unit(:,i))
               fragments%v_r_mag(i) = dot_product(fragments%vc(:,i), fragments%r_unit(:,i))
            end do
            allocate(v_t_initial, source=fragments%v_t_mag)
            do while(tol < TOL_MIN)

               ! ! Find the local kinetic energy minimum for the system that conserves linear and angular momentum
               objective_function = lambda_obj(tangential_objective_function, lfailure)
               fval = tangential_objective_function(v_t_output(:), lfailure)

               call minimize_bfgs(objective_function, nfrag-6, v_t_initial(7:nfrag), tol, MAXLOOP, lfailure, v_t_output)
               fval = tangential_objective_function(v_t_initial(:), lfailure)
               fragments%v_t_mag(7:nfrag) = v_t_output(:)
               ! Now that the KE-minimized values of the i>6 fragments are found, calculate the momentum-conserving solution for tangential velociteis
               v_t_initial(7:nfrag) = fragments%v_t_mag(7:nfrag)

               fragments%v_t_mag(1:nfrag) = solve_fragment_tan_vel(v_t_mag_input=v_t_initial(7:nfrag), lfailure=lfailure)

               ! Perform one final shift of the radial velocity vectors to align with the center of mass of the collisional system (the origin)
               fragments%vb(:,1:nfrag) = fraggle_util_vmag_to_vb(fragments%v_r_mag(1:nfrag), fragments%r_unit(:,1:nfrag), fragments%v_t_mag(1:nfrag), &
                                                            fragments%t_unit(:,1:nfrag), fragments%mass(1:nfrag), impactors%vbcom(:)) 
               do concurrent (i = 1:nfrag)
                  fragments%vc(:,i) = fragments%vb(:,i) - impactors%vbcom(:)
               end do

               ! Now do a kinetic energy budget check to make sure we are still within the budget.
               kefrag = 0.0_DP
               do concurrent(i = 1:nfrag)
                  kefrag(i) = fragments%mass(i) * dot_product(fragments%vc(:,i), fragments%vc(:,i))
               end do
               fragments%ke_orbit = 0.5_DP * sum(kefrag(:))

               ! If we are over the energy budget, flag this as a failure so we can try again
               call fragments%get_angular_momentum()
               lfailure = ((fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit) < 0.0_DP)
               if (.not.lfailure) exit
               tol = tol * 2_DP ! Keep increasing the tolerance until we converge on a solution
               ! Reduce fragment spins to try to get a better solution
               fragments%rot(:,2:nfrag) = fragments%rot(:,2:nfrag) * 0.9_DP
            end do
            if (lfailure) then
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, " ")
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Tangential velocity failure diagnostics")
               call fragments%get_angular_momentum()
               L_frag_tot = fragments%Lspin(:) + fragments%Lorbit(:)
               write(message, *) .mag.(L_frag_tot(:)) / (.mag.fragments%L_budget(:))
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "|L_remainder| : " // trim(adjustl(message)))
               write(message, *) fragments%ke_budget
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_budget     : " // trim(adjustl(message)))
               write(message, *) fragments%ke_spin
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_spin       : " // trim(adjustl(message)))
               write(message, *) fragments%ke_orbit
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_tangential : " // trim(adjustl(message)))
               write(message, *) fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_radial     : " // trim(adjustl(message)))
            end if
         end select
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
      
            associate(impactors => collider%impactors,  nfrag => collider%fragments%nbody)
               select type(fragments => collider%fragments)
               class is (fraggle_fragments(*))
                  lfailure = .false.
                  ! We have 6 constraint equations (2 vector constraints in 3 dimensions each)
                  ! The first 3 are that the linear momentum of the fragments is zero with respect to the collisional barycenter
                  ! The second 3 are that the sum of the angular momentum of the fragments is conserved from the pre-impact state
                  L_lin_others(:) = 0.0_DP
                  L_orb_others(:) = 0.0_DP
                  do i = 1, nfrag
                     if (i <= 2 * NDIM) then ! The tangential velocities of the first set of bodies will be the unknowns we will solve for to satisfy the constraints
                        A(1:3, i) = fragments%mass(i) * fragments%t_unit(:, i) 
                        A(4:6, i) = fragments%mass(i) * fragments%rmag(i) * (fragments%r_unit(:, i) .cross. fragments%t_unit(:, i))
                     else if (present(v_t_mag_input)) then
                        vtmp(:) = v_t_mag_input(i - 6) * fragments%t_unit(:, i)
                        L_lin_others(:) = L_lin_others(:) + fragments%mass(i) * vtmp(:)
                        L(:) = fragments%mass(i) * (fragments%rc(:, i) .cross. vtmp(:)) 
                        L_orb_others(:) = L_orb_others(:) + L(:)
                     end if
                  end do
                  b(1:3) = -L_lin_others(:)
                  b(4:6) = fragments%L_budget(:) - fragments%Lspin(:) - L_orb_others(:)
                  allocate(v_t_mag_output(nfrag))
                  v_t_mag_output(1:6) = solve_linear_system(A, b, 6, lfailure)
                  if (present(v_t_mag_input)) v_t_mag_output(7:nfrag) = v_t_mag_input(:)
               end select
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
            real(DP), dimension(NDIM,collider%fragments%nbody) :: vc, vb
            real(DP), dimension(collider%fragments%nbody) :: v_t_new, kearr
            real(DP) :: keo
            real(DP), save :: fval_scale = 1.0_DP
      
            associate(impactors => collider%impactors,  nfrag => collider%fragments%nbody)
               select type(fragments => collider%fragments)
               class is (fraggle_fragments(*))
                  lfailure = .false.
            
                  v_t_new(:) = solve_fragment_tan_vel(v_t_mag_input=v_t_mag_input(:), lfailure=lfailure)
                  vb(:,1:nfrag) = fraggle_util_vmag_to_vb(fragments%v_r_mag(1:nfrag), fragments%r_unit(:,1:nfrag), v_t_new(1:nfrag), &
                                                          fragments%t_unit(:,1:nfrag), fragments%mass(1:nfrag), impactors%vbcom(:)) 
                  do concurrent (i = 1:nfrag)
                     vc(:,i) = vb(:,i) - impactors%vbcom(:)
                  end do
                  kearr = 0.0_DP
                  do concurrent(i = 1:nfrag)
                     kearr(i) = fragments%mass(i) * dot_product(vc(:,i), vc(:,i))
                  end do
                  keo = 0.5_DP * sum(kearr(:))
                  fval = keo 
                  if (lfirst_func) then
                     fval_scale = keo
                     lfirst_func = .false.
                  end if
                  fval = keo / fval_scale
                  lfailure = .false.
               end select
            end associate

            return
         end function tangential_objective_function

   end subroutine fraggle_generate_tan_vel


   subroutine fraggle_generate_rad_vel(collider, lfailure)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! 
      !! Adjust the fragment velocities to set the fragment orbital kinetic energy. This will minimize the difference between the fragment kinetic energy and the energy budget
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: collider !! Fraggle collision system object
      logical,                  intent(out)   :: lfailure  !! Logical flag indicating whether this step fails or succeeds! 
      ! Internals
      real(DP), parameter                   :: TOL_MIN = FRAGGLE_ETOL   ! This needs to be more accurate than the tangential step, as we are trying to minimize the total residual energy
      real(DP), parameter                   :: TOL_INIT = 1e-14_DP
      real(DP), parameter                   :: VNOISE_MAG = 1e-10_DP !! Magnitude of the noise to apply to initial conditions to help minimizer find a solution in case of failure
      integer(I4B), parameter               :: MAXLOOP = 100
      real(DP)                              :: ke_radial, tol, fval
      integer(I4B)                          :: i
      real(DP), dimension(:), allocatable   :: v_r_initial
      real(DP), dimension(collider%fragments%nbody)       :: vnoise
      type(lambda_obj)                      :: objective_function
      character(len=STRMAX)                 :: message
      real(DP), dimension(:), allocatable :: v_r_output

      associate(impactors => collider%impactors,  nfrag => collider%fragments%nbody)
         select type(fragments => collider%fragments)
         class is (fraggle_fragments(*))
            ! Set the "target" ke for the radial component
            ke_radial = fragments%ke_budget - fragments%ke_spin - fragments%ke_orbit

            do i = 1, nfrag
               fragments%v_t_mag(i) = dot_product(fragments%vc(:,i), fragments%t_unit(:,i))
               fragments%v_r_mag(i) = dot_product(fragments%vc(:,i), fragments%r_unit(:,i))
            end do
            
            allocate(v_r_initial, source=fragments%v_r_mag)

            ! Initialize the lambda function using a structure constructor that calls the init method
            ! Minimize the ke objective function using the BFGS optimizer
            objective_function = lambda_obj(radial_objective_function)
            tol = TOL_INIT
            do while(tol < TOL_MIN)
               fval = radial_objective_function(v_r_initial)
               call minimize_bfgs(objective_function, nfrag, v_r_initial, tol, MAXLOOP, lfailure, v_r_output)
               fragments%v_r_mag = v_r_output
               if (.not.lfailure) exit
               tol = tol * 2 ! Keep increasing the tolerance until we converge on a solution
               v_r_initial(:) = fragments%v_r_mag(:)
               call random_number(vnoise(1:nfrag)) ! Adding a bit of noise to the initial conditions helps it find a solution more often
               vnoise(:) = 1.0_DP + VNOISE_MAG * (2 * vnoise(:) - 1._DP)
               v_r_initial(:) = v_r_initial(:) * vnoise(:)
            end do
            
            ! Shift the radial velocity vectors to align with the center of mass of the collisional system (the origin)
            fragments%ke_orbit = 0.0_DP
            fragments%vb(:,1:nfrag) = fraggle_util_vmag_to_vb(fragments%v_r_mag(1:nfrag), fragments%r_unit(:,1:nfrag), &
                                 fragments%v_t_mag(1:nfrag), fragments%t_unit(:,1:nfrag), fragments%mass(1:nfrag), impactors%vbcom(:)) 
            do i = 1, nfrag
               fragments%vc(:, i) = fragments%vb(:, i) - impactors%vbcom(:)
               fragments%ke_orbit = fragments%ke_orbit + fragments%mass(i) * dot_product(fragments%vc(:, i), fragments%vc(:, i))
            end do
            fragments%ke_orbit = 0.5_DP * fragments%ke_orbit

            lfailure = abs((fragments%ke_budget - (fragments%ke_orbit + fragments%ke_spin)) / fragments%ke_budget) > FRAGGLE_ETOL
            if (lfailure) then
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, " ")
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Radial velocity failure diagnostics")
               write(message, *) fragments%ke_budget
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_budget     : " // trim(adjustl(message)))
               write(message, *) fragments%ke_spin
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_spin       : " // trim(adjustl(message)))
               write(message, *) fragments%ke_orbit
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_orbit : " // trim(adjustl(message)))
               write(message, *) fragments%ke_budget - (fragments%ke_orbit + fragments%ke_spin)
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "ke_remainder  : " // trim(adjustl(message)))
            end if
         end select
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
            real(DP), dimension(collider%fragments%nbody)       :: kearr
            real(DP)                              :: keo, ke_radial, rotmag2, vmag2
     
            associate(impactors => collider%impactors,  nfrag => collider%fragments%nbody)
               select type(fragments => collider%fragments)
               class is (fraggle_fragments(*))
                  allocate(v_shift, mold=fragments%vb)
                  v_shift(:,:) = fraggle_util_vmag_to_vb(v_r_mag_input, fragments%r_unit, fragments%v_t_mag, fragments%t_unit, fragments%mass, impactors%vbcom) 
                  do i = 1,fragments%nbody
                     v_shift(:,i) = v_shift(:,i) - impactors%vbcom(:)
                     rotmag2 = fragments%rot(1,i)**2 + fragments%rot(2,i)**2 + fragments%rot(3,i)**2
                     vmag2 = v_shift(1,i)**2 + v_shift(2,i)**2 + v_shift(3,i)**2
                     kearr(i) = fragments%mass(i) * (fragments%Ip(3, i) * fragments%radius(i)**2 * rotmag2 + vmag2) 
                  end do
                  keo = 2 * fragments%ke_budget - sum(kearr(:))
                  ke_radial = fragments%ke_budget - fragments%ke_orbit - fragments%ke_spin
                  ! The following ensures that fval = 0 is a local minimum, which is what the BFGS method is searching for
                  fval = (keo / (2 * ke_radial))**2
               end select
            end associate
      
            return
         end function radial_objective_function

   end subroutine fraggle_generate_rad_vel

end submodule s_fraggle_generate
