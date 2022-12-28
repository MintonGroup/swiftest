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

            ! Use the simple collision model to generate initial conditions
            ! Compute the "before" energy/momentum and the budgets
            call self%get_energy_and_momentum(nbody_system, param, lbefore=.true.)
            call self%collision_simple_disruption%disrupt(nbody_system, param, t)
            call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)
            call self%set_budgets()

            ! Minimize difference between energy/momentum and budgets
            call fraggle_generate_minimize(self, lfailure_local)

            call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)
            dEtot = self%Etot(2) - self%Etot(1)
            dLmag = .mag. (self%Ltot(:,2) - self%Ltot(:,1))

            lfailure_local = (dEtot > 0.0_DP) 
            if (lfailure_local) then
               write(message, *) dEtot, abs(dEtot + impactors%Qloss) / FRAGGLE_ETOL
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
      real(DP), parameter :: TOL_MIN = 1.0e-5_DP
      real(DP), parameter :: TOL_INIT = 1e-6_DP
      integer(I4B), parameter :: MAXLOOP = 50
      real(DP), dimension(collider%fragments%nbody) :: input_v
      real(DP), dimension(:), allocatable :: output_v
      type(lambda_obj) :: Efunc
      real(DP) :: tol, fval
      integer(I4B) :: i, nelem

      associate(impactors => collider%impactors,  nfrag => collider%fragments%nbody)
         select type(fragments => collider%fragments)
         class is (fraggle_fragments(*))

            nelem = nfrag 
            lfailure = .false.
               ! Find the local kinetic energy minimum for the nbody_system that conserves linear and angular momentum
            Efunc = lambda_obj(E_objective_function)
            tol = TOL_INIT

            fragments%v_r_unit(:,:) = .unit. fragments%vc(:,:)
            fragments%vmag(:) = .mag. fragments%vc(:,1:nfrag) 
            fragments%rot(:,1:nfrag) = fragments%rot(:,1:nfrag) * 1e-12_DP
            do while(tol < TOL_MIN)

               input_v(:) = fragments%vmag(1:nfrag)
               fval = E_objective_function(input_v)
               call minimize_bfgs(Efunc, nelem, input_v, tol, MAXLOOP, lfailure, output_v)
               fval = E_objective_function(output_v)
               input_v(:) = output_v(:)

               fragments%vmag(1:nfrag) = output_v(1:nfrag)

               do concurrent(i=1:nfrag)
                  fragments%vc(:,i) = abs(fragments%vmag(i)) * fragments%v_r_unit(:,i)
               end do

               ! Set spins in order to conserve angular momentum
               call collision_util_set_spins(fragments)

               if (.not.lfailure) exit
               tol = tol * 2 ! Keep increasing the tolerance until we converge on a solution
            end do


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
   
            associate(impactors => collider%impactors, nfrag => collider%fragments%nbody)
               select type(fragments => collider%fragments)
               class is (fraggle_fragments(*))
                  allocate(tmp_frag, source=fragments)
                  tmp_frag%vmag(1:nfrag) = val_input(1:nfrag)
                  do concurrent(i=1:nfrag)
                     tmp_frag%vc(:,i) = abs(tmp_frag%vmag(i)) * tmp_frag%v_r_unit(:,i)
                  end do

                  call collision_util_shift_vector_to_origin(tmp_frag%mass, tmp_frag%vc)
                  ! Set spins in order to conserve angular momentum
                  call collision_util_set_spins(tmp_frag)

                  ! Get the current kinetic energy of the system
                  call tmp_frag%get_kinetic_energy()
                  deltaE = tmp_frag%ke_budget - (tmp_frag%ke_orbit + tmp_frag%ke_spin)
            
                  ! Use the deltaE as the basis of our objective function, with a higher penalty for having excess kinetic energy compared with having a deficit
                  if (deltaE < 0.0_DP) then
                     fval = deltaE**8
                  else
                     fval = deltaE**2
                  end if
                  deallocate(tmp_frag)
                  
               end select
            end associate
         
            return
         end function E_objective_function


   end subroutine fraggle_generate_minimize

end submodule s_fraggle_generate
