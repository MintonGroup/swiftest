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

   real(DP), parameter     :: FRAGGLE_LTOL = 1e-4_DP !10 * epsilon(1.0_DP)
   real(DP), parameter     :: FRAGGLE_ETOL = 1e-12_DP

contains

   module subroutine fraggle_generate(self, nbody_system, param, t)
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic disruption collision
      !! 
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: self
      class(base_nbody_system),           intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),             intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
      real(DP),                           intent(in)    :: t            !! Time of collision
      ! Internals
      integer(I4B)          :: i, ibiggest, nfrag
      character(len=STRMAX) :: message 
      real(DP) :: dpe

      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(impactors => self%impactors, status => self%status)

            select case (impactors%regime) 
            case (COLLRESOLVE_REGIME_HIT_AND_RUN)
               call self%hitandrun(nbody_system, param, t)
               return
            case (COLLRESOLVE_REGIME_MERGE, COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
               call self%merge(nbody_system, param, t) ! Use the default collision model, which is merge
               return
            case(COLLRESOLVE_REGIME_DISRUPTION)
               message = "Disruption between"
            case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
               message = "Supercatastrophic disruption between"
            case default 
               write(*,*) "Error in swiftest_collision, unrecognized collision regime"
               call util_exit(FAILURE)
            end select
            call self%set_mass_dist(param) 
            call self%disrupt(nbody_system, param, t)

            dpe = self%pe(2) - self%pe(1) 
            nbody_system%Ecollisions = nbody_system%Ecollisions - dpe 
            nbody_system%Euntracked  = nbody_system%Euntracked + dpe 

            associate (fragments => self%fragments)
               ! Populate the list of new bodies
               nfrag = fragments%nbody
               write(message, *) nfrag
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")
               select case(impactors%regime)
               case(COLLRESOLVE_REGIME_DISRUPTION)
                  status = DISRUPTED
                  ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
                  fragments%id(1) = pl%id(ibiggest)
                  fragments%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
                  param%maxid = fragments%id(nfrag)
               case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
                  status = SUPERCATASTROPHIC
                  fragments%id(1:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag)]
                  param%maxid = fragments%id(nfrag)
               end select

               call collision_resolve_mergeaddsub(nbody_system, param, t, status)
            end associate
         end associate
      end select
      end select
      return
   end subroutine fraggle_generate


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
      integer(I4B), parameter              :: MAXTRY = 100
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

         if (param%lflatten_interactions) then
            lk_plpl = allocated(pl%k_plpl)
            if (lk_plpl) deallocate(pl%k_plpl)
         else 
            lk_plpl = .false.
         end if
         call ieee_set_flag(ieee_all, .false.) ! Set all fpe flags to quiet

         !call self%set_natural_scale()

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
            call self%set_budgets()
            call fraggle_generate_pos_vec(self)
            call fraggle_generate_rot_vec(self)
            call fraggle_generate_vel_vec(self)
            call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)
            exit
            dEtot = self%Etot(2) - self%Etot(1)
            dLmag = .mag. (self%Ltot(:,2) - self%Ltot(:,1))

            lfailure_local = (dEtot > FRAGGLE_ETOL) 
            if (lfailure_local) then
               write(message, *) "dEtot: ",dEtot, "dEtot/Qloss", dEtot / impactors%Qloss
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Fraggle failed due to energy gain: " // &
                                                        trim(adjustl(message)))
               cycle
               lfailure_local = .false.
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

         !call self%set_original_scale()

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


   module subroutine fraggle_generate_hitandrun(self, nbody_system, param, t) 
      !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Create the fragments resulting from a non-catastrophic hit-and-run collision
      !! 
      implicit none
      ! Arguments
      class(collision_fraggle),   intent(inout) :: self         !! Collision system object
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(base_parameters),   intent(inout) :: param        !! Current run configuration parameters with SyMBA additions
      real(DP),                 intent(in)    :: t            !! Time of collision
      ! Result
      integer(I4B)                            :: status           !! Status flag assigned to this outcome
      ! Internals
      integer(I4B)                            :: i, ibiggest, nfrag, jtarg, jproj
      logical                                 :: lpure 
      character(len=STRMAX) :: message
      real(DP) :: dpe


      select type(nbody_system)
      class is (swiftest_nbody_system)
      select type(pl => nbody_system%pl)
      class is (swiftest_pl)
         associate(impactors => self%impactors)
            message = "Hit and run between"
            call collision_io_collider_message(nbody_system%pl, impactors%id, message)
            call swiftest_io_log_one_message(COLLISION_LOG_OUT, trim(adjustl(message)))

            if (impactors%mass(1) > impactors%mass(2)) then
               jtarg = 1
               jproj = 2
            else
               jtarg = 2
               jproj = 1
            end if

            ! The frag disruption model (and its extended types allow for non-pure hit and run. 
            if (impactors%mass_dist(2) > 0.9_DP * impactors%mass(jproj)) then ! Pure hit and run, so we'll just keep the two bodies untouched
               call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Pure hit and run. No new fragments generated.")
               nfrag = 0
               call self%collision_basic%hitandrun(nbody_system, param, t)
               lpure = .true.
               return
            else ! Imperfect hit and run, so we'll keep the largest body and destroy the other
               lpure = .false.
               call self%set_mass_dist(param)

               ! Generate the position and velocity distributions of the fragments
               call self%get_energy_and_momentum(nbody_system, param, lbefore=.true.)
               call self%disrupt(nbody_system, param, t, lpure)
               call self%get_energy_and_momentum(nbody_system, param, lbefore=.false.)

               dpe = self%pe(2) - self%pe(1) 
               nbody_system%Ecollisions = nbody_system%Ecollisions - dpe 
               nbody_system%Euntracked  = nbody_system%Euntracked + dpe 

               if (lpure) then
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Should have been a pure hit and run instead")
                  nfrag = 0
               else
                  nfrag = self%fragments%nbody
                  write(message, *) nfrag
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, "Generating " // trim(adjustl(message)) // " fragments")
               end if
            end if

            ibiggest = impactors%id(maxloc(pl%Gmass(impactors%id(:)), dim=1))
            self%fragments%id(1) = pl%id(ibiggest)
            self%fragments%id(2:nfrag) = [(i, i = param%maxid + 1, param%maxid + nfrag - 1)]
            param%maxid = self%fragments%id(nfrag)
            status = HIT_AND_RUN_DISRUPT
            call collision_resolve_mergeaddsub(nbody_system, param, t, status)
         end associate
      end select
      end select

      return
   end subroutine fraggle_generate_hitandrun


   module subroutine fraggle_generate_pos_vec(collider)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Initializes the position vectors of the fragments around the center of mass based on the collision style.
      !! For hit and run with disruption, the fragments are generated in a random cloud around the smallest of the two colliders (body 2)
      !! For disruptive collisions, the fragments are generated in a random cloud around the impact point. Bodies are checked for overlap and
      !! regenerated if they overlap.
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: collider !! Fraggle collision system object
      ! Internals
      real(DP)  :: dis
      real(DP), dimension(NDIM,2) :: fragment_cloud_center
      real(DP), dimension(2) :: fragment_cloud_radius
      logical, dimension(collider%fragments%nbody) :: loverlap
      integer(I4B) :: i, j, loop
      logical :: lcat, lhitandrun
      integer(I4B), parameter :: MAXLOOP = 10000
      real(DP) :: rdistance
      real(DP), parameter :: fail_scale = 1.1_DP ! Scale factor to apply to cloud radius and distance if cloud generation fails


      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)
         lcat = (impactors%regime == COLLRESOLVE_REGIME_SUPERCATASTROPHIC) 
         lhitandrun = (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) 

         ! We will treat the first two fragments of the list as special cases. 
         ! Place the first two bodies at the centers of the two fragment clouds, but be sure they are sufficiently far apart to avoid overlap
         rdistance = .mag. (impactors%rc(:,2) - impactors%rc(:,1)) - sum(fragments%radius(1:2))
         rdistance = min(0.5_DP*rdistance, 1e-6_DP*impactors%radius(2))

         fragment_cloud_radius(:) = impactors%radius(:)

         loverlap(:) = .true.
         do loop = 1, MAXLOOP
            if (.not.any(loverlap(:))) exit
            fragment_cloud_center(:,1) = impactors%rc(:,1) + rdistance * impactors%bounce_unit(:)
            fragment_cloud_center(:,2) = impactors%rc(:,2) - rdistance * impactors%bounce_unit(:)
            do concurrent(i = 1:nfrag, loverlap(i))
               if (i < 3) then
                  fragments%rc(:,i) = fragment_cloud_center(:,i)
               else
                  ! Make a random cloud
                  call random_number(fragments%rc(:,i))
   
                  ! Make the fragment cloud symmertic about 0
                  fragments%rc(:,i) = 2 *(fragments%rc(:,i) - 0.5_DP)
   
                  j = fragments%origin_body(i)
   
                  ! Scale the cloud size
                  fragments%rc(:,i) = fragment_cloud_radius(j) * fragments%rc(:,i)
   
                  ! Shift to the cloud center coordinates
                  fragments%rc(:,i) = fragments%rc(:,i) + fragment_cloud_center(:,j)
               end if
            end do

            ! Check for any overlapping bodies.
            loverlap(:) = .false.
            do j = 1, nfrag
               do i = j + 1, nfrag
                  dis = .mag.(fragments%rc(:,j) - fragments%rc(:,i))
                  loverlap(i) = loverlap(i) .or. (dis <= (fragments%radius(i) + fragments%radius(j))) 
               end do
            end do
            rdistance = rdistance * fail_scale
            fragment_cloud_radius(:) = fragment_cloud_radius(:) * fail_scale
         end do

         call collision_util_shift_vector_to_origin(fragments%mass, fragments%rc)
         call collider%set_coordinate_system()

         do concurrent(i = 1:nfrag)
            fragments%rb(:,i) = fragments%rc(:,i) + impactors%rbcom(:)
         end do

         impactors%rbcom(:) = 0.0_DP
         do concurrent(i = 1:nfrag)
            impactors%rbcom(:) = impactors%rbcom(:) + fragments%mass(i) * fragments%rb(:,i) 
         end do
         impactors%rbcom(:) = impactors%rbcom(:) / fragments%mtot
      end associate

      return
   end subroutine fraggle_generate_pos_vec


   module subroutine fraggle_generate_rot_vec(collider)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Calculates the spins of a collection of fragments such that they conserve angular momentum 
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: collider !! Fraggle collision system object
      ! Internals
      real(DP), dimension(NDIM) :: Lbefore, Lafter, Lspin, rotdir
      real(DP) :: v_init, v_final, mass_init, mass_final, rotmag
      integer(I4B) :: i

      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)

         ! Torque the first body based on the change in angular momentum betwen the pre- and post-impact system assuming a simple bounce
         mass_init = impactors%mass(2)
         mass_final = sum(fragments%mass(2:nfrag))
         v_init = .mag.(impactors%vb(:,2) - impactors%vb(:,1))
         v_final = sqrt(mass_init / mass_final * v_init**2 - 2 * impactors%Qloss / mass_final)

         Lbefore(:) = mass_init * (impactors%rb(:,2) - impactors%rb(:,1)) .cross. (impactors%vb(:,2) - impactors%vb(:,1))
          
         Lafter(:) = mass_final * (impactors%rb(:,2) - impactors%rb(:,1)) .cross. (v_final * impactors%bounce_unit(:))
         Lspin(:) = impactors%Lspin(:,1) + (Lbefore(:) - Lafter(:))
         fragments%rot(:,1) = Lspin(:) / (fragments%mass(1) * fragments%radius(1)**2 * fragments%Ip(3,1))

         ! Add in some random spin noise. The magnitude will be scaled by the before-after amount and the direction will be random
         do concurrent(i = 2:nfrag)
            call random_number(rotdir)
            rotdir = rotdir - 0.5_DP
            rotdir = .unit. rotdir
            call random_number(rotmag)
            rotmag = rotmag * .mag. (Lbefore(:) - Lafter(:)) / ((nfrag - 1) * fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i))
            fragments%rot(:,i) = rotmag * rotdir
         end do

      end associate

      return
   end subroutine fraggle_generate_rot_vec


   module subroutine fraggle_generate_vel_vec(collider)
      !! Author:  David A. Minton
      !!
      !! Generates an initial velocity distribution. For disruptions, the velocity magnitude is set to be
      !! 2x the escape velocity of the colliding pair. For hit and runs the velocity magnitude is set to be
      !! 2x the escape velocity of the smallest of the two bodies.
      implicit none
      ! Arguments
      class(collision_fraggle), intent(inout) :: collider !! Fraggle collision system object
      ! Internals
      integer(I4B) :: i,j, loop, istart
      logical :: lhitandrun, lnoncat
      real(DP), dimension(NDIM) :: vimp_unit, rimp, vrot, Lresidual
      real(DP) :: vmag, vesc, rotmag, ke_residual, ke_per_dof
      integer(I4B), dimension(collider%fragments%nbody) :: vsign
      real(DP), dimension(collider%fragments%nbody) :: vscale, mass_vscale, ke_avail
      integer(I4B), parameter :: MAXLOOP = 100

      associate(fragments => collider%fragments, impactors => collider%impactors, nfrag => collider%fragments%nbody)
         lhitandrun = (impactors%regime == COLLRESOLVE_REGIME_HIT_AND_RUN) 
         lnoncat = (impactors%regime /= COLLRESOLVE_REGIME_SUPERCATASTROPHIC) ! For non-catastrophic impacts, make the fragments act like ejecta and point away from the impact point

         ! The fragments will be divided into two "clouds" based on identified origin body. 
         ! These clouds will collectively travel like two impactors bouncing off of each other. 
         where(fragments%origin_body(:) == 1)
            vsign(:) = -1
         elsewhere
            vsign(:) = 1
         end where

         ! The minimum fragment velocity will be set by the escape velocity
         if (lhitandrun) then
            vesc = sqrt(2 * impactors%Gmass(2) / impactors%radius(2))
         else
            vesc = sqrt(2 * sum(impactors%Gmass(:)) / sum(impactors%radius(:)))
         end if

         ! Scale the magnitude of the velocity by the distance from the impact point
         ! This will reduce the chances of fragments colliding with each other immediately, and is more physically correct  
         do concurrent(i = 1:nfrag)
            rimp(:) = fragments%rc(:,i) - impactors%rbimp(:) 
            vscale(i) = .mag. rimp(:) / (.mag. (impactors%rb(:,2) - impactors%rb(:,1)))
         end do
         vscale(:) = vscale(:)/minval(vscale(:))

         ! Give the fragment velocities a random value that is scaled with fragment mass
         call random_number(mass_vscale)
         mass_vscale(:) = (mass_vscale(:) + 1.0_DP) / 2
         mass_vscale(:) = mass_vscale(:) * (fragments%mtot / fragments%mass(:))**(0.125_DP) ! The 1/8 power is arbitrary. It just gives the velocity a small mass dependence
         mass_vscale(:) = mass_vscale(:) / minval(mass_vscale(:))

         ! Set the velocities of all fragments using all of the scale factors determined above
         do concurrent(i = 1:nfrag)
            j = fragments%origin_body(i)
            vrot(:) = impactors%rot(:,j) .cross. (fragments%rc(:,i) - impactors%rc(:,j))
            if (lhitandrun) then
               if (i == 1) then
                  fragments%vc(:,1) = impactors%vc(:,1)
               else
                  vmag = .mag.impactors%vc(:,2) / (maxval(mass_vscale(:) * maxval(vscale(:))))
                  fragments%vc(:,i) = vmag * mass_vscale(i) * vscale(i) * impactors%bounce_unit(:) * vsign(i) + vrot(:)
               end if
            else
               ! Add more velocity dispersion to disruptions vs hit and runs.
               vmag = vesc * vscale(i) * mass_vscale(i)
               rimp(:) = fragments%rc(:,i) - impactors%rbimp(:)
               vimp_unit(:) = .unit. rimp(:)
               fragments%vc(:,i) = vmag *  (impactors%bounce_unit(:) + vimp_unit(:)) * vsign(i) + vrot(:)
            end if
         end do

         if (lhitandrun) then
            istart = 2
         else 
            istart = 1
         end if
         do loop = 1, MAXLOOP
            call collider%set_coordinate_system()
            call fragments%get_kinetic_energy()
            ke_residual = fragments%ke_budget - (fragments%ke_orbit + fragments%ke_spin)
            ! Make sure we don't take away too much orbital kinetic energy, otherwise the fragment can't escape
            ke_avail(:) = fragments%ke_orbit_frag(:) - impactors%Gmass(1)*impactors%mass(2)/fragments%vmag(:)
            ke_per_dof = -ke_residual/((nfrag - istart + 1))
            do concurrent(i = istart:nfrag, ke_avail(i) > ke_per_dof)
               fragments%vmag(i) = sqrt(2 * (fragments%ke_orbit_frag(i) - ke_per_dof)/fragments%mass(i))
               fragments%vc(:,i) = fragments%vmag(i) * fragments%v_unit(:,i)
            end do
            ! do concurrent(i = istart:nfrag, fragments%ke_spin_frag(i) > ke_per_dof)
            !    fragments%rotmag(i) = sqrt(2 * (fragments%ke_spin_frag(i) - ke_per_dof)/(fragments%mass(i) * fragments%radius(i)**2 * fragments%Ip(3,i)))
            !    fragments%rot(:,i) = fragments%rotmag(i) * .unit.fragments%rot(:,i)
            ! end do
            call fragments%get_kinetic_energy()
            ke_residual = fragments%ke_budget - (fragments%ke_orbit + fragments%ke_spin)

            ! Check for any residual angular momentum, and if there is any, put it into spin of the biggest body
            call collider%set_coordinate_system()
            call fragments%get_angular_momentum()
            Lresidual(:) = fragments%L_budget(:) - (fragments%Lorbit(:) + fragments%Lspin(:)) 
            rotmag = .mag. fragments%rot(:,1)
            fragments%rot(:,1) = fragments%rot(:,1) + Lresidual(:) / (fragments%mass(1) * fragments%radius(1)**2 * fragments%Ip(:,1))
            rotmag = .mag. fragments%rot(:,1)

            if (ke_residual >= 0.0_DP) exit

         end do

         do concurrent(i = 1:nfrag)
            fragments%vb(:,i) = fragments%vc(:,i) + impactors%vbcom(:)
         end do

         impactors%vbcom(:) = 0.0_DP
         do concurrent(i = 1:nfrag)
            impactors%vbcom(:) = impactors%vbcom(:) + fragments%mass(i) * fragments%vb(:,i) 
         end do
         impactors%vbcom(:) = impactors%vbcom(:) / fragments%mtot

      end associate
      return
   end subroutine fraggle_generate_vel_vec


end submodule s_fraggle_generate
