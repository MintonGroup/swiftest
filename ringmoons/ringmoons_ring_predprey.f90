!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_ring_predprey
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Evolves the ring aggregate mass and velocity dispersion according to the predator/prey model of 
!                 Esposito et al. (2012)
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_ring_predprey(dt,ring,ring)
!
!  Notes       : 
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_ring_predprey(swifter_pl1P,ring,seeds,dt,stepfail)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_predprey
   implicit none

! Arguments
   type(swifter_pl),pointer                  :: swifter_pl1P
   type(ringmoons_ring), intent(inout)       :: ring
   type(ringmoons_seeds), intent(in)         :: seeds
   real(DP), intent(in)                      :: dt
   logical(lgt), intent(out)                 :: stepfail

! Internals
   integer(I4B)                              :: rkn,i
   real(DP),dimension(0:ring%N+1)            :: kGm,kv2,Gm_pdisk,v2_pdisk,Gmf_pdisk,v2f_pdisk,tau,r_pdisk,r_hstar,Q,nu
   real(DP),dimension(2:4),parameter         :: rkh = (/0.5_DP, 0.5_DP, 1._DP/)
   integer(I4B),dimension(4),parameter       :: rkmult = (/1, 2, 2, 1/)

! Executable code

   stepfail = .false.

   v2_pdisk(:) = (ring%vrel_pdisk(:))**2
   Gm_pdisk(:) = ring%Gm_pdisk(:)
   kGm(:) = 0._DP
   kv2(:) = 0._DP
   v2f_pdisk(:) = 0._DP 
   Gmf_pdisk(:) = 0._DP

   !write(*,*) 'pred prey input'
   !write(*,*) 'dt = ',dt
   !write(*,*) maxval(ring%r_pdisk(:)) * DU2CM

   do rkn = 1,4 ! Runge-Kutta steps 
      if (rkn > 1) then
         v2_pdisk(:) = (ring%vrel_pdisk(:))**2 + rkh(rkn) * kv2(:)
         Gm_pdisk(:) = ring%Gm_pdisk(:)        + rkh(rkn) * kGm(:)

         if (any(v2_pdisk < 0.0_DP) .or. (any(Gm_pdisk < 0.0_DP))) then
            stepfail = .true.
            exit
         end if
      end if
      r_pdisk(:) = (3 * Gm_pdisk(:) / (4 * PI * ring%rho_pdisk(:)))**(1._DP / 3._DP) 
      tau(:) = PI * r_pdisk(:)**2 * ring%Gsigma(:) / ring%Gm_pdisk(:)
      r_hstar(:) = ring%r(:) * (2 * Gm_pdisk(:) /(3._DP * swifter_pl1P%mass))**(1._DP/3._DP) / (2 * r_pdisk(:)) 
      Q(:) = ring%w(:) * sqrt(v2_pdisk(:)) / (3.36_DP * ring%Gsigma(:))

      where ((ring%Gm(:) > N_DISK_FACTOR * ring%Gm_pdisk(:)) .and. (Q(:) < 1._DP))
         nu(:) = ringmoons_viscosity(ring%Gsigma(:), Gm_pdisk(:), v2_pdisk(:), r_pdisk(:), r_hstar(:), Q(:), tau(:), ring%w(:))
         kv2(:) = dt * ringmoons_ring_dvdt(Gm_pdisk(:),v2_pdisk(:),tau(:),nu(:),ring%w(:)) 
         kGm(:) = dt * ringmoons_ring_dMdt(Gm_pdisk(:),v2_pdisk(:),tau(:),ring%w(:))
      elsewhere
         kv2(:) = 0.0_DP
         kGm(:) = 0.0_DP
      end where
      v2f_pdisk(:) = v2f_pdisk(:) + rkmult(rkn) * kv2(:)
      Gmf_pdisk(:) = Gmf_pdisk(:) + rkmult(rkn) * kGm(:)

   end do
   if (stepfail) return


   v2f_pdisk(:) = (ring%vrel_pdisk(:))**2 + v2f_pdisk(:) / 6._DP
   Gmf_pdisk(:) = ring%Gm_pdisk(:) + Gmf_pdisk(:) / 6._DP


   if (any(v2f_pdisk < 0.0_DP)) then
      stepfail = .true.
      return
   end if
   if (any(Gmf_pdisk < 0.0_DP)) then
      stepfail = .true.
      return
   end if

   ring%vrel_pdisk(:) = sqrt(v2f_pdisk(:))
   ring%Gm_pdisk(:) = Gmf_pdisk(:)
   ring%r_pdisk = (3 * ring%Gm_pdisk(:) / (4 * PI * ring%rho_pdisk(:)))**(1._DP / 3._DP)

   call ringmoons_update_ring(swifter_pl1P,ring)



   !write(*,*) 'pred prey output'
   !write(*,*) maxval(ring%r_pdisk(:)) * DU2CM
   !read(*,*)
   
   return
end subroutine ringmoons_ring_predprey
