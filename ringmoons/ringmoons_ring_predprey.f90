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
subroutine ringmoons_ring_predprey(swifter_pl1P,ring,seeds,dtin,stepfail,dtnew)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_predprey
   implicit none

! Arguments
   type(swifter_pl),pointer             :: swifter_pl1P
   type(ringmoons_ring), intent(inout)  :: ring
   type(ringmoons_seeds), intent(in)    :: seeds
   real(DP), intent(in)                 :: dtin
   logical(lgt), intent(out)            :: stepfail   
   real(DP),intent(out)                 :: dtnew

! Internals
   integer(I4B)                         :: rkn,i,rki,loop
   real(DP),dimension(0:ring%N+1)       :: Gmi, v2i, Gm,v2,tau,r,r_hstar,Q,nu
   real(DP),dimension(0:ring%N+1,rkfo)   :: kGm,kv2
   real(DP),dimension(0:ring%N+1)     :: v2f,Gmf,dv2,dGm
   integer(I4B),dimension(0:ring%N+1)     :: loopcounter
   real(DP),parameter :: TOL = 1e-9_DP
   real(DP),dimension(0:ring%N+1)       :: Ev2,EGm,sarr,dt,dtleft,dtmin
   real(DP),parameter                   :: DTMIN_FAC = 1.0e-5_DP
   logical(lgt),dimension(0:ring%N+1)   :: ringmask,goodbin
   real(DP)                             :: mass_limit,rad_limit

! Executable code
   dt(:) = min(1e-1_DP / ring%w(:),dtin)
   dtmin(:) = DTMIN_FAC * dt(:)
   
   dtnew = dtin
   v2i(:) = (ring%vrel_pdisk(:))**2
   Gmi(:) = ring%Gm_pdisk(:)

   v2f(:) = v2i(:)
   Gmf(:) = Gmi(:)
   rad_limit = 1e-1_DP / DU2CM
   mass_limit = 4._DP / 3._DP * PI * rad_limit**3 * maxval(ring%rho_pdisk(:))
   loopcounter = 0


   where (ring%Gm(:) > epsilon(1._DP) * maxval(ring%Gm(:)))
      ringmask(:) = .true.
      dtleft(:) = dt(:)
   elsewhere
      ringmask(:) = .false.
      dtleft(:) = 0.0_DP
   end where 
  
   do loop = 1, LOOPMAX 
      if (loop == LOOPMAX) then
         write(*,*) 'max loop reached in preprey!'
         stepfail = .true.
         dtnew = 0.5_DP * dtin
         return
      end if
      kGm(:,:) = 0._DP
      kv2(:,:) = 0._DP
      goodbin(:) = ringmask(:)
      where (ringmask(:)) loopcounter(:) = loopcounter(:) + 1

      do rkn = 1,rkfo ! Runge-Kutta-Fehlberg steps 
         where (goodbin(:))
            v2(:) = v2i(:) + matmul(kv2(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
            Gm(:) = Gmi(:) + matmul(kGm(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))

            where((v2(:) < 0.0_DP).or.(GM(:) < 0.0_DP))
               goodbin(:) = .false.
            elsewhere 
               Q(:) = ring%w(:) * sqrt(v2(:)) / (3.36_DP * ring%Gsigma(:))
               r(:) = (3 * Gm(:) / (4 * PI * ring%rho_pdisk(:)))**(1._DP / 3._DP) 
               r_hstar(:) = ring%r(:) * (2 * Gm(:) /(3._DP * swifter_pl1P%mass))**(1._DP/3._DP) / (2 * r(:)) 
               tau(:) = PI * r(:)**2 * ring%Gsigma(:) / Gm(:)
               nu(:) = ringmoons_viscosity(ring%Gsigma(:), Gm(:), v2(:), r(:), r_hstar(:), Q(:), tau(:), ring%w(:))
               kv2(:,rkn) = dt(:) * ringmoons_ring_dvdt(Gm(:),v2(:),tau(:),nu(:),ring%w(:)) 
               kGm(:,rkn) = dt(:) * ringmoons_ring_dMdt(Gm(:),v2(:),r(:),tau(:),ring%w(:))
            end where
         end where
      end do
      where (goodbin(:))
         dv2(:) = matmul(kv2(:,1:rkfo-1), rkf4_coeff(1:rkfo-1))
         dGm(:) = matmul(kGm(:,1:rkfo-1), rkf4_coeff(1:rkfo-1))
         v2f(:) = v2i(:) + dv2(:)
         Gmf(:) = Gmi(:) + dGm(:)
         where((v2f(:) < 0.0_DP).or.(GMf(:) < 0.0_DP))
            goodbin(:) = .false.
         end where
      end where
      where (goodbin(:))
         Ev2(:) = abs(matmul(kv2(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
         EGm(:) = abs(matmul(kGm(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
         sarr(:) = (TOL / (2 * max(Ev2(:),EGm(:))))**(0.25_DP)
         where((sarr(:) < 1._DP).and.(dt(:) > dtmin(:)))
            goodbin(:) = .false.
            dt(:) = 0.5_DP * sarr(:) * dt(:)
         elsewhere 
            dv2(:) = v2f(:) - v2i(:)
            dGm(:) = Gmf(:) - Gmi(:)
            v2i(:) = v2f(:)
            Gmi(:) = Gmf(:)
            
            dtleft(:) = dtleft(:) - dt(:)
            where (dtleft(:) <= 0.0_DP) 
               ringmask(:) = .false.
            elsewhere(Gmf(:) < mass_limit) 
               ringmask(:) = .false.
            elsewhere
               dt(:) = min(0.9_DP * sarr(:) * dt(:),dtleft(:))
            endwhere
         end where
      elsewhere (ringmask(:))
         dt(:) = 0.5_DP * dt(:)
         sarr(:) = 1._DP
      end where

      if (all(.not.ringmask(:))) exit

   end do

   ring%vrel_pdisk(:) = sqrt(v2f(:))
   ring%Gm_pdisk(:) = Gmf(:)
   ring%r_pdisk(:) = (3 * ring%Gm_pdisk(:) / (4 * PI * ring%rho_pdisk(:)))**(1._DP / 3._DP)
   call ringmoons_update_ring(swifter_pl1P,ring)

   stepfail = .false. 
   dtnew = dtin 
   return
end subroutine ringmoons_ring_predprey
