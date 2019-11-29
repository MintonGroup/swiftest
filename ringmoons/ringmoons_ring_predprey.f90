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
   real(DP),dimension(6,5),parameter    :: rkf45_btab = reshape( & ! Butcher tableau for Runge-Kutta-Fehlberg method
      (/        1./4.,       1./4.,          0.,            0.,           0.,           0.,&
                3./8.,      3./32.,      9./32.,            0.,           0.,           0.,&
              12./13., 1932./2197., -7200./2197.,  7296./2197.,           0.,           0.,&
                   1.,   439./216.,          -8.,   3680./513.,   -845./4104.,          0.,&
                1./2.,     -8./27.,           2., -3544./2565.,   1859./4104.,    -11./40./), shape(rkf45_btab))
   real(DP),dimension(6),parameter      :: rkf5_coeff =  (/ 16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55. /)
   real(DP),dimension(6),parameter      :: rkf4_coeff =  (/ 25./216., 0., 1408./2565. ,  2197./4104. , -1./5. ,     0. /)
   real(DP),dimension(0:ring%N+1,6)     :: kGm,kv2
   real(DP),dimension(0:ring%N+1)     :: v2f,Gmf
   real(DP),parameter :: TOL = 1e-8_DP
   real(DP),dimension(0:ring%N+1)       :: Ev2,EGm,sarr,dt,dtleft
   real(DP),parameter                   :: DTMIN = 1.0e-12_DP
   logical(lgt),dimension(0:ring%N+1)   :: ringmask,goodbin

! Executable code
   dt(:) = 2 * PI / ring%w(:)  !dtin
   
   dtnew = dtin
   v2i(:) = (ring%vrel_pdisk(:))**2
   Gmi(:) = ring%Gm_pdisk(:)

   v2f(:) = v2i(:)
   Gmf(:) = Gmi(:)


   where ((ring%Gm(:) > N_DISK_FACTOR * ring%Gm_pdisk(:))) 
      ringmask(:) = .true.
      dtleft(:) = dt(:)
   elsewhere
      ringmask(:) = .false.
      dtleft(:) = 0.0_DP
   end where 
  
   do loop = 1, LOOPMAX 
      if (loop == LOOPMAX) then
         stepfail = .true.
         dtnew = 0.5_DP * dtin
         return
      end if
      stepfail = .false.
      kGm(:,:) = 0._DP
      kv2(:,:) = 0._DP
      goodbin(:) = ringmask(:)

      !write(*,*) 'pred prey input'
      !write(*,*) 'dt: ',dt, 'dtleft: ',dtleft
      !write(*,*) maxval(ring%r_pdisk(:)) * DU2CM,minval(ring%r_pdisk(:)) * DU2CM
      !write(*,*) maxval(ring%vrel_pdisk(:)) * DU2CM / TU2S,minval(ring%vrel_pdisk(:)) * DU2CM / TU2S

      do rkn = 1,6 ! Runge-Kutta-Fehlberg steps 
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
         v2f(:) = v2i(:) + matmul(kv2(:,1:5), rkf4_coeff(1:5))
         Gmf(:) = Gmi(:) + matmul(kGm(:,1:5), rkf4_coeff(1:5))
         where((v2f(:) < 0.0_DP).or.(GMf(:) < 0.0_DP))
            goodbin(:) = .false.
         end where
      end where
      where (goodbin(:))
         Ev2(:) = abs(matmul(kv2(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
         EGm(:) = abs(matmul(kGm(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
         sarr(:) = (TOL / (2 * max(Ev2(:),EGm(:))))**(0.25_DP)
         where((sarr(:) < 1._DP).and.(dt(:) > DTMIN * dtin))
            goodbin(:) = .false.
            dt(:) = 0.5_DP * sarr(:) * dt(:)
         elsewhere 
            v2i(:) = v2f(:)
            Gmi(:) = Gmf(:)
            dtleft(:) = dtleft(:) - dt(:)
            where (dtleft(:) <= 0.0_DP) 
               ringmask(:) = .false.
            elsewhere
               dt(:) = min(0.9_DP * sarr(:) * dt(:),dtleft(:))
            endwhere
         end where
      elsewhere (ringmask(:))
         dt(:) = 0.5_DP * dt(:)
         sarr(:) = 1._DP
      end where
      !write(*,*) 'loop ',loop
      !write(*,*) count(ringmask),count(goodbin),minval(dt,mask=goodbin),maxval(dtleft,mask=ringmask)
      !do i = 0,ring%N+1
      !    write(*,*) i,goodbin(i),sarr(i),dt(i),dtleft(i),v2f(i),Gmf(i)
      !end do 
      !read(*,*)
      if (all(.not.ringmask(:))) exit
   end do

   ring%vrel_pdisk(:) = sqrt(v2f(:))
   ring%Gm_pdisk(:) = Gmf(:)
   ring%r_pdisk(:) = (3 * ring%Gm_pdisk(:) / (4 * PI * ring%rho_pdisk(:)))**(1._DP / 3._DP)
   call ringmoons_update_ring(swifter_pl1P,ring)

  
   dtnew = dtin 
   return
end subroutine ringmoons_ring_predprey
