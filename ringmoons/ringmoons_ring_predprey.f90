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
   integer(I4B)                         :: rkn,i,j,rki,loop,goodval
   real(DP),dimension(0:ring%N+1)       :: Gmi, v2i, Gm,v2,tau,r,r_hstar,Q,nu
   real(DP),dimension(0:ring%N+1,rkfo)   :: kGm,kv2
   real(DP),dimension(0:ring%N+1)     :: v2f,Gmf,dv2,dGm
   integer(I4B),dimension(0:ring%N+1)     :: loopcounter
   real(DP),parameter :: TOL = 1e-16_DP
   integer(I4B),parameter               :: NRKFLOOP = 10000
   real(DP),dimension(0:ring%N+1)       :: Ev2,EGm,sarr,dt,dtleft,dtmin
   real(DP),parameter                   :: DTMIN_FAC = 1.0e-5_DP
   logical(lgt),dimension(0:ring%N+1)   :: ringmask,goodbin
   real(DP)                             :: mass_limit,rad_limit,upper_rad_limit,upper_mass_limit,rtmp,x1,x2
   interface
      function nufunc(Gm_pdisk,rho_pdisk,GMP,Gsigma,ringr,ringw) result(ans)
      use module_parameters 
      implicit none
      real(DP),intent(in) :: Gm_pdisk,rho_pdisk,GMP,Gsigma,ringr,ringw
      real(DP) :: ans
      end function nufunc

      integer function BrentRoots(AFunction,x1,Tolerance,d1,d2,d3,d4,d5)
      use module_parameters
      implicit none
      real(DP),external :: AFunction
      real(DP),intent(inout) :: x1
      real(DP),intent(in) :: Tolerance
      real(DP),intent(in) :: d1,d2,d3,d4,d5
      end function BrentRoots

      elemental function v2sticking_threshold(Gm_pdisk,r_pdisk) result(vth2)
      use module_parameters
      implicit none
      real(DP),intent(in) :: Gm_pdisk,r_pdisk
      real(DP) :: vth2
      end function v2sticking_threshold
   end interface


! Executable code
   dt(:) = min(1.0e0_DP / ring%w(:),dtin)
   dtmin(:) = DTMIN_FAC * dt(:)
   
   dtnew = dtin
   v2i(:) = (ring%vrel_pdisk(:))**2
   Gmi(:) = ring%Gm_pdisk(:)

   v2f(:) = v2i(:)
   Gmf(:) = Gmi(:)
   rad_limit = 1e-1_DP / DU2CM
   upper_rad_limit = 1e6_DP / DU2cm
   mass_limit = 4._DP / 3._DP * PI * rad_limit**3 * maxval(ring%rho_pdisk(:))
   upper_mass_limit = 4._DP / 3._DP * PI * upper_rad_limit**3 * maxval(ring%rho_pdisk(:))
!   loopcounter = 0


!   where (ring%Gm(:) > epsilon(1._DP) * maxval(ring%Gm(:)))
!      ringmask(:) = .true.
!      dtleft(:) = dt(:)
!   elsewhere
!      ringmask(:) = .false.
!      dtleft(:) = 0.0_DP
!   end where 
!  
!   do loop = 1, NRKFLOOP
!      !if (loop == LOOPMAX) then
!      !   write(*,*) 'max loop reached in preprey!'
!      !   stepfail = .true.
!      !   dtnew = 0.5_DP * dtin
!      !   return
!      !end if
!      kGm(:,:) = 0._DP
!      kv2(:,:) = 0._DP
!      goodbin(:) = ringmask(:)
!      where (ringmask(:)) loopcounter(:) = loopcounter(:) + 1
!
!      do rkn = 1,rkfo ! Runge-Kutta-Fehlberg steps 
!         where (goodbin(:))
!            v2(:) = v2i(:) + matmul(kv2(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
!            Gm(:) = Gmi(:) + matmul(kGm(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
!
!            where((v2(:) < 0.0_DP).or.(GM(:) < 0.0_DP))
!               goodbin(:) = .false.
!            elsewhere 
!               Q(:) = ring%w(:) * sqrt(v2(:)) / (3.36_DP * ring%Gsigma(:))
!               r(:) = (3 * Gm(:) / (4 * PI * ring%rho_pdisk(:)))**(1._DP / 3._DP) 
!               r_hstar(:) = ring%r(:) * (2 * Gm(:) /(3._DP * swifter_pl1P%mass))**(1._DP/3._DP) / (2 * r(:)) 
!               tau(:) = PI * r(:)**2 * ring%Gsigma(:) / Gm(:)
!               nu(:) = ringmoons_viscosity(ring%Gsigma(:), Gm(:), v2(:), r(:), r_hstar(:), Q(:), tau(:), ring%w(:))
!               kv2(:,rkn) = dt(:) * ringmoons_ring_dvdt(Gm(:),v2(:),tau(:),nu(:),ring%w(:)) 
!               kGm(:,rkn) = dt(:) * ringmoons_ring_dMdt(Gm(:),v2(:),r(:),tau(:),ring%w(:))
!            end where
!         end where
!      end do
!      where (goodbin(:))
!         dv2(:) = matmul(kv2(:,1:rkfo-1), rkf4_coeff(1:rkfo-1))
!         dGm(:) = matmul(kGm(:,1:rkfo-1), rkf4_coeff(1:rkfo-1))
!         v2f(:) = v2i(:) + dv2(:)
!         Gmf(:) = Gmi(:) + dGm(:)
!         where((v2f(:) < 0.0_DP).or.(GMf(:) < 0.0_DP))
!            goodbin(:) = .false.
!         end where
!      end where
!      where (goodbin(:))
!         Ev2(:) = abs(matmul(kv2(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
!         EGm(:) = abs(matmul(kGm(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
!         sarr(:) = (TOL / (2 * max(Ev2(:),EGm(:))))**(0.25_DP)
!         where((sarr(:) < 1._DP).and.(dt(:) > dtmin(:)))
!            goodbin(:) = .false.
!            dt(:) = 0.5_DP * sarr(:) * dt(:)
!         elsewhere 
!            !dv2(:) = v2f(:) - v2i(:)
!            !dGm(:) = Gmf(:) - Gmi(:)
!            v2i(:) = v2f(:)
!            Gmi(:) = Gmf(:)
!            
!            dtleft(:) = dtleft(:) - dt(:)
!            where (dtleft(:) <= 0.0_DP) 
!               ringmask(:) = .false.
!            elsewhere(Gmf(:) < mass_limit) 
!               ringmask(:) = .false.
!            elsewhere
!               dt(:) = min(0.9_DP * sarr(:) * dt(:),dtleft(:))
!            endwhere
!         end where
!      elsewhere (ringmask(:))
!         dt(:) = 0.5_DP * dt(:)
!         sarr(:) = 1._DP
!      end where
!
!      if (all(.not.ringmask(:))) exit
!
!   end do
!

   do i = 1,ring%N
      if (ring%Gm(i) > epsilon(1._DP) * maxval(ring%Gm(:))) then
         Gmf(i)  = ring%Gm_pdisk(i)
         goodval = BrentRoots(nufunc,Gmf(i),TOL*Gmf(i),ring%rho_pdisk(i),swifter_pl1P%mass,ring%Gsigma(i),ring%r(i),ring%w(i))
         if (goodval /= 0) Gmf(i) = ring%Gm_pdisk(i)
      end if
   end do
   
   ring%Gm_pdisk(:) = max(Gmf(:),mass_limit)
   ring%r_pdisk(:) = (3 * ring%Gm_pdisk(:) / (4 * PI * ring%rho_pdisk(:)))**(1._DP / 3._DP)
   v2f(:) = v2sticking_threshold(ring%Gm_pdisk(:),ring%r_pdisk(:))
   ring%vrel_pdisk(:) = sqrt(v2f(:))
   call ringmoons_update_ring(swifter_pl1P,ring)

   stepfail = .false. 
   dtnew = dtin 
   return
end subroutine ringmoons_ring_predprey

elemental function v2sticking_threshold(Gm_pdisk,r_pdisk) result(vth2)
use module_parameters
implicit none

real(DP),intent(in) :: Gm_pdisk,r_pdisk
real(DP) :: vth2

real(DP),parameter                     :: MBOUNCE_G = 2.1e-10_DP !See Weidling et al. (2012) eq. 8

!Threshold velocity for sticking
!See Esposito et al. (2012) and Blum (2006)
!vth2 = 1.0e0_DP

!See Weidling et al. (2012)
vth2 = ((Gm_pdisk * MU2GM / GU) / MBOUNCE_G)**(-5._DP / 18._DP) 
! Convert from m/s to system units
vth2 = vth2 * 100._DP * TU2S / DU2CM
vth2 = max(vth2**2,2 * Gm_pdisk / r_pdisk)

return
end function v2sticking_threshold



function nufunc(Gm_pdisk,rho_pdisk,GMP,Gsigma,ringr,ringw) result(ans)
use module_parameters 
use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_predprey
implicit none

real(DP),intent(in) :: Gm_pdisk,rho_pdisk,GMP,Gsigma,ringr,ringw
real(DP) :: ans

real(DP) :: nu,Tcoll,Torb,r_pdisk,r_hstar,Q,tau,v2_pdisk,S2
real(DP)                               :: eps2 ! coefficient of restitution
real(DP), parameter                    :: Vc = 0.0077 ! See Schmidt et al. (2006) eqn. 14.14
real(DP), parameter                    :: eps_exponent = -0.234_DP
real(DP), parameter                    :: eps_constant = 0.89_DP ! See Brisset et al. (2019)
interface
   elemental function v2sticking_threshold(Gm_pdisk,r_pdisk) result(vth2)
   use module_parameters
   implicit none
   real(DP),intent(in) :: Gm_pdisk,r_pdisk
   real(DP) :: vth2
   end function v2sticking_threshold
end interface



!eps2 = (v2_pdisk * (DU2CM / TU2S)**2 / Vc**2)**(eps_exponent)
!eps2 = 0.9_DP * exp(-0.22_DP * (sqrt(v2_pdisk) * DU2CM / TU2S)) + 0.01_DP * (sqrt(v2_pdisk) * DU2CM / TU2S)**(-0.6_DP)
eps2 = eps_constant**2

r_pdisk = ((3 * Gm_pdisk) / (4 * PI * rho_pdisk))**(1._DP / 3._DP) 

v2_pdisk = v2sticking_threshold(Gm_pdisk,r_pdisk)

Q = ringw * sqrt(v2_pdisk) / (3.36_DP * Gsigma)

r_hstar = ringr * (2 * Gm_pdisk /(3._DP * GMP))**(1._DP/3._DP) / (2 * r_pdisk) 
tau = PI * r_pdisk**2 * Gsigma / Gm_pdisk

nu = ringmoons_viscosity(Gsigma, Gm_pdisk, v2_pdisk, r_pdisk, r_hstar, Q, tau, ringw)
Torb = 2 * PI / ringw
Tcoll = Torb / (4 * tau) 

S2 = 2.25_DP * ringw**2  ! Shear rate squared

ans = v2_pdisk / Tcoll * (1._DP - eps2) - nu * S2


end function nufunc


!*****************************************************
!*      Program to demonstrate the real domain       *
!*               Zbrent subroutine                   *
!* ------------------------------------------------- *
!* Reference:  BORLAND MATHEMATICAL LIBRARY          *
!*                                                   *
!*                F90 version by J-P Moreau, Paris.  *
!*                Modified by David Minton           *
!* ------------------------------------------------- *
!* Example:    Find a real root of f(x)=(x+1)^5      *
!*                                                   *
!* SAMPLE RUN:                                       *
!*                                                   *
!*  Input interval (X1,X2):                          *
!*                                                   *
!*        X1 = -2                                    *
!*        X2 =  0                                    *
!*                                                   *
!*  Convergence criterion: 1e-10                     *
!*  Maximum number of iterations: 10                 *
!*                                                   *
!*  The estimated root is:                           *
!*                                                   *
!*        X = -1.000000                              *
!*                                                   *
!*  The associated Y value is Y = 0.000000           *
!*                                                   *
!*  The number of iterations was: 2                  *
!* Returns: 0 converged
!*         -1 root not bracketed
!*         -2 maxed the loops
!*                                                   *
!*****************************************************
 !real(DP) Function BrentRoots_double( x1, x2, Tolerance,  &
 !                           maxIterations,      &
 !                           valueAtRoot,        &
 !                           niter, error )  
  integer function BrentRoots(AFunction,x1,Tolerance,d1,d2,d3,d4,d5)
   use module_parameters
   implicit none


  real(DP),external :: AFunction
  real(DP),intent(inout) :: x1
  real(DP),intent(in) :: Tolerance
  real(DP),intent(in) :: d1,d2,d3,d4,d5
  integer,parameter :: maxIterations=100
  real(DP) :: valueAtRoot

  real(DP),parameter :: FPP = 1.d-11
  real(DP),parameter :: nearzero = 1.d-20

  integer, external :: RootBracketed
  real(DP),external :: Minimum
  integer :: niter,error

  real(DP) :: resultat,AA,BB,CC,DD,EE,FA,FB,FC,Tol1,PP,QQ,RR,SS,xm
  integer :: i, done
   integer,parameter :: NTRY = 50
   integer,parameter :: NBRACKET = 5
   real(DP),parameter :: FIRSTFACTOR = 1.6_DP
   integer :: br,j
   real(DP) :: startx,x2,factor,f1,f2
   
   factor = FIRSTFACTOR

   startx = x1
      bracket: do br = 1,NBRACKET
      !if (br == NBRACKEt)  write(*,*) 'could not bracket root ',x1,x2,f1,f2,factor
         x1 = startx
         x2 = startx ** factor

         ! First bracket the root
         f1 = AFunction(x1,d1,d2,d3,d4,d5)
         f2 = AFunction(x2,d1,d2,d3,d4,d5)

         do j=1,NTRY
            if (f1 * f2 < 0._DP) exit bracket
            if (abs(f1) < abs(f2)) then
               x1 = max(exp(log(x1) + factor * (log(x1) - log(x2))),tiny(1._DP))
               f1 = AFunction(x1,d1,d2,d3,d4,d5)
            else
               x2 = max(exp(log(x2) + factor * (log(x2) - log(x1))),tiny(1._DP))
               f2 = AFunction(x2,d1,d2,d3,d4,d5)
            end if
         end do
         factor = 0.5_DP * (factor + 1._DP)
      end do bracket



  i = 0
  done = 0
  error = 0
  AA = x1
  BB = x2
  FA = AFunction(AA,d1,d2,d3,d4,d5)
  FB = AFunction(BB,d1,d2,d3,d4,d5)
  if (RootBracketed(FA,FB).eq.0) then 
    error = -1
    resultat=x1
  else 
    FC = FB;
    do while (done.eq.0.and.i < maxIterations)
      if (RootBracketed(FC,FB).eq.0) then
        CC = AA; FC = FA; DD = BB - AA; EE = DD
      endif
      if (abs(FC) < abs(FB)) then
        AA = BB; BB = CC; CC = AA
        FA = FB; FB = FC; FC = FA
      endif
      Tol1 = 2.d0 * FPP * abs(BB) + 0.5d0 * Tolerance
      xm = 0.5d0 * (CC-BB)
      if ((abs(xm) <= Tol1).or.(abs(FA) < nearzero)) then
        ! A root has been found
        resultat = BB;
        done = 1
        valueAtRoot = AFunction(resultat,d1,d2,d3,d4,d5)
      else 
        if ((abs(EE) >= Tol1).and.(abs(FA) > abs(FB))) then
          SS = FB/ FA;
          if (abs(AA - CC) < nearzero) then
            PP = 2.d0 * xm * SS;
            QQ = 1.d0 - SS;
          else 
            QQ = FA/FC;
            RR = FB /FC;
            PP = SS * (2.d0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.d0));
            QQ = (QQ - 1.d0) * (RR - 1.d0) * (SS - 1.d0);
          endif
          if (PP > nearzero) QQ = -QQ;
          PP = abs(PP);
          if ((2.d0*PP)<Minimum(3.d0*xm *QQ-abs(Tol1*QQ),abs(EE*QQ))) then
            EE = DD;  DD = PP/QQ;
          else 
            DD = xm;   EE = DD;
          endif
        else 
          DD = xm;
          EE = DD;
        endif
        AA = BB;
        FA = FB;
        if (abs(DD) > Tol1) then 
          BB = BB + DD;
        else 
          if (xm > 0) then 
            BB = BB + abs(Tol1)
          else 
            BB = BB - abs(Tol1)
          endif
        endif
        FB = AFunction(BB,d1,d2,d3,d4,d5)
        i=i+1
      endif
   end do
    if (i >= maxIterations) error = -2
  endif
  niter = i
  x1 = resultat
  BrentRoots = error
  return
end function BrentRoots

! TRUE if x1*x2 negative
integer Function RootBracketed(x1,x2)
  use module_parameters
  implicit none
  real(DP) x1,x2 
  integer resultat
  if ((x1 > 0.and.x2 > 0).or.(x1 < 0.and.x2 < 0)) then 
    resultat = 0
  else
    resultat = 1
  endif
  RootBracketed = resultat
end function RootBracketed

! returns the minimum of two real numbers
real(DP) Function Minimum(x1,x2) 
  use module_parameters
  implicit none
  real(DP) x1,x2,resultat
  if (x1 < x2) then
    resultat = x1
  else 
    resultat = x2
  endif
  Minimum = resultat
end function Minimum
