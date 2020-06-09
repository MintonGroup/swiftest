!**********************************************************************************************************************************
!
!  Unit Name   : util_regime
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : util
!  Language    : Fortran 90/95
!
!  Description : Index input real array into ascending numerical order using Quicksort algorithm
!
!  Input
!    Arguments : arr   : array to index
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : index : index table for sorted array
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL util_regime(symba_plA, index1, index2)
!
!  Notes       : Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
!                Vetterling, and Flannery, 2nd ed., pp. 1173-4
!
!**********************************************************************************************************************************
SUBROUTINE util_regime(Mcenter, m1, m2, rad1, rad2, xh1, xh2, vh1, vh2, den1, den2, regime, Mlr, Mslr)

! Modules
     USE swiftest
     USE module_symba
     USE module_swiftest
     USE module_helio
     USE module_nrutil
     USE module_swiftestalloc
     USE module_interfaces, EXCEPT_THIS_ONE => util_regime
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(OUT)                 :: regime
     REAL(DP), INTENT(OUT)                     :: Mlr, Mslr 
     REAL(DP), INTENT(IN)                      :: Mcenter, m1, m2, rad1, rad2, den1, den2  
     REAL(DP), DIMENSION(NDIM), INTENT(IN)     :: xh1, xh2, vh1, vh2

! Internals
     REAL(DP)                      :: alpha, Aint,b,bcrit,fgamma,l,Lint, mu, phi, theta
     REAL(DP)                      :: QRD_pstar, QR, QR_supercat, QRD_lr
     REAL(DP)                      :: vcr, verosion, vescp, vhill, vimp, vsupercat
     REAL(DP)                      :: mint, mtot
     REAL(DP)                      :: Rp, Rhill, a1, E
! Constants
     INTEGER(I4B)                  :: N1 = 1  !number of objects with mass equal to the largest remnant from LS12
     INTEGER(I4B)                  :: N2 = 2  !number of objects with mass larger than second largest remnant from LS12
     INTEGER(I4B)                  :: N1g = 2  !number of objects with mass equal to the largest remnant from LS12 if Mp = Mtarg
     INTEGER(I4B)                  :: N2g = 4  !number of objects with mass larger than second largest remnant from LS12 if Mp = Mtarg
     REAL(DP)                      :: density1 = 1000.0_DP !standard density parameter from LS12 [kg/m3]
     REAL(DP)                      :: G = 6.674e-11 !Gravitational constant [Nm2/kg2]
     REAL(DP)                      :: c_star = 1.8_DP !3.0 #3.0# #5#1.8 #1.8 #Measure of dissipation of energy within the target (Chambers frag.f90)
     REAL(DP)                      :: mu_bar = 0.37_DP !0.385#0.37#0.3333# 3.978 # 1/3 material parameter for hydrodynamic planet-size bodies (LS12)
     REAL(DP)                      :: beta = 2.85_DP !slope of SFD for remnants from LS12 2.85
     REAL(DP)                      :: c1 = 2.43_DP !Ls12 constants
     REAL(DP)                      :: c2 = -0.0408_DP !Ls12 constants
     REAL(DP)                      :: c3 = 1.86_DP !Ls12 constants
     REAL(DP)                      :: c4 = 1.08_DP !Ls12 constants
     REAL(DP)                      :: c5 = 2.5_DP !Ls12 constants
     REAL(DP)                      :: crufu = (2.0_DP-3.0_DP*0.36_DP) ! central potential variable from Rufu et al. 2019

! Executable code

     Vimp = NORM2(vh2(:) - vh1(:))
     b = calc_b(xh2, vh2, rad2, xh1, vh1, rad1)
     l = (rad1 + rad2)*(1-b)
     IF (l < 2*rad2) THEN
          alpha = (l**2.0_DP)*(3*rad2-l)/(4*(rad2**3.0_DP))
     ELSE
          alpha = 1.0_DP
     END IF 

! shouldn't this be mass of index1 + index1 children?
      Vimp = NORM2(vh2(:) - vh1(:))
      b = calc_b(xh2, vh2, rad2, xh1, vh1, rad1)
      l = (rad1 + rad2)*(1-b)
      E = (NORM2(vh1)**2)/2 - G*Mcenter/NORM2(xh1)
      a1 = - G*M_Planet/2/E
      mu = (m1*m2)/mtot
      IF (l < 2*rad2) THEN
           alpha = (l**2.0_DP)*(3*rad2-l)/(4*(rad2**3.0_DP))
      ELSE
           alpha = 1.0_DP
      END IF 
      mtot = m1 + m2 
      Rp = (3*(m1/den1+alpha*m2/den2)/(4.0_DP * PI))**(1.0_DP/3.0_DP) ! (Mustill et al. 2019)
     !Calculate Vescp
      vescp = SQRT(2*GC*(m1+alpha*m2)/(Rp))
     !Calculate Rhill
      Rhill = a1*(m1/3/(Mcenter+m1))**(1.0_DP/3.0_DP)
     !Calculate Vhill
      if ((rad2 + rad1) < Rhill) then 
        vhill = sqrt(2.0_DP * G * m1 * ((Rhill ** 2.0_DP - Rhill * (rad1 + rad2)) / (Rhill ** 2.0_DP - 0.5_DP * (rad1+rad2) ** 2.0_DP)) / (rad1+rad2))
      else
        vhill = vesc_p
      end if 
     !Calculate QR_pstar
      QRD_pstar = calc_erosion(m1, m2, alpha)*(vhill/vesc_p)**crufu !rufu et al. eq (3)
     Write(*,*) "QRD_pstar", QRD_pstar
     !Calculate verosion
      QR_erosion = 2.0_DP * ((1.0_DP - m1) / mtot) * QRD_pstar
      verosion = (2.0_DP * QR_erosion * mtot / mu)** (1.0_DP / 2.0_DP)
     Write(*,*) "verosion", verosion 
      QR = 0.5*mu*(vimp**2.0_DP)/mtot
     Write(*,*) "QR", QR
     !Calculate Mass largest remnant Mlr 
      Mlr = (1.0_DP - 0.5 * QR / QRD_pstar) * (mtot)  ! [kg] #(Eq 5)
     Write(*,*) "Mlr", Mlr 
     !Calculate vsupercat
      QR_supercat = 1.8_DP * QRD_pstar
      vsupercat = ( 2.0_DP * QR_supercat * mtot / mu ) ** (1.0_DP / 2.0_DP)
     !Calculate Vcr
      fgamma = (m1 - m2) / mtot
      theta = 1 - b
      vcr = vescp * (c1 * fgamma * theta ** c5 + c2 * fgamma + c3 * theta ** c5 + c4)
      bcrit = rad1/(rad1+rad2)
     !Calculate mint
      Phi = 2.0_DP * ACOS((l - rad2) / rad2)
      Aint = (rad2 ** 2.0_DP) * (PI - (Phi - sin(Phi)) / 2.0_DP)
      Lint = 2.0_DP * (rad2 ** 2.0_DP - (rad2 - l / 2.0_DP) ** 2.0_DP) ** (1.0_DP/2.0_DP)
      mint = Aint * Lint  ![kg]
      IF( vimp < vescp) THEN
        regime = COLRESOLVE_REGIME_MERGE !perfect merging regime
      ELSE IF (vimp < verosion) THEN 
        IF (b<bcrit) THEN 
          regime = COLRESOLVE_REGIME_MERGE !partial accretion regime"
        ELSE IF ((b>bcrit) .AND. (Vimp < Vcr)) THEN 
          regime = COLRESOLVE_REGIME_MERGE ! graze and merge
        ELSE 
          Mlr = m1
          Mslr = (m2 + mint) * (1.0_DP - 0.5_DP * QR / QRD_lr)
          regime = COLLRESOLVE_REGIME_HIT_AND_RUN !hit and run
        END IF 
      ELSE IF (vimp > verosion .AND. vimp < vsupercat) THEN 
        IF ((m2 < 1e-3 * m1)) THEN 
          regime = COLRESOLVE_REGIME_MERGE !cratering regime"
        ELSE 
          Mslr = (mtot * ((3.0_DP - beta) * (1.0_DP - (N1 * Mlr / mtot)))) / (N2 * beta)  ! (Eq 37)
          regime = COLLRESOLVE_REGIME_DISRUPTION !disruption
        END IF 
      ELSE IF (vimp > vsupercat) THEN 
        Mlr = mtot * (0.1_DP * ((QR / (QRD_pstar * 1.8_DP)) ** (-1.5_DP)))     !Eq (44)
        Mslr = (mtot * ((3.0_DP - beta) * (1.0_DP - (N1 * Mlr / mtot)))) / (N2 * beta)  ! (Eq 37)
        regime = COLLRESOLVE_REGIME_SUPERCATASTROPHIC ! supercatastrophic
      ELSE 
        WRITE(*,*) "Error no regime found in util_regime"
      END IF 
    RETURN 


! Functions
contains
function calc_erosion(Mtarg,Mp,alpha) result(ans)
   implicit none
   real(DP),intent(in) :: Mtarg, Mp, alpha
   real(DP)            :: QRD_star1, mu_alpha, mu, QRD_star, RD_pstar, QR
   real(DP)            :: ans
   ! calc mu, mu_alpha
   mu = (Mtarg * Mp) / (Mtarg + Mp)  ! [kg]
   mu_alpha = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)  ! [kg]
   ! calc QRD_star1
   QRD_star1 = (c_star * 4.0_DP * PI * density1 * G * (RC1 ** 2.0_DP)) / 5.0_DP
   ! calc QRD_star
   QRD_star = QRD_star1 * (((Mp / Mtarg + 1.0_DP) ** 2.0_DP) / (4.0_DP * Mp / Mtarg)) ** (2.0_DP / (3.0_DP * mu_bar) - 1.0_DP)  !(Eq 23)
   ! calc QRD_pstar, V_pstar
   QRD_pstar = ((mu / mu_alpha) ** (2.0_DP - 3.0_DP * mu_bar / 2.0_DP)) * QRD_star  ! (Eq 15)
   
   ans = QRD_pstar
   return
end function calc_erosion

function calc_QRD_lr(Mp,Mtarg,Mint) result(ans)
   implicit none
   real(DP),intent(in) :: Mp, Mtarg, Mint
   real(DP) :: ans, Mtlr, mu, gammalr, QRD_star1, c_star, G, V_star1, QRD_lr, mu_bar, QR, verosion
   c_star = 1.8
   G = 6.674e-11
   mu_bar = 0.37
   ! calc Mtlr, RC1, mu, gammalr
   Mtlr =  Mint + Mp
   RC1 = ((3.0_DP * Mtlr) / (4.0_DP * PI * density1)) ** (1.0_DP / 3.0_DP) ! [m]
   mu = (Mint * Mp) / (Mint + Mp) ! [kg]
   gammalr = Mint / Mp
   ! calc QRD_star1, V_star1
   QRD_star1 = (c_star * 4.0_DP * PI * density1 * G * (RC1 ** 2.0_DP)) / 5.0_DP
   V_star1 = ((2.0_DP * QRD_star1 * (Mint + Mp)) / mu) ** (1.0_DP / 2.0_DP)
   ! calc QRD_lr, QR, verosion
   QRD_lr = QRD_star1 * (((gammalr + 1.0_DP) ** 2.0_DP) / (4.0_DP * gammalr)) ** (2.0_DP / (3.0_DP * mu_bar) - 1.0_DP) !(Eq 52)
   QR = 2.0_DP * (1.0_DP - (Mint / Mtlr)) * QRD_lr
   verosion = ((2.0_DP * QR * Mtlr) / mu) ** (1.0_DP / 2.0_DP)     !(Eq 54)

   ans = verosion
   return
end function calc_QRD_lr

function calc_b(Mp_pos, Mp_vel, Mp_r, Mtarg_pos, Mtarg_vel, Mtarg_r) result(b)
   implicit none
   real(DP), intent(in), DIMENSION(3) :: Mp_pos, Mp_vel, Mtarg_pos, Mtarg_vel
   real(DP), intent(in) :: Mp_r, Mtarg_r
   real(DP) :: ans, h_sq, b, dvel_sq
   real(DP), DIMENSION(3) :: dpos, dvel, h

   dpos(1) = mtarg_pos(1) - mp_pos(1)
   dpos(2) = mtarg_pos(2) - mp_pos(2)
   dpos(3) = mtarg_pos(3) - mp_pos(3)

   dvel(1) = mtarg_vel(1) - mp_vel(1)
   dvel(2) = mtarg_vel(2) - mp_vel(2)
   dvel(3) = mtarg_vel(3) - mp_vel(3)

   h(1) = (dpos(2) * dvel(3)) - (dpos(3) * dvel(2))
   h(2) = (dpos(3) * dvel(1)) - (dpos(1) * dvel(3))
   h(3) = (dpos(1) * dvel(2)) - (dpos(2) * dvel(1))

   h_sq = (h(1) * h(1)) + (h(2) * h(2)) + (h(3) * h(3))
   dvel_sq = (dvel(1) * dvel(1)) + (dvel(2) * dvel(2)) + (dvel(3) * dvel(3))

   b = (h_sq / (((Mp_r + Mtarg_r) ** 2.0_DP) * dvel_sq)) ** (1.0_DP / 2.0_DP)

   return
end function calc_b

END SUBROUTINE util_regime
!**********************************************************************************************************************************
!
!  Author(s)   : C.Wishard and J.Pouplin
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
