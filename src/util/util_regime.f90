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
SUBROUTINE util_regime(symba_plA, m1, m2, rad1, rad2, xh1, xh2, vh1, vh2, index1, index2, regime, Mlr, Mslr)

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
     TYPE(symba_pl), INTENT(INOUT) :: symba_plA
     INTEGER(I4B), INTENT(IN)      :: index1, index2
     INTEGER(I4B), INTENT(OUT)     :: regime
     REAL(DP), INTENT(OUT)         :: Mlr, Mslr, m1, m2, rad1, rad2
     REAL(DP), DIMENSION(NDIM)     :: xh1, xh2, vh1, vh2

! Internals
     REAL(DP)                      :: b,l,mu,Vescp,V_pstar, Rp, mtot, RC1
     REAL(DP)                      :: alpha, QRD_pstar, QR, QR_supercat, QRD_lr, V_lr, vimp, bcrit
     REAL(DP)                      :: Vcr, V_supercat, Mint, Lint, Aint, fgamma, theta, rtarg, n2g, n2, phi
     REAL(DP)                      :: c1, c2,c3,c4,c5, rho1, rho2, beta, c_star, density1, g, mp, mtarg, mu_bar, n1, n1g
     REAL(DP), DIMENSION(3)        :: ans
! Constants
     density1 = 1000 ![kg/m3]
     G = 6.674e-11 !Gravitational constant [Nm2/kg2]
     c_star = 1.8!3.0 #3.0# #5#1.8 #1.8 #Measure of dissipation of energy within the target (Chambers frag.f90)
     mu_bar = 0.37!0.385#0.37#0.3333# 3.978 # 1/3 material parameter for hydrodynamic planet-size bodies (LS12)
     beta = 2.85 !slope of SFD for remnants from LS12 2.85
     N1 = 1  !number of objects with mass equal to the largest remnant from LS12
     N2 = 2  !number of objects with mass larger than second largest remnant from LS12
     N1g = 2  !number of objects with mass equal to the largest remnant from LS12 if Mp = Mtarg
     N2g = 4  !number of objects with mass larger than second largest remnant from LS12 if Mp = Mtarg
     c1 = 2.43
     c2 = -0.0408
     c3 = 1.86
     c4 = 1.08
     c5 = 5.0_DP/2.0_DP



! Executable code


! shouldn't this be mass of index1 + index1 children?
     !m1 = symba_plA%helio%swiftest%mass(index1)/GC
     !m2 = symba_plA%helio%swiftest%mass(index2)/GC
     !xh1(:) = symba_plA%helio%swiftest%xh(:,index1)
     !xh2(:) = symba_plA%helio%swiftest%xh(:,index2)
     !vh1(:) = symba_plA%helio%swiftest%vh(:,index1)
     !vh2(:) = symba_plA%helio%swiftest%vh(:,index2)
     Vimp = NORM2(vh2(:) - vh1(:))
     !rad1 = symba_plA%helio%swiftest%radius(index1)
     !rad2 = symba_plA%helio%swiftest%radius(index2)
     b = calc_b(xh2, vh2, rad2, xh1, vh1, rad1)
     l = (rad1 + rad2)*(1-b)
     IF (l < 2*rad2) THEN
          alpha = (l**2.0_DP)*(3*rad2-l)/(4*(rad2**3.0_DP))
     ELSE
          alpha = 1.0_DP
     END IF 

     mtot = m1 + m2 
     Rp = (3*(m1+alpha*m2)/(4.0_DP * PI * rho1)**(1.0_DP/3.0_DP))
     
     Vescp = SQRT(2*GC*mtot/(rad1+rad2))

     ans(:) = calc_erosion(m1, m2, alpha)
     QRD_pstar = ans(1)
     V_pstar = ans(2)
     V_lr = ans(3)

     mu = (m1*m2)/mtot
     QR = 0.5*mu*(vimp**2.0_DP)/mtot
     Mlr = (1.0_DP - 0.5 * QR / QRD_pstar) * (mtot)  ! [kg] #(Eq 5)
     Phi = 2.0_DP * ACOS((l - Rp) / Rp)
     Aint = (Rp ** 2.0_DP) * (PI - (Phi - sin(Phi)) / 2.0_DP)
     Lint = 2.0_DP * (Rtarg ** 2.0_DP - (Rtarg - l / 2.0_DP) ** 2.0_DP) ** (1.0_DP/2.0_DP)
     Mint = Aint * Lint  ![kg]
     QRD_lr = calc_QRD_lr(m2, m1, Mint)
     QR_supercat = 1.8_DP * QRD_pstar
     V_supercat = ((QR_supercat * mtot)/ (0.5 * mu)) ** (1.0_DP / 2.0_DP)
     fgamma = (m1 - m2) / mtot
     theta = 1 - b
     Vcr = Vescp * (c1 * fgamma * theta ** c5 + c2 * fgamma + c3 * theta ** c5 + c4)
     bcrit = rad1/(rad1+rad2)

     IF( Vimp < Vescp) THEN
          WRITE(*,*) "regime perfect merging regime"
          regime = MERGED
     ELSE IF (b < bcrit) THEN
          WRITE(*,*) "regime non grazing"
          IF (Vimp < V_lr) THEN
               WRITE(*,*) "regime partial accretion regime"
               regime = MERGED
          ELSE IF (Vimp > V_lr .AND. Vimp < V_supercat) THEN
               IF (m2 < 1e-3 * m1) THEN
                    WRITE(*,*) "regime cratering"
                    regime = MERGED
               ELSE 
                    WRITE(*,*) "regime disruption"
                    Mslr = ((m1 + m2) * ((3.0_DP - beta) * (1.0_DP - (N1 * Mlr / Mtot)))) / (N2 * beta)  ! (Eq 37)
                    regime = DISRUPTION
               END IF
          ELSE IF (Vimp > V_supercat) THEN 
               WRITE(*,*) "regime supercatastrophic"
               Mlr = Mtot * (0.1_DP * ((QR / (QRD_pstar * 1.8_DP)) ** (-1.5_DP)))     !Eq (44)
               Mslr = ((m1 + m2) * ((3.0_DP - beta) * (1.0_DP - (N1 * Mlr / Mtot)))) / (N2 * beta)  ! (Eq 37)
               regime = SUPERCATASTROPHIC
          ELSE 
               WRITE(*,*) "error"
          END IF
    ELSE IF  (b > bcrit) THEN
          WRITE(*,*) "regime grazing"
          IF (Vimp < V_lr) THEN
               IF (Vimp < Vcr) THEN
                    WRITE(*,*) "regime graze and merge"
                    regime = MERGED
               ELSE
                    WRITE(*,*) "regime hit and run"
                    Mlr = Mtarg
                    Mslr = (Mp + Mint) * (1.0_DP - 0.5_DP * QR / QRD_lr)
                    regime = HIT_AND_RUN
               END IF
               
          ELSE IF (Vimp > V_lr .AND. Vimp < V_supercat) THEN
               IF (m2 < 1e-3 * m1) THEN 
                    WRITE(*,*) "regime cratering"
                    regime = MERGED
               ELSE
                    WRITE(*,*) "regime disruption"
                     Mslr = ((m1 + m2) * ((3.0_DP - beta) * (1.0_DP - (N1 * Mlr / Mtot)))) / (N2 * beta)  ! (Eq 37)
                     regime = DISRUPTION
               END IF 
          ELSE IF (Vimp > V_supercat) THEN 
               WRITE(*,*) "regime supercatastrophic"
               Mlr = Mtot * (0.1_DP * ((QR / (QRD_pstar * 1.8_DP)) ** (-1.5_DP)))     !Eq (44)
               Mslr = ((m1 + m2) * ((3.0_DP - beta) * (1.0_DP - (N1 * Mlr / Mtot)))) / (N2 * beta)  ! (Eq 37)
               regime = SUPERCATASTROPHIC
          END IF 
     END IF 

    RETURN 
! Functions
contains
function calc_erosion(Mtarg,Mp,alpha) result(ans)
   implicit none
   real(DP),intent(in) :: Mtarg, Mp, alpha
   real(DP) :: QRD_star1, V_star1, RC1, mu_alpha, mu, QRD_star, mu_bar, V_star, QRD_pstar, V_pstar, QR, Verosion, G
   real(DP), DIMENSION(3) :: ans
   mu_bar = 0.37
   G = 6.674e-11
   ! calc RC1
   RC1 = (3.0_DP * (Mp + Mtarg)) /  ((4.0_DP * PI * density1) ** (1.0_DP / 3.0_DP))  ! [m]
   ! calc mu, mu_alpha
   mu = (Mtarg * Mp) / (Mtarg + Mp)  ! [kg]
   mu_alpha = (Mtarg * alpha * Mp) / (Mtarg + alpha * Mp)  ! [kg]
   ! calc QRD_star1, V_star1
   QRD_star1 = (c_star * 4.0_DP * PI * density1 * G * (RC1 ** 2.0_DP)) / 5.0_DP
   V_star1 = ((2.0_DP * QRD_star1 * (Mtarg + Mp)) /  mu) ** (1.0_DP / 2.0_DP)
   ! calc QRD_star, V_star
   QRD_star = QRD_star1 * (((Mp / Mtarg + 1.0_DP) ** 2.0_DP) / (4.0_DP * Mp / Mtarg)) ** (2.0_DP / (3.0_DP * mu_bar) - 1.0_DP)  !(Eq 23)
   V_star = V_star1 * (((Mp / Mtarg + 1.0_DP) ** 2.0_DP) / (4.0_DP * Mp / Mtarg)) ** (1.0_DP / (3.0_DP * mu_bar))  ! (Eq 22)
   ! calc QRD_pstar, V_pstar
   QRD_pstar = ((mu / mu_alpha) ** (2.0_DP - 3.0_DP * mu_bar / 2.0_DP)) * QRD_star  ! (Eq 15)
   V_pstar = (2.0_DP * QRD_pstar * (Mtarg + Mp) / mu) ** (1.0_DP / 2.0_DP)  ! (Eq 16)
   ! calc QR, Verosion
   QR = 2.0_DP * ((1.0_DP - Mtarg) / (Mtarg + Mp)) * QRD_pstar
   Verosion = ((2.0_DP * QR * (Mtarg + Mp)) / (Mtarg * Mp / (Mtarg + Mp))) ** (1.0_DP / 2.0_DP)

   ans(1) = QRD_pstar
   ans(2) = V_pstar
   ans(3) = Verosion
   return
end function calc_erosion

function calc_QRD_lr(Mp,Mtarg,Mint) result(ans)
   implicit none
   real(DP),intent(in) :: Mp, Mtarg, Mint
   real(DP) :: ans, Mtlr, mu, gammalr, QRD_star1, c_star, G, V_star1, QRD_lr, mu_bar, QR, V_lr
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
   ! calc QRD_lr, QR, V_lr
   QRD_lr = QRD_star1 * (((gammalr + 1.0_DP) ** 2.0_DP) / (4.0_DP * gammalr)) ** (2.0_DP / (3.0_DP * mu_bar) - 1.0_DP) !(Eq 52)
   QR = 2.0_DP * (1.0_DP - (Mint / Mtlr)) * QRD_lr
   V_lr = ((2.0_DP * QR * Mtlr) / mu) ** (1.0_DP / 2.0_DP)     !(Eq 54)

   ans = V_lr
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
