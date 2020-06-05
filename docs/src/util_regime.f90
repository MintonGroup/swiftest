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
SUBROUTINE util_regime(symba_plA, index1, index2)

! Modules
     USE module_parameters
     USE module_symba
     USE module_swiftest
     USE module_helio
     USE module_nrutil
     USE module_swiftestalloc
     USE module_interfaces, EXCEPT_THIS_ONE => util_resize_pl
     IMPLICIT NONE

! Arguments
     TYPE(symba_pl), INTENT(INOUT) :: symba_plA
     INTEGER(I4B), INTENT(IN)      :: index1, index2

! Internals
     TYPE(symba_pl)                :: 
     REAL(DP)                      :: GU, m1,m2,rad1,rad2,b,l,mu,Vescp,V_pstar, Rp, mtot
     REAL(DP)                      :: alpha, QRD_pstar, QR, QR_supercat, QRD_lr, V_lr, vimp, bcrit
     REAL(DP)                      :: Vcr, V_supercat, Vescp, Mlr, Mint, Lint, Aint, fgamma, theta
     REAL(DP)                      :: c1, c2,c3,c4,c5, rho1, rho2 
     REAL(DP), DIMENSION(NDIM)     :: xh1,xh2,vh1,vh2
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
     c5 = 5/2
! Executable code
     GU = GC * 
     m1 = symba_pl%helio%swiftest%mass(index1)/GU
     m2 = symba_pl%helio%swiftest%mass(index2)/GU
     xh1(:) = symba_pl%helio%swiftest%xh(:,index1)
     xh2(:) = symba_pl%helio%swiftest%xh(:,index2)
     vh1(:) = symba_pl%helio%swiftest%vh(:,index1)
     vh2(:) = symba_pl%helio%swiftest%vh(:,index2)
     vimp = NORM2(vh2(:) - vh1(:))
     rad1 = symba_pl%helio%swiftest%radius(index1)
     rad2 = symba_pl%helio%swiftest%radius(index2)
     b = 
     l = (rad1 + rad2)*(1-b)
     IF (l < 2*rad2) THEN
          alpha = (l**2.0_DP)*(3*rad2-l)/(4*(rad2**3.0_DP))
     ELSE
          alpha = 1.0_DP
     END IF 

     mtot = m1 + m2 
     Rp = (3*(m1+alpha*m2)/(4.0_DP * PI * rho1)**(1.0_DP/3.0_DP))
     
     Vescp = SQRT(2*GU*mtot/(rad1+rad2))
     QRD_pstar = 
     V_pstar = 
     mu = (m1*m2)/mtot
     QR = 0.5*mu*(vimp**2.0_DP)/mtot
     Mlr = (1.0_DP - 0.5 * QR / QRD_pstar) * (mtot)  ! [kg] #(Eq 5)
     Phi = 2.o_DP * ACOS((l - Rp) / Rp)
     Aint = (Rp ** 2.0_DP) * (PI - (Phi - sin(Phi)) / 2.0_DP)
     Lint = 2.0_DP * (Rtarg ** 2.0_DP - (Rtarg - l / 2.0_DP) ** 2.0_DP) ** (1.0_DP/2.0_DP)
     Mint = Aint * Lint  # [kg]
     V_lr = calc_erosion(Mtarg, Mp, alpha)
     QRD_lr, V_lr = calc_QRD_lr(Mp, Mtarg, Mint)
     QR_supercat = 1.8_DP * QRD_pstar
     V_supercat = ((QR_supercat * mtot)/ (0.5 * mu)) ** (1.0_DP / 2.0_DP)
     fgamma = (m1 - m2) / mtot
     theta = 1 - b
     Vcr = Vescp * (c1 * fgamma * theta ** c5 + c2 * fgamma + c3 * theta ** c5 + c4)
     bcrit = rad1/(rad1+rad2)

     IF( Vi < Vescp) THEN
          WRITE(*,*) "regime perfect merging regime"
     ELSE IF (b < bcrit) THEN
          WRITE(*,*) "regime non grazing"
          IF (Vi < V_lr) THEN
               WRITE(*,*) "regime partial accretion regime"
          ELSE IF (vimp > V_lr .AND. Vi < V_supercat) THEN
               IF (m2 < 1e-3 * m1) THEN
                    WRITE(*,*) "regime cratering"
               ELSE 
                    WRITE(*,*) "regime disruption"
               END IF
          ELSE IF (Vi > V_supercat) THEN 
               WRITE(*,*) "regime supercatastrophic"
          ELSE 
               WRITE(*,*) "error"
          END IF
    ELSE IF  (b > bcrit) THEN
          WRITE(*,*) "regime grazing"
          IF (vimp < V_lr) THEN
               IF (Vi < Vcr) THEN
                    WRITE(*,*) "regime graze and merge"
               ELSE
                    WRITE(*,*) "regime hit and run"
               END IF
               
          ELSE IF (vimp > V_lr .AND. Vi < V_supercat) THEN
               IF (m2 < 1e-3 * m1) THEN 
                    WRITE(*,*) "regime cratering"
               ELSE
                    WRITE(*,*) "regime disruption"
               END IF 
          ELSE IF (vimp > V_supercat) THEN 
               WRITE(*,*) "regime supercatastrophic"
          END IF 
     END IF 

     RETURN

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
