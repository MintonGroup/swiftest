submodule (util) s_util_regime
contains
   module procedure util_regime
   !! author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and  David A. Minton
   !!
   !! Determines the collision regime defined by Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. 
   !!                                               I. Outcome Regimes and Scaling Laws 745, 79. 
   !!                                               https://doi.org/10.1088/0004-637X/745/1/79
   !! 
   use swiftest
   real(DP)              :: a1, alpha, aint, b, bcrit, e, fgamma, l, lint, mu, phi, theta
   real(DP)              :: qr, qrd_pstar, qr_erosion, qr_supercat
   real(DP)              :: vcr, verosion, vescp, vhill, vimp, vsupercat
   real(DP)              :: mint, mtot
   real(DP)              :: rp, rhill 
! constants
   integer(I4B)            :: n1 = 1  !number of objects with mass equal to the largest remnant from LS12
   integer(I4B)            :: n2 = 2  !number of objects with mass larger than second largest remnant from LS12
   !integer(I4B)            :: n1g = 2  !number of objects with mass equal to the largest remnant from LS12 if mp = mtarg
   !integer(I4B)            :: n2g = 4  !number of objects with mass larger than second largest remnant from LS12 if mp = mtarg
   real(DP)              :: density1 = 1000.0_DP !standard density parameter from LS12 [kg/m3]
   real(DP)              :: g = 6.674e-11 !gravitational constant [nm2/kg2]
   real(DP)              :: c_star = 1.8_DP !3.0 #3.0# #5#1.8 #1.8 #measure of dissipation of energy within the target (chambers frag.f90)
   real(DP)              :: mu_bar = 0.37_DP !0.385#0.37#0.3333# 3.978 # 1/3 material parameter for hydrodynamic massive body-size bodies (LS12)
   real(DP)              :: beta = 2.85_DP !slope of sfd for remnants from LS12 2.85
   real(DP)              :: c1 = 2.43_DP !LS12 constants
   real(DP)              :: c2 = -0.0408_DP !LS12 constants
   real(DP)              :: c3 = 1.86_DP !LS12 constants
   real(DP)              :: c4 = 1.08_DP !LS12 constants
   real(DP)              :: c5 = 2.5_DP !LS12 constants
   real(DP)              :: crufu = (2.0_DP-3.0_DP*0.36_DP) ! central potential variable from rufu et al. 2019

! executable code
    vimp = norm2(vh2(:) - vh1(:))
    b = calc_b(xh2, vh2, rad2, xh1, vh1, rad1)
    l = (rad1 + rad2)*(1-b)
    e = (norm2(vh1)**2.0_DP)/2.0_DP - g*mcenter/norm2(xh1)
    a1 = - g*mcenter/2.0_DP/e
    mtot = m1 + m2 
    mu = (m1*m2)/mtot
    if (l < 2*rad2) then
        !calculate mint
        phi = 2.0_DP * acos((l - rad2) / rad2)
        aint = (rad2 ** 2.0_DP) * (pi - (phi - sin(phi)) / 2.0_DP)
        lint = 2.0_DP * (rad2 ** 2.0_DP - (rad2 - l / 2.0_DP) ** 2.0_DP) ** (1.0_DP/2.0_DP)
        mint = aint * lint  ![kg]
        alpha = (l**2.0_DP)*(3.0_DP*rad2-l)/(4.0_DP*(rad2**3.0_DP))
    else
       alpha = 1.0_DP
       mint = m2
    end if 
    rp = (3.0_DP*(m1/den1+alpha*m2/den2)/(4.0_DP * pi))**(1.0_DP/3.0_DP) ! (mustill et al. 2019)
   !calculate vescp
    vescp = sqrt(2.0_DP*g*(mtot)/(rp)) !mustill et al. 2018 eq 6 
   !calculate rhill
    rhill = a1*(m1/3.0_DP/(mcenter+m1))**(1.0_DP/3.0_DP)
   !calculate vhill
    if ((rad2 + rad1) < rhill) then 
      vhill = sqrt(2.0_DP * g * m1 * ((rhill ** 2.0_DP - rhill * (rad1 + rad2)) / &
      (rhill ** 2.0_DP - 0.5_DP * (rad1+rad2) ** 2.0_DP)) / (rad1+rad2))
    else
      vhill = vescp
    end if 
   !calculate qr_pstar
    qrd_pstar = calc_qrd_pstar(m1, m2, alpha)*(vhill/vescp)**crufu !rufu et al. eq (3)
   !calculate verosion
    qr_erosion = 2.0_DP * (1.0_DP - m1 / mtot) * qrd_pstar
    verosion = (2.0_DP * qr_erosion * mtot / mu)** (1.0_DP / 2.0_DP)
    qr = mu*(vimp**2.0_DP)/mtot/2.0_DP
    !qrd_lr = calc_qrd_lr(m2, m1, mint)
   !calculate mass largest remnant mlr 
    mlr = (1.0_DP - qr / qrd_pstar / 2.0_DP) * (mtot)  ! [kg] #(eq 5)
   !calculate vsupercat
    qr_supercat = 1.8_DP * qrd_pstar
    vsupercat = ( 2.0_DP * qr_supercat * mtot / mu ) ** (1.0_DP / 2.0_DP)
   !calculate vcr
    fgamma = (m1 - m2) / mtot
    theta = 1.0_DP - b
    vcr = vescp * (c1 * fgamma * theta ** c5 + c2 * fgamma + c3 * theta ** c5 + c4)
    bcrit = rad1/(rad1+rad2)

    if( vimp < vescp) then
      regime = collresolve_regime_merge !perfect merging regime
       mlr = mtot
       mslr = 0.0_DP
    else if (vimp < verosion) then 
      if (b<bcrit) then
      regime = collresolve_regime_merge !partial accretion regime"
       mlr = mtot
       mslr = 0.0_DP
      else if ((b>bcrit) .and. (vimp < vcr)) then
      regime = collresolve_regime_merge ! graze and merge
       mlr = mtot
       mslr = 0.0_DP
      else
       mlr = m1
       mslr = calc_qrd_rev(m2,m1,mint,den1,den2,vimp)
       regime = collresolve_regime_hit_and_run !hit and run
      end if 
    else if (vimp > verosion .and. vimp < vsupercat) then
      if ((m2 < 0.001_DP * m1)) then 
      regime = collresolve_regime_merge !cratering regime"
       mlr = mtot
       mslr = 0.0_DP
      else 
       mslr = (mtot * ((3.0_DP - beta) * (1.0_DP - (n1 * mlr / mtot)))) / (n2 * beta)  ! (eq 37)
       regime = collresolve_regime_disruption !disruption
      end if 
    else if (vimp > vsupercat) then 
       mlr = mtot * (0.1_DP * ((qr / (qrd_pstar * 1.8_DP)) ** (-1.5_DP)))   !eq (44)
       mslr = (mtot * ((3.0_DP - beta) * (1.0_DP - (n1 * mlr / mtot)))) / (n2 * beta)  ! (eq 37)
      regime = collresolve_regime_supercatastrophic ! supercatastrophic
    else 
      write(*,*) "error no regime found in util_regime"
    end if 
    return 

   end procedure util_regime

   module procedure calc_qrd_pstar
      use swiftest
      implicit none
      real(DP)        :: qrd_star1, mu_alpha, mu, qrd_star, qrd_pstar
      ! calc mu, mu_alpha
      mu = (mtarg * mp) / (mtarg + mp)  ! [kg]
      mu_alpha = (mtarg * alpha * mp) / (mtarg + alpha * mp)  ! [kg]
      ! calc qrd_star1
      qrd_star1 = (c_star * 4.0_DP * pi * density1 * g * (rp ** 2.0_DP)) / 5.0_DP
      ! calc qrd_star
      qrd_star = qrd_star1 * (((mp / mtarg + 1.0_DP) ** 2.0_DP) / (4.0_DP * mp / mtarg)) ** (2.0_DP / (3.0_DP * mu_bar) - 1.0_DP)  !(eq 23)
      ! calc qrd_pstar, v_pstar
      qrd_pstar = ((mu / mu_alpha) ** (2.0_DP - 3.0_DP * mu_bar / 2.0_DP)) * qrd_star  ! (eq 15)
      
      ans = qrd_pstar
      return
   end procedure calc_qrd_pstar

   module procedure calc_qrd_rev
      use swiftest
      implicit none
      real(DP) :: mtot_rev, mu_rev, gamma_rev, qrd_star1, qrd_star, mu_alpha_rev
      real(DP) :: qrd_pstar, rc1, qr_rev, qrd_pstar_rev, mslr, qr_supercat_rev
      ! calc mtlr, rc1, mu, gammalr
      mtot_rev =  mint + mp
      rc1 = (3.0_DP*(mint/den1+mp/den2)/(4.0_DP * pi))**(1.0_DP/3.0_DP) ! [m] mustill et al 2018
      mu_rev = (mint * mp) / mtot_rev ! [kg] eq 49 LS12
      mu_alpha_rev = (mtarg * alpha * mp) / (mtarg + alpha * mp)
      gamma_rev = mint / mp ! eq 50 LS12
      !calc qr_rev
      qr_rev = mu_rev * (vimp ** 2.0_DP) / (2.0_DP * mtot_rev)
      ! calc qrd_star1, v_star1
      qrd_star1 = (c_star * 4.0_DP * pi * mtot_rev * g ) / rc1 / 5.0_DP
      ! calc qrd_pstar_rev
      qrd_star = qrd_star1 * (((gamma_rev + 1.0_DP) ** 2.0_DP) / (4.0_DP * gamma_rev)) ** (2.0_DP / (3.0_DP * mu_bar) - 1.0_DP) !(eq 52)
      qrd_pstar = qrd_star * ((mu_rev / mu_alpha_rev) ** (2.0_DP - 3.0_DP * mu_bar / 2.0_DP))
      qrd_pstar_rev = qrd_pstar *(vhill/vescp)**crufu !rufu et al. eq (3)
      !calc qr_supercat_rev
      qr_supercat_rev = 1.8_DP * qrd_pstar_rev 
      !v_supercat_rev = ( 2.0_DP * qr_supercat_rev * mtot_rev / mu_rev ) ** (1.0_DP / 2.0_DP)
      if (qr_rev > qr_supercat_rev ) then 
      mslr = mtot_rev * (0.1_DP * ((qr_rev / (qrd_pstar_rev * 1.8_DP)) ** (-1.5_DP)))   !eq (44)
      else if ( qr_rev < qrd_pstar_rev ) then 
      mslr = mp 
      else 
      mslr = (1.0_DP - qr_rev / qrd_pstar_rev / 2.0_DP) * (mtot_rev)  ! [kg] #(eq 5)
      end if 

      if ( mslr > mp ) mslr = mp !check conservation of mass
      ans = mslr

      return
   end procedure calc_qrd_rev

   module procedure calc_b
      use swiftest
      implicit none
      real(DP) :: h_sq, dvel_sq
      real(DP), dimension(3) :: DPos, dvel, h

      DPos(1) = mtarg_pos(1) - mp_pos(1)
      DPos(2) = mtarg_pos(2) - mp_pos(2)
      DPos(3) = mtarg_pos(3) - mp_pos(3)

      dvel(1) = mtarg_vel(1) - mp_vel(1)
      dvel(2) = mtarg_vel(2) - mp_vel(2)
      dvel(3) = mtarg_vel(3) - mp_vel(3)

      h(1) = (DPos(2) * dvel(3)) - (DPos(3) * dvel(2))
      h(2) = (DPos(3) * dvel(1)) - (DPos(1) * dvel(3))
      h(3) = (DPos(1) * dvel(2)) - (DPos(2) * dvel(1))

      h_sq = (h(1) * h(1)) + (h(2) * h(2)) + (h(3) * h(3))
      dvel_sq = (dvel(1) * dvel(1)) + (dvel(2) * dvel(2)) + (dvel(3) * dvel(3))

      b = (h_sq / (((mp_r + mtarg_r) ** 2.0_DP) * dvel_sq)) ** (1.0_DP / 2.0_DP)

      return
   end procedure calc_b

end submodule s_util_regime
