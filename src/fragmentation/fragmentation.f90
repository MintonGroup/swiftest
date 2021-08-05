submodule(swiftest_classes) s_fragmentation
   use swiftest
contains
   module subroutine fragmentation_regime(Mcb, m1, m2, rad1, rad2, xh1, xh2, vb1, vb2, den1, den2, regime, Mlr, Mslr, mtiny, Qloss)
      !! Author: Jennifer L.L. Pouplin, Carlisle A. Wishard, and David A. Minton
      !!
      !! Determine the collisional regime of two colliding bodies. 
      !! Current version requires all values to be converted to SI units prior to calling the function
      !!       References:
      !!       Kokubo, E., Genda, H., 2010. Formation of Terrestrial Planets from Protoplanets Under a Realistic Accretion 
      !!          Condition. ApJL 714, L21. https://doi.org/10.1088/2041-8205/714/1/L21
      !!       Leinhardt, Z.M., Stewart, S.T., 2012. Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling 
      !!          Laws 745, 79. https://doi.org/10.1088/0004-637X/745/1/79
      !!       Mustill, A.J., Davies, M.B., Johansen, A., 2018. The dynamical evolution of transiting planetary systems including 
      !!          a realistic collision prescription. Mon Not R Astron Soc 478, 2896–2908. https://doi.org/10.1093/mnras/sty1273
      !!       Rufu, R., Aharonson, O., 2019. Impact Dynamics of Moons Within a Planetary Potential. J. Geophys. Res. Planets 124, 
      !!          1008–1019. https://doi.org/10.1029/2018JE005798
      !!       Stewart, S.T., Leinhardt, Z.M., 2012. Collisions between Gravity-dominated Bodies. II. The Diversity of Impact 
      !!          Outcomes during the End Stage of Planet Formation. ApJ 751, 32. https://doi.org/10.1088/0004-637X/751/1/32
      !!
      implicit none
   ! Arguments
      integer(I4B), intent(out)         :: regime
      real(DP), intent(out)          :: Mlr, Mslr
      real(DP), intent(in)           :: Mcb, m1, m2, rad1, rad2, den1, den2, mtiny 
      real(DP), dimension(:), intent(in)   :: xh1, xh2, vb1, vb2
      real(DP), intent(out)          :: Qloss !! The residual energy after the collision 

      return
   end subroutine fragmentation_regime

end submodule s_fragmentation