submodule (swiftest_classes) s_orbel_el2xv
contains
   module procedure orbel_el2xv
   !! author: David A. Minton
   !!
   !! Compute osculating orbital elements from relative Cartesian position and velocity
   !!  All angular measures are returned in radians
   !!      If inclination < TINY, longitude of the ascending node is arbitrarily set to 0
   !!
   !!      If eccentricity < sqrt(TINY), argument of pericenter is arbitrarily set to 0
   !!
   !!      ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
   !!
   !! Adapted from Martin Duncan's el2xv.f
   !! DATE WRITTEN:  May 11, 1992.
   !! REVISIONS: May 26 - now use better Kepler solver for ellipses
   !!  and hyperbolae called EHYBRID.F and FHYBRID.F
   use swiftest
   implicit none
   real(DP), dimension(NDIM) :: x, v

  end procedure orbel_el2xv
end submodule s_orbel_el2xv
