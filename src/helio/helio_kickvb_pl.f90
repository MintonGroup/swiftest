submodule (helio) s_helio_kickvb
contains
   module procedure helio_kickvb
   !! author: David A. Minton
   !!
   !! Kick barycentric velocities of plAnets
   !! Includes vectorized version
   !!
   !! Adapted from David E. Kaufmann's Swifter routine helio_kickvb.f90
   !! Adapted from Hal Levison's Swift routine kickvh.f 
   use swiftest
   integer(I4B)          :: i, npl

   helio_plA%vb(:,2:npl) = helio_plA%vb(:,2:npl) + helio_plA%ah(:,2:npl) * dt

   return

   end procedure helio_kickvb
end submodule s_helio_kickvb