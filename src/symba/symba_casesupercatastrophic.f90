!**********************************************************************************************************************************
!
!  Unit Name   : symba_casesupercatastrophic
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge planets
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-planet encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_casesupercatastrophic(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
   symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

! Modules
   USE swiftest
   USE module_helio
   USE module_symba
   USE module_interfaces, EXCEPT_THIS_ONE => symba_casesupercatastrophic
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
   INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, nplplenc, fragmax
   REAL(DP), INTENT(IN)                             :: t, dt
   REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: mres, rres
   REAL(DP), DIMENSION(:), INTENT(IN)            :: vbs
   REAL(DP), DIMENSION(:), INTENT(INOUT)         :: x1, x2, v1, v2
   TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
   TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
   INTEGER(I4B)                                     :: nfrag, i, k, index1, index2, frags_added
   INTEGER(I4B)                                     :: index1_parent, index2_parent
   INTEGER(I4B)                                     :: name1, name2, nstart
   REAL(DP)                                         :: mtot, msun, avg_d, d_p1, d_p2, semimajor_encounter, e, q, semimajor_inward
   REAL(DP)                                         :: rhill_p1, rhill_p2, r_circle, theta, radius1, radius2, r_smallestcircle
   REAL(DP)                                         :: m_rem, m_test, mass1, mass2, enew, eold, A, B, v_col
   REAL(DP)                                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com
   REAL(DP)                                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag, m1m2_10
   REAL(DP), DIMENSION(NDIM)                        :: vnew, xr, mv, l, kk, p

   !TEMPORARY
   interface 
      function cross_product_supercatastrophic(ar1,ar2) result(ans)
         use swiftest
         implicit none
         real(DP),dimension(3),intent(in) :: ar1,ar2
         real(DP),dimension(3)             :: ans
      end function cross_product_supercatastrophic
   end interface

! Executable code
     
   WRITE(*,*) "ENTERING CASESUPERCATASTROPHIC"
   ! Set the maximum number of fragments to be added in a Supercatastrophic Disruption collision (nfrag)
   nfrag = 10
   ! Pull in the information about the two particles involved in the collision  
   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)
   index1_parent = symba_plA%index_parent(index1)
   index2_parent = symba_plA%index_parent(index2)
   name1 = symba_plA%helio%swiftest%name(index1)
   name2 = symba_plA%helio%swiftest%name(index2)
   mass1 = symba_plA%helio%swiftest%mass(index1) ! the mass of the first particle in the collision NOT INCLUDING all it's children
   mass2 = symba_plA%helio%swiftest%mass(index2)
   radius1 = symba_plA%helio%swiftest%radius(index1)
   radius2 = symba_plA%helio%swiftest%radius(index2)
   msun = symba_plA%helio%swiftest%mass(1)

   ! Find COM
   x_com = ((x1(1) * m1) + (x2(1) * m2)) / (m1 + m2)
   y_com = ((x1(2) * m1) + (x2(2) * m2)) / (m1 + m2)
   z_com = ((x1(3) * m1) + (x2(3) * m2)) / (m1 + m2)

   vx_com = ((v1(1) * m1) + (v2(1) * m2)) / (m1 + m2)
   vy_com = ((v1(2) * m1) + (v2(2) * m2)) / (m1 + m2)
   vz_com = ((v1(3) * m1) + (v2(3) * m2)) / (m1 + m2)

   ! Find Collision velocity
   v_col = NORM2(v2(:) - v1(:))

   ! Find energy pre-frag
   eold = 0.5_DP*(m1*DOT_PRODUCT(v1(:), v1(:)) + m2*DOT_PRODUCT(v2(:), v2(:)))
   xr(:) = x2(:) - x1(:)
   eold = eold - (m1*m2/(SQRT(DOT_PRODUCT(xr(:), xr(:)))))

   WRITE(*, *) "Supercatastrophic disruption between particles ", name1, " and ", name2, " at time t = ",t
     
   ! Add both particles involved in the collision to mergesub_list
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = SUPERCATASTROPHIC
   mergesub_list%xh(:,nmergesub) = x1(:)
   mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = radius1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = SUPERCATASTROPHIC
   mergesub_list%xh(:,nmergesub) = x2(:)
   mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = radius2

   ! Go through the encounter list and look for particles actively encoutering in this timestep
   ! Prevent them from having further encounters in this timestep by setting status in plplenc_list to MERGED
   DO k = 1, nplplenc 
      IF ((plplenc_list%status(k) == ACTIVE) .AND. &
         ((index1 == plplenc_list%index1(k) .OR. index2 == plplenc_list%index2(k)) .OR. &
         (index2 == plplenc_list%index1(k) .OR. index1 == plplenc_list%index2(k)))) THEN
            plplenc_list%status(k) = MERGED
      END IF
   END DO

   ! Set the status of the particles in symba_plA to DISRUPTION
   symba_plA%helio%swiftest%status(index1) = SUPERCATASTROPHIC
   symba_plA%helio%swiftest%status(index2) = SUPERCATASTROPHIC

   l(:) = (v2(:) - v1(:)) / NORM2(v2(:)-v1(:))
   p(:) = cross_product_supercatastrophic(xr(:) / NORM2(xr(:)), l(:))
   kk(:) = cross_product_supercatastrophic(l(:),p(:))

   ! Calculate the positions of the new fragments in a circle with a radius large enough to space
   ! all fragments apart by a distance of rhill_p1 + rhill_p2
   rhill_p1 = symba_plA%helio%swiftest%rhill(index1_parent)
   rhill_p2 = symba_plA%helio%swiftest%rhill(index2_parent)
   r_smallestcircle = (RHSCALE * rhill_p1 + RHSCALE * rhill_p2) / (2.0_DP * sin(PI /2.0_DP))

   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
   ! If they are going to be added interior to this orbit, give a warning
   IF (semimajor_inward > (semimajor_encounter - r_smallestcircle)) THEN
      WRITE(*,*) "WARNING in symba_casesupercatastrophic: Timestep is too large to resolve fragments."
   END IF
   ! If not, continue through all possible fragments to be added
      mtot = 0.0_DP ! running total mass of new fragments
      mv = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
      frags_added = 0 ! running total number of new fragments
      m1m2_10 = 0.1_DP * (m1 + m2) ! one tenth the total initial mass of the system used to check the size of the fragments
      nstart = nmergeadd

      d_p1 = (3.0_DP * m1) / (4.0_DP * PI * (rad1 ** 3.0_DP))
      d_p2 = (3.0_DP * m2) / (4.0_DP * PI * (rad2 ** 3.0_DP))
      avg_d = ((m1 * d_p1) + (m2 * d_p2)) / (m1 + m2)

         ! If we are adding the first and largest fragment (lr), check to see if its mass is SMALLER than one tenth the total
         ! mass of the system aka if it is too small to resolve. If so, add a fragment with a mass of one tenth the total mass 
         ! of the system and calculate its radius.
         IF ((mres(1) < m1m2_10)) THEN
            DO i = 1, nfrag
               frags_added = frags_added + 1
               nmergeadd = nmergeadd + 1
               mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
               mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
               mergeadd_list%ncomp(nmergeadd) = 2
               mergeadd_list%mass(nmergeadd) = m1m2_10
               mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * avg_d))  & 
                  ** (1.0_DP / 3.0_DP)
               mtot = mtot + mergeadd_list%mass(nmergeadd) 
            END DO 
         END  IF
         ! If we are adding the first and largest fragment (lr), check to see if its mass is LARGER than one tenth the total 
         ! mass of the system aka if it is large enough to resolve. If so, its mass and radius should be taken from 
         ! util_regime.
         IF ((mres(1) > m1m2_10)) THEN
            frags_added = frags_added + 1
            nmergeadd = nmergeadd + 1
            mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
            mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%mass(nmergeadd) = mres(1)
            mergeadd_list%radius(nmergeadd) = rres(1)
            mtot = mtot + mergeadd_list%mass(nmergeadd) 
            DO i = 2, nfrag
               frags_added = frags_added + 1
               nmergeadd = nmergeadd + 1
               mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
               mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
               mergeadd_list%ncomp(nmergeadd) = 2
               mergeadd_list%mass(nmergeadd) = (m1 + m2 - mres(1)) / (nfrag - 1.0_DP)
               mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * avg_d))  & 
                  ** (1.0_DP / 3.0_DP)  
               mtot = mtot + mergeadd_list%mass(nmergeadd)
            END DO
         END IF 


         !END IF
         ! If we are adding more than one fragment
         !IF ((i > 1) .AND. (mres(1) > m1m2_10)) THEN
            ! m_rem is the mass needed to be "made up for" in fragments, mres(1) is the mass of the largest fragments 
            ! that has already been added, and m1 and m2 are the masses of the original particles involved in the collision.
         !   m_rem = (m1 + m2) - (mergeadd_list%mass(nmergeadd))
            ! Check if these fragments will be large enough to be resolved
         !   IF (m_rem > (1.0_DP / 10.0_DP)*mres(1))) THEN

               ! If yes, add a fragment using Durda et al 2007 Figure 2 Supercatastrophic: N = (1.5e5)e(-1.3*D) for the mass
         !      frags_added = frags_added + 1
         !      nmergeadd = nmergeadd + 1
         !      mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
         !      mergeadd_list%status(nmergeadd) = SUPERCATASTROPHIC
         !      mergeadd_list%ncomp(nmergeadd) = 2
         !      mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 

               ! Create a "test mass" using Durda et al 2007 Figure 2 Supercatastrophic: N = (1.5e5)e(-1.3*D)
               !m_test = (((- 1.0_DP / 2.6_DP) * log(i / (1.5_DP * 10.0_DP ** 5))) ** 3.0_DP) * ((4.0_DP / 3.0_DP) &
               !   * PI * avg_d)
               ! If the test mass is smaller than the mass that needs to be "made up for", add it.
               !IF (m_test < m_rem) THEN
               !   mergeadd_list%mass(nmergeadd) = m_test
               ! If not, aka if the test mass is too large, then add a fragment with a mass equal to the difference between
               ! the sum of the mass of the parents and the total mass already added. 
               !ELSE
               !   mergeadd_list%mass(nmergeadd) = (m1 + m2) - mtot 
               !END IF

               ! Calculate the radius of the fragments using the weighted average density of the parents. 
         !      mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * avg_d))  & 
         !         ** (1.0_DP / 3.0_DP) 
         !      mtot = mtot + mergeadd_list%mass(nmergeadd)
         !   ELSE 
         !      mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
         !      mergeadd_list%radius(nmergeadd) = (((3.0_DP/4.0_DP) * PI) * (mergeadd_list%mass(nmergeadd) / avg_d)) &
         !         ** (1.0_DP / 3.0_DP)
         !   END IF
         !END IF

   r_circle = (RHSCALE * rhill_p1 + RHSCALE * rhill_p2) / (2.0_DP * sin(PI / frags_added))
   theta = (2.0_DP * PI) / frags_added

   DO i=1, frags_added

      ! Increment around the circle for positions of fragments
      x_frag = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1) + x_com
      y_frag = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
      z_frag = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

      A = v_col * (m1 + m2) * (1.0_DP / mergeadd_list%mass(nstart + i))

      vx_frag = ((A * cos(theta * i))*l(1)) + ((A * sin(theta * i))*p(1)) + vx_com
      vy_frag = ((A * cos(theta * i))*l(2)) + ((A * sin(theta * i))*p(2)) + vy_com
      vz_frag = ((A * cos(theta * i))*l(3)) + ((A * sin(theta * i))*p(3)) + vz_com

         !Conservation of Angular Momentum for velocities of fragments
         !A = (((x1(2) * v1(3) * m1) - (x1(3) * v1(2) * m1)) + ((x2(2) * v2(3) * m2) - (x2(3) * v2(2) * m2))) &
         !   / mergeadd_list%mass(nmergeadd)
         !B = (((x1(3) * v1(1) * m1) - (x1(1) * v1(3) * m1)) + ((x2(3) * v2(1) * m2) - (x2(1) * v2(3) * m2))) &
         !   / mergeadd_list%mass(nmergeadd)
         !vx_frag = ((1.0_DP / frags_added) * (B / z_frag)) - vbs(1)
         !vy_frag = ((1.0_DP / frags_added) * (-A / z_frag)) - vbs(2)
         !vz_frag = vz_com - vbs(3)

      mergeadd_list%xh(1,nstart + i) = x_frag
      mergeadd_list%xh(2,nstart + i) = y_frag 
      mergeadd_list%xh(3,nstart + i) = z_frag                                                    
      mergeadd_list%vh(1,nstart + i) = vx_frag
      mergeadd_list%vh(2,nstart + i) = vy_frag
      mergeadd_list%vh(3,nstart + i) = vz_frag
      ! Tracking linear momentum.
      mv = mv + (mergeadd_list%mass(nstart + i) * mergeadd_list%vh(:,nstart + i))
   END DO

   WRITE(*, *) "Number of fragments added: ", frags_added
   ! Calculate energy after frag                                                                           
   vnew(:) = mv / mtot    ! COM of new fragments                               
   enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
   eoffset = eoffset + eold - enew

   ! Update fragmax to account for new fragments
   fragmax = fragmax + frags_added
   WRITE(*,*) "LEAVING CASESUPERCATASTROPHIC"

   RETURN 
END SUBROUTINE symba_casesupercatastrophic

function cross_product_supercatastrophic(ar1,ar2) result(ans)
   use swiftest
   implicit none
   
   real(DP),dimension(3),intent(in) :: ar1,ar2
   real(DP),dimension(3)             :: ans

   ans(1) = ar1(2) * ar2(3) - ar1(3) * ar2(2)
   ans(2) = ar1(3) * ar2(1) - ar1(1) * ar2(3)
   ans(3) = ar1(1) * ar2(2) - ar1(2) * ar2(1)

   return 
end function cross_product_supercatastrophic
