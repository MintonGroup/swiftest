!**********************************************************************************************************************************
!
!  Unit Name   : symba_casedisruption
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
!  Invocation  : CALL symba_casedisruption(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
   symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

! Modules
   USE swiftest
   USE module_helio
   USE module_symba
   USE module_interfaces, EXCEPT_THIS_ONE => symba_casedisruption
   IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
   INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, nplplenc, fragmax
   REAL(DP), INTENT(IN)                             :: t, dt
   REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
   REAL(DP), DIMENSION(3), INTENT(INOUT)            :: mres, rres
   REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
   REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
   TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
   TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
   INTEGER(I4B)                                     :: nfrag, i, k, index1, index2, frags_added
   INTEGER(I4B)                                     :: index1_parent, index2_parent
   INTEGER(I4B)                                     :: name1, name2
   REAL(DP)                                         :: mtot, msun, avg_d, d_p1, d_p2, semimajor_encounter, e, q, semimajor_inward
   REAL(DP)                                         :: rhill_p1, rhill_p2, r_circle, theta, radius1, radius2
   REAL(DP)                                         :: m_rem, m_test, mass1, mass2, enew, eold, A, B
   REAL(DP)                                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com
   REAL(DP)                                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag
   REAL(DP), DIMENSION(NDIM)                        :: vnew, xr, mv, l, k, p

! Executable code

   WRITE(*,*) "ENTERING CASEDISRUPTION"
   ! Set the maximum number of fragments to be added in a Disruption collision (nfrag)
   nfrag = 3 
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

   ! Find energy pre-frag
   eold = 0.5_DP*(m1*DOT_PRODUCT(v1(:), v1(:)) + m2*DOT_PRODUCT(v2(:), v2(:)))
   xr(:) = x2(:) - x1(:)
   eold = eold - (m1*m2/(SQRT(DOT_PRODUCT(xr(:), xr(:)))))
   
   WRITE(*, *) "Disruption between particles ", name1, " and ", name2, " at time t = ",t
     
   ! Add both particles involved in the collision to mergesub_list
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = DISRUPTION
   mergesub_list%xh(:,nmergesub) = x1(:)
   mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = radius1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = DISRUPTION
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
   symba_plA%helio%swiftest%status(index1) = DISRUPTION
   symba_plA%helio%swiftest%status(index2) = DISRUPTION

   ! Calculate the positions of the new fragments in a circle with a radius large enough to space
   ! all fragments apart by a distance of rhill_p1 + rhill_p2
   rhill_p1 = symba_plA%helio%swiftest%rhill(index1_parent)
   rhill_p2 = symba_plA%helio%swiftest%rhill(index2_parent)
   r_circle = (rhill_p1 + rhill_p2) / (2.0_DP * sin(PI / nfrag))
   theta = (2.0_DP * PI) / nfrag
   l(:) = (v2(:)-v1(:))/NORM2((v2(:)-v1(:))
   p(:) = CROSS_PRODUCT(xr(:)/NORM2(xr(:), l(:)))
   k(:) = CROSS_PRODUCT(l(:),p(:))
   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
   ! If they are going to be added interior to this orbit, give a warning
   IF (semimajor_inward > (semimajor_encounter - r_circle)) THEN
      WRITE(*,*) "WARNING in symba_casedisruption: Timestep is too large to resolve fragments."
   ELSE
      ! If not, continue through all possible fragments to be added
      mtot = 0.0_DP ! running total mass of new fragments
      mv = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
      frags_added = 0 ! running total number of new fragments
      DO i = 1, nfrag
         ! If we are adding the first and largest fragment (lr), it's mass and radius should be taken from 
         ! util_regime while it's position and velocity should be calculated on the circle of radius 
         ! r_circle as described above.
         IF (i == 1) THEN
            frags_added = frags_added + 1
            nmergeadd = nmergeadd + 1
            mergeadd_list%status(nmergeadd) = DISRUPTION
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
            mergeadd_list%mass(nmergeadd) = mres(1)
            mergeadd_list%radius(nmergeadd) = rres(1)
            mtot = mtot + mergeadd_list%mass(nmergeadd)                             
         END IF
         ! If we are adding the second fragment (slr), it's mass and radius should be taken from 
         ! util_regime while it's position and velocity should be calculated on the circle of 
         ! radius r_circle as described above.
         IF (i == 2) THEN
            ! frags_added is the actual number of fragments added to the simulation vs nfrag which is the total possible
            frags_added = frags_added + 1
            nmergeadd = nmergeadd + 1
            mergeadd_list%status(nmergeadd) = DISRUPTION
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
            mergeadd_list%mass(nmergeadd) = mres(2)!*GC/(DU2CM**3 / (MU2GM * TU2S**2))/MU2GM
            mergeadd_list%radius(nmergeadd) = rres(2)!/DU2CM
            mtot = mtot + mergeadd_list%mass(nmergeadd)                            
         END IF
         ! If we are doing more fragments
         IF (i > 2) THEN
            ! m_rem is the mass needed to be "made up for" in fragments, mres(1) and mres(2) are the mass of the largest 
            ! and second largest fragments that have already been added, and m1 and m2 are the masses of the original 
            ! particles involved in the collision.
            d_p1 = (3.0_DP * m1) / (4.0_DP * PI * (rad1 ** 3.0_DP))
            d_p2 = (3.0_DP * m2) / (4.0_DP * PI * (rad2 ** 3.0_DP))
            avg_d = ((m1 * d_p1) + (m2 * d_p2)) / (m1 + m2)
            m_rem = (m1 + m2) - (mres(1) + mres(2))!*GC/(DU2CM**3 / (MU2GM * TU2S**2))/MU2GM

            ! Check if these fragments will be large enough to be resolved    
            IF (m_rem > (m2) / 1000.0_DP) THEN
               ! If yes, add a fragment using Durda et al 2007 Figure 2 Supercatastrophic: N = (1.5e5)e(-1.3*D) for the mass
               frags_added = frags_added + 1
               nmergeadd = nmergeadd + 1
               mergeadd_list%status(nmergeadd) = DISRUPTION
               mergeadd_list%ncomp(nmergeadd) = 2
               mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i
               ! Create a "test mass" using Durda et al 2007 Figure 2 Supercatastrophic: N = (1.5e5)e(-1.3*D)
               m_test = (((- 1.0_DP / 2.6_DP) * log(i / (1.5_DP * 10.0_DP ** 5.0_DP))) ** 3.0_DP) * ((4.0_DP / 3.0_DP) &
                  * PI * avg_d)
               ! If the test mass is smaller than the mass that needs to be "made up for", add it. 
               IF (m_test < m_rem) THEN
                  mergeadd_list%mass(nmergeadd) = m_test
               ! If not, aka if the test mass is too large, then add a fragment with a mass equal to the difference between
               ! the sum of the mass of the parents and the total mass already added.   
               ELSE
                  mergeadd_list%mass(nmergeadd) = (m1 + m2) - mtot 
               END IF 

               ! Calculate the radius of the fragments using the weighted average density of the parents. 
               mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * avg_d))  & 
                  ** (1.0_DP / 3.0_DP) 
               mtot = mtot + mergeadd_list%mass(nmergeadd) 

            ! If these fragments will NOT be large enough to be resolved, add the remaining mass that we need to 
            ! "make up for" to the mass of the most recent fragment and recalculate the radius. 
            ELSE 
               mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
               mergeadd_list%radius(nmergeadd) = (((3.0_DP/4.0_DP) * PI) * (mergeadd_list%mass(nmergeadd) / avg_d)) &
                  ** (1.0_DP / 3.0_DP)                                                            
            END IF
         END IF

         ! Increment around the circle for positions of fragments
         x_frag = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1)+ x_com
         y_frag = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
         z_frag = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

         !Conservation of Angular Momentum for velocities of fragments
         A = (((x1(2) * v1(3) * m1) - (x1(3) * v1(2) * m1)) + ((x2(2) * v2(3) * m2) - (x2(3) * v2(2) * m2))) &
            / mergeadd_list%mass(nmergeadd)
         B = (((x1(3) * v1(1) * m1) - (x1(1) * v1(3) * m1)) + ((x2(3) * v2(1) * m2) - (x2(1) * v2(3) * m2))) &
            / mergeadd_list%mass(nmergeadd)
         vx_frag = ((1.0_DP / frags_added) * (B / z_frag)) - vbs(1)
         vy_frag = ((1.0_DP / frags_added) * (-A / z_frag)) - vbs(2)
         vz_frag = vz_com - vbs(3)

         mergeadd_list%xh(1,nmergeadd) = x_frag
         mergeadd_list%xh(2,nmergeadd) = y_frag 
         mergeadd_list%xh(3,nmergeadd) = z_frag                                                    
         mergeadd_list%vh(1,nmergeadd) = vx_frag
         mergeadd_list%vh(2,nmergeadd) = vy_frag
         mergeadd_list%vh(3,nmergeadd) = vz_frag
         ! Tracking linear momentum. 
         mv = mv + (mergeadd_list%mass(nmergeadd) * mergeadd_list%vh(:,nmergeadd))
      END DO
   END IF
   WRITE(*, *) "Number of fragments added: ", frags_added
   ! Calculate energy after frag                                                                           
   vnew(:) = mv / mtot    ! COM of new fragments                               
   enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
   eoffset = eoffset + eold - enew

   ! Update fragmax to account for new fragments
   fragmax = fragmax + frags_added
   WRITE(*,*) "LEAVING CASEDISRUPTION"
   RETURN 
END SUBROUTINE symba_casedisruption

