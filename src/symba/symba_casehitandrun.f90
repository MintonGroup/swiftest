!**********************************************************************************************************************************
!
!  Unit Name   : symba_casehitandrun
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
!  Invocation  : CALL symba_casehitandrun(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
   symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

! Modules
   USE swiftest
   USE module_helio
   USE module_symba
   USE module_interfaces, EXCEPT_THIS_ONE => symba_casehitandrun
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
   INTEGER(I4B)                                     :: index1_parent, index2_parent, index_keep_parent, index_rm_parent
   INTEGER(I4B)                                     :: name1, name2, index_keep, index_rm, name_keep, name_rm
   REAL(DP)                                         :: mtot, msun, d_rm, m_rm, r_rm, x_rm, y_rm, z_rm, vx_rm, vy_rm, vz_rm 
   REAL(DP)                                         :: rhill_keep, r_circle, theta, radius1, radius2, e, q, semimajor_encounter
   REAL(DP)                                         :: m_rem, m_test, mass1, mass2, enew, eold, semimajor_inward, A, B
   REAL(DP)                                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com, mass_keep, mass_rm
   REAL(DP)                                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag, rad_keep, rad_rm
   REAL(DP), DIMENSION(NDIM)                        :: vnew, xr, mv, xh_keep, xh_rm, vh_keep, vh_rm


! Executable code

   WRITE(*,*) "ENTERING CASEHITANDRUN"

   ! Set the maximum number of fragments to be added in a Hit and Run collision (nfrag)
   nfrag = 4
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

   ! Determine which of the two particles in the collision is larger where mass INCLUDES the mass of all their children
   IF (m2 > m1) THEN
      index_keep = index2
      index_rm = index1
      mass_keep = m2
      mass_rm = m1
      rad_keep = rad2
      rad_rm = rad1
      xh_keep = x2
      xh_rm = x1
      vh_keep = v2
      vh_rm = v1
      index_keep_parent = symba_plA%index_parent(index_keep)
      index_rm_parent = symba_plA%index_parent(index_rm)
      name_keep = symba_plA%helio%swiftest%name(index_keep)
      name_rm = symba_plA%helio%swiftest%name(index_rm)
   ELSE
      index_keep = index1
      index_rm = index2
      mass_keep = m1
      mass_rm = m2
      rad_keep = rad1
      rad_rm = rad2
      xh_keep = x1
      xh_rm = x2
      vh_keep = v1
      vh_rm = v2
      index_keep_parent = symba_plA%index_parent(index_keep)
      index_rm_parent = symba_plA%index_parent(index_rm)
      name_keep = symba_plA%helio%swiftest%name(index_keep)
      name_rm = symba_plA%helio%swiftest%name(index_rm)
   END IF

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

   WRITE(*, *) "Hit and run between particles ", name1, " and ", name2, " at time t = ",t
   WRITE(*, *) "Particle ", name_keep, " survives; Particle ", name_rm, " is fragmented."

   ! Add both particles involved in the collision to mergesub_list
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = HIT_AND_RUN 
   mergesub_list%xh(:,nmergesub) = x1(:)
   mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = rad1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = HIT_AND_RUN
   mergesub_list%xh(:,nmergesub) = x2(:)
   mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = rad2

   ! Go through the encounter list and look for particles actively encoutering in this timestep
   ! Prevent them from having further encounters in this timestep by setting status in plplenc_list to MERGED
   DO k = 1, nplplenc 
      IF ((plplenc_list%status(k) == ACTIVE) .AND. &
         ((index1 == plplenc_list%index1(k) .OR. index2 == plplenc_list%index2(k)) .OR. &
         (index2 == plplenc_list%index1(k) .OR. index1 == plplenc_list%index2(k)))) THEN
            plplenc_list%status(k) = MERGED
      END IF
   END DO

   ! Set the status of the particles in symba_plA to HIT_AND_RUN
   symba_plA%helio%swiftest%status(index1) = HIT_AND_RUN
   symba_plA%helio%swiftest%status(index2) = HIT_AND_RUN

   ! Calculate the positions of the new fragments in a circle of radius rhill_keep
   rhill_keep = symba_plA%helio%swiftest%rhill(index_keep_parent)
   r_circle = rhill_keep
   theta = (2.0_DP * PI) / (nfrag - 1)

   mtot = 0.0_DP ! running total mass of new fragments
   mv = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
   frags_added = 0 ! running total number of new fragments

   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
   ! If they are going to be added interior to this orbit, give a warning
   IF (semimajor_inward > (semimajor_encounter - r_circle)) THEN
      WRITE(*,*) "WARNING in symba_casehitandrun: Timestep is too large to resolve fragments."
   ELSE
   ! If not, continue through all possible fragments to be added
      DO i = 1, nfrag
         m_rm = mass_rm
         r_rm = rad_rm
         x_rm = xh_rm(1)
         y_rm = xh_rm(2)
         z_rm = xh_rm(3)
         vx_rm = vh_rm(1)
         vy_rm = vh_rm(2)
         vz_rm = vh_rm(3)
         d_rm = (3.0_DP * m_rm) / (4.0_DP * PI * (r_rm ** 3.0_DP))
         ! If we are adding the first and largest fragment (lr), it's characteristics should be equal to those of the kept parent
         IF (i == 1) THEN
            nmergeadd = nmergeadd + 1
            mergeadd_list%status(nmergeadd) = HIT_AND_RUN
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_keep)
            mergeadd_list%mass(nmergeadd) = mass_keep
            mergeadd_list%radius(nmergeadd) = rad_keep
            mergeadd_list%xh(:,nmergeadd) = xh_keep
            mergeadd_list%vh(:,nmergeadd) = vh_keep
            mtot = mtot + mergeadd_list%mass(nmergeadd)                       
         END IF
         ! If we are adding the second fragment, it's mass and radius should be taken from util_regime while it's position
         ! and velocity should be calculated on the circle of radius rhill_keep as described above.
         IF (i == 2) THEN
            ! frags_added is the actual number of fragments added to the simulation vs nfrag which is the total possible
            frags_added = frags_added + 1
            nmergeadd = nmergeadd + 1
            mergeadd_list%status(nmergeadd) = HIT_AND_RUN
            mergeadd_list%ncomp(nmergeadd) = 2
            mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i - 1
            mergeadd_list%mass(nmergeadd) = mres(2)
            mergeadd_list%radius(nmergeadd) = rres(2) 
            mtot = mtot + mergeadd_list%mass(nmergeadd) 

            ! Increment around the circle for positions of fragments
            x_frag = (r_circle * cos(theta * i)) + x_com
            y_frag = (r_circle * sin(theta * i)) + y_com
            z_frag = z_com

            ! Conservation of Angular Momentum for velocities of fragments
            A = ((y_rm * vz_rm * m_rm) - (z_rm * vy_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
            B = ((z_rm * vx_rm * m_rm) - (x_rm * vz_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
            vx_frag = ((1.0_DP / frags_added) * (B / z_frag)) - vbs(1)
            vy_frag = ((1.0_DP / frags_added) * (-A / z_frag)) - vbs(2)
            vz_frag = vz_com - vbs(3)

            mergeadd_list%xh(1,nmergeadd) = x_frag
            mergeadd_list%xh(2,nmergeadd) = y_frag 
            mergeadd_list%xh(3,nmergeadd) = z_frag                                                    
            mergeadd_list%vh(1,nmergeadd) = vx_frag
            mergeadd_list%vh(2,nmergeadd) = vy_frag
            mergeadd_list%vh(3,nmergeadd) = vz_frag        
         END IF
         ! If we are doing more fragments
         IF (i > 2) THEN
            ! m_rem is the mass needed to be "made up for" in fragments, m_rm is the mass of the removed particle, 
            ! and mres(2) is the mass of the second largest fragment that has already been added
            m_rem = m_rm - mres(2)
            ! Check if these fragments will be large enough to be resolved 
            IF (m_rem > (m_rm) / 1000.0_DP) THEN
               ! If yes, add a fragment using Durda et al 2007 Figure 2 Supercatastrophic: N = (1.5e5)e(-1.3*D) for the mass
               frags_added = frags_added + 1
               nmergeadd = nmergeadd + 1
               mergeadd_list%status(nmergeadd) = HIT_AND_RUN
               mergeadd_list%ncomp(nmergeadd) = 2
               mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i - 1
               ! Create a "test mass" using Durda et al 2007 Figure 2 Supercatastrophic: N = (1.5e5)e(-1.3*D)
               m_test = (((- 1.0_DP / 2.6_DP) * log(i / (1.5_DP * 10.0_DP ** 5.0_DP))) ** 3.0_DP) * ((4.0_DP / 3.0_DP) &
                  * PI * d_rm)

               ! If the test mass is smaller than the mass that needs to be "made up for", add it. 
               IF (m_test < m_rem) THEN
                  mergeadd_list%mass(nmergeadd) = m_test
               ! If not, aka if the test mass is too large, then add a fragment with a mass equal to the difference between
               ! the sum of the mass of the parents and the total mass already added.   
               ELSE
                  mergeadd_list%mass(nmergeadd) = (m1 + m2) - mtot
               END IF 

               ! Calculate the radius of the fragments using the density of the removed particle. 
               mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * d_rm))  & 
                  ** (1.0_DP / 3.0_DP) 

               ! Increment around the circle for positions of fragments
               mtot = mtot + mergeadd_list%mass(nmergeadd)
               x_frag = (r_circle * cos(theta * i)) + x_com
               y_frag = (r_circle * sin(theta * i)) + y_com
               z_frag = z_com

               ! Conservation of Angular Momentum for velocities of fragments
               A = ((y_rm * vz_rm * m_rm) - (z_rm * vy_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
               B = ((z_rm * vx_rm * m_rm) - (x_rm * vz_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
               vx_frag = ((1.0_DP / frags_added) * (B / z_frag)) - vbs(1)
               vy_frag = ((1.0_DP / frags_added) * (-A / z_frag)) - vbs(2)
               vz_frag = vz_com - vbs(3)

               mergeadd_list%xh(1,nmergeadd) = x_frag
               mergeadd_list%xh(2,nmergeadd) = y_frag 
               mergeadd_list%xh(3,nmergeadd) = z_frag                                                    
               mergeadd_list%vh(1,nmergeadd) = vx_frag
               mergeadd_list%vh(2,nmergeadd) = vy_frag
               mergeadd_list%vh(3,nmergeadd) = vz_frag 
            ! Check if these fragments will NOT be large enough to be resolved AND we have only added one fragment 
            ! previously (aka the slr). This is the perfect hit and run case.   
            ELSE IF ((m_rem < (m_rm) / 1000.0_DP) .AND. frags_added == 1) THEN
               ! If yes, update the mass of the slr to be the mass of the removed particle and give it all the
               ! characteristics of the removed particle
               mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_rm)
               mergeadd_list%mass(nmergeadd) = mass_rm
               mergeadd_list%radius(nmergeadd) = rad_rm
               mergeadd_list%xh(:,nmergeadd) = xh_rm
               mergeadd_list%vh(:,nmergeadd) = vh_rm
               mtot = mtot - mres(2) + mass_rm
            ! If these fragments will NOT be large enough to be resolved but we have added more than one fragment
            ! previously, add the remaining mass that we need to "make up for" to the mass of the most recent
            ! fragment and recalculate the radius. 
            ELSE 
               mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
               mergeadd_list%radius(nmergeadd) = (((3.0_DP/4.0_DP) * PI) * (mergeadd_list%mass(nmergeadd) / d_rm)) &
                  ** (1.0_DP / 3.0_DP)                                                            
            END IF  
         END IF
         ! Tracking linear momentum.                                            
         mv = mv + (mergeadd_list%mass(nmergeadd) * mergeadd_list%vh(:,nmergeadd))
      END DO
   END IF
   WRITE(*, *) "Number of fragments added: ", (frags_added)
   ! Calculate energy after frag                                                                           
   vnew(:) = mv / mtot    ! COM of new fragments                               
   enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
   eoffset = eoffset + eold - enew
   ! Update fragmax to account for new fragments
   fragmax = fragmax + frags_added
   WRITE(*,*) "LEAVING CASEHITANDRUN"
   RETURN 
END SUBROUTINE symba_casehitandrun
