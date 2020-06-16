!**********************************************************************************************************************************
!
!  Unit Name   : symba_casehitandrun
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge massive bodies
!
!  Input
!    Arguments : t            : time
!                npl          : number of massive bodies
!                nsppl        : number of spilled massive bodies
!                symba_pl1P   : pointer to head of SyMBA massive body structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA massive body structure linked-list
!                nplplenc     : number of massive body-massive body encounters
!                plplenc_list : array of massive body-massive body encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of massive bodies
!                nsppl        : number of spilled massive bodies
!                symba_pl1P   : pointer to head of SyMBA massive body structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA massive body structure linked-list
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
     use swiftest, EXCEPT_THIS_ONE => symba_casehitandrun
     IMPLICIT NONE

! Arguments
   INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
   INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, nplplenc, fragmax
   REAL(DP), INTENT(IN)                             :: t, dt
   REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: mres, rres
   REAL(DP), DIMENSION(:), INTENT(IN)               :: vbs
   REAL(DP), DIMENSION(:), INTENT(INOUT)            :: x1, x2, v1, v2
   TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
   TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
   TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
   INTEGER(I4B)                                     :: nfrag, i, k, index1, index2, frags_added
   INTEGER(I4B)                                     :: index1_parent, index2_parent, index_keep_parent, index_rm_parent
   INTEGER(I4B)                                     :: name1, name2, index_keep, index_rm, name_keep, name_rm, nstart
   REAL(DP)                                         :: mtot, msun, d_rm, m_rm, r_rm, x_rm, y_rm, z_rm, vx_rm, vy_rm, vz_rm 
   REAL(DP)                                         :: rhill_keep, r_circle, theta, radius1, radius2, e, q, semimajor_encounter
   REAL(DP)                                         :: m_rem, m_test, mass1, mass2, enew, eold, semimajor_inward, A, B, v_col
   REAL(DP)                                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com, mass_keep, mass_rm
   REAL(DP)                                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag, rad_keep, rad_rm
   REAL(DP), DIMENSION(NDIM)                        :: vnew, xr, mv, xh_keep, xh_rm, vh_keep, vh_rm, l, kk, p

   !TEMPORARY
   interface 
      function cross_product_hitandrun(ar1,ar2) result(ans)
         use swiftest
         implicit none
         real(DP),dimension(3),intent(in) :: ar1,ar2
         real(DP),dimension(3)             :: ans
      end function cross_product_hitandrun
   end interface

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

   ! Find Collision velocity
   v_col = NORM2(v2(:) - v1(:))

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

   l(:) = (v2(:) - v1(:)) / NORM2(v2(:)-v1(:))
   p(:) = cross_product_hitandrun(xr(:) / NORM2(xr(:)), l(:))
   kk(:) = cross_product_hitandrun(l(:),p(:))

   mtot = 0.0_DP ! running total mass of new fragments
   mv = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation
   frags_added = 0 ! running total number of new fragments
   nstart = nmergeadd + 1 ! start of new fragments in mergeadd_list
   ! Increment around the circle for positions of fragments
   ! Calculate the positions of the new fragments in a circle of radius rhill_keep
   rhill_keep = symba_plA%helio%swiftest%rhill(index_keep_parent)
   r_circle = rhill_keep

   ! Check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   CALL orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
   ! If they are going to be added interior to this orbit, give a warning
   IF (semimajor_inward > (semimajor_encounter - r_circle)) THEN
      WRITE(*,*) "WARNING in symba_casehitandrun: Timestep is too large to resolve fragments."
   END IF

   ! The largest fragment = the kept parent particle
   nmergeadd = nmergeadd + 1
   mergeadd_list%status(nmergeadd) = HIT_AND_RUN
   mergeadd_list%ncomp(nmergeadd) = 2
   mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_keep)
   mergeadd_list%mass(nmergeadd) = mass_keep
   mergeadd_list%radius(nmergeadd) = rad_keep
   mergeadd_list%xh(:,nmergeadd) = xh_keep
   mergeadd_list%vh(:,nmergeadd) = vh_keep
   mtot = mtot + mergeadd_list%mass(nmergeadd) 


   ! Pure Hit & Run
   IF (mres(2) > m2 * 0.9_DP) THEN
      frags_added = frags_added + 1
      nmergeadd = nmergeadd + 1
      mergeadd_list%status(nmergeadd) = HIT_AND_RUN
      mergeadd_list%ncomp(nmergeadd) = 2
      mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i - 1
      mergeadd_list%mass(nmergeadd) = mass_rm
      mergeadd_list%radius(nmergeadd) = rad_rm
      mergeadd_list%xh(:,nmergeadd) = xh_rm(:)
      mergeadd_list%vh(:,nmergeadd) = vh_rm(:)
      mtot = mtot + mergeadd_list%mass(nmergeadd)
   ELSE       
      DO i = 1, nfrag
         m_rm = mass_rm
         r_rm = rad_rm
         !x_rm = xh_rm(1)
         !y_rm = xh_rm(2)
         !z_rm = xh_rm(3)
         !vx_rm = vh_rm(1)
         !vy_rm = vh_rm(2)
         !vz_rm = vh_rm(3)
         d_rm = (3.0_DP * m_rm) / (4.0_DP * PI * (r_rm ** 3.0_DP))

         m_rem = m_rm - mres(2)
         frags_added = frags_added + 1
         nmergeadd = nmergeadd + 1
         mergeadd_list%status(nmergeadd) = HIT_AND_RUN
         mergeadd_list%ncomp(nmergeadd) = 2
         mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i - 1
         mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 
         mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * d_rm))  & 
            ** (1.0_DP / 3.0_DP) 

            ! Check if these fragments will NOT be large enough to be resolved AND we have only added one fragment 
            ! previously (aka the slr). This is the perfect hit and run case.   
         !ELSE IF ((i > 2) .AND. (mres(2) > m2 * 0.9_DP) .AND. frags_added == 1) THEN
            ! If yes, update the mass of the slr to be the mass of the removed particle and give it all the
            ! characteristics of the removed particle
         
         !   mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_rm)
         !   mergeadd_list%mass(nmergeadd) = mass_rm
         !   mergeadd_list%radius(nmergeadd) = rad_rm
         !   mergeadd_list%xh(:,nmergeadd) = xh_rm
         !   mergeadd_list%vh(:,nmergeadd) = vh_rm
         !   mtot = mtot - mres(2) + mass_rm

            ! If these fragments will NOT be large enough to be resolved but we have added more than one fragment
            ! previously, add the remaining mass that we need to "make up for" to the mass of the most recent
            ! fragment and recalculate the radius. 
            !ELSE 
            !   mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
            !   mergeadd_list%radius(nmergeadd) = (((3.0_DP/4.0_DP) * PI) * (mergeadd_list%mass(nmergeadd) / d_rm)) &
            !      ** (1.0_DP / 3.0_DP)                                                            
            !END IF  
      END DO
   END IF

   IF (frags_added > 1) THEN
         theta = (2.0_DP * PI) / (frags_added)
         DO i=1, frags_added
            ! Increment around the circle for positions of fragments
            x_frag = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1) + x_com
            y_frag = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
            z_frag = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

            !vx_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m2 * v2(1)))) !- vbs(1)
            !vy_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m2 * v2(2)))) !- vbs(2)
            !vz_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m2 * v2(3)))) !- vbs(3)

            A = v_col * m2 * (1.0_DP / mergeadd_list%mass(nstart + i))

            vx_frag = ((A * cos(theta * i))*l(1)) + ((A * sin(theta * i))*p(1)) + vh_rm(1) !+ vx_com
            vy_frag = ((A * cos(theta * i))*l(2)) + ((A * sin(theta * i))*p(2)) + vh_rm(2) !+ vy_com
            vz_frag = ((A * cos(theta * i))*l(3)) + ((A * sin(theta * i))*p(3)) + vh_rm(3) !+ vz_com

            ! Conservation of Angular Momentum for velocities of fragments
            !A = ((y_rm * vz_rm * m_rm) - (z_rm * vy_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
            !B = ((z_rm * vx_rm * m_rm) - (x_rm * vz_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
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


function cross_product_hitandrun(ar1,ar2) result(ans)
   use swiftest
   implicit none
   
   real(DP),dimension(3),intent(in) :: ar1,ar2
   real(DP),dimension(3)             :: ans

   ans(1) = ar1(2) * ar2(3) - ar1(3) * ar2(2)
   ans(2) = ar1(3) * ar2(1) - ar1(1) * ar2(3)
   ans(3) = ar1(1) * ar2(2) - ar1(2) * ar2(1)

   return 
end function cross_product_hitandrun
