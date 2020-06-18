submodule (util) s_symba_casehitandrun 
contains
   module procedure symba_casehitandrun 
   !! author: Jennifer L. L. Pouplin and Carlisle A. Wishard
   !!
   !! Hit-and-run collision
   use swiftest
implicit none
   integer(I4B)                       :: nfrag, i, k, index1, index2, frags_added
   integer(I4B)                       :: index1_parent, index2_parent, index_keep_parent, index_rm_parent
   integer(I4B)                       :: name1, name2, index_keep, index_rm, name_keep, name_rm, nstart
   real(DP)                         :: mtot, msun, d_rm, m_rm, r_rm, x_rm, y_rm, z_rm, vx_rm, vy_rm, vz_rm 
   real(DP)                         :: rhill_keep, r_circle, theta, radius1, radius2, e, q, semimajor_encounter
   real(DP)                         :: m_rem, m_test, mass1, mass2, enew, eold, semimajor_inward, a, b, v_col
   real(DP)                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com, mass_keep, mass_rm, rhill_rm
   real(DP)                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag, rad_keep, rad_rm
   real(DP)                         :: r_smallestcircle
   real(DP), dimension(ndim)                :: vnew, xr, mv, xh_keep, xh_rm, vh_keep, vh_rm, l, kk, p

   !temporary
   interface 
    function cross_product_hitandrun(ar1,ar2) result(ans)
       use swiftest
       implicit none
       real(DP),dimension(3),intent(in) :: ar1,ar2
       real(DP),dimension(3)         :: ans
    end function cross_product_hitandrun
   end interface

! executable code

   write(*,*) "entering casehitandrun"

   ! set the maximum number of fragments to be added in a hit and run collision (nfrag)
   nfrag = 4
   ! pull in the information about the two particles involved in the collision 
   index1 = plplenc_list%index1(index_enc)
   index2 = plplenc_list%index2(index_enc)
   index1_parent = symba_plA%index_parent(index1)
   index2_parent = symba_plA%index_parent(index2)
   name1 = symba_plA%name(index1)
   name2 = symba_plA%name(index2)
   mass1 = symba_plA%mass(index1) ! the mass of the first particle in the collision not including all it's children
   mass2 = symba_plA%mass(index2)
   radius1 = symba_plA%radius(index1)
   radius2 = symba_plA%radius(index2)
   msun = symba_plA%mass(1)

   ! determine which of the two particles in the collision is larger where mass includes the mass of all their children
   if (m2 > m1) then
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
    index_keep_parent = index2_parent
    index_rm_parent = index1_parent
    name_keep = name2
    name_rm = name1
   else
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
    index_keep_parent = index1_parent
    index_rm_parent = index2_parent
    name_keep = name1
    name_rm = name2
   end if

   ! find com
   x_com = ((x1(1) * m1) + (x2(1) * m2)) / (m1 + m2)
   y_com = ((x1(2) * m1) + (x2(2) * m2)) / (m1 + m2)
   z_com = ((x1(3) * m1) + (x2(3) * m2)) / (m1 + m2)

   vx_com = ((v1(1) * m1) + (v2(1) * m2)) / (m1 + m2)
   vy_com = ((v1(2) * m1) + (v2(2) * m2)) / (m1 + m2)
   vz_com = ((v1(3) * m1) + (v2(3) * m2)) / (m1 + m2)

   ! find collision velocity
   v_col = norm2(v2(:) - v1(:))

   ! find energy pre-frag
   eold = 0.5_DP*(m1*dot_product(v1(:), v1(:)) + m2*dot_product(v2(:), v2(:)))
   xr(:) = x2(:) - x1(:)
   eold = eold - (m1*m2/(sqrt(dot_product(xr(:), xr(:)))))

   write(*, *) "hit and run between particles ", name1, " and ", name2, " at time t = ",t
   write(*, *) "particle ", name_keep, " survives; particle ", name_rm, " is fragmented."

   ! add both particles involved in the collision to mergesub_list
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name1
   mergesub_list%status(nmergesub) = hit_and_run 
   mergesub_list%xh(:,nmergesub) = x1(:)
   mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass1
   mergesub_list%radius(nmergesub) = rad1
   nmergesub = nmergesub + 1
   mergesub_list%name(nmergesub) = name2
   mergesub_list%status(nmergesub) = hit_and_run
   mergesub_list%xh(:,nmergesub) = x2(:)
   mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
   mergesub_list%mass(nmergesub) = mass2
   mergesub_list%radius(nmergesub) = rad2

   ! go through the encounter list and look for particles actively encoutering in this timestep
   ! prevent them from having further encounters in this timestep by setting status in plplenc_list to merged
   do k = 1, nplplenc 
    if ((plplenc_list%status(k) == active) .and. &
       ((index1 == plplenc_list%index1(k) .or. index2 == plplenc_list%index2(k)) .or. &
       (index2 == plplenc_list%index1(k) .or. index1 == plplenc_list%index2(k)))) then
        plplenc_list%status(k) = merged
    end if
   end do

   ! set the status of the particles in symba_plA to hit_and_run
   symba_plA%status(index1) = hit_and_run
   symba_plA%status(index2) = hit_and_run

   l(:) = (v2(:) - v1(:)) / norm2(v2(:)-v1(:))
   p(:) = cross_product_hitandrun(xr(:) / norm2(xr(:)), l(:))
   kk(:) = cross_product_hitandrun(l(:),p(:))

   mtot = 0.0_DP ! running total mass of new fragments
   mv = 0.0_DP   ! running sum of m*v of new fragments to be used in com calculation
   frags_added = 0 ! running total number of new fragments
   nstart = nmergeadd + 1 ! start of new fragments in mergeadd_list
   ! increment around the circle for positions of fragments
   ! calculate the positions of the new fragments in a circle of radius rhill_keep
   rhill_keep = symba_plA%rhill(index_keep_parent)
   rhill_rm = symba_plA%rhill(index_rm_parent)
   r_smallestcircle = (rhscale * rhill_rm + rhscale * rhill_keep) / (2.0_DP * sin(pi / 2.0_DP))

   ! check that no fragments will be added interior of the smallest orbit that the timestep can reliably resolve
   semimajor_inward = ((dt * 32.0_DP) ** 2.0_DP) ** (1.0_DP / 3.0_DP)
   call orbel_xv2aeq(x1, v1, msun, semimajor_encounter, e, q)
   ! if they are going to be added interior to this orbit, give a warning
   if (semimajor_inward > (semimajor_encounter - r_smallestcircle)) then
    write(*,*) "warning in symba_casehitandrun: timestep is too large to resolve fragments."
   end if

   ! the largest fragment = the kept parent particle
   nmergeadd = nmergeadd + 1
   mergeadd_list%status(nmergeadd) = hit_and_run
   mergeadd_list%ncomp(nmergeadd) = 2
   mergeadd_list%name(nmergeadd) = symba_plA%name(index_keep)
   mergeadd_list%mass(nmergeadd) = mass_keep
   mergeadd_list%radius(nmergeadd) = rad_keep
   mergeadd_list%xh(:,nmergeadd) = xh_keep
   mergeadd_list%vh(:,nmergeadd) = vh_keep
   mtot = mtot + mergeadd_list%mass(nmergeadd) 


   ! pure hit & run
   if (mres(2) > m2 * 0.9_DP) then
    frags_added = frags_added + 1
    nmergeadd = nmergeadd + 1
    mergeadd_list%status(nmergeadd) = hit_and_run
    mergeadd_list%ncomp(nmergeadd) = 2
    mergeadd_list%name(nmergeadd) = config%nplmax + config%ntpmax + fragmax + i - 1
    mergeadd_list%mass(nmergeadd) = mass_rm
    mergeadd_list%radius(nmergeadd) = rad_rm
    mergeadd_list%xh(:,nmergeadd) = xh_rm(:)
    mergeadd_list%vh(:,nmergeadd) = vh_rm(:)
    mtot = mtot + mergeadd_list%mass(nmergeadd)
   else     
    do i = 1, nfrag
       m_rm = mass_rm
       r_rm = rad_rm
       !x_rm = xh_rm(1)
       !y_rm = xh_rm(2)
       !z_rm = xh_rm(3)
       !vx_rm = vh_rm(1)
       !vy_rm = vh_rm(2)
       !vz_rm = vh_rm(3)
       d_rm = (3.0_DP * m_rm) / (4.0_DP * pi * (r_rm ** 3.0_DP))

       m_rem = m_rm - mres(2)
       frags_added = frags_added + 1
       nmergeadd = nmergeadd + 1
       mergeadd_list%status(nmergeadd) = hit_and_run
       mergeadd_list%ncomp(nmergeadd) = 2
       mergeadd_list%name(nmergeadd) = config%nplmax + config%ntpmax + fragmax + i - 1
       mergeadd_list%mass(nmergeadd) = m_rem / (nfrag - 1) 
       mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * pi * d_rm))  & 
        ** (1.0_DP / 3.0_DP) 

        ! check if these fragments will not be large enough to be resolved and we have only added one fragment 
        ! previously (aka the slr). this is the perfect hit and run case.   
       !else if ((i > 2) .and. (mres(2) > m2 * 0.9_DP) .and. frags_added == 1) then
        ! if yes, update the mass of the slr to be the mass of the removed particle and give it all the
        ! characteristics of the removed particle
       
       !   mergeadd_list%name(nmergeadd) = symba_plA%name(index_rm)
       !   mergeadd_list%mass(nmergeadd) = mass_rm
       !   mergeadd_list%radius(nmergeadd) = rad_rm
       !   mergeadd_list%xh(:,nmergeadd) = xh_rm
       !   mergeadd_list%vh(:,nmergeadd) = vh_rm
       !   mtot = mtot - mres(2) + mass_rm

        ! if these fragments will not be large enough to be resolved but we have added more than one fragment
        ! previously, add the remaining mass that we need to "make up for" to the mass of the most recent
        ! fragment and recalculate the radius. 
        !else 
        !   mergeadd_list%mass(nmergeadd) = mergeadd_list%mass(nmergeadd) + m_rem
        !   mergeadd_list%radius(nmergeadd) = (((3.0_DP/4.0_DP) * pi) * (mergeadd_list%mass(nmergeadd) / d_rm)) &
        !    ** (1.0_DP / 3.0_DP)                                    
        !end if  
    end do
   end if

   if (frags_added > 1) then
       r_circle = (rhscale * rhill_keep + rhscale * rhill_rm) / (2.0_DP * sin(pi / frags_added))
       theta = (2.0_DP * pi) / (frags_added)
       do i=1, frags_added
        ! increment around the circle for positions of fragments
        x_frag = (r_circle * cos(theta * i))*l(1) + (r_circle * sin(theta * i))*p(1) + x_com
        y_frag = (r_circle * cos(theta * i))*l(2) + (r_circle * sin(theta * i))*p(2) + y_com
        z_frag = (r_circle * cos(theta * i))*l(3) + (r_circle * sin(theta * i))*p(3) + z_com

        !vx_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m2 * v2(1)))) !- vbs(1)
        !vy_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m2 * v2(2)))) !- vbs(2)
        !vz_frag = ((1.0_DP / frags_added) * (1.0_DP / mergeadd_list%mass(nstart + i)) * ((m2 * v2(3)))) !- vbs(3)

        a = v_col * m2 * (1.0_DP / mergeadd_list%mass(nstart + i))

        vx_frag = ((a * cos(theta * i))*l(1)) + ((a * sin(theta * i))*p(1)) + vh_rm(1) !+ vx_com
        vy_frag = ((a * cos(theta * i))*l(2)) + ((a * sin(theta * i))*p(2)) + vh_rm(2) !+ vy_com
        vz_frag = ((a * cos(theta * i))*l(3)) + ((a * sin(theta * i))*p(3)) + vh_rm(3) !+ vz_com

        ! conservation of angular momentum for velocities of fragments
        !a = ((y_rm * vz_rm * m_rm) - (z_rm * vy_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
        !b = ((z_rm * vx_rm * m_rm) - (x_rm * vz_rm * m_rm)) / mergeadd_list%mass(nmergeadd)
        !vx_frag = ((1.0_DP / frags_added) * (b / z_frag)) - vbs(1)
        !vy_frag = ((1.0_DP / frags_added) * (-a / z_frag)) - vbs(2)
        !vz_frag = vz_com - vbs(3)

        mergeadd_list%xh(1,nstart + i) = x_frag
        mergeadd_list%xh(2,nstart + i) = y_frag 
        mergeadd_list%xh(3,nstart + i) = z_frag                                
        mergeadd_list%vh(1,nstart + i) = vx_frag
        mergeadd_list%vh(2,nstart + i) = vy_frag
        mergeadd_list%vh(3,nstart + i) = vz_frag 

        ! tracking linear momentum.                            
        mv = mv + (mergeadd_list%mass(nmergeadd) * mergeadd_list%vh(:,nmergeadd))
       end do 
   end if
   write(*, *) "number of fragments added: ", (frags_added)
   ! calculate energy after frag                                             
   vnew(:) = mv / mtot    ! com of new fragments                   
   enew = 0.5_DP*mtot*dot_product(vnew(:), vnew(:))
   eoffset = eoffset + eold - enew
   ! update fragmax to account for new fragments
   fragmax = fragmax + frags_added
   write(*,*) "leaving casehitandrun"
   return 
 end procedure symba_casehitandrun
end submodule s_symba_casehitandrun
