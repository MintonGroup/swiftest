! Copyright 2026 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(ringmoons) s_ringmoons_torque
    use swiftest
contains
    
    module function ringmoons_torque_lindblad_ring(self,cb,seed,param) result(Torque)
        !! author: David A. Minton
        !!
        !! Calculates the lindblad torques between each ring element and given  given satellite. 
        !! 
        !! The function returns total torque acting on the satellite, and stores the torques acting on each ring element ! in the 
        !! ring. 
        !! 
        !! References:
        !! Tajeddine, R., Nicholson, P.D., Longaretti, P.-Y., Moutamid, M.E., Burns, J.A., 2017. What Confines the Rings of Saturn? 
        !!     ApJS 232, 28. https://doi.org/10.3847/1538-4365/aa8c09

        !!
        !! Adapted from Andy Hesselbrock's ringmoons Python scripts.
        implicit none
        ! Arguments
        class(ringmoons_ring),      intent(inout) :: self
        class(swiftest_cb),         intent(in)    :: cb 
        class(ringmoons_seed),      intent(inout) :: seed
        class(swiftest_parameters), intent(in)    :: param
        real(DP),dimension(0:self%nbins+1)        :: Torque
        ! Internals
        real(DP)                               :: asat,esat,isat,msat
        integer(I4B)                           :: i,j,m,inner_outer_sign,il,w,w1,w2,js, mshep
        real(DP)                               :: beta, Amk, nw,lap,dlap,da3,Xs,Xlo,Xhi,rlo,rhi,mratio,lind_factor,awidth,Xw2
        real(DP), parameter                    :: g = 2.24_DP !! see Tajeddine et al. (2017) eq. 7
        integer(I4B),parameter :: M_MAX = 100     
            !! Maximum number of Lindblad modes to compute
        real(DP),dimension(2:M_MAX)            :: Xring, aring
        logical,dimension(0:self%nbins+1)      :: T_mask
        integer(I4B),dimension(2:M_MAX)        :: w1_arr,w2_arr
        real(DP),dimension(0:self%nbins+1)     :: iTorque
        real(DP),dimension(2:M_MAX,-1:1),save :: lapm,dlapm,marr 
            !! Laplace coefficients and mode array
        real(DP),dimension(2:M_MAX),save :: mfac
            !! Mode factor for computation
        real(DP),parameter         :: RAD_LIMIT_M = 0.001_DP 
            !! Lower limit on disk particle radius in meters
        logical, save :: lfirst = .true.

        ! For performance reasons, we compute a table of Laplace coefficient terms in Am,k from Tajeddine et al. (2017) eq. 27
        ! lapm and dlapm are the two Laplace coefficient terms 
        ! marr is an array used to locate the resonance in X space given a satellite location Xs. It is based on the mode and 
        ! Kepler's 3rd law (which is where the 1/3 power comes from)
        ! mfac is the leading term of the equation for the torque. See Tajeddine et al. (2017) eq. 4
        if (lfirst) then
            do m = 2, M_MAX
                do inner_outer_sign = -1,1,2
                    beta =  (1._DP + inner_outer_sign * 1.0_DP / real(m, kind=DP))**(-inner_outer_sign * 2._DP / 3._DP)
                    lapm(m,inner_outer_sign)  = 2 * m * compute_laplace_coefficient(beta,m,0.5_DP,0) 
                    dlapm(m,inner_outer_sign) = beta * compute_laplace_coefficient(beta,m,0.5_DP,1) 
                    marr(m,inner_outer_sign) = (1._DP + inner_outer_sign / real(m, kind=DP))**(1._DP / 3._DP)
                end do
                mfac(m) = 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) 
            end do
            lfirst = .false.
        end if
        Torque(:) = 0.0_DP
        esat = 0.0_DP
        isat = 0.0_DP

        associate(ring => self, nbins => self%nbins)
            do i = 1, seed%nbody
                associate(asat => seed%a(i), msat => seed%mass(i))
                    mratio = abs(msat / cb%mass)
                    ! Mask out any ring bins that don't have enough mass in them
                    where (ring%sigma(0:nbins+1) > 1000 * VSMALL)
                        T_mask(0:nbins+1) = .true.
                    elsewhere
                        T_mask(0:nbins+1) = .false. 
                    end where

                    Xs = 2 * sqrt(asat)
                    ! Resonance width: See Longaretti sec. 5.3.1
                    awidth = asat * sqrt(mratio)
                    Xw2 = 0.5_DP * Xs**2 * sqrt(mratio)
                    Xlo = ring%X_inner + ring%deltaX * ring%inside
                    Xhi = ring%X_outer
                    rlo = 0.25_DP * Xlo**2
                    rhi = 0.25_DP * Xhi**2

                    ! Just do the first order resonances for now. The full suite of resonances will come later
                    iTorque(:) = 0.0_DP
                    ! Calculate the number of modes that are separated by at least 1 bin width in X space
                    mshep = max(2,min(M_MAX - 1,ceiling(0.5_DP * (sqrt(1._DP + 4._DP / 3._DP * Xs / ring%deltaX) - 1._DP))))
                    
                    ! Inner then outer Lindblads
                    do il = -1,1,2
                        ! Calculate resonance location space
                        Xring(2:mshep+1) = Xs * marr(2:mshep+1,il)
                        aring(2:mshep+1) = 0.25_DP * Xring(2:mshep+1)**2

                        ! Calculate bin boundaries of resonance using its width in X space
                        where((Xring(2:mshep) > Xlo).and.(Xring(2:mshep) <= Xhi))
                            w1_arr(2:mshep) = min(max(ceiling((sqrt(Xring(2:mshep)**2-Xw2)-ring%X_inner)/ring%deltaX),0),nbins+1)
                            w2_arr(2:mshep) = min(max(ceiling((sqrt(Xring(2:mshep)**2+Xw2)-ring%X_inner)/ring%deltaX),0),nbins+1)
                        elsewhere
                            w1_arr(2:mshep) = 0
                            w2_arr(2:mshep) = 0
                        end where
                        
                        do m  = 2, mshep
                            if ((Xring(m) > Xlo).and.(Xring(m) < Xhi)) then
                                beta = (Xs / Xring(m))**(il * 2)
                                lap  =  lapm(m,il)
                                dlap = dlapm(m,il)
                                Amk = 0.5_DP * (lap + dlap)
                                w1 = w1_arr(m)
                                w2 = w2_arr(m)
                                nw = real(abs(w2 - w1) + 1,kind=DP)
                                ! Calculate the 1st order Lindblad torques and distribute them over the bins that include the resonance
                                lind_factor = il * mfac(m) / nw * aring(m)**4 * (beta * mratio * Amk)**2 
                                where(T_mask(w1:w2)) 
                                    iTorque(w1:w2) = iTorque(w1:w2) + lind_factor * ring%sigma(w1:w2) * (ring%nkep(w1:w2))**2
                                endwhere
                            end if
                        end do

                        ! Add in shepherding torque
                        if ((aring(mshep+1) > rhi).or.(aring(mshep+1) < rlo)) cycle
                        j = ring%find_bin(aring(mshep+1)) !ring location of resonance
                        da3 = il * max(abs((aring(mshep+1) - asat)**3),epsilon(asat))
                        if (T_mask(j)) then 
                            iTorque(j) = iTorque(j) + g**2 / 6._DP * aring(mshep+1)**3 / da3 * mratio**2  &
                                                * ring%sigma(j) * (ring%nkep(j))**2 * aring(mshep+1)**4
                        end if
                    end do

                    ! Apply the torques to the ring and the seed. The ring gets a torque equal and opposite to the satellite, but 
                    ! distributed over the bins that include the resonance.
                    seed%Torque(i) = seed%Torque(i) - sum(iTorque(:))
                    Torque(:) = Torque(:) + iTorque(:)
                end associate
            end do
        end associate
        return
    end function ringmoons_torque_lindblad_ring

    module subroutine ringmoons_torque_tidal_seed(self,cb,param)
        !! author: David A. Minton
        !!
        !! Calculates the tidal torque acting on the seed by the central body. 
        !! 
        !! Constant Q tidal model. See, for example, Cheng et al. (2014) eq. (16) for zero obliquity and eccentricity.
        !!
        !! References:
        !! Cheng, W.H., Lee, M.H., Peale, S.J., 2014. Complete tidal evolution of Pluto-Charon. Icarus 233, 242–258. 
        !!      https://doi.org/10.1016/j.icarus.2014.01.046

        implicit none
        class(ringmoons_seed),      intent(inout) :: self
        class(swiftest_cb),         intent(in)    :: cb
        class(swiftest_parameters), intent(in)    :: param
        ! Internals
        real(DP),dimension(self%nbody) :: n
        associate(seed => self, Ns => self%nbody)
            n(1:Ns) = sqrt((seed%mu(1:Ns)) / seed%a(1:Ns)**3)
            seed%Ttide(1:Ns) = sign(1._DP,cb%rot(3) * DEG2RAD - n(1:Ns))                  &
                                * 1.5_DP * seed%a(1:Ns) * n(1:Ns) * (cb%k2/cb%Q) &
                                * (seed%mass(1:Ns) / cb%mass)                &
                                * (cb%radius / seed%a(1:Ns))**5              &
                                * (seed%mass(1:Ns) * sqrt(cb%Gmass / seed%a(1:Ns)))
            seed%Torque(1:Ns) = seed%Torque(1:Ns) + seed%Ttide(1:Ns)
        end associate
        return
        end subroutine ringmoons_torque_tidal_seed

    module subroutine ringmoons_torque_yarkovsky_schach_ring(self, cb, param, Torque)
        !! author: Kaustub P. Anand
        !! 
        !! Calculates the Yarkovsky-Schach torque acting on the ring bin. Torque is averaged over 1 orbit around the planet for a given angular tilt.
        !!
        !! References:
        !! Ferich, et al, 2022 (ADD doi)
        !! Veras, et al 2015 (ADD doi)
        implicit none
        ! Arguments
        class(ringmoons_ring),      intent(inout)               :: self
        class(swiftest_cb),         intent(in)                  :: cb 
        class(swiftest_parameters), intent(in)                  :: param
        real(DP),dimension(0:self%nbins+1), intent(out)         :: Torque
        ! Internals
        real(DP), dimension(0:self%nbins+1)                     :: YS_Torque
        real(DP), dimension(0:self%nbins+1)                     :: a_ys_mag ! YS acceleration magnitude

        associate(ring => self, nbins => self%nbins)
            a_ys_mag(1:nbins) = ring%rot_k * (1 - ring%albedo) * param%L_SUN_sys * sqrt(param%inv_c2) &
                                    / (16.0_DP * PI * (ring%a_pl)**2 * ring%sigma(1:nbins))

            YS_Torque(:) = -1.0_DP * a_ys_mag(1:nbins) * ring%r(1:nbins) * sin(ring%delta(1:nbins) * DEG2RAD / 2.0_DP) * ring%Y_21(1:nbins) / PI 
            Torque(:) = Torque(:) + YS_Torque(:)
        end associate

    end subroutine ringmoons_torque_yarkovsky_schach_ring

end submodule s_ringmoons_torque
