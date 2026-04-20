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
    
    module function ringmoons_torque_lindblad_ring(self,cb,asat,esat,isat,msat,param) result(Torque)
        !! author: David A. Minton
        !!
        !! Calculates the lindblad torques between each ring element and given  given satellite. 
        !! 
        !! The function unction returns total torque acting on the satellite, and stores the torques acting on each ring element 
        !! in the ring.
        !!
        !! Adapted from Andy Hesselbrock's ringmoons Python scripts.
        implicit none
        ! Arguments
        class(ringmoons_ring),      intent(inout) :: self
        class(swiftest_cb),         intent(in)    :: cb 
        real(DP),                   intent(in)    :: asat,esat,isat,msat
        class(swiftest_parameters), intent(in)    :: param
        real(DP),dimension(0:self%nbins+1)        :: Torque
        ! Internals
        integer(I4B)                           :: i,j,m,inner_outer_sign,il,w,w1,w2,js, mshep
        real(DP)                               :: dTorque, beta, Amk, nw,lap,dlap,da3,ar,Xr,Xs,Xlo,Xhi,rlo,rhi,Gfac,lind_factor,awidth
        real(DP), parameter                    :: g = 2.24_DP
        integer(I4B),parameter :: M_MAX = 100     
            !! Maximum number of Lindblad modes to compute
        real(DP),dimension(2:M_MAX)            :: aring
        logical,dimension(0:self%nbins+1)      :: T_mask
        integer(I4B),dimension(2:M_MAX)        :: w1_arr,w2_arr
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


        associate(ring => self)
            ! Mask out any ring bins that don't have enough mass in them
            where (ring%sigma(0:ring%nbins+1) > 1000 * VSMALL)
                T_mask(0:ring%nbins+1) = .true.
            elsewhere
                T_mask(0:ring%nbins+1) = .false. 
            end where

            Xs = 2 * sqrt(asat)
            ! Resonance width: See Longaretti sec. 5.3.1
            awidth = asat * sqrt(msat/cb%mass)
            Xlo = ring%X_inner + ring%deltaX * ring%inside
            Xhi = ring%X_outer
            rlo = 0.25_DP * Xlo**2
            rhi = 0.25_DP * Xhi**2

            ! Just do the first order resonances for now. The full suite of resonances will come later
            Torque(:) = 0.0_DP

            ! Calculate the number of modes that are separated by at least 1 bin width in X space
            mshep = max(2,min(M_MAX - 1,ceiling(0.5_DP * (sqrt(1._DP + 4._DP / 3._DP * Xs / ring%deltaX) - 1._DP))))
            
            ! Inner then outer Lindblads
            do il = -1,1,2
                ! Calculate resonance location space
                aring(2:mshep+1) = 0.25_DP * (Xs * marr(2:mshep+1,il))**2

                ! Calculate bin boundaries of resonance using its width in X space
                w1_arr(:) = 0
                w2_arr(:) = 0
                do m=2,mshep
                    if((aring(m) > rlo).and.(aring(m) < rhi)) then
                        w1_arr(m)=ring%find_bin(aring(m)-awidth)
                        w2_arr(m)=ring%find_bin(aring(m)+awidth)
                    end if
                end do
                
                do m  = 2, mshep
                    if ((aring(m) > rlo).and.(aring(m) < rhi)) then
                        Xr = 2 * sqrt(aring(m))
                        beta = (Xs / Xr)**(il * 2)
                        lap  =  lapm(m,il)
                        dlap = dlapm(m,il)
                        Amk = 0.5_DP * (lap + dlap)
                        w1 = w1_arr(m)
                        w2 = w2_arr(m)
                        nw = real(abs(w2 - w1) + 1,kind=DP)
                        ! Calculate the 1st order Lindblad torques and distribute them over the bins that include the resonance
                        lind_factor = il * mfac(m) / nw * aring(m)**4 * (beta * Amk)**2 
                        where(T_mask(w1:w2)) 
                            Torque(w1:w2) = Torque(w1:w2) + lind_factor * ring%sigma(w1:w2) * (ring%wkep(w1:w2))**2
                        endwhere
                    end if
                end do

                ! Add in shepherding torque
                if ((aring(mshep+1) > rhi).or.(aring(mshep+1) < rlo)) cycle
                j = ring%find_bin(aring(mshep+1)) !ring location of resonance
                da3 = il * max(abs((aring(mshep+1) - asat)**3),epsilon(asat))
                if (T_mask(j)) then 
                    Torque(j) = Torque(j) + g**2 / 6._DP * aring(mshep+1)**3 &
                                          / da3 * ring%sigma(j) * (ring%wkep(j))**2 * aring(mshep+1)**4
                end if
            end do
        end associate
        return
    end function ringmoons_torque_lindblad_ring

    module subroutine ringmoons_torque_tidal_seed(self,cb,param)
        !! author: David A. Minton
        !!
        !! Calculates the tidal torque acting on the seed by the central body
        !! Constant Q tidal model, see eq. 16 in [Cheng, Lee, and Peale 2014](https://doi.org/10.1016/j.icarus.2014.01.046) 
        implicit none
        class(ringmoons_seed),      intent(inout) :: self
        class(swiftest_cb),         intent(in)    :: cb
        class(swiftest_parameters), intent(in)    :: param
        ! Internals
        real(DP),dimension(self%nbody) :: n
        associate(seed => self)
            n(:) = sqrt((seed%mu(:)) / seed%a(:)**3)
            seed%Ttide(:) = sign(1._DP,cb%rot(3) - n(:)) * &
                    3 * seed%a(:) * n * (cb%k2 / cb%Q) * &
                    (seed%mass(:) / cb%mass) * &
                    (cb%radius / seed%a(:))**5 
        end associate
        return
        end subroutine ringmoons_torque_tidal_seed


end submodule s_ringmoons_torque
