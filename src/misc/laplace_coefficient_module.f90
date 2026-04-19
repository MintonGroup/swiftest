! Copyright 2026 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module laplace_coefficients
    !! author: David A. Minton
    !! 
    !! Contains a 4th order Runge-Kutta-Fehlberg ODE solver and a linear system of equations solver
    use globals
    use, intrinsic :: ieee_exceptions
    private
    public :: compute_laplace_coefficient

    contains

        pure recursive function compute_laplace_coefficient(alpha,j,s,n) result(ans)
            !! author: David A. Minton
            !! 
            !! Computes the Laplace coefficient b_s^(j)(alpha) and its derivatives using the recursive relations given in 
            !! Murray & Dermott (1999), including the case where alpha is close to 1 where the G(y) series solution is used.
            implicit none
            ! Arguments
            real(DP),     intent(in)  :: alpha
                !! alpha is the ratio of the inner to outer semimajor axis (a_inner/a_outer) for the coefficient being computed
            real(DP),     intent(in)  :: s
                !! s is the order of the Laplace coefficient being computed. 
            integer(I4B), intent(in)  :: j
                !! j is the index of the Laplace coefficient being computed.
            integer(I4B), intent(in)  :: n
                !! n is the order of the derivative. n=0 gives the coefficient itself, n=1 gives the first derivative, etc.
            real(DP)                  :: ans         
                !! The value of the Laplace coefficient or its derivative, depending on n.
            ! Internals
            real(DP) :: T1,T2,T3,T4
            
            if (n == 0) then
                ans = laplace(alpha,j,s)
                return
            else if (n==1) then
                T1 = laplace(alpha,j - 1,s + 1._DP)
                T2 = -2 * alpha * laplace(alpha,j,s + 1._DP)
                T3 = laplace(alpha, j + 1,s + 1._DP)
                ans = s * (T1 + T2 + T3)
                return
            else 
                T1 = compute_laplace_coefficient(alpha,j - 1,s + 1._DP,n - 1)
                T2 = -2 * alpha * compute_laplace_coefficient(alpha,j,s + 1._DP,n - 1)
                T3 = compute_laplace_coefficient(alpha,j + 1,s + 1._DP,n - 1)
                T4 = -2 * (n - 1) * compute_laplace_coefficient(alpha,j,s + 1._DP,n - 2)
                ans = s * (T1 + T2 + T3 + T4)
                return
            end if
        end function compute_laplace_coefficient


        pure function laplace(alpha,jp,s,ver) result(ans)
            implicit none
            real(DP),intent(in) :: alpha,s
            integer(I4B), intent(in) :: jp
            integer(I4B),intent(in),optional :: ver
            real(DP)           :: ans
            real(DP) :: num,denom,T1,tmp,F,F0
            real(DP),parameter :: tol=1e-25_DP
            real(DP),parameter :: alphaG=0.9_DP ! Switchover alpha to go to G series instead of F series
            integer(I4B) :: i,j,k,v
            real(DP),parameter :: x1=0.4_DP ! For G series near alpha=1. Intermediate value that converges fast using F series
            real(DP),parameter :: x2=0.8_DP ! For G series near alpha=1. Intermediate value that converges fast using F series
            real(DP),dimension(2,2) :: G,Ginv
            real(DP) :: A0,C2,F1,F2,det

            j = abs(jp)
            if (.not.present(ver)) then
                if (alpha < alphaG) then 
                    v = 1
                else
                    v = 2
                end if
            else 
                v = ver
            end if
            select case(v)
            case(1) ! F(x) series
                T1 = 1._DP
                do i = 0, j - 1
                    num = s + real(i,kind=DP)
                    denom = real(j - i, kind=DP)
                    T1 = T1 * num / denom
                end do
                T1 = T1 * alpha**j
                F = Fseries(alpha**2,s,j,tol)
                ans = 2._DP * T1 * F
            case(2) ! G(y) series: based on solution to problem 6.2 in Murray & Dermott (1999)

                ! Use an "easy" problem for the F(x) series to bootstrap a solution to the G(y) series coefficients A0 and C2
                G(1,1) = Gseries(1._DP-x1,1._DP,0._DP,s,jp,tol)
                G(1,2) = Gseries(1._DP-x1,0._DP,1._DP,s,jp,tol)
                G(2,1) = Gseries(1._DP-x2,1._DP,0._DP,s,jp,tol)
                G(2,2) = Gseries(1._DP-x2,0._DP,1._DP,s,jp,tol)
                det = (G(1,1) * G(2,2) - G(1,2) * G(2,1))
                Ginv(1,1) = G(2,2)
                Ginv(1,2) = -G(1,2)
                Ginv(2,1) = -G(2,1)
                Ginv(2,2) = G(1,1)
                Ginv = Ginv / det
                F1 = Fseries(x1,s,jp,tol)
                F2 = Fseries(x2,s,jp,tol)

                A0 = Ginv(1,1) * F1 + Ginv(1,2) * F2
                C2 = Ginv(2,1)*F1+Ginv(2,2)*F2
                F = Gseries(1._DP - alpha**2,A0,C2,s,j,tol)

                T1 = 1._DP
                do i = 0,j - 1
                    num = s + real(i, kind=DP)
                    denom = real(j-i, kind=DP)
                    T1 = T1 * num / denom
                end do
                T1 = T1 * alpha**j
                ans = 2._DP * T1 * F

            end select

            return
        end function laplace

        pure function Fseries(x,s,jp,tol) result(ans)
            implicit none
            real(DP),intent(in) :: x,s,tol
            integer(I4B),intent(in) :: jp
            real(DP) :: ans
            real(DP) :: F0,F,tmp,num,denom
            integer(I4B) :: i,k

            F0 = 0._DP
            F = 1._DP
            k = 1
            do 
                tmp = 1.0_DP
                do i=1,k
                    num = (s + real(i - 1, kind=DP)) * (s + real(jp + i - 1,kind=DP))
                    denom = real(i * (jp + i),kind=DP)
                    tmp = tmp * num / denom * x
                end do
                F0 = F
                F = F + tmp
                k = k + 1
                if (abs(1._DP - F / F0) < tol) exit 
            end do
            ans = F
            return
        end function Fseries

        pure function Gseries(y,A0,C2,s,jp,tol) result(ans)
            ! implements G(y) series (based on answer to problem 6.2 in Murray &  Dermott (1999))
            ! C2 is either B0 (for s=1/2) or A_{2s-1} (for s>1/2)
            implicit none
            real(DP),intent(in) :: y,A0,C2,s,tol
            integer(I4B),intent(in) :: jp
            real(DP) :: ans
            integer(I4B) :: s2,l,tic,toc
            integer(I4B),parameter :: lmax=10000000
            real(DP) :: G0,Bl2s1,Bl2s2,Bl,Al
            real(DP),dimension(:),allocatable :: Bold

            s2 = nint(2 * s)
            allocate(Bold(s2))
            ans = 1e-8_DP
            G0 = -1._DP
            Al = A0
            l = 0
            if (s2 == 1) then
                Bl = C2
                do l = 0,lmax
                    G0 = ans
                    ans = ans + Al * y**(l - s2 + 1) + log(y) * Bl * y**l
                    if (abs(1._DP - ans / G0) < tol) exit 
                    Bold(1) = Bl
                    Bl = Blp1(Bl,s,jp,l)
                    Al = Alp1(Al,Bold(1),Bl,s,jp,l)
                end do
            else
                Bl = B0func(A0,s,jp)
                Bl2s1 = 0._DP
                Bl2s2 = 0._DP
                tic = 1
                Bold(tic) = Bl
                do l =0, lmax
                    G0 = ans
                    ans = ans + Al*y**(l - s2 + 1) + log(y) * Bl * y**l
                    if ((l > s2 - 2) .and. (abs(1._DP - ans / G0) < tol)) exit 
                    if (l /= s2 - 2) then 
                        if (l - s2 + 1 >= 0) then
                        toc = tic - s2 + 1
                        if (toc <= 0) toc = size(Bold) + toc
                        Bl2s1 = Bold(toc)
                        end if
                        if (l - s2 + 2 >= 0) then
                        toc = tic - s2 + 2
                        if (toc <= 0) toc = size(Bold) + toc
                        Bl2s2 = Bold(toc)
                        end if
                        Al = Alp1(Al,Bl2s1,Bl2s2,s,jp,l)
                    else
                        Al = C2
                    end if
                    Bl = Blp1(Bl,s,jp,l)
                    tic = tic + 1
                    if (tic > size(Bold)) tic = 1
                    Bold(tic) = Bl ! tic keeps a looping array of old Bl values
                end do

            end if


            deallocate(Bold)

            return
        end function Gseries

        pure function Alp1(Al,Bl2s1,Bl2s2,s,jp,l) result(ans)
            ! evaluates B_l coefficients (based on answer to problem 6.2 in Murray &  Dermott (1999))
            implicit none
            real(DP),intent(in) :: Al,Bl2s1,Bl2s2,s
            integer,intent(in) :: jp,l
            real(DP) :: ans
            integer :: s2,n,start

            s2 = nint(2 * s)
            ans = Al * ((l - s2 + 1) * (l + jp + 1) + s * (s + real(jp, kind=DP))) + &
                    Bl2s1 * (2 * l - s2 + jp + 2) - Bl2s2 * (2 * l - s2 + 3)
            ans = ans / real(((l + 1) * (l - s2 + 2)),kind=DP)

            return
        end function Alp1

        pure function Blp1(Bl,s,jp,l) result(ans)
            ! evaluates B_l coefficients (based on answer to problem 6.2 in Murray & Dermott (1999) )
            implicit none
            real(DP),intent(in) :: Bl,s
            integer(I4B),intent(in) :: jp,l
            real(DP) :: ans
            integer(I4B) :: n

            ans = Bl * (l * (2 * s + jp + l) + s * (s + real(jp,kind=DP))) / ((l + 1) * (2 * s + real(l,kind=DP)))

            return
        end function Blp1

        pure function B0func(A0,s,jp) result(ans)
            implicit none
            real(DP),intent(in) :: A0,s
            integer,intent(in) :: jp
            real(DP) :: ans
            integer :: n,s2

            s2 = nint(2 * s)
            ans = A0
            do n = 1, s2 - 2
                ans = ans * (real((n - s2) * (n + jp),kind=DP) + s * (s + real(jp,kind=DP))) &
                            / real(n * (n - s2 + 1),kind=DP)
            end do
            ans = ans * (real(1 - s2 - jp,kind=DP) + s * (s + real(jp,kind=DP))) / real(s2 - 1,kind=DP)
            
            return
        end function B0func

end module laplace_coefficients