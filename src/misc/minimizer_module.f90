!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module minimizer
   !! author: David A. Minton
   !! 
   !! Includes the Broyden-Fletcher-Goldfarb-Shanno minimizer used by Fraggle
   use globals
   use lambda_function
   use solver
   use, intrinsic :: ieee_exceptions
   private
   public :: minimize_bfgs

   type(ieee_status_type) :: original_fpe_status
   logical, dimension(:), allocatable :: fpe_flag 

   interface 
      module subroutine minimize_bfgs(f, N, x0, eps, maxloop, lerr, x1)
         implicit none
         integer(I4B),           intent(in)               :: N
         class(lambda_obj),      intent(inout)            :: f
         real(DP), dimension(:), intent(in)               :: x0
         real(DP),               intent(in)               :: eps
         integer(I4B),           intent(in)               :: maxloop
         logical,                intent(out)              :: lerr
         real(DP), dimension(:), intent(out), allocatable :: x1
      end subroutine minimize_bfgs

      module function gradf(f, N, x1, dx, lerr) result(grad)
         implicit none
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x1
         real(DP),               intent(in)    :: dx
         logical,                intent(out)   :: lerr
         real(DP), dimension(N)                :: grad
      end function gradf

      module function minimize1D(f, x0, S, N, eps, lerr) result(astar)
         implicit none
         integer(I4B),           intent(in)  :: N
         class(lambda_obj),      intent(inout)  :: f
         real(DP), dimension(:), intent(in)  :: x0, S
         real(DP),               intent(in)  :: eps
         logical,                intent(out) :: lerr
         real(DP)                            :: astar
      end function minimize1D

      module function n2one(f, x0, S, N, a, lerr) result(fnew)
         implicit none
         integer(I4B),           intent(in) :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in) :: x0, S
         real(DP),               intent(in) :: a
         logical,                intent(out) :: lerr
         real(DP) :: fnew
      end function n2one

      module subroutine bracket(f, x0, S, N, gam, step, lo, hi, lerr)
         implicit none
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(in)    :: gam, step
         real(DP),               intent(inout) :: lo
         real(DP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
      end subroutine bracket

      module subroutine golden(f, x0, S, N, eps, lo, hi, lerr) 
         implicit none
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(in)    :: eps
         real(DP),               intent(inout) :: lo
         real(DP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
      end subroutine golden

      module subroutine quadfit(f, x0, S, N, eps, lo, hi, lerr) 
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(in)    :: eps
         real(DP),               intent(inout) :: lo
         real(DP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
      end subroutine quadfit
   end interface

   contains

      module subroutine minimize_bfgs(f, N, x0, eps, maxloop, lerr, x1)
         !! author: David A. Minton
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function implements the Broyden-Fletcher-Goldfarb-Shanno method to determine the minimum of a function of N variables.  
         !! It recieves as input:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   N       : Number of variables of function f
         !!   x0      : Initial starting value of x
         !!   eps     : Accuracy of 1 - dimensional minimization at each step
         !!   maxloop : Maximum number of loops to attempt to find a solution
         !! The outputs include
         !!   lerr :  Returns .true. if it could not find the minimum
         !! Returns
         !!   x1   :  Final minimum (all 0 if none found)
         !!   0 = No miniumum found
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0
         real(DP),               intent(in)    :: eps
         integer(I4B),           intent(in)    :: maxloop
         logical,                intent(out)   :: lerr
         ! Result
         real(DP), dimension(:), intent(out), allocatable :: x1
         ! Internals
         integer(I4B) ::  i, j, k, l, conv
         real(DP), parameter     :: graddelta = 1e-4_DP !! Delta x for gradient calculations
         real(DP), dimension(N) :: S               !! Direction vectors 
         real(DP), dimension(N,N) :: H             !! Approximated inverse Hessian matrix 
         real(DP), dimension(N) :: grad1           !! gradient of f 
         real(DP), dimension(N) :: grad0           !! old value of gradient 
         real(DP) :: astar                         !! 1D minimized value 
         real(DP), dimension(N) :: y, P
         real(DP), dimension(N,N) :: PP, PyH, HyP
         real(DP), save :: yHy, Py
         real(DP) :: Hnorm

         call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
         call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
         if (.not.allocated(fpe_flag)) allocate(fpe_flag(size(ieee_usual)))

         lerr = .false.
         if (allocated(x1)) deallocate(x1)
         allocate(x1, source=x0)
         ! Initialize approximate Hessian with the identity matrix (i.e. begin with method of steepest descent) 
         ! Get initial gradient and initialize arrays for updated values of gradient and x
         H(:,:) = reshape([((0._DP, i=1, j-1), 1._DP, (0._DP, i=j+1, N), j=1, N)], [N,N])  
         grad0 = gradf(f, N, x0(:), graddelta, lerr)
         if (lerr) then
            call ieee_set_status(original_fpe_status)
            return
         end if
         grad1(:) = grad0(:)
         do i = 1, maxloop 
            !check for convergence
            conv = count(abs(grad1(:)) > eps)
            if (conv == 0) exit 
            S(:) = -matmul(H(:,:), grad1(:))
            astar = minimize1D(f, x1, S, N, graddelta, lerr)
            if (lerr) exit
            ! Get new x values 
            P(:) = astar * S(:) 
            x1(:) = x1(:) + P(:)
            ! Calculate new gradient
            grad0(:) = grad1(:)
            grad1 = gradf(f, N, x1, graddelta, lerr)
            y(:) = grad1(:) - grad0(:)
            Py = sum(P(:) * y(:))
            ! set up factors for H matrix update 
            yHy = 0._DP
            !$omp do simd schedule(static)&
            !$omp firstprivate(N, y, H) &
            !$omp reduction(+:yHy)
            do k = 1, N 
               do j = 1, N
                  yHy = yHy + y(j) * H(j,k) * y(k)
               end do
            end do
            !$omp end do simd
            ! prevent divide by zero (convergence) 
            if (abs(Py) < N**2 * tiny(Py)) exit
            ! set up update 
            PyH(:,:) = 0._DP
            HyP(:,:) = 0._DP
            !$omp parallel do default(private) schedule(static)&
            !$omp shared(N, PP, P, y, H) &
            !$omp reduction(+:PyH, HyP)
            do k = 1, N 
               do j = 1, N
                  PP(j, k) = P(j) * P(k)
                  do l = 1, N
                     PyH(j, k) = PyH(j, k) + P(j) * y(l) * H(l,k)
                     HyP(j, k) = HyP(j, k) + P(k) * y(l) * H(j,l)
                  end do
               end do
            end do
            !$omp end parallel do 
            ! update H matrix 
            H(:,:) = H(:,:) + ((1._DP - yHy / Py) * PP(:,:) - PyH(:,:) - HyP(:,:)) / Py
            ! Normalize to prevent it from blowing up if it takes many iterations to find a solution
            Hnorm = 0.0_DP
            do concurrent (j = 1:N,k=1:N,abs(H(j,k))>sqrt(10*tiny(1.0_DP)))
               Hnorm = Hnorm + H(j,k)**2
            end do
            Hnorm = sqrt(Hnorm) 
            ! Stop everything if there are any exceptions to allow the routine to fail gracefully
            call ieee_get_flag(ieee_usual, fpe_flag)
            if (any(fpe_flag)) exit 
            if (i == maxloop) then
               lerr = .true.
            end if
            where(abs(H(:,:)) < sqrt(10*tiny(1.0_DP)) )
               H(:,:) = 0.0_DP
            endwhere
            H(:,:) = H(:,:) / Hnorm
         end do
         call ieee_get_flag(ieee_usual, fpe_flag)
         lerr = lerr .or. any(fpe_flag)  
         call ieee_set_status(original_fpe_status)

         return 

      end subroutine minimize_bfgs


      module function gradf(f, N, x1, dx, lerr) result(grad)
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! Purpose:  Estimates the gradient of a function using a central difference
         !! approximation
         !! Inputs:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   N    :  number of variables N
         !!   x1   :  x value array
         !!   dx   :  step size to use when calculating derivatives
         !! Outputs: 
         !!   lerr : .true. if an error occurred. Otherwise returns .false.
         !! Returns
         !!   grad :  N sized array containing estimated gradient of f at x1
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x1
         real(DP),               intent(in)    :: dx
         logical,                intent(out)   :: lerr
         ! Result
         real(DP), dimension(N)                :: grad
         ! Internals
         integer(I4B) :: i, j
         real(DP), dimension(N) :: xp, xm
         real(DP) :: fp, fm
         logical :: lerrp, lerrm

         do i = 1, N
            do j = 1, N
               if (j == i) then
                  xp(j) = x1(j) + dx
                  xm(j) = x1(j) - dx
               else
                  xp(j) = x1(j)
                  xm(j) = x1(j)
               end if
            end do
            select type (f)
            class is (lambda_obj_err)
               fp = f%eval(xp)
               lerrp = f%lerr
               fm = f%eval(xm)
               lerrm = f%lerr
               lerr = lerrp .or. lerrm
            class is (lambda_obj)
               fp = f%eval(xp)
               fm = f%eval(xm)
               lerr = .false.
            end select
            grad(i) = (fp - fm) / (2 * dx)
            if (lerr) return
         end do
         return 
      end function gradf


      module function minimize1D(f, x0, S, N, eps, lerr) result(astar)
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This program find the minimum of a function of N variables in a single direction
         !! S using in sequence:
         !!    1.  A Bracketing method
         !!    2.  The golden section method
         !!    3.  A quadratic polynomial fit
         !! Inputs
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   N    :  Number of variables of function f
         !!   eps  :  Accuracy of 1 - dimensional minimization at each step
         !! Output
         !!   lerr : .true. if an error occurred. Otherwise returns .false.
         !! Returns
         !!   astar      :  Final minimum along direction S
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)  :: N
         class(lambda_obj),      intent(inout)  :: f
         real(DP), dimension(:), intent(in)  :: x0, S
         real(DP),               intent(in)  :: eps
         logical,                intent(out) :: lerr
         ! Result
         real(DP)                            :: astar
         ! Internals
         integer(I4B) :: num = 0
         real(DP), parameter :: step = 0.7_DP     !! Bracketing method step size   
         real(DP), parameter :: gam = 1.2_DP      !! Bracketing method expansion parameter   
         real(DP), parameter :: greduce = 0.2_DP  !! Golden section method reduction factor   
         real(DP), parameter :: greduce2 = 0.1_DP ! Secondary golden section method reduction factor   
         real(DP) :: alo, ahi                     !! High and low values for 1 - D minimization routines   
         real(DP), parameter :: a0 = epsilon(1.0_DP)       !! Initial guess of alpha   
      
         alo = a0
         call bracket(f, x0, S, N, gam, step, alo, ahi, lerr)
         if (lerr) then
            !write(*,*) "BFGS bracketing step failed!"
            !write(*,*) "alo: ",alo, "ahi: ", ahi
            return 
         end if
         if (abs(alo - ahi) < eps) then
            astar = alo
            lerr = .false.
            return 
         end if
         call golden(f, x0, S, N, greduce, alo, ahi, lerr)
         if (lerr) then
            !write(*,*) "BFGS golden section step failed!"
            return 
         end if
         if (abs(alo - ahi) < eps) then
            astar = alo
            lerr = .false.
            return 
         end if
         call quadfit(f, x0, S, N, eps, alo, ahi, lerr)
         if (lerr) then
            !write(*,*) "BFGS quadfit failed!"
            return 
         end if
         if (abs(alo - ahi) < eps) then
            astar = alo
            lerr = .false.
            return 
         end if 
         ! Quadratic fit method won't converge, so finish off with another golden section   
         call golden(f, x0, S, N, greduce2, alo, ahi, lerr)
         if (.not. lerr) astar = (alo + ahi) / 2.0_DP
         return 
      end function minimize1D


      module function n2one(f, x0, S, N, a, lerr) result(fnew)
         implicit none
         ! Arguments
         integer(I4B),           intent(in) :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in) :: x0, S
         real(DP),               intent(in) :: a
         logical,                intent(out) :: lerr
         ! Return
         real(DP) :: fnew
         ! Internals
         real(DP), dimension(N) :: xnew
         integer(I4B) :: i
         
         xnew(:) = x0(:) + a * S(:)
         fnew = f%eval(xnew(:))
         select type(f)
         class is (lambda_obj_err)
            lerr = f%lerr
         class is (lambda_obj)
            lerr = .false.
         end select
         return 
      end function n2one


      module subroutine bracket(f, x0, S, N, gam, step, lo, hi, lerr)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This subroutine brackets the minimum.  It recieves as input:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   gam  :  expansion parameter
         !!   step :  step size
         !!   lo   :  initial guess of lo bracket value
         !! The outputs include
         !!   lo   :  lo bracket
         !!   hi   :  hi bracket
         !!   lerr : .true. if an error occurred. Otherwise returns .false.
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(in)    :: gam, step
         real(DP),               intent(inout) :: lo
         real(DP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
         ! Internals
         real(DP) :: a0, a1, a2, atmp, da
         real(DP) :: f0, f1, f2
         integer(I4B) :: i, j
         integer(I4B), parameter :: MAXLOOP = 100 ! maximum number of loops before method is determined to have failed   
         real(DP), parameter :: eps = epsilon(lo) ! small number precision to test floating point equality   

         ! set up initial bracket points   
         a0 =  lo
         da = step
         a1 = a0 + da
         a2 = a0 + 2 * da
         f0 = n2one(f, x0, S, N, a0, lerr)
         if (lerr) return
         f1 = n2one(f, x0, S, N, a1, lerr)
         if (lerr) return
         f2 = n2one(f, x0, S, N, a2, lerr)
         if (lerr) return
         ! loop over bracket method until either min is bracketed method fails   
         do i = 1, MAXLOOP 
            if ((f0 > f1) .and. (f1 < f2)) then  ! Minimum was found   
               lo = a0
               hi = a2
               return 
            else if ((f0 >= f1) .and. (f1 > f2)) then ! Function appears to decrease   
               da = da * gam
               atmp = a2 + da
               a0 = a1
               a1 = a2
               a2 = atmp
               f0 = f1
               f1 = f2
               f2 = n2one(f, x0, S, N, a2, lerr)
            else if ((f0 < f1) .and. (f1 <= f2)) then ! Function appears to increase   
               da = da * gam
               atmp = a0 - da
               a2 = a1
               a1 = a0
               a0 = atmp
               f2 = f1
               f0 = n2one(f, x0, S, N, a0, lerr)
            else if ((f0 < f1) .and. (f1 > f2)) then ! We are at a peak. Pick the direction that descends the fastest
               da = da * gam
               if (f2 > f0) then ! LHS is lower than RHS
                  atmp = a2 + da
                  a0 = a1
                  a1 = a2
                  a2 = atmp
                  f0 = f1
                  f1 = f2
                  f2 = n2one(f, x0, S, N, a2, lerr)
               else ! RHS is lower than LHS
                  atmp = a0 - da
                  a2 = a1
                  a1 = a0
                  a0 = atmp
                  f2 = f1
                  f1 = f2
                  f0 = n2one(f, x0, S, N, a0, lerr)
               end if
            else if ((f0 > f1) .and. (abs(f2 - f1) <= eps)) then ! Decrasging but RHS equal   
               da = da * gam
               atmp = a2 + da
               a2 = atmp
               f2 = n2one(f, x0, S, N, a2, lerr)
            else if ((abs(f0 - f1) < eps) .and. (f1 < f2)) then ! Increasing but LHS equal   
               da = da * gam
               atmp = a0 - da
               a0 = atmp
               f0 = n2one(f, x0, S, N, a0, lerr)
            else  ! all values equal. Expand in either direction and try again
               a0 = a0 - da
               a2 = a2 + da
               f0 = n2one(f, x0, S, N, a0, lerr)
               if (lerr) exit ! An error occurred while evaluating the function
               f2 = n2one(f, x0, S, N, a2, lerr)
            end if
            if (lerr) exit ! An error occurred while evaluating the function
         end do
         lerr = .true.
         return ! no minimum found   
      end subroutine bracket


      module subroutine golden(f, x0, S, N, eps, lo, hi, lerr) 
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function uses the golden section method to reduce the starting interval lo, hi by some amount sigma.  
         !! It recieves as input:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   x0   :  Array of size N of initial x values
         !!   S    :  Array of size N that determines the direction of minimization
         !!   gam  :  expansion parameter
         !!   eps  :  reduction interval in range (0 < sigma < 1) such that:
         !!             hi(new) - lo(new) = eps * (hi(old) - lo(old))
         !!   lo   :  initial guess of lo bracket value
         !! The outputs include
         !!   lo   :  lo bracket
         !!   hi   :  hi bracket
         !!   lerr : .true. if an error occurred. Otherwise returns .false.
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(in)    :: eps
         real(DP),               intent(inout) :: lo
         real(DP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
         ! Internals 
         real(DP), parameter :: tau = 0.5_DP * (sqrt(5.0_DP) - 1.0_DP)  ! Golden section constant   
         integer(I4B), parameter :: MAXLOOP = 40 ! maximum number of loops before method is determined to have failed (unlikely, but could occur if no minimum exists between lo and hi)   
         real(DP) :: i0 ! Initial interval value   
         real(DP) :: a1, a2
         real(DP) :: f1, f2
         integer(I4B) :: i, j

         i0 =  hi - lo
         a1 =  hi - tau * i0
         a2 =  lo + tau * i0
         f1 = n2one(f, x0, S, N, a1, lerr)
         if (lerr) return
         f2 = n2one(f, x0, S, N, a2, lerr)
         if (lerr) return
         do i = 1, MAXLOOP 
            if (abs((hi - lo) / i0) <= eps) return ! interval reduced to input amount   
            if (f2 > f1) then
               hi = a2
               a2 = a1
               f2 = f1
               a1 = hi - tau * (hi - lo)
               f1 = n2one(f, x0, S, N, a1, lerr)
            else 
               lo = a1
               a1 = a2
               f2 = f1
               a2 = hi - (1.0_DP - tau) * (hi - lo)
               f2 = n2one(f, x0, S, N, a2, lerr)
            end if
            if (lerr) exit
         end do
         lerr = .true.
         return ! search took too many iterations - no minimum found   
      end subroutine golden


      module subroutine quadfit(f, x0, S, N, eps, lo, hi, lerr) 
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - 
         !! This function uses a quadratic polynomial fit to locate the minimum of a function
         !! to some accuracy eps.  It recieves as input:
         !!   f%eval(x) : lambda function object containing the objective function as the eval metho
         !!   lo    :  low bracket value
         !!   hi    :  high bracket value
         !!   eps   :  desired accuracy of final minimum location
         !! The outputs include
         !!   lo   :  final minimum location
         !!   hi   :  final minimum location
         !! Notes: Uses the ieee_exceptions intrinsic module to allow for graceful failure due to floating point exceptions, which won't terminate the run.
         !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -   
         implicit none
         ! Arguments
         integer(I4B),           intent(in)    :: N
         class(lambda_obj),      intent(inout) :: f
         real(DP), dimension(:), intent(in)    :: x0, S
         real(DP),               intent(in)    :: eps
         real(DP),               intent(inout) :: lo
         real(DP),               intent(out)   :: hi
         logical,                intent(out)   :: lerr
         ! Internals 
         integer(I4B), parameter :: MAXLOOP = 20 ! maximum number of loops before method is determined to have failed.   
         real(DP) :: a1, a2, a3, astar   ! three points for the polynomial fit and polynomial minimum   
         real(DP) :: f1, f2, f3, fstar   ! three function values for the polynomial and polynomial minimum   
         real(DP), dimension(3) :: row_1, row_2, row_3, rhs, soln        ! matrix for 3 equation solver (gaussian elimination)   
         real(DP), dimension(3,3) :: lhs
         real(DP) :: d1, d2, d3, aold, denom, errval
         integer(I4B) :: i

         lerr = .false.
         ! Get initial a1, a2, a3 values   
         a1 =  lo
         a2 =  lo + 0.5_DP * (hi - lo)
         a3 =  hi
         aold = a1
         astar = a2
         f1 = n2one(f, x0, S, N, a1, lerr)
         if (lerr) return
         f2 = n2one(f, x0, S, N, a2, lerr)
         if (lerr) return
         f3 = n2one(f, x0, S, N, a3, lerr)
         if (lerr) return
         do i = 1, MAXLOOP 
            ! check to see if convergence is reached and exit   
            errval = abs((astar - aold) / astar)
            call ieee_get_flag(ieee_usual, fpe_flag)
            if (any(fpe_flag)) then
               !write(*,*) 'quadfit fpe'
               !write(*,*) 'aold : ',aold
               !write(*,*) 'astar: ',astar
               lerr = .true.
               exit
            end if
            if (errval < eps) then
               lo = astar
               hi = astar
               exit
            end if
            ! Set up nbody_system for gaussian elimination equation solver   
            row_1 = [1.0_DP, a1, a1**2]
            row_2 = [1.0_DP, a2, a2**2]
            row_3 = [1.0_DP, a3, a3**2]
            rhs = [f1, f2, f3]
            lhs(1, :) = row_1
            lhs(2, :) = row_2
            lhs(3, :) = row_3
            ! Solve nbody_system of equations   
            soln(:) = solve_linear_system(lhs, rhs, 3, lerr)
            call ieee_set_flag(ieee_all, .false.) ! Set all flags back to quiet
            call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
            if (lerr) then
               !write(*,*) 'quadfit fpe:'
               !write(*,*) 'util_solve_linear_system failed'
               exit
            end if
            aold = astar
            if (soln(2) == soln(3)) then ! Handles the case where they are both 0. 0/0 is an unhandled exception
               astar = -0.5_DP
            else
               astar =  -soln(2) / (2 * soln(3))
            end if
            call ieee_get_flag(ieee_usual, fpe_flag)
            if (any(fpe_flag)) then
               !write(*,*) 'quadfit fpe'
               !write(*,*) 'soln(2:3): ',soln(2:3)
               !write(*,*) 'a1, a2, a3'
               !write(*,*) a1, a2, a3
               !write(*,*) 'f1, f2, f3'
               !write(*,*) f1, f2, f3
               lerr = .true.
               exit
            end if
            fstar = n2one(f, x0, S, N, astar, lerr)
            if (lerr) exit
            ! keep the three closest a values to astar and discard the fourth  
            d1 = abs(a1 - astar)
            d2 = abs(a2 - astar)
            d3 = abs(a3 - astar)

            if (d1 > d2) then
               if (d1 > d3) then
                  f1 = fstar
                  a1 = astar
               else if (d3 > d2) then
                  f3 = fstar
                  a3 = astar
               end if
            else 
               if (d2 > d3) then
                  f2 = fstar
                  a2 = astar
               else if (d3 > d1) then
                  f3 = fstar
                  a3 = astar
               end if
            end if
         end do
         if (lerr) return
         lo = a1
         hi = a3
         return 
      end subroutine quadfit


end module minimizer