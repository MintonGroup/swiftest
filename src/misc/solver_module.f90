!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module solver
   !! author: David A. Minton
   !! 
   !! Contains a 4th order Runge-Kutta-Fehlberg ODE solver and a linear system of equations solver
   use globals
   use base
   use lambda_function
   use, intrinsic :: ieee_exceptions
   private
   public :: solve_linear_system, solve_rkf45

   interface solve_linear_system
      module procedure solve_linear_system_dp
      module procedure solve_linear_system_qp
   end interface

   interface
      module function solve_rkf45(f, y0in, t1, dt0, tol) result(y1)
         implicit none
         class(lambda_obj),      intent(inout) :: f    !! lambda function object that has been initialized to be a function of derivatives. The object will return with components lastarg and lasteval set
         real(DP), dimension(:), intent(in)    :: y0in !! Initial value at t=0
         real(DP),               intent(in)    :: t1   !! Final time
         real(DP),               intent(in)    :: dt0  !! Initial step size guess
         real(DP),               intent(in)    :: tol  !! Tolerance on solution
         real(DP), dimension(:), allocatable   :: y1   !! Final result
      end function solve_rkf45
   end interface


   contains

      function solve_linear_system_dp(A,b,n,lerr) result(x)
         !! Author: David A. Minton
         !!
         !! Solves the linear equation of the form A*x = b for x. 
         !!   A is an (n,n) arrays
         !!   x and b are (n) arrays
         !! Uses Gaussian elimination, so will have issues if nbody_system is ill-conditioned.
         !! Uses quad precision intermidiate values, so works best on small arrays.
         implicit none
         ! Arguments
         integer(I4B),             intent(in)  :: n
         real(DP), dimension(:,:), intent(in)  :: A
         real(DP), dimension(:),   intent(in)  :: b
         logical,                  intent(out) :: lerr
         ! Result
         real(DP), dimension(n)                :: x
         ! Internals
         real(QP), dimension(:), allocatable :: qx
         type(ieee_status_type) :: original_fpe_status
         logical, dimension(:), allocatable :: fpe_flag 

         call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
         call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
         allocate(fpe_flag(size(ieee_usual)))

         qx = solve_wbs(ge_wpp(real(A, kind=QP), real(b, kind=QP)))

         call ieee_get_flag(ieee_usual, fpe_flag)
         lerr = any(fpe_flag) 
         if (lerr .or. (any(abs(qx) > huge(x))) .or. (any(abs(qx) < tiny(x)))) then
            x = 0.0_DP
         else
            x = real(qx, kind=DP)
         end if
         call ieee_set_status(original_fpe_status)

         return
      end function solve_linear_system_dp


      function solve_linear_system_qp(A,b,n,lerr) result(x)
         !! Author: David A. Minton
         !!
         !! Solves the linear equation of the form A*x = b for x. 
         !!   A is an (n,n) arrays
         !!   x and b are (n) arrays
         !! Uses Gaussian elimination, so will have issues if nbody_system is ill-conditioned.
         !! Uses quad precision intermidiate values, so works best on small arrays.
         implicit none
         ! Arguments
         integer(I4B),             intent(in) :: n
         real(QP), dimension(:,:), intent(in) :: A
         real(QP), dimension(:),   intent(in) :: b
         logical,                  intent(out) :: lerr
         ! Result
         real(QP), dimension(n)  :: x
         ! Internals
         type(ieee_status_type) :: original_fpe_status
         logical, dimension(:), allocatable :: fpe_flag 

         call ieee_get_status(original_fpe_status) ! Save the original floating point exception status
         call ieee_set_flag(ieee_all, .false.) ! Set all flags to quiet
         allocate(fpe_flag(size(ieee_usual)))

         x = solve_wbs(ge_wpp(A, b))

         call ieee_get_flag(ieee_usual, fpe_flag)
         lerr = any(fpe_flag) 
         if (lerr) x = 0.0_DP
         call ieee_set_status(original_fpe_status) 

         return
      end function solve_linear_system_qp


      function solve_wbs(u) result(x) ! solve with backward substitution
         !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
         implicit none
         ! Arguments
         real(QP), intent(in), dimension(:,:), allocatable  :: u
         ! Result
         real(QP), dimension(:), allocatable :: x
         ! Internals
         integer(I4B)             :: i,n

         n = size(u, 1)
         if (allocated(x)) deallocate(x)
         if (.not.allocated(x)) allocate(x(n))
         if (any(abs(u) < tiny(1._DP)) .or. any(abs(u) > huge(1._DP))) then 
            x(:) = 0._DP
            return
         end if
         call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
         do i = n, 1, -1 
            x(i) = (u(i, n + 1) - sum(u(i, i + 1:n) * x(i + 1:n))) / u(i, i)
         end do   
         return
      end function solve_wbs


      function ge_wpp(A, b) result(u) ! gaussian eliminate with partial pivoting
         !! Solve  Ax=b  using Gaussian elimination then backwards substitution.
         !!   A being an n by n matrix.
         !!   x and b are n by 1 vectors. 
         !! Based on code available on Rosetta Code: https://rosettacode.org/wiki/Gaussian_elimination#Fortran
         implicit none
         ! Arguments
         real(QP), dimension(:,:), intent(in) :: A
         real(QP), dimension(:),   intent(in) :: b
         ! Result
         real(QP), dimension(:,:), allocatable :: u
         ! Internals
         integer(I4B) :: i,j,n,p
         real(QP)     ::  upi

         n = size(a, 1)
         allocate(u(n, (n + 1)))
         u = reshape([A, b], [n, n + 1])
         call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
         do j = 1, n
            p = maxloc(abs(u(j:n, j)), 1) + j - 1 ! maxloc returns indices between (1, n - j + 1)
            if (p /= j) u([p, j], j) = u([j, p], j)
            u(j + 1:, j) = u(j + 1:, j) / u(j, j)
            do i = j + 1, n + 1
               upi = u(p, i)
               if (p /= j) u([p, j], i) = u([j, p], i)
               u(j + 1:n, i) = u(j + 1:n, i) - upi * u(j + 1:n, j)
            end do
         end do
         return
      end function ge_wpp


      module function solve_rkf45(f, y0in, t1, dt0, tol) result(y1)
         !! author: David A. Minton
         !!
         !! Implements the 4th order Runge-Kutta-Fehlberg ODE solver for initial value problems of the form f=dy/dt, y0 = y(t=0), solving for y1 = y(t=t1). Uses a 5th order adaptive step size control.
         !! Uses a lambda function object as defined in the lambda_function module
         implicit none
         ! Arguments   
         class(lambda_obj),      intent(inout) :: f    !! lambda function object that has been initialized to be a function of derivatives. The object will return with components lastarg and lasteval set
         real(DP), dimension(:), intent(in)    :: y0in !! Initial value at t=0
         real(DP),               intent(in)    :: t1   !! Final time
         real(DP),               intent(in)    :: dt0  !! Initial step size guess
         real(DP),               intent(in)    :: tol  !! Tolerance on solution
         ! Result
         real(DP), dimension(:), allocatable   :: y1  !! Final result
         ! Internals
         integer(I4B),                          parameter :: MAXREDUX = 1000 !! Maximum number of times step size can be reduced
         real(DP),                              parameter :: DTFAC = 0.95_DP !! Step size reduction safety factor (Value just under 1.0 to prevent adaptive step size control from discarding steps too aggressively)
         integer(I4B),                          parameter :: RKS = 6         !! Number of RK stages
         real(DP),     dimension(RKS, RKS - 1), parameter :: rkf45_btab = reshape( & !! Butcher tableau for Runge-Kutta-Fehlberg method
            (/       1./4.,       1./4.,          0.,            0.,           0.,           0.,&
                     3./8.,      3./32.,      9./32.,            0.,           0.,           0.,&
                  12./13., 1932./2197., -7200./2197.,  7296./2197.,           0.,           0.,&
                        1.,   439./216.,          -8.,   3680./513.,   -845./4104.,          0.,&
                     1./2.,     -8./27.,           2., -3544./2565.,   1859./4104.,    -11./40./), shape(rkf45_btab))
         real(DP), dimension(RKS),  parameter   :: rkf4_coeff =  (/ 25./216., 0., 1408./2565. ,  2197./4104. , -1./5.,      0. /)
         real(DP), dimension(RKS),  parameter   :: rkf5_coeff =  (/ 16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55. /)
         real(DP), dimension(:, :), allocatable :: k                !! Runge-Kutta coefficient vector
         real(DP), dimension(:),   allocatable  :: ynorm            !! Normalized y value used for adaptive step size control
         real(DP), dimension(:),   allocatable  :: y0           !! Value of y at the beginning of each substep
         integer(I4B)                           :: Nvar             !! Number of variables in problem
         integer(I4B)                           :: rkn              !! Runge-Kutta loop index
         real(DP)                               :: t, x1, dt, trem      !! Current time, step size and total time remaining
         real(DP)                               :: s, yerr, yscale  !!  Step size reduction factor, error in dependent variable, and error scale factor
         integer(I4B)                           :: i

         allocate(y0, source=y0in)
         allocate(y1, mold=y0)
         allocate(ynorm, mold=y0)
         Nvar = size(y0)
         allocate(k(Nvar, RKS))

         dt = dt0

         trem = t1
         t = 0._DP
         do
            yscale = norm2(y0(:))
            do i = 1, MAXREDUX
               select type(f)
               class is (lambda_obj_tvar)
                  do rkn = 1, RKS
                     y1(:) = y0(:) + matmul(k(:, 1:rkn - 1), rkf45_btab(2:rkn, rkn - 1))
                     if (rkn == 1) then
                        x1 = t
                     else
                        x1 = t + rkf45_btab(1,rkn-1)
                     end if
                     k(:, rkn) = dt * f%evalt(y1(:), t)
                  end do
               class is (lambda_obj)
                  do rkn = 1, RKS
                     y1(:) = y0(:) + matmul(k(:, 1:rkn - 1), rkf45_btab(2:rkn, rkn - 1))
                     k(:, rkn) = dt * f%eval(y1(:))
                  end do
               end select
               ! Now determine if the step size needs adjusting
               ynorm(:) = matmul(k(:,:), (rkf5_coeff(:) - rkf4_coeff(:))) / yscale
               yerr = norm2(ynorm(:)) 
               s = (tol / (2 * yerr))**(0.25_DP)
               dt = min(s * DTFAC * dt, trem) ! Alter step size either up or down, but never bigger than the remaining time
               if (s >= 1.0_DP) exit ! Good step!
               if (i == MAXREDUX) then
                  write(*,*) "Something has gone wrong in util_solve_rkf45!! Step size reduction has gone too far this time!"
                  call util_exit(FAILURE)
               end if
            end do
         
            ! Compute new value then step ahead in time
            y1(:) = y0(:) + matmul(k(:, :), rkf4_coeff(:))
            trem = trem - dt
            t = t + dt
            if (trem <= 0._DP) exit
            y0(:) = y1(:)
         end do

         return
      end function solve_rkf45


end module solver