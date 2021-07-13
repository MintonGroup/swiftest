submodule(swiftest_classes) s_util_solve
   use swiftest
contains
   function util_solve_rkf45(f, y0in, t1, dt0, tol) result(y1)
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
         (/        1./4.,       1./4.,          0.,            0.,           0.,           0.,&
                  3./8.,      3./32.,      9./32.,            0.,           0.,           0.,&
               12./13., 1932./2197., -7200./2197.,  7296./2197.,           0.,           0.,&
                     1.,   439./216.,          -8.,   3680./513.,   -845./4104.,          0.,&
                  1./2.,     -8./27.,           2., -3544./2565.,   1859./4104.,    -11./40./), shape(rkf45_btab))
      real(DP), dimension(RKS),  parameter   :: rkf4_coeff =  (/ 25./216., 0., 1408./2565. ,  2197./4104. , -1./5.,      0. /)
      real(DP), dimension(RKS),  parameter   :: rkf5_coeff =  (/ 16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55. /)
      real(DP), dimension(:, :), allocatable :: k                !! Runge-Kutta coefficient vector
      real(DP), dimension(:),   allocatable  :: ynorm            !! Normalized y value used for adaptive step size control
      real(DP), dimension(:),   allocatable  :: y0               !! Value of y at the beginning of each substep
      integer(I4B)                           :: Nvar             !! Number of variables in problem
      integer(I4B)                           :: rkn              !! Runge-Kutta loop index
      real(DP)                               :: dt, trem         !! Current step size and total time remaining
      real(DP)                               :: s, yerr, yscale  !!  Step size reduction factor, error in dependent variable, and error scale factor
      integer(I4B)                           :: i, n     

      allocate(y0, source=y0in)
      allocate(y1, mold=y0)
      allocate(ynorm, mold=y0)
      Nvar = size(y0)
      allocate(k(Nvar, RKS))

      dt = dt0

      trem = t1
      do
         yscale = norm2(y0(:))
         do i = 1, MAXREDUX
            do rkn = 1, RKS
               y1(:) = y0(:) + matmul(k(:, 1:rkn - 1), rkf45_btab(2:rkn, rkn - 1))
               k(:, rkn) = dt * f%eval(y1(:))
            end do
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
         if (trem <= 0._DP) exit
         y0(:) = y1(:)
      end do

      return
   end function util_solve_rkf45

end submodule s_util_solve