! submodule(swiftest_classes) s_util_solve
!    use swiftest
! contains
!    subroutine util_solve_rkf45(xv, dt0, tol)
!       !! author: David A. Minton
!       !!
!       !! Implements the 4th order Runge-Kutta-Fehlberg ODE solver for initial value problemswith 5th order adaptive step size control.
!       implicit none
!       ! Arguments
!       real(DP), dimension(:), intent(in) :: xv !! The dependent variable 
!       real(DP), intent(in) :: dt0, tol ! Output cadence time step (also used as initial step size guess) and error tolerance
!       integer :: i, n, nsteps ! The number of steps to generate output
!       integer, parameter :: maxredux = 1000 ! Maximum number of times step size can be reduced
!       real(DP),dimension(:), allocatable :: y,y0,ynorm !  Internal temporary variable used to store intermediate results until total number of steps is known
!       integer, parameter :: rks = 6 ! Number of RK stages
!       real(DP),dimension(rks, rks - 1),parameter :: rkf45_btab = reshape( & ! Butcher tableau for Runge-Kutta-Fehlberg method
!          (/        1./4.,       1./4.,          0.,            0.,           0.,           0.,&
!                   3./8.,      3./32.,      9./32.,            0.,           0.,           0.,&
!                12./13., 1932./2197., -7200./2197.,  7296./2197.,           0.,           0.,&
!                      1.,   439./216.,          -8.,   3680./513.,   -845./4104.,          0.,&
!                   1./2.,     -8./27.,           2., -3544./2565.,   1859./4104.,    -11./40./), shape(rkf45_btab))
!       real(DP),dimension(rks),parameter :: rkf4_coeff =  (/ 25./216., 0., 1408./2565. ,  2197./4104. , -1./5.,      0. /)
!       real(DP),dimension(rks),parameter :: rkf5_coeff =  (/ 16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55. /)
!       real(DP), dimension(:, :), allocatable :: k ! Runge-Kutta coefficient vector
!       integer :: rkn  ! Runge-Kutta loop index
!       integer :: ndim ! Number of dimensions of the problem
!       real(DP) :: dt, trem ! Current step size and total time remaining
!       real(DP) :: s, yerr, yscale  !  Step size reduction factor, error in dependent variable, and error scale factor
!       real(DP), parameter :: dtfac = 0.95_DP ! Step size reduction safety factor (Value just under 1.0 to prevent adaptive step size control from discarding steps too aggressively)
!       real(DP) :: dtmean  ! Mean step size 
!       integer :: ntot ! Total number of steps (used in mean step size calculation)
!       real(DP) :: xscale, vscale

!       ndim = size(xv, 1)
!       nsteps = size(xv, 2)
!       allocate(k(ndim, rks))
!       allocate(y(ndim))
!       allocate(y0(ndim))
!       allocate(ynorm(ndim))

!       dt = dt0
!       dtmean = 0.0_DP
!       ntot = 0

!       do n = 2, nsteps
!          y0(:) = xv(:, n - 1)
!          trem = dt0
!          do
!             yscale = norm2(y0(:))
!             xscale = norm2(y0(1:2))
!             vscale = norm2(y0(3:4))
!             do i = 1, maxredux
!                do rkn = 1, rks
!                   y(:) = y0(:) + matmul(k(:, 1:rkn - 1), rkf45_btab(2:rkn, rkn - 1))
!                   k(:, rkn) = dt * derivs(y(:))
!                end do
!                ! Now determine if the step size needs adjusting
!                ynorm(:) = matmul(k(:,:), (rkf5_coeff(:) - rkf4_coeff(:)))
!                ynorm(1:2) = ynorm(1:2) / xscale
!                ynorm(3:4) = ynorm(3:4) / vscale
!                !ynorm(:) = ynorm(:) / yscale
!                yerr = norm2(ynorm(:)) 
!                s = (tol / (2 * yerr))**(0.25_DP)
!                dt = min(s * dtfac * dt, trem) ! Alter step size either up or down
!                if (s >= 1.0_DP) exit ! Good step!
!                if (i == maxredux) then
!                   write(*,*) 'Something has gone wrong!!'
!                   stop
!                end if
!             end do
         
!             ! Compute new value
!             y(:) = y0(:) + matmul(k(:, :), rkf4_coeff(:))
!             trem = trem - dt
!             ntot = ntot + 1
!             dtmean = dtmean + dt
!             if (trem <= 0._DP) exit
!             y0(:) = y(:)
!          end do

!          xv(:,n) = y(:)
!       end do

!       dtmean = dtmean / ntot
!       write(*,*) 'Total number of steps taken: ',ntot
!       write(*,*) 'Mean step size: ', dtmean / (2 * pi)

!       deallocate(k,y,y0)

!       return

!    end subroutine util_solve_rkf45

! end submodule s_util_solve