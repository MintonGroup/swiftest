submodule(tides) s_tides_step_spin
   use swiftest

contains

   module subroutine tides_step_spin_system(self, param, t, dt)
      !! author: Jennifer L.L. Pouplin and David A. Minton
      !!
      !! Integrates the spin equations for central and massive bodies of the nbody_system subjected to tides.
      implicit none
      ! Arguments
      class(base_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(base_parameters),   intent(in)    :: param  !! Current run configuration parameters  
      real(DP),                     intent(in)    :: t     !! Simulation time
      real(DP),                     intent(in)    :: dt    !! Current stepsize
      ! Internals
      real(DP), dimension(:), allocatable       :: rot0, rot1
      real(DP)                                  :: subt
      real(DP), parameter                       :: tol=1e-6_DP !! Just a guess at the moment
      real(DP)                                  :: subdt 

      select type(self)
      class is (swiftest_nbody_system)
         associate(pl => self%pl, npl => self%pl%nbody, cb => self%cb)
            allocate(rot0(NDIM*(npl+1)))
            ! rot0 = [pack(pl%rot(:,1:npl),.true.), pack(cb%rot(:),.true.)]
            ! Use this space call the ode_solver, passing tides_spin_derivs as the function:
            ! subdt = dt / 20._DP
            ! rot1(:) = swiftest_util_solve_rkf45(lambda_obj(tides_spin_derivs, subdt, pl%rbeg, pl%rend), rot0, dt, subdt,tol)
            ! ! Recover with unpack
            ! pl%rot(:,1:npl) = unpack(rot1...
            ! cb%rot(:) = unpack(rot1...
         end associate
      end select

      return
   end subroutine tides_step_spin_system


   module function tides_spin_derivs(rot_pl_cb, t, dt, rbeg, rend) result(drot) !! Need to add more arguments so we can pull in mass, radius, Ip, J2, etc...
      !! author: Jennifer L.L. Pouplin and David A. Minton
      !!
      !! function used to calculate the derivatives that are fed to the ODE solver
      implicit none
      ! Arguments
      real(DP), dimension(:,:),     intent(in) :: rot_pl_cb !! Array of rotations. The last element is the central body, and all others are massive bodies
      real(DP),                     intent(in) :: t         !! Current time, which is used to interpolate the massive body positions
      real(DP),                     intent(in) :: dt        !! Total step size
      real(DP), dimension(:,:),     intent(in) :: rbeg
      real(DP), dimension(:,:),     intent(in) :: rend
      ! Internals
      real(DP), dimension(:,:), allocatable    :: drot
      real(DP), dimension(:), allocatable      :: flatrot
      real(DP), dimension(NDIM)                :: N_Tcb, N_Rcb, N_Tpl, N_Rpl, xinterp
      real(DP)                                 :: C_cb, C_pl, r_dot_rot_cb, r_dot_rot_pl, rmag
      integer(I4B)                             :: i, n


      n = size(rot_pl_cb,2)
      if (allocated(drot)) deallocate(drot) 
      allocate(drot, mold=rot_pl_cb)
      drot(:,:) = 0.0_DP
      do i = 1,n-1
         xinterp(:) = rbeg(:,i) + t / dt * (rend(:,i) - rbeg(:,i))
         ! Calculate Ncb and Npl as a function of xinterp
         !drot(:,i) = -Mcb / (Mcb + Mpl(i)) * (N_Tpl + N_Rpl)
         !drot(:,n) = drot(:,n) - Mcb / (Mcb + Mpl(i) * (N_Tcb + N_Rcb)
         !
      end do

      return
   end function tides_spin_derivs

   module function tides_derivs_eval(self, x, t) result(y)
      implicit none
      ! Arguments
      class(tides_derivs_func), intent(inout) :: self
      real(DP), dimension(:), intent(in) :: x
      real(DP),                 intent(in) :: t
      ! Result
      real(DP), dimension(:), allocatable  :: y
      if (associated(self%lambdaptr_tides_deriv)) then
         y = self%lambdaptr_tides_deriv(x, t, self%dt, self%rbeg, self%rend)
      else
         stop "Lambda function was not initialized"
      end if

      return
   end function tides_derivs_eval

   module function tides_derivs_init(lambda, dt, rbeg, rend) result(f)
      implicit none
      ! Arguments
      procedure(tidederiv)                     :: lambda
      real(DP),                     intent(in) :: dt
      real(DP), dimension(:,:),     intent(in) :: rbeg
      real(DP), dimension(:,:),     intent(in) :: rend
      ! Result
      type(tides_derivs_func)                  :: f
      f%lambdaptr_tides_deriv => lambda
      f%dt = dt
      allocate(f%rbeg, source = rbeg)
      allocate(f%rend, source = rend)

      return
   end function tides_derivs_init

end submodule s_tides_step_spin