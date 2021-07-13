submodule(swiftest_classes) s_tides_step_spin
   use swiftest

   type, extends(lambda_obj_tvar) :: tides_derivs_func 
      !! Base class for an lambda function object. This object takes no additional arguments other than the dependent variable x, an array of real numbers
      procedure(tidederiv), pointer, nopass :: lambdaptr_tides_deriv 
      real(DP), dimension(:,:), allocatable :: xbeg
      real(DP), dimension(:,:), allocatable :: xend
      real(DP)                              :: dt
   contains
      generic   :: init => tides_derivs_init
      procedure :: evalt => tides_derivs_eval
      procedure, nopass :: tides_derivs_init
   end type
   interface lambda_obj
      module procedure tides_derivs_init
   end interface
   abstract interface
      function tidederiv(x, t, dt, xbeg, xend) result(y)
         ! Template for a 0 argument function
         import DP, swiftest_nbody_system
         real(DP), dimension(:),     intent(in) :: x
         real(DP),                     intent(in) :: t
         real(DP),                     intent(in) :: dt
         real(DP), dimension(:,:),     intent(in) :: xbeg
         real(DP), dimension(:,:),     intent(in) :: xend
         real(DP), dimension(:), allocatable    :: y
      end function
   end interface

contains
   module subroutine tides_step_spin_system(self, param, t, dt)
      !! author: Jennifer L.L. Pouplin and David A. Minton
      !!
      !! Integrates the spin equations for central and massive bodies of the system subjected to tides.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters  
      real(DP),                     intent(in)    :: t     !! Simulation time
      real(DP),                     intent(in)    :: dt    !! Current stepsize
      ! Internals
      real(DP), dimension(:), allocatable       :: rot0, rot1
      real(DP)                                  :: subt
      real(DP), parameter                       :: tol=1e-6_DP !! Just a guess at the moment
      real(DP)                                  :: subdt 

      associate(pl => self%pl, npl => self%pl%nbody, cb => self%cb)
         allocate(rot0(NDIM*(npl+1)))
         rot0 = [pack(pl%rot(:,1:npl),.true.), pack(cb%rot(:),.true.)]
         ! Use this space call the ode_solver, passing tides_spin_derivs as the function:
         subdt = dt / 20._DP
         !rot1(:) = util_solve_rkf45(lambda_obj(tides_spin_derivs, subdt, pl%xbeg, pl%xend), rot0, dt, subdt tol)
         ! Recover with unpack
         !pl%rot(:,1:npl) = unpack(rot1...
         !cb%rot(:) = unpack(rot1...
      end associate

      return
   end subroutine tides_step_spin_system

   function tides_spin_derivs(rot_pl_cb, t, dt, xbeg, xend) result(drot) !! Need to add more arguments so we can pull in mass, radius, Ip, J2, etc...
      !! author: Jennifer L.L. Pouplin and David A. Minton
      !!
      !! function used to calculate the derivatives that are fed to the ODE solver
      implicit none
      ! Arguments
      real(DP), dimension(:,:),     intent(in) :: rot_pl_cb !! Array of rotations. The last element is the central body, and all others are massive bodies
      real(DP),                     intent(in) :: t         !! Current time, which is used to interpolate the massive body positions
      real(DP),                     intent(in) :: dt        !! Total step size
      real(DP), dimension(:,:),     intent(in) :: xbeg
      real(DP), dimension(:,:),     intent(in) :: xend
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
         xinterp(:) = xbeg(:,i) + t / dt * (xend(:,i) - xbeg(:,i))
         ! Calculate Ncb and Npl as a function of xinterp
         !drot(:,i) = -Mcb / (Mcb + Mpl(i)) * (N_Tpl + N_Rpl)
         !drot(:,n) = drot(:,n) - Mcb / (Mcb + Mpl(i) * (N_Tcb + N_Rcb)
         !
      end do

   end function tides_spin_derivs

   function tides_derivs_eval(self, x, t) result(y)
      implicit none
      ! Arguments
      class(tides_derivs_func), intent(inout) :: self
      real(DP), dimension(:), intent(in) :: x
      real(DP),                 intent(in) :: t
      ! Result
      real(DP), dimension(:), allocatable  :: y
      if (associated(self%lambdaptr_tides_deriv)) then
         y = self%lambdaptr_tides_deriv(x, t, self%dt, self%xbeg, self%xend)
      else
         error stop "Lambda function was not initialized"
      end if
   end function tides_derivs_eval

   function tides_derivs_init(lambda, dt, xbeg, xend) result(f)
      implicit none
      ! Arguments
      procedure(tidederiv)                     :: lambda
      real(DP),                     intent(in) :: dt
      real(DP), dimension(:,:),     intent(in) :: xbeg
      real(DP), dimension(:,:),     intent(in) :: xend
      ! Result
      type(tides_derivs_func)                  :: f
      f%lambdaptr_tides_deriv => lambda
      f%dt = dt
      allocate(f%xbeg, source = xbeg)
      allocate(f%xend, source = xend)
      return
   end function tides_derivs_init
end submodule s_tides_step_spin