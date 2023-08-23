!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module lambda_function
   !! author: David A. Minton
   !!
   !! Defines a class that can enable objects that behave like lambda functions.
   !! 
   !! To use this class, define a type of either lambda_obj or lambda_obj_err, or extend the lambda_obj class as necessary, such that an interface that matches the function you wish to lambdafy.
   !! Once defined, the lambda object can evaluate itself by calling the type-bound procedure eval. e.g. f%eval(x) (or f%eval(x, lerr), f%eval(x, [argument list], etc)) 
   !! 
   !! ********************************************************************************************************************************************************************************************
   !! Example - Defining a lambda function f(x,rval,ival) where rval and ival are a real and integer argument, respectively. This implementation uses an abstract interface, though this is not
   !! strictly necessary unless you want to bind more than one function with the same interface.
   !! ********************************************************************************************************************************************************************************************
   !!
   !! module lambda_new 
   !!    use swiftest ! This will bring in the lambda_function module
   !!    ! Define types in a module
   !!
   !!    type, extends(lambda_obj) :: lambda_obj_ri_args
   !!       procedure(abstract_lambda_ri_args), pointer, nopass :: lambdaptr_ri_args => null()
   !!       real(DP)                                            :: rval     !! Real parameter
   !!       integer(I4B)                                        :: ival     !! Integer paramete
   !!    contains
   !!       generic           :: init => lambda_ri_args_init
   !!       procedure         :: eval => lambda_ri_args_eval
   !!       procedure, nopass :: lambda_ri_args_init
   !!       final             :: lambda_ri_args_destroy
   !!    end type
   !!    interface lambda_obj
   !!       module procedure lambda_ri_args_init
   !!    end interface
   !! 
   !!    abstract interface
   !!       function abstract_lambda_ri_args(x, rval, ival) result(y)
   !!          !Template for the lambda function
   !!          import DP, I4B
   !!          real(DP), dimension(:), intent(in) :: x        !! Dependent variable
   !!          real(DP),               intent(in) :: rval     !! Real parameter
   !!          integer(I4B),           intent(in) :: ival     !! Integer parameter
   !!          real(DP)                           :: y        !! Real result
   !!       end function
   !!    end interface
   !! 
   !! contains
   !!    type(lambda_obj_ri_args) function lambda_ri_args_init(lambda, rval, ival)
   !!       !! Initializes the lambda function parameters (can be used as a structure constructor)
   !!       implicit none
   !!       procedure(abstract_lambda_ri_args)  :: lambda !! The lambda function that will be passed
   !!       real(DP),     intent(in)            :: rval     !! Real parameter
   !!       integer(I4B), intent(in)            :: ival     !! Integer parameter     
   !!    
   !!       ! Assign the procedure passed to this function to the procedure pointer 
   !!       lambda_ri_args_init%lambdaptr_ri_args => lambda 
   !! 
   !!       ! Assign the argument values
   !!       lambda_ri_args_init%rval = rval 
   !!       lambda_ri_args_init%ival = ival
   !!       return
   !!    end function lambda_ri_args_init
   !!
   !!    function lambda_ri_args_eval(self, x) result(y)
   !!       !! Defines the evaluation method, allowing the lambda function to be called with a single argument
   !!       implicit none
   !!       class(lambda_obj_ri_args),      intent(inout) :: self
   !!       real(DP), dimension(:),         intent(in) :: x
   !!       real(DP)                                   :: y
   !!    
   !!      if (associated(self%lambdaptr_ri_args)) then
   !!         y = self%lambdaptr_ri_args(x, self%rval, self%ival)
   !!          self%lastval = y
   !!          if (allocated(self%lastarg)) deallocate(self%lastarg)
   !!          allocate(self%lastarg, source=x)
   !!       else
   !!          stop "Lambda function was not initialized"
   !!       end if
   !!    end function lambda_ri_args_eval
   !!
   !!    subroutine lambda_ri_args_destroy(self)
   !!       !! Finalizer method. Use this as a template for cleaning up the object upon destruction, such as nullifying pointers 
   !!       implicit none
   !!       type(lambda_obj_ri_args) :: self
   !!       if (associated(self%lambdaptr_ri_args)) nullify(self%lambdaptr_ri_args)
   !!    end subroutine lambda_ri_args_destroy
   !!
   !!    function example_function(x, rval, ival) result(y)
   !!       !This is the actual function you are going to use as the lambda function. Its interface must match the abstract interface previously defined
   !!       implicit none
   !!       ! Arguments
   !!       real(DP), dimension(:), intent(in) :: x
   !!       real(DP),               intent(in) :: rval
   !!       integer(I4B),           intent(in) :: ival 
   !!       ! Result
   !!       real(DP)                           :: y
   !!       ! Internals
   !!       integer(I4B) :: i, n
   !!       n = size(x)
   !!       y = 42._DP * ival
   !!       do i = 1, n
   !!          y = y + x(i)**2
   !!       end do
   !!       return
   !!    end function example_function
   !! end module lambda_new
   !!
   !! program usage
   !!    use swiftest
   !!    use lambda_new
   !!    implicit none
   !!    type(lambda_obj_ri_args) :: f
   !!    real(DP) :: sigma_par
   !!    integer(I4B) :: iwonky, i,j
   !!    real(DP), dimension(12) :: xarr
   !!   
   !!    sigma_par = 3.14_DP
   !!    iwonky = 13
   !!
   !!    f = lambda_obj(example_function, sigma_par, iwonky)
   !!    do i = 1, 10
   !!       xarr(:) = [(j * 0.25_DP / i, j=1, 12)]
   !!       write(*,*) i,f%eval(xarr)
   !!    end do
   !! end program usage
   !! ********************************************************************************************************************************************************************************************
   
   use globals
   implicit none
   public

   type :: lambda_obj 
      !! Base class for an lambda function object. This object takes no additional arguments other than the dependent variable x, an array of real numbers
      procedure(lambda0), pointer, nopass :: lambdaptr => null()
      real(DP) :: lastval
      real(DP),dimension(:), allocatable :: lastarg
   contains
      generic   :: init => lambda_init_0
      procedure :: eval => lambda_eval_0
      procedure, nopass :: lambda_init_0
      final     :: lambda_destroy
   end type

   type, extends(lambda_obj) :: lambda_obj_err
      !! Extended class for an lambda function object. This object takes allows for the return of a logical error flag during evaluation of the function.
      procedure(lambda0err), pointer, nopass :: lambdaptr_err => null()
      logical   :: lerr     
   contains
      generic   :: init => lambda_init_0_err
      procedure :: eval => lambda_eval_0_err
      procedure, nopass :: lambda_init_0_err
   end type

   type, extends(lambda_obj) :: lambda_obj_tvar
      !! Base class for an lambda function object. This object takes no additional arguments other than the dependent variable x, an array of real numbers
      procedure(lambda0tvar), pointer, nopass :: lambdaptr_tvar => null()
      real(DP) :: t 
   contains
      generic   :: init => lambda_init_tvar
      procedure :: evalt => lambda_eval_tvar
      procedure, nopass :: lambda_init_tvar
   end type
   interface lambda_obj
      module procedure lambda_init_0
      module procedure lambda_init_0_err
      module procedure lambda_init_tvar
   end interface

   abstract interface
      function lambda0(x) result(y)
         ! Template for a 0 argument function
         import DP
         real(DP), dimension(:), intent(in) :: x
         real(DP)                           :: y
      end function

      function lambda0err(x, lerr) result(y)
         ! Template for a 0 argument function that returns an error value
         import DP
         real(DP), dimension(:), intent(in)  :: x
         logical,                intent(out) :: lerr
         real(DP)                            :: y
      end function

      function lambda0tvar(x, t) result(y)
         ! Template for a 0 argument function that returns an error value
         import DP
         real(DP), dimension(:), intent(in)  :: x
         real(DP),               intent(in)  :: t
         real(DP), dimension(:), allocatable :: y
      end function
   end interface

   contains
      type(lambda_obj) function lambda_init_0(lambda)
         implicit none
         ! Arguments
         procedure(lambda0)             :: lambda
         lambda_init_0%lambdaptr => lambda
         return
      end function lambda_init_0

      type(lambda_obj_err) function lambda_init_0_err(lambda, lerr)
         implicit none
         ! Arguments
         procedure(lambda0err)  :: lambda
         logical, intent(in) :: lerr
         lambda_init_0_err%lambdaptr_err => lambda
         lambda_init_0_err%lerr = lerr
         return
      end function lambda_init_0_err

      type(lambda_obj_tvar) function lambda_init_tvar(lambda, t)
         implicit none
         ! Arguments
         procedure(lambda0tvar)             :: lambda
         real(DP), intent(in)               :: t
         lambda_init_tvar%lambdaptr_tvar => lambda
         lambda_init_tvar%t = t
         return
      end function lambda_init_tvar
   
      function lambda_eval_0(self, x) result(y)
         implicit none
         ! Arguments
         class(lambda_obj),      intent(inout) :: self
         real(DP), dimension(:), intent(in) :: x
         ! Result
         real(DP)                      :: y
         if (associated(self%lambdaptr)) then
            y = self%lambdaptr(x)
            self%lastval = y
            if (allocated(self%lastarg)) deallocate(self%lastarg)
            allocate(self%lastarg, source=x)
         else
            stop "Lambda function was not initialized"
         end if
      end function lambda_eval_0

      function lambda_eval_0_err(self, x) result(y)
         implicit none
         ! Arguments
         class(lambda_obj_err),  intent(inout) :: self
         real(DP), dimension(:), intent(in) :: x
         ! Result
         real(DP)                      :: y
         if (associated(self%lambdaptr_err)) then
            y = self%lambdaptr_err(x, self%lerr)
            self%lastval = y
            if (allocated(self%lastarg)) deallocate(self%lastarg)
            allocate(self%lastarg, source=x)
         else
            stop "Lambda function was not initialized"
         end if
      end function lambda_eval_0_err

      function lambda_eval_tvar(self, x, t) result(y)
         implicit none
         ! Arguments
         class(lambda_obj_tvar), intent(inout) :: self
         real(DP), dimension(:), intent(in) :: x
         real(DP),               intent(in) :: t
         ! Result
         real(DP), dimension(:), allocatable :: y
         if (associated(self%lambdaptr_tvar)) then
            y = self%lambdaptr_tvar(x,t)
         else
            stop "Lambda function was not initialized"
         end if
      end function lambda_eval_tvar

      subroutine lambda_destroy(self)
         implicit none
         type(lambda_obj) :: self
         if (associated(self%lambdaptr)) nullify(self%lambdaptr)
      end subroutine lambda_destroy

end module lambda_function

