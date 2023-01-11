!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module tides
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods used to determine close encounters
   use base
   use lambda_function
   implicit none
   public


   type, extends(lambda_obj_tvar) :: tides_derivs_func 
      !! Base class for an lambda function object. This object takes no additional arguments other than the dependent variable x, an array of real numbers
      procedure(tidederiv), pointer, nopass :: lambdaptr_tides_deriv 
      real(DP), dimension(:,:), allocatable :: rbeg
      real(DP), dimension(:,:), allocatable :: rend
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
      function tidederiv(x, t, dt, rbeg, rend) result(y)
         ! Template for a 0 argument function
         import DP, base_nbody_system
         real(DP), dimension(:),     intent(in) :: x
         real(DP),                     intent(in) :: t
         real(DP),                     intent(in) :: dt
         real(DP), dimension(:,:),     intent(in) :: rbeg
         real(DP), dimension(:,:),     intent(in) :: rend
         real(DP), dimension(:), allocatable    :: y
      end function
   end interface


   interface
      module subroutine tides_kick_getacch_pl(self, nbody_system)
         implicit none
         class(base_object),         intent(inout) :: self   !! Swiftest massive body object
         class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      end subroutine tides_kick_getacch_pl

      module function tides_derivs_init(lambda, dt, rbeg, rend) result(f)
         implicit none
         procedure(tidederiv)                     :: lambda
         real(DP),                     intent(in) :: dt
         real(DP), dimension(:,:),     intent(in) :: rbeg
         real(DP), dimension(:,:),     intent(in) :: rend
         type(tides_derivs_func)                  :: f
      end function tides_derivs_init

      module function tides_derivs_eval(self, x, t) result(y)
         class(tides_derivs_func), intent(inout) :: self
         real(DP), dimension(:), intent(in) :: x
         real(DP),                 intent(in) :: t
         real(DP), dimension(:), allocatable  :: y
      end function tides_derivs_eval

      module function tides_spin_derivs(rot_pl_cb, t, dt, rbeg, rend) result(drot) !! Need to add more arguments so we can pull in mass, radius, Ip, J2, etc...
         real(DP), dimension(:,:),     intent(in) :: rot_pl_cb !! Array of rotations. The last element is the central body, and all others are massive bodies
         real(DP),                     intent(in) :: t         !! Current time, which is used to interpolate the massive body positions
         real(DP),                     intent(in) :: dt        !! Total step size
         real(DP), dimension(:,:),     intent(in) :: rbeg
         real(DP), dimension(:,:),     intent(in) :: rend
         real(DP), dimension(:,:), allocatable    :: drot
      end function tides_spin_derivs

      module subroutine tides_step_spin_system(self, param, t, dt)
         implicit none
         class(base_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
         class(base_parameters),   intent(in)    :: param  !! Current run configuration parameters  
         real(DP),                     intent(in)    :: t     !! Simulation time
         real(DP),                     intent(in)    :: dt    !! Current stepsize
      end subroutine tides_step_spin_system
   
   end interface


end module 