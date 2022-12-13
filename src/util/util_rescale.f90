!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_rescale
   use swiftest
contains
   module subroutine util_rescale_system(self, param, mscale, dscale, tscale)
      !! author: David A. Minton
      !!
      !! Rescales an nbody system to a new set of units. Inputs are the multipliers on the mass (mscale), distance (dscale), and time units (tscale). 
      !! Rescales all united quantities in the system, as well as the mass conversion factors, gravitational constant, and Einstein's constant in the parameter object.
      implicit none
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters. Returns with new values of the scale vactors and GU
      real(DP),                     intent(in)    :: mscale, dscale, tscale !! Scale factors for mass, distance, and time units, respectively. 
      ! Internals
      real(DP) :: vscale

      param%MU2KG = param%MU2KG * mscale
      param%DU2M = param%DU2M * dscale
      param%TU2S = param%TU2S * tscale

      ! Calculate the G for the system units
      param%GU = GC / (param%DU2M**3 / (param%MU2KG * param%TU2S**2))

      if (param%lgr) then
         ! Calculate the inverse speed of light in the system units
         param%inv_c2 = einsteinC * param%TU2S / param%DU2M
         param%inv_c2 = (param%inv_c2)**(-2)
      end if

      vscale = dscale / tscale

      associate(cb => self%cb, pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody)

         cb%mass = cb%mass / mscale
         cb%Gmass = param%GU * cb%mass 
         cb%radius = cb%radius / dscale
         cb%rb(:) = cb%rb(:) / dscale
         cb%vb(:) = cb%vb(:) / vscale
         cb%rot(:) = cb%rot(:) * tscale
         pl%mass(1:npl) = pl%mass(1:npl) / mscale
         pl%Gmass(1:npl) = param%GU * pl%mass(1:npl) 
         pl%radius(1:npl) = pl%radius(1:npl) / dscale
         pl%rh(:,1:npl) = pl%rh(:,1:npl) / dscale
         pl%vh(:,1:npl) = pl%vh(:,1:npl) / vscale
         pl%rb(:,1:npl) = pl%rb(:,1:npl) / dscale
         pl%vb(:,1:npl) = pl%vb(:,1:npl) / vscale
         pl%rot(:,1:npl) = pl%rot(:,1:npl) * tscale

      end associate


      return
   end subroutine util_rescale_system 

end submodule s_util_rescale