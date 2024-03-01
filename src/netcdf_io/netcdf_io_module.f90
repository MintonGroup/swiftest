! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module netcdf_io
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Base type definitions. This allows the collision and encounter modules to be defined before the swiftest module.
   use globals
   use base
   implicit none
   public


   !! This derived datatype stores the NetCDF ID values for each of the variables included in the NetCDF data file. This is used as
   !! the base class defined in base
   type, abstract :: netcdf_parameters
      character(STRMAX) :: file_name
         !! Name of the output file
      logical :: lfile_is_open = .false.
         !! Flag indicating that the linked file is currently open
      integer(I4B) :: out_type
         !! output type (will be assigned either NF90_DOUBLE or NF90_FLOAT, depending on the user parameter)
      integer(I4B) :: id
         !! ID for the output file
      integer(I4B) :: tslot = 1
         !! The current time slot that gets passed to the NetCDF reader/writer
      integer(I4B) :: max_tslot = 0
         !! Records the last index value of time in the NetCDF file
      integer(I4B), dimension(:), allocatable :: idvals
         !! Array of id values in this NetCDF file
      integer(I4B) :: idslot = 1
         !! The current id slot that gets passed to the NetCDF reader/writer
      integer(I4B) :: max_idslot = 0
         !! Records the last index value of id in the NetCDF file

      ! Dimension ids and variable names
      character(NAMELEN) :: str_dimname = "string32"
         !! name of the character string dimension
      integer(I4B) :: str_dimid
         !! ID for the character string dimension
      character(NAMELEN) :: time_dimname = "time"
         !! name of the time dimension 
      integer(I4B) :: time_dimid
         !! ID for the time dimension 
      integer(I4B) :: time_varid
         !! ID for the time variable
      character(NAMELEN) :: name_dimname = "name"
         !! name of the particle name dimension
      integer(I4B) :: name_dimid
         !! ID for the particle name dimension
      integer(I4B) :: name_varid
         !! ID for the particle name variable
      character(NAMELEN) :: space_dimname = "space"
         !! name of the space dimension
      integer(I4B) :: space_dimid
         !! ID for the space dimension
      integer(I4B) :: space_varid
         !! ID for the space variable
      character(len=1), dimension(3) :: space_coords = ["x","y","z"]
         !! The space dimension coordinate labels

      ! Non-dimension ids and variable names
      character(NAMELEN) :: id_varname = "id"
         !! name of the particle id variable
      integer(I4B) :: id_varid
         !! ID for the id variable 
      character(NAMELEN) :: status_varname = "status"
         !! name of the particle status variable
      integer(I4B) :: status_varid
         !! ID for the status variable 
      character(NAMELEN) :: ptype_varname = "particle_type"
         !! name of the particle type variable
      integer(I4B) :: ptype_varid
         !! ID for the particle type variable
      character(NAMELEN) :: npl_varname = "npl"
         !! name of the number of active massive bodies variable
      integer(I4B) :: npl_varid
         !! ID for the number of active massive bodies variable
      character(NAMELEN) :: ntp_varname = "ntp"
         !! name of the number of active test particles variable
      integer(I4B) :: ntp_varid
         !! ID for the number of active test particles variable
      character(NAMELEN) :: nplm_varname = "nplm"
         !! name of the number of active fully interacting massive bodies variable (SyMBA)
      integer(I4B) :: nplm_varid 
         !! ID for the number of active fully interacting massive bodies variable (SyMBA)
      character(NAMELEN) :: a_varname = "a"
         !! name of the semimajor axis variable 
      integer(I4B) :: a_varid 
         !! ID for the semimajor axis variable 
      character(NAMELEN) :: e_varname = "e"
         !! name of the eccentricity variable 
      integer(I4B) :: e_varid 
         !! ID for the eccentricity variable 
      character(NAMELEN) :: inc_varname = "inc"
         !! name of the inclination variable 
      integer(I4B) :: inc_varid 
         !! ID for the inclination variable 
      character(NAMELEN) :: capom_varname = "capom"
         !! name of the long. asc. node variable 
      integer(I4B) :: capom_varid 
         !! ID for the long. asc. node variable 
      character(NAMELEN) :: omega_varname = "omega"
         !! name of the arg. of periapsis variable 
      integer(I4B) :: omega_varid 
         !! ID for the arg. of periapsis variable 
      character(NAMELEN) :: capm_varname = "capm"
         !! name of the mean anomaly variable 
      integer(I4B) :: capm_varid 
         !! ID for the mean anomaly variable 
      character(NAMELEN) :: varpi_varname = "varpi"
         !! name of the long. of periapsis variable 
      integer(I4B) :: varpi_varid 
         !! ID for the long. of periapsis variable 
      character(NAMELEN) :: lam_varname = "lam"
         !! name of the mean longitude variable 
      integer(I4B) :: lam_varid 
         !! ID for the mean longitude variable 
      character(NAMELEN) :: f_varname = "f"
         !! name of the true anomaly variable 
      integer(I4B) :: f_varid 
         !! ID for the true anomaly variable 
      character(NAMELEN) :: cape_varname = "cape"
         !! name of the eccentric anomaly variable 
      integer(I4B) :: cape_varid 
         !! ID for the eccentric anomaly variable 
      character(NAMELEN) :: rh_varname = "rh"
         !! name of the heliocentric position vector variable
      integer(I4B) :: rh_varid 
         !! ID for the heliocentric position vector variable 
      character(NAMELEN) :: vh_varname = "vh"
         !! name of the heliocentric velocity vector variable
      integer(I4B) :: vh_varid 
         !! ID for the heliocentric velocity vector variable 
      character(NAMELEN) :: gr_pseudo_vh_varname = "gr_pseudo_vh"
         !! name of the heliocentric pseudovelocity vector variable (GR)
      integer(I4B) :: gr_pseudo_vh_varid
         !! ID for the heliocentric pseudovelocity vector variable (used in GR)
      character(NAMELEN) :: Gmass_varname = "Gmass"
         !! name of the G*mass variable
      integer(I4B) :: Gmass_varid 
         !! ID for the G*mass variable
      character(NAMELEN) :: mass_varname = "mass"
         !! name of the mass variable
      integer(I4B) :: mass_varid 
         !! ID for the mass variable
      character(NAMELEN) :: rhill_varname = "rhill"
         !! name of the hill radius variable
      integer(I4B) :: rhill_varid 
         !! ID for the hill radius variable
      character(NAMELEN) :: radius_varname = "radius"
         !! name of the radius variable
      integer(I4B) :: radius_varid
         !! ID for the radius variable
      character(NAMELEN) :: Ip_varname = "Ip"
         !! name of the principal moment of inertial variable
      integer(I4B) :: Ip_varid 
         !! ID for the axis principal moment of inertia variable
      character(NAMELEN) :: rot_varname = "rot"
         !! name of the rotation vector variable
      integer(I4B) :: rot_varid 
         !! ID for the rotation vector variable
      character(NAMELEN) :: rotphase_varname = "rotphase"
         !! name of the rotation phase variable
      integer(I4B) :: rotphase_varid
         !! ID for the rotation phase variable
      character(NAMELEN) :: j2rp2_varname = "j2rp2"
         !! name of the j2rp2 variable
      integer(I4B) :: j2rp2_varid 
         !! ID for the j2 variable
      character(NAMELEN) :: j4rp4_varname = "j4rp4"
         !! name of the j4pr4 variable
      integer(I4B) :: j4rp4_varid 
         !! ID for the j4 variable
      character(NAMELEN) :: c_lm_varname = "c_lm"
         !! name for the c_lm array
      integer(I4B) :: c_lm_varid 
         !! ID for the c_lm aqrray
      character(NAMELEN) :: k2_varname = "k2"
         !! name of the Love number variable
      integer(I4B) :: k2_varid 
         !! ID for the Love number variable
      character(NAMELEN) :: q_varname = "Q"
         !! name of the energy dissipation variable
      integer(I4B) :: Q_varid 
         !! ID for the energy dissipation variable
      character(NAMELEN) :: ke_orb_varname = "KE_orb"
         !! name of the system orbital kinetic energy variable
      integer(I4B) :: KE_orb_varid
         !! ID for the system orbital kinetic energy variable
      character(NAMELEN) :: ke_spin_varname = "KE_spin"
         !! name of the system spin kinetic energy variable
      integer(I4B) :: KE_spin_varid
         !! ID for the system spin kinetic energy variable
      character(NAMELEN) :: pe_varname = "PE"
         !! name of the system potential energy variable
      integer(I4B) :: PE_varid 
         !! ID for the system potential energy variable
      character(NAMELEN) :: be_varname = "BE"
         !! name of the system binding energy variable
      integer(I4B) :: BE_varid 
         !! ID for the system binding energy variable
      character(NAMELEN) :: te_varname = "TE"
         !! name of the system binding energy variable
      integer(I4B) :: TE_varid 
         !! ID for the system binding energy variable
      character(NAMELEN) :: L_orbit_varname = "L_orbit"
         !! name of the orbital angular momentum vector variable
      integer(I4B) :: L_orbit_varid 
         !! ID for the system orbital angular momentum vector variable
      character(NAMELEN) :: L_spin_varname = "L_spin"
         !! name of the spin angular momentum vector variable
      integer(I4B) :: L_spin_varid 
         !! ID for the system spin angular momentum vector variable
      character(NAMELEN) :: L_escape_varname = "L_escape"
         !! name of the escaped angular momentum vector variable
      integer(I4B) :: L_escape_varid
         !! ID for the escaped angular momentum vector variable
      character(NAMELEN) :: E_collisions_varname = "E_collisions"
         !! name of the escaped angular momentum y variable 
      integer(I4B) :: E_collisions_varid
         !! ID for the energy lost in collisions variable
      character(NAMELEN) :: E_untracked_varname = "E_untracked" 
         !! name of the energy that is untracked due to loss (due to mergers and body energy for escaped bodies)
      integer(I4B) :: E_untracked_varid 
         !! ID for the energy that is untracked due to loss (due to mergers and body energy for escaped bodies)
      character(NAMELEN) :: GMescape_varname = "GMescape"
         !! name of the G*Mass of bodies that escape the system
      integer(I4B) :: GMescape_varid
         !! ID for the G*Mass of bodies that escape the system
      character(NAMELEN) :: origin_type_varname = "origin_type"
         !! name of the origin type variable 
      integer(I4B) :: origin_type_varid
         !! ID for the origin type
      character(NAMELEN) :: origin_time_varname = "origin_time"
         !! name of the time of origin variable
      integer(I4B) :: origin_time_varid
         !! ID for the origin time
      character(NAMELEN) :: collision_id_varname = "collision_id"
         !! name of the collision id variable
      integer(I4B) :: collision_id_varid
         !! Netcdf ID for the origin collision ID
      character(NAMELEN) :: origin_rh_varname = "origin_rh" 
         !! name of the heliocentric position vector of the body at the time of origin variable
      integer(I4B) :: origin_rh_varid
         !! ID for the origin position vector variable
      character(NAMELEN) :: origin_vh_varname = "origin_vh" 
         !! name of the heliocentric velocity vector of the body at the time of origin variable
      integer(I4B) :: origin_vh_varid
         !! ID for the origin velocity vector component
      character(NAMELEN) :: discard_time_varname = "discard_time"
         !! name of the time of discard variable
      integer(I4B) :: discard_time_varid
         !! ID for the time of discard variable
      character(NAMELEN) :: discard_rh_varname = "discard_rh" 
         !! name of the heliocentric position vector of the body at the time of discard variable
      integer(I4B) :: discard_rh_varid
         !! ID for the heliocentric position vector of the body at the time of discard variable
      character(NAMELEN) :: discard_vh_varname = "discard_vh" 
         !! name of the heliocentric velocity vector of the body at the time of discard variable
      integer(I4B) :: discard_vh_varid 
         !! ID for the heliocentric velocity vector of the body at the time of discard variable
      character(NAMELEN) :: discard_body_id_varname = "discard_body_id"
         !! name of the id of the other body involved in the discard
      integer(I4B) :: discard_body_id_varid
         !! ID for the id of the other body involved in the discard
      logical :: lpseudo_vel_exists = .false. 
         !! Logical flag to indicate whether or not the pseudovelocity vectors were present in an old file.

      ! Gravitational harmonics ids and variable names
      logical :: lc_lm_exists = .false.
         !! Logical flag to indicate whether or not the c_lm array was present in an old file.
      character(NAMELEN) :: sign_dimname = "sign"
         !! name of the sign dimension for c_lm
      integer(I4B) :: sign_dimid
         !! ID for sign dimension
      integer(I4B) :: sign_varid
         !! ID for sign variable
      integer(I4B), dimension(2) :: sign_coords = [-1,1]
         !! The sign dimension coordinate labels
      character(NAMELEN) :: l_dimname = "l"
         !! name of l dimension for c_lm
      integer(I4B) :: l_dimid
         !! ID for the l dimension for c_lm
      integer(I4B) :: l_varid
         !! ID for the l variable
      character(NAMELEN) :: m_dimname = "m"
         !! name of m dimension for c_lm
      integer(I4B) :: m_dimid
         !! ID for the m dimension for c_lm
      integer(I4B) :: m_varid
         !! ID for the m variable
      integer(I4B) :: m_dim_max
         !! Maximum value of the m dimension
      integer(I4B) :: l_dim_max
         !! Maximum value of the l dimension

   contains
      procedure :: close       => netcdf_io_close      
         !! Closes an open NetCDF file
      procedure :: find_tslot  => netcdf_io_find_tslot 
         !! Finds the time dimension index for a given value of t
      procedure :: find_idslot => netcdf_io_find_idslot
         !! Finds the id dimension index for a given value of id
      procedure :: get_idvals  => netcdf_io_get_idvals 
         !! Gets the valid id numbers currently stored in this dataset
      procedure :: sync        => netcdf_io_sync        
         !! Syncrhonize the disk and memory buffer of the NetCDF file (e.g. commit the frame files stored in memory to disk) 
   end type netcdf_parameters

   interface
      module subroutine netcdf_io_check(status, call_identifier)
         implicit none
         integer, intent (in) :: status 
            !! The status code returned by a NetCDF function
         character(len=*), intent(in), optional :: call_identifier 
            !! String that indicates which calling function caused the error for diagnostic purposes
      end subroutine netcdf_io_check
   
      module subroutine netcdf_io_close(self)
         implicit none
         class(netcdf_parameters),intent(inout) :: self
            !! Parameters used to identify a particular NetCDF dataset
      end subroutine netcdf_io_close

      module subroutine netcdf_io_get_idvals(self)
         implicit none
         class(netcdf_parameters), intent(inout) :: self
            !! Parameters used to identify a particular NetCDF dataset
      end subroutine netcdf_io_get_idvals

      module subroutine netcdf_io_find_tslot(self, t, tslot)
         implicit none
         class(netcdf_parameters), intent(inout) :: self  
            !! Parameters used to identify a particular NetCDF dataset
         real(DP),intent(in) :: t 
            !! The value of time to search for
         integer(I4B), intent(out) :: tslot
            !! The index of the time slot where this data belongs
      end subroutine netcdf_io_find_tslot

      module subroutine netcdf_io_find_idslot(self, id, idslot)
         implicit none
         class(netcdf_parameters), intent(inout) :: self 
            !! Parameters used to identify a particular NetCDF dataset
         integer(I4B), intent(in) :: id 
            !! The value of id to search for
         integer(I4B), intent(out) :: idslot 
            !! The index of the id slot where this data belongs
      end subroutine netcdf_io_find_idslot
   
      module subroutine netcdf_io_sync(self)
         implicit none
         class(netcdf_parameters), intent(inout) :: self 
            !! Parameters used to identify a particular NetCDF dataset
      end subroutine netcdf_io_sync
   end interface 


end module netcdf_io
