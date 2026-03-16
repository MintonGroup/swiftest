##########################################################
A guide to adding a new kick to Swiftest
##########################################################

.. rubric:: by Kaustub Anand

Here, we provide a guide on how to add a kick or additional feature, including new variables and particle features, to Swiftest.
This guide is written in the context of adding a new kick or acceleration feature that requires new properties for the particles.
However, the general steps can be extended to adding other features as well. The guide is split up into two sections, the Fortran 
and Python sides of Swiftest. The former handles the calculations and data, while the latter helps with user interface.

Important nomenclature:
=======================

- Nbody system = ``nbody_system`` = nbody system object
- Central Body = ``cb`` = ``nbody_system%cb``
- Particles are of type ``body`` with subclasses:
    - Massive Particles = ``pl`` = type ``body`` of class ``swiftest_pl``
    - Test particles = ``tp`` = type ``body`` of class ``swiftest_tp``
- Input parameters = ``param`` = simulation paramters object
- NetCDF file object = ``nc``

- *file.f90* is a file or folder.
- ``code`` is a block of code.

Note: small details such as ensuring that internal variables are defined for a given code block and correct syntax are not mentioned and left to the user.

Fortran
========
The Fortran side of Swiftest is where all the calculations for the simulation are performed. Going into *swiftest/src/* To add a kick:

- Check and add flags for the kick to the parameter object ``param`` in *base_module.f90*:

    - define and add a flag for the kick as an attribute for the ``base_parameters`` class in *base_module.f90*.
        - For example, a flag for the Yarkovsky effect would be ``lyarkovsky``

        .. code-block:: Fortran

            type, abstract :: base_parameters
                .
                .
                logical :: lyarkovsky = .false. 
                    !! Turn on Yarkovsky effect
                .
                .   

- Write the kick in a new file *swiftest_\*.f90*:

    - Create a new fortran file in *swiftest/src/swiftest* following the naming convention of *swiftest_\*.f90* where * is the name of the kick. 
        - You can also add a subroutine for a kick inside another related file.
        - For example, a Yarkovsky related kick could have a new file called *swiftest_yarkovsky.f90* or be added to a file with other radiation kicks such as *swiftest_radiation.f90*.
    - The file has to be set up as a ``submodule`` that contains various ``module subroutines`` for individual kicks. Below is an example of a yarkovsky kick for massive particles.
        
        .. code-block:: Fortran

            submodule (swiftest) s_swiftest_yarkovsky

            contains

                module subroutine swiftest_getacch_yarkovsky_pl(self, nbody_system, param)
                    !!
                    !! Calculate the Yarkovsky effect on massive bodies. 
                    !!
                    implicit none
                    ! Arguments
                    class(swiftest_pl),           intent(inout) :: self
                        !! Swiftest body object
                    class(swiftest_nbody_system), intent(inout) :: nbody_system
                        !! Swiftest nbody system object
                    class(swiftest_parameters),   intent(in)    :: param
                        !! Current run configuration parameters
                    .
                    .
                end subroutine swiftest_getacch_yarkovsky_pl
                    .
                    .
                    .
            end submodule s_swiftest_yarkovsky
    - ``subroutines`` follow the naming convention of ``swiftest_getacch_NAME_PARTICLETYPE``.
    - Write the formulas for calculating the kick in the appropriate ``subroutine`` add the kick to the ``i``-th particle's acceleration vector.
        .. code-block:: Fortran
            
            .
            .
            module subroutine swiftest_getacch_yarkovsky_pl(self, nbody_system, param)
                implicit none
                .
                .
                a_kick(:) = a_kick_mag * a_kick_direction(:)
                pl%ah(:, i) = pl%ah(:, i) + a_kick(:)
                .
                .
            end subroutine swiftest_getacch_yarkovsky_pl
            .
            .

- Add the kick and necessary particle parmeters to *swiftest_module.f90*:
    
    - Define and add any new particle parameters required by the kick to the appropriate particle type(s) (``swiftest_pl``, ``swiftest_tp``, or ``swiftest_body`` for both particles). 
    - The variables should be defiend as ``allocatable`` with the right dimensions. The dimension variables in Swiftest are ``time``, ``name``, and ``space``. ``space`` is split up into ``x``, ``y``, and ``z`` and is for vector variables in 3D.

        - For example, we can add emissivity and albedo to massive particles (``pl``) for the Yarkovsky kicks. They do not change with time or space, and only vary per particle . Their dimension is thus only ``name``. 

        .. code-block:: Fortran

            type, abstract, extends(swiftest_body) :: swiftest_pl
                .
                .
                .
                real(DP),                dimension(:), allocatable   :: albedo
                    !! Bond albedo for radiation acceleration calculations
                real(DP),                dimension(:), allocatable   :: emissivity
                    !! Emissivity for Yarkovsky acceleration calculations
                .
                .
                .
            contains
                .
                .

    - Add the kick with a simple name as a procedure attribute for the appropriate particle type(s). 
        - For example, we can add the Yarkovsky effect to ``swiftest_pl`` only.

        .. code-block:: Fortran

            type, abstract, extends(swiftest_body) :: swiftest_pl
                .
                .
                .
            contains
                .
                .
            procedure :: accel_yarkovsky => swiftest_yarkovsky_getacch_pl 
                !! Compute the heliocentric accelerations of bodies due to the Yarkovsky effect
    
    - Add the kick submodule definition in an appropriate interface.

        .. code-block:: fortran
            
                .
                .
                .
            interface
                .
                .
                .
                module subroutine swiftest_yarkovsky_getacch_pl(self, nbody_system, param)
                    !!
                    !! Calculate the Yarkovsky effect on massive bodies. 
                    !!
                    implicit none
                    ! Arguments
                    class(swiftest_pl),           intent(inout) :: self
                        !! Swiftest body object
                    class(swiftest_nbody_system), intent(inout) :: nbody_system
                        !! Swiftest nbody system object
                    class(swiftest_parameters),   intent(in)    :: param
                        !! Current run configuration parameters
                end subroutine swiftest_yarkovsky_getacch_pl
                .
                .
                .
            end interface
                .
                .
                .

- Add the kicks to the code flow in *whm_kick.f90* and *helio_kick.f90*:

    - All the current integrators (RMVS, WHM, SyMBA, HELIO) run though *whm_kick.f90* and *helio_kick.f90*. 
    - Add the appropriate flag check and kick function per particle type in ``whm_kick_getacch_pl(...)``, ``whm_kick_getacch_tp(...)``, ``helio_kick_getacch_pl(...)``, and/or ``helio_kick_getacch_tp(...)``.
        
        .. code-block:: fortran

            module subroutine helio_kick_getacch_pl(self, nbody_system, param, t, lbeg)
                .
                .
                if (param%lyarkovsky) call pl%accel_yarkovsky(nbody_system, param)
                .
                .
                return
            end subroutine helio_kick_getacch_pl

- In the home directory of Swiftest, add the new files to *src/CMakeLists.txt* under the right section for fortran compilation:

    - A kick will typically go under ``STRICT_MATH_FILES``

    .. code-block:: Fortran

        # Add the source files
        SET(STRICT_MATH_FILES
            .
            .
            ${SRC}/swiftest/swiftest_yarkovsky.f90
            .
            .
        ) 

- General ``param`` I/O checks in *swiftest_io.f90*:

    - Here we will add in checks for inputting and outputting the data variables. This will ensure that variables are defined, read in, and printed out correctly.
    - We will start with reading the ``param`` object in ``swiftest_io_param_reader()``.
    - Define a simple name for the effect, ex: ``"YARKOVSKY"``. This will be used for reading the ``param`` file from Python and defined again later on the Python side.

    .. code-block:: fortran

        module subroutine swiftest_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
            .
            .
            associate(param => self) 
                .
                .
                case("YARKOVSKY")
                    call swiftest_io_toupper(param_value)
                    if (param_value == "YES" .or. param_value == 'T') param%lyarkovsky = .true.
                .
                .
            .
            .
        end subroutine swiftest_io_param_reader
    
    - In the same ``subroutine``, you can add other ``param`` flag checks. For example, the Yarkovsky kick also requires rotation of particles. 
    - We can add a check to ensure that they are both turned on outside the ``do`` loop and/or do any other required unit coversions for constants.

    .. code-block:: fortran

        module subroutine swiftest_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
            .
            .
            associate(param => self) 
                .
                .
                case("YARKOVSKY")
                    call swiftest_io_toupper(param_value)
                    if (param_value == "YES" .or. param_value == 'T') param%lyarkovsky = .true.
                .
                .
            .
            .
            if (param%lyarkovsky .and. .not. param%lrotation) then
                write(iomsg,*) 'Yarkovsky forces require rotation to be turned on'
                iostat = -1
                return
            end if
            .
            .
            ! Calculate Solar Luminosity in system units and turn on gr for inv_c2 calculation if radiation forces are enabled
            if (param%lradiation .or. param%lyarkovsky) then
                param%L_SUN_sys = L_SUN / param%MU2KG / param%DU2M**2 * param%TU2S**3
                param%sigma_sys = SIGMA /param%MU2KG * param%TU2S**3 ! system units / K^4
                param%lgr = .true.
            end if
            .
            .
        end subroutine swiftest_io_param_reader
                
    - We will ensure the flag is correct when writing an updated ``param`` object in ``swiftest_io_param_writer()``.

    .. code-block:: fortran

        module subroutine swiftest_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
            .
            .
            associate(param => self)
                .
                .
                call io_param_writer_one("YARKOVSKY", param%lyarkovsky, unit)
                .
                .
            .
            .   
        end subroutine swiftest_io_param_writer


Now we can tackle the data and particle parameter I/O. If you have not added any new particle parameters/features (in our example, ``albedo`` and ``emissivity``) or new variables to the data, you can ignore this section.
Swiftest uses hdf5 and NetCDF files for data handling. This is handled by Xarray in Python. We will split the I/O section into general Swiftest data checks and then NetCDF I/O checks.

- NetCDF data handling I/O checks in *netcdf_io_module.f90*

    - We assume the additions are not a new dimension variable.
    - Add new parameters' and variables' to the ``netcdf_parameters`` class. We will define a name of type ``character(NAMELEN)`` and a variable ID of type ``integer(I4B)``.

    .. code-block:: fortran

        !! This derived datatype stores the NetCDF ID values for each of the variables included in the NetCDF data file. This is used as
        !! the base class defined in base
        type, abstract :: netcdf_parameters
            .
            .
            character(NAMELEN) :: albedo_varname = "albedo"
                !! name of the albedo variable
            integer(I4B) :: albedo_varid 
                !! ID for the albedo variable
            character(NAMELEN) :: emissivity_varname = "emissivity"
                !! name of the emissivity variable
            integer(I4B) :: emissivity_varid 
                !! ID for the emissivity variable
            .
            .

- NetCDF data handling checks for new parameter variables in *swiftest_io.f90*:

    - A good resource for the NetCDF functions is here: https://docs.unidata.ucar.edu/netcdf-fortran/current/f90-variables.html

    - We will define the NetCDF variables and add checks when initalizing the output file in ``swiftest_io_netcdf_initialize_output()``.
    - The general format for adding the parameter definition and NetCDF check is as below. Other NetCDF functions (ex: ``nf90_def_var``, ``nf90_inq_varid``, etc.) follow a similar structure.

    .. code-block:: fortran

        call netcdf_io_check( nf90_def_var(NetCDF ID, new parameter variable name, NetCDF output type, [new parmeter dimension IDs], new parameter variable ID), &
                                "error message"  )
    
    - The dimension IDs are ``nc%space_dimid``, ``nc%name_dimid``, and ``nc%time_dimid``. The ``[]`` are only necessary with variables of more than 1 dimension and should follow the given order.
    - For a new parameter called ``newparameter`` that varies with ``space``, ``name``, and ``time``, this addition would look like this:
                            
    .. code-block:: fortran

        call netcdf_io_check( nf90_def_var(nc%id, nc%newparameter_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%newparameter_varid), &
                            "netcdf_io_initialize_output nf90_def_var newparameter_varid"  )
    
    - For a new parameter called ``newparameter2`` that varies with ``name``, and ``time``, this addition would look like this:
                            
    .. code-block:: fortran

        call netcdf_io_check( nf90_def_var(nc%id, nc%newparameter2_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%newparameter2_varid), &
                            "netcdf_io_initialize_output nf90_def_var newparameter2_varid"  )


    - It is also important to only define the new variables if the user sets the appropriate ``param`` flag. 
    - Continuing with the Yarkovsky code example, we add ``albedo`` and ``emissivity`` that only vary with the ``name`` dimension:

    .. code-block:: fortran

        module subroutine swiftest_io_netcdf_initialize_output(self, param)
            .
            .
            associate(nc => self)
                .
                .
                if (param%lyarkovsky) then
                    call netcdf_io_check( nf90_def_var(nc%id, nc%albedo_varname, nc%out_type, nc%name_dimid, nc%albedo_varid), &
                                        "netcdf_io_initialize_output nf90_def_var albedo_varid"  )
                    call netcdf_io_check( nf90_def_var(nc%id, nc%emissivity_varname, nc%out_type, nc%name_dimid, nc%emissivity_varid), &
                                        "netcdf_io_initialize_output nf90_def_var emissivity_varid"  )
                .
                .
            .
            .
        end subroutine swiftest_io_netcdf_initialize_output

    - Next, we handle the checks when opening the NetCDF file in ``swiftest_io_netcdf_open()``. Here we are reading in the variable IDs.

    .. code-block:: fortran

        module subroutine swiftest_io_netcdf_open(self, param, readonly)
            .
            .
            associate(nc => self)
                .
                .
                if (param%lyarkovsky) then
                    call netcdf_io_check( nf90_inq_varid(nc%id, nc%albedo_varname, nc%albedo_varid), &
                                        "swiftest_io_netcdf_open nf90_inq_varid albedo_varid" )
                    call netcdf_io_check( nf90_inq_varid(nc%id, nc%emissivity_varname, nc%emissivity_varid), &
                                        "swiftest_io_netcdf_open nf90_inq_varid emissivity_varid" )
                .
                .
        end subroutine swiftest_io_netcdf_open
    
    - After initialization, we can now read in the data in ``swiftest_read_frame_system()``. 
    - We will use ``nf90_get_var()`` that follows a similar but slightly different format with the ``start`` and ``count`` arguments. 
    - ``start`` defines the index of value to read per dimension and ``count`` defines the amount of indices along each dimension.
    - The general format is as so:

    .. code-block:: fortran

        call netcdf_io_check( nf90_get_var(NetCDF ID, new parameter variable ID, temporary_variable, start=[...], count=[...]), &
                                        "error message"  )
    
    - The temporary variable in Swiftest is ``rtemp`` for scalar variables and ``vectemp`` for vector variables. The actual values are read in later automatically.
    - For vector variables that vary with ``space``, ``name``, and ``time``, you can define ``start = [1, 1, tslot], count = [NDIM, idmax, 1]``.
        - ``NDIM`` is the pre-defined umber of space dimensions = 3.
        - ``idmax`` is the highest id of the bodies, i.e., the total number of bodies in the data.
    - For scale variables that vary with ``name`` or ``time`` or both, you can define ``start = [1, tslot], count = [idmax, 1]``.

    - Lastly, we also check if the variable is allocated or not and then ``pack`` the values into the parameter for the appropriate bodies using ``plmask`` or ``tpmask``.

    .. code-block:: fortran

        module function swiftest_io_netcdf_read_frame_system(self, nc, param) result(ierr)
            .
            .
            associate(cb => self%cb, pl => self%pl, tp => self%tp, tslot => nc%tslot)
                .
                .
                if (param%lyarkovsky) then
                    call netcdf_io_check( nf90_get_var(nc%id, nc%albedo_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                        "netcdf_io_read_frame_system nf90_getvar albedo_varid"  )
                    if (.not.allocated(pl%albedo)) allocate(pl%albedo(npl))
                    if (npl > 0) pl%albedo(:) = pack(rtemp, plmask)

                    call netcdf_io_check( nf90_get_var(nc%id, nc%emissivity_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                        "netcdf_io_read_frame_system nf90_getvar emissivity_varid"  )
                    if (.not.allocated(pl%emissivity)) allocate(pl%emissivity(npl))
                    if (npl > 0) pl%emissivity(:) = pack(rtemp, plmask)
                .
                .
            end subroutine swiftest_io_read_frame_system

    - We then write out the values to the NetCDF data file in ``swiftest_io_write_frame_body()``:

        - Looping through all the bodies, we put in the necessary data for the ``j``-th using ``nf90_put_var()``. This function is akin to ``nf90_get_var()`` and uses ``start`` and ``count``.
        - For a vector variable, you can define ``start = [1, idslot, tslot], count = [NDIM, 1, 1]``
        - For a scalar variable, you can define ``start = [idslot, tslot]`` (``count`` is not needed).

        - For writing out our the ``albedo`` and ``emissivity`` variables in the Yarkovsky kick for only ``pl`` particles:

        .. code-block:: fortran

            module subroutine swiftest_io_netcdf_write_frame_body(self, nc, param)
                .
                .
                associate(n => self%nbody, tslot => nc%tslot)
                    .
                    .
                    do i = 1, n
                        j = ind(i)
                        .
                        .
                        select type(self)  
                        class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
                            .
                            .
                            if (param%lyarkovsky) then
                                call netcdf_io_check( nf90_put_var(nc%id, nc%albedo_varid, self%albedo(j), start=[idslot, tslot]), &
                                            "netcdf_io_write_frame_body nf90_put_var body albedo_varid"  )
                                call netcdf_io_check( nf90_put_var(nc%id, nc%emissivity_varid, self%emissivity(j), start=[idslot, tslot]), &
                                            "netcdf_io_write_frame_body nf90_put_var body emissivity_varid"  )
                            .
                            .
            end subroutine swiftest_io_netcdf_write_frame_body
    
    - Depending on the properties and dimensions of the new variables, you may have to define it in other subroutines in *swiftest_io.f90* too. Please read the descriptions of each function and use test runs to check.
    - Some potential subroutines include ``swiftest_io_netcdf_read_in_system()``, ``swiftest_io_netcdf_write_frame_cb()``, ``swiftest_io_netcdf_read_hdr_system()``, etc.

Python
=======

Now we will navigate to the Python side of Swiftest. The relevant directory is *swiftest/swiftest/*

- Add the new feature as a valid parameter flag in *io.py*

    - We will add the previously defined name of this to ``newfeaturelist``. Also add it to ``bool_param`` if it is a boolean argument.
    - We defined this as "YARKOVSKY" earlier.

    .. code-block:: python

        newfeaturelist = (
            .
            .
            "YARKOVSKY",
            .
            .
            )
        
        bool_param = [
            .
            .
            "YARKVOSKY",
            .
            .
            ]

- Add and check param flags in *simulation.py*.

- Add variables to *simulation.py*

