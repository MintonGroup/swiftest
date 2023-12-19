.. currentmodule:: swiftest

.. _api:

#########################
Swiftest API reference
#########################

This section of the documentation provides a detailed reference for the Production classes in the Swiftest project.

Simulation
==========

The Simulation class is the main class for the swiftest project. It is used to create a simulation of a crater on a given target 
body. The Simulation class is used to generate craters of a given size and morphology based on the production function, morphology
function, and crater scaling relationship model. The surface of the target body is represented by a Surface attribute called
`surf`, which is derived from a UxDataset object. This is an unstructured grid dataset that contains data for the target body surface.

Creating a Simulation
---------------------

.. autosummary::
    :toctree: generated/

    Simulation

Running a Simulation
-------

.. autosummary::
    :toctree: generated/

    Simulation.run

