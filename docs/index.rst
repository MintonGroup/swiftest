.. module:: swiftest

Swiftest
========

Swiftest is a software packaged designed to model the dynamical evolution of gravitational systems. Swiftest is a re-write of the
`Swifter <https://www.boulder.swri.edu/swifter/>`__ software package that incorporates modern programming techniques and performance
improvements. 
Swiftest contains the following numerical integrators:

-  **Wisdom-Holman Mapping (WHM)** - A symplectic n-body mapping method.
   See `Wisdom & Holman (1991) <https://ui.adsabs.harvard.edu/abs/1991AJ....102.1528W/abstract>`__.
-  **Regularized Mixed Variable Symplectic (RMVS)** - An extension of WHM that is capable of handling close encounters between test
   particles and massive bodies. See `Levison & Duncan (1994) <https://www.sciencedirect.com/science/article/pii/S0019103584710396?via%3Dihub>`__.
-  **Democratic Heliocentric (HELIO)** - A symplectic integrator that uses the democratic heliocentric coordinate frame. See 
-  `Duncan, Levison, & Lee (1998) <https://iopscience.iop.org/article/10.1086/300541>`__.
-  **Symplectic Massive Body Algorithm (SyMBA)** - An extension of HELIO that is capable of handling close encounters between massive bodies.
   See `Duncan, Levison, & Lee (1998) <https://iopscience.iop.org/article/10.1086/300541>`__.

Swiftest also includes the collisional fragmentation algorithm **Fraggle**, an addition to the SyMBA integrator. Fraggle is designed to
resolve collisions between massive bodies by determining the collisional regime, derived from the work of `Leinhardt & Stewart
(2012) <https://iopscience.iop.org/article/10.1088/0004-637X/745/1/79>`__, and generating the appropriate mass distribution of fragments. Swiftest
fully incorporates collisional fragments into the gravitational system, evolving these new bodies along with pre-existing bodies, including
their growth and any future fragmentation events in which they are involved.


**Useful links**:
`Home <https://swiftest.readthedocs.io/en/latest/>`__ |
`Code Repository <https://github.com/profminton/swiftest>`__ |
`Issues <https://github.com/profminton/swiftest/issues>`__ |
`Discussions <https://github.com/profminton/swiftest/discussions>`__ |
`Releases <https://github.com/profminton/swiftest/releases>`__ |


.. grid:: 1 1 2 2
    :gutter: 2

    .. grid-item-card:: Getting started
        :img-top: _static/index_getting_started.svg
        :link: getting-started-guide/index
        :link-type: doc

        New to *swiftest*? Check out the getting started guides. They contain an
        introduction to *Swiftest's* main concepts and links to additional tutorials.

    .. grid-item-card::  User guide
        :img-top: _static/index_user_guide.svg
        :link: user-guide/index
        :link-type: doc

        The user guide provides in-depth information on the
        key concepts of Swiftest with useful background information and explanation.

    .. grid-item-card::  API reference
        :img-top: _static/index_api_reference.svg
        :link: api
        :link-type: doc

        The reference guide contains a detailed description of the Swiftest API.
        The reference describes how the methods work and which parameters can
        be used. It assumes that you have an understanding of the key concepts.

    .. grid-item-card::  Developer guide
        :img-top: _static/index_contribute.svg
        :link: contributing
        :link-type: doc

        Saw a typo in the documentation? Want to improve existing functionalities?
        The contributing guidelines will guide you through the process of improving
        Swiftest.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For users

   Getting Started <getting-started-guide/index>
   User Guide <user-guide/index>
   API Reference <api>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For developers/contributors

   Contributing Guide <contributing>
   GitHub repository <https://github.com/profminton/swiftest>

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Community

   GitHub discussions <https://github.com/profminton/swiftest/discussions>
