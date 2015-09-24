.. _chap-api:

**********
Public API
**********


Open object-oriented architecture
=================================

As well as extensibility through sequences, which was described in :ref:`chap-projects`, in the section :ref:`sec-sequential-integration`, the open architecture of xSPDE allows functional extensions.

This further type of extensibility permits user definable functions to be specified in the ``in`` metadata structures. These functions have default values in xSPDE, and can simply be used as is. It is also possible to have user defined functions that satisfy the interface definitions instead. This is achieved simply by including the relevant function handles and parameters in the input metadata.

This input metadata includes both data and methods acting on the data, in the tradition of object-oriented programs. Yet there is no strict class typing. Users are encouraged to adapt the xSPDE program by adding more input parameters and methods to the input structures. Internal parameters and function handles stored in the :ref:`lattice structure <sec-lattice-structure>` ``r`` are available to all user-defined simulation functions. Note that use of pre-existing reserved names is not advisable.

*Such unorthodox object orientation is deliberate*.

The xSPDE software architecture is intended to be easily extended, and users are encouraged to develop their own libraries. Because this may require new functions and parameters, the internal data architecture is as open as possible.

For example, to define your own integration function, include in the xSPDE input the line:

::

    in.step = @Mystep;

Next, include anywhere on your Matlab path, the function definition, for example:

::

    function a = Mystep(a,xi,dt,r)
        % a = Mystep(a,xi,dt,r) propagates a step my way.
        ...
        a = ...;
    end


Lattice, coordinates and time
=============================

Time and space
--------------

The default lattice for plotted output data is rectangular, with periodic boundary conditions in space, and

::

    dx(i) = in.ranges(i) / (in.points(i) - 1)

The time index is ``1``, and the space index ``i`` ranges from ``2`` to :attr:`in.dimension`. The maximum space-time dimension is ``in.dimension = 4``, while ``in.ranges(i)`` is the time and space duration of the simulation, and ``in.points(i)`` is the total number of plotted points in the ``i``-th direction.

Time is advanced in basic integration steps that are equal to or smaller than ``dx(1)``, for purposes of controlling and reducing errors:

::

    dt = dx(1) / (in.steps * nc)

Here, :attr:`in.steps` is the minimum number of steps used per plotted point, and ``nc = 1, 2`` is the check number. If ``nc = 1``, the run uses coarse time-divisions. If ``nc = 2`` the steps are halved in size for error-checking. Error-checking can be turned off if not required.

Functions
---------

The xSPDE program is function oriented: user specified functions define initial conditions, equations and observables, amongst other things.

Function arguments are always in the following order:

-  the field ``a`` or initial random variable ``v``;
-  the stochastic noise ``z`` or other fields;
-  non-field arguments;
-  the grid structure ``r`` and any previous grid structure needed.

The first argument, ``a``, is a real or complex vector field. This is a matrix whose first dimension is the field index. The second dimension is the lattice index.

The second argument, ``z``, if needed, is a real random noise, corresponding to :math:`\zeta` in the mathematical notation. This is a matrix whose first dimension is the noise index. The second dimension is the lattice index.

The last function argument is the current :ref:`lattice structure <sec-lattice-structure>`, ``r``. This contains data about the current integration lattice. The most important constants are :attr:`r.t`, the current time, and the space coordinates, :attr:`r.x`, :attr:`r.y`, :attr:`r.z`. Other data stored in the lattice structure is explained in later chapters.

Arrays
------

In all function calls, the variables used are matrices. The most important first dimension used is the field length :attr:`in.fields`. The second dimension in all arrays is the lattice index, with a length ``n.lattice = ensembles(1) * points(2) * ... * points(dimension)``. Here ``ensembles(1)`` is the number of stochastic samples integrated as an array.

For reference, the field dimensions are:

- ``a, da, L = [in.fields, r.n.lattice]``;
- ``v = [r.n.random, r.n.lattice]``;
- ``z = [r.n.noise, r.n.lattice]``;
- ``D.x, r.x, r.kx = [1, r.n.lattice]``;
- ``o = [1, r.n.lattice]``.

Each observable is defined by a function in a cell array with length :attr:`in.graphs`.


xSPDE functions
===============

.. function:: xspde(input)

    This is the combined xSPDE function. It accepts a simulation sequence, ``input``. As well as generating graphs, it returns an array ``[e, data]``, where ``e`` are error estimates. It calls the functions :func:`xsim` and :func:`xgraph`.


.. function:: xsim(input)

    This is the xSPDE simulation function. It accepts a simulation sequence in the ``input`` cell array and returns ``[e, data]``, where ``e`` are error estimates, together with a cell array of simulated ``data``. The data are: mean values of functions, error bars and sampling errors. This can be run as a stand-alone function.

.. function:: xgraph(data,input)

    This is the xSPDE graphics function. It takes computed simulation ``data`` and ``input`` parameter specifications. It plots graphs, and returns the maximum difference ``ec`` from comparisons. The ``data`` should have as many cells as ``input`` cells, for sequences. If ``data = ''``, then an HDF5 data file will be read using the file-name specified in ``input``. If ``data = 'filename.h5'`` or ``data = 'filename.mat'`` then the specified file is read as data. Note that ``.h5`` indicates an HDF5 file format, and ``.mat`` indicates a Matlab internal file format.


Simulation parameters
---------------------

For each simulation in the ``input`` sequence, the input and functions are specified as a data structure, ``in``. These can be entered either interactively or as part of a simulation function file. The function file approach allows recycling and editing, so it is better for a large project.

There are extensive default preferences to simplify the inputs. If any inputs are omitted, there are default values which are set by inpreferences in all cases. These defaults are changed by editing the inpreferences function. The :func:`xgrpreferences` function is used to supply graphics default values.

**For vector or cell inputs, an input shorter than required is padded to the right using default values.**


.. _sec-input:

Input parameters and user functions
===================================

A sequence of simulations is obtained from inputs in a cell array, as ``input = {in1, in2, ...}``. The input parameters of each simulation in the sequence are specified in a Matlab structure. If there is one simulation, just one structure can be input, without the braces. This data is also passed to the :func:`xgraph` function. The inputs are numbers, vectors, strings, functions and cell arrays. All xSPDE metadata has preferred values, so only changes from the preferences need to be input. The resulting data is stored internally as a sequence of structures in a cell array, to describe the simulation sequence.

The standard form of each parameter value is:

::

    in.label = parameter

The inputs are scalar or vector parameters or function handles. Quantities relating to graphed averages are cell arrays, indexed by the graph number. The available inputs, with their default values in brackets, are as follows.

Simulation metadata, including all preferred default values that were used in a particular simulation, is also stored for reference in any xSPDE output files. This is done in both the ``.mat`` and the ``.h5`` output files, so the entire simulation can be easily reconstructed or changed.

Note that inputs can be numbers, vectors, strings or cells arrays. To simplify the inputs, some conventions are used, as follows:

- All input data has default values
- Vector inputs of numbers are enclosed in square brackets, ``[...]``.
- Where multiple inputs of strings, functions or vectors are needed they should be enclosed in curly brackets, ``{...}``, to create a cell array.
- Vector or cell array inputs with only one member don’t require brackets.
- Incomplete or partial vector or cell array inputs are filled in with the last applicable default value.
- New function definitions can be just handles pointing elsewhere.


Parameters
----------

.. attribute:: in.name

    *Default:* ``' '``

    Name used to label simulation, usually corresponding to the equation or problem solved. This can be added or removed from graphs using the :attr:`in.headers` Boolean variable, as explained in the section on graphics parameters.

    ::

        in.name = 'your project name'

.. attribute:: in.dimension

    *Default:* ``1``

    The total space-time dimension is labelled, unsurprisingly,

    ::

        in.dimension = 1...4

.. attribute:: in.fields

    *Default:* ``1``

    These are real or complex variables stored at each lattice point, and are the independent variables for integration. The fields are vectors that can have any dimension.

    ::

        in.fields = 1, 2, ...

.. attribute:: in.randoms

    *Default:* :attr:`in.fields`

    This gives the number of random fields generated per lattice point for the initial noise, in coordinate and momentum space. Set to zero (``in.randoms = 0``) for no random fields. Random fields can be correlated either in ordinary or momentum spaces. The second input is the dimension of random fields in momentum space. It can be left out if zero. Note that ``in.randoms = in.randoms(1) + in.randoms(2)``:

    ::

        in.randoms = [in.randoms(1), in.randoms(2)] >= 0

.. attribute:: in.noises

    *Default:* :attr:`in.fields`

    This gives the number of stochastic noises generated per lattice point for both the initial noise and the integration noise, in coordinate and momentum space. Set to zero (``in.noises = 0``) for no noises. This is the number of *rows* in the noise-vector. Noises can be correlated either in ordinary or momentum spaces. The second input is the dimension of noises in k-space. It can be left out if zero. Note that ``in.noise = noises(1) + noises(2)``:

    ::

        in.noises = [in.noises(1), in.noises(2)] >= 0.

.. attribute:: in.ranges

    *Default:* ``[10, 10, ...]``

    Each lattice dimension has a coordinate range, given by:

    ::

        in.ranges = [in.ranges(1), ..., in.ranges(dimension)]

    In the temporal graphs, the first coordinate is plotted over ``0:in.ranges(1)``. All other coordinates are plotted over ``-in.ranges(n)/2:in.ranges(n)/2``. The default value is ``10`` in each dimension.

.. attribute:: in.points

    *Default:* ``[49, 35, ..., 35]``

    The rectangular lattice of points plotted for each dimension are defined by a vector giving the number of points in each dimension:

    ::

        in.points = [in.points(1), ..., in.points(in.dimension)]

    The default values are simply given as a rough guide for initial calculations. Large, high dimensional lattices take more time to integrate. Increasing :attr:`in.points` improves graphics resolution, and gives better accuracy in each relevant dimension as well, but requires more memory. Speed is improved when the lattice points are a product of small prime factors.

.. attribute:: in.steps

    *Default:* ``1``

    Number of time-steps per plotted point. The total number of integration steps in a simulation is therefore ``in.steps * (in.points(1)-1)``. Thus, :attr:`in.steps` can be increased to improve the accuracy, but gives no change in graphics resolution. **Increase** steps to give a **lower** time-discretization error:

    ::

        in.steps = 1, 2, ...

.. attribute:: in.ensembles

    *Default:* ``[1, 1, 1]``

    Number of independent stochastic trajectories simulated. This is specified in three levels to allow maximum parallelism. The first gives within-thread parallelism, allowing vector instructions. The second gives a number of independent trajectories calculated serially. The third gives multi-core parallelism, and requires the Matlab parallel toolbox. Either ``in.ensembles(2)`` or ``in.ensembles(3)`` are required if sampling error-bars are to be calculated.

    ::

        in.ensembles = [in.ensembles(1), in.ensembles(2), in.ensembles(3)] >= 1

    The *total* number of stochastic trajectories or samples is ``ensembles(1) * ensembles(2) * ensembles(3)``.

.. attribute:: in.transforms

    *Default:* ``{0}``

    **Cell array** that defines the different transform spaces used to calculate field observables. This has the structure

    ::

        in.transforms{n} = [t(1), ..., t(4)] >= 0

    There is one transform vector per observable. The ``j``-th index, ``t(j)``, indicates a Fourier transform on the ``j``-th axis. The normalization of the Fourier transform is such that the :math:`k=0` value in momentum space corresponds to the integral over space, with an additional factor of :math:`1/\sqrt{2\pi}`. This gives a Fourier integral which is symmetrically normalized in ordinary and momentum space. The Fourier transform is such that
    :math:`k=0` is the *central* value.

.. attribute:: in.olabels

    *Default:* ``{'a_1', ...}``

    **Cell array** of labels for the graph axis observable functions. These are text labels that are used on the graph axes. The default value is ``'a_1'`` if the default observable is used, otherwise it is blank. This is overwritten by any subsequent label input when the graphics program is run:

    ::

        in.olabels{n} = 'string'

.. attribute:: in.c

    This starting letter is always reserved to store user-specified constants and parameters. All inputs --- including ``c`` data --- are copied into the data files and also the lattice structure ``r``. It is passed to user functions, and can be any data.

    ::

        in.c = anything


Invariant inputs
----------------

The following can’t be changed during a sequence in the current xSPDE version --- the specified values for the first simulation will be used:

#. The extrapolation order

#. The number of ensembles (2)

#. The number of ensembles (3)


Input functions
---------------

A stochastic equation solver requires the definition of an initial distribution and a time derivative. In xSPDE, the time derivatives is divided up into a linear term including space derivatives, used to define an interaction picture, and the remaining derivatives. In addition, one must define quantities to be averaged over during the simulation, called graphs in xSPDE. These are all defined as functions, specified below.

.. attribute:: in.initial

    *Default:* :func:`xinitial`

    Initializes the fields :math:`a` for the first simulation in a sequence. The initial Gaussian random field variable, ``v``, has unit variance if :attr:`in.dimension` is ``1`` or else is delta-correlated in space, with variance ``1/r.dV`` (:math:`\equiv 1/(dx_2...dx_d)`) for :math:`d` space-time dimensions. If specified in the input, ``ra`` has a first dimension of :attr:`r.n.random`, otherwise the default is :attr:`in.fields`. The default set by :func:`xinitial` is ``a = 0``.

.. attribute:: in.transfer

    *Default:* :func:`xtransfer`

    Initializes the fields :math:`a` for subsequent calculations in a sequence. Otherwise, this function behaves in a similar way to :attr:`in.initial`. The function includes the previous field ``a0`` and lattice ``r0``. The default set by :func:`xtransfer` is ``a = a0``.

.. attribute:: in.da

    *Default:* :func:`xda`

    Calculates derivatives :math:`da` of the equation. The noise vector, ``z``, has variance :math:`1/(dx_{1}..dx_{d})`, for dimension :math:`d \le 4`, and a first dimension of :attr:`r.n.noise` whose default value is :attr:`in.fields`. If specified by the two elements of the :attr:`in.noises` vector, ``z`` can have a different first dimension from :attr:`in.fields`. This can also include noise correlated in momentum space.

.. attribute:: in.linear

    *Default:* :func:`xlinear`

    A user-definable function which returns the linear coefficients :math:`L` in Fourier space. This is a function of the differential operator ``D``. The default is zero. Here ``D`` is a structure with components ``D.x``, ``D.y``, ``D.z``. Each component has an array dimension the same as the coordinate lattice.

.. attribute:: in.observe

    *Default:* cell array of :func:`xobserve`

    **Cell array** of function handles that take the current field and returns a real observable ``o`` with dimension of ``[1, n.lattice]``. The default observable is the first real field amplitude. Note the use of braces for cell arrays! One can also input these individually as ``in.observe{1} = @(a,r) f(a,r)``, using an inline anonymous function. The total number of observe functions is stored internally as :attr:`in.graphs`. The fields ``a`` passed in the input are transformed according to the :attr:`in.transforms` metadata.

.. attribute:: in.rfilter(r)

    *Default:* :func:`xrfilter`

    Returns the momentum-space filters for the input random terms. Each component has an array dimension the same as the coordinate lattice, that is, the return dimension is ``[r.randoms(2), r.n.lattice]``.

.. attribute:: in.nfilter(r)

    *Default:* :func:`xnfilter`

    Returns the momentum-space filters for the propagation noise terms. Each component has an array dimension the same as the coordinate lattice, that is, the return dimension is ``[r.noises(2), r.n.lattice]``.


Advanced input parameters
-------------------------

More advanced input parameters, which don’t usually need to be changed from default values, are as follows:

.. attribute:: in.iterations

    *Default:* ``4``

    For iterative algorithms like the implicit midpoint method, the iteration count is set here, typically around 3-4. Will increase the integration accuracy if set higher, but it may be better to increase :attr:`in.steps` if this is needed. With non-iterated algorithms, this input is not used:

    ::

        in.iterations = 1, 2, ...

.. attribute:: in.errorchecks

    *Default:* ``2``

    This defines how many times the integration is carried out for error-checking purposes. If :attr:`in.errorchecks` is `1`, there is one integration, but no checking at smaller time-steps. For error checking, set ``in.errorchecks = 2``, which repeats the calculation at a shorter time-step --- but with identical noise --- to obtain the error bars, taking three times longer overall:

    ::

        in.errorchecks = 1, 2

.. attribute:: in.order

    *Default:* ``1``

    This is the extrapolation order, which is **only** used if ``in.errorchecks = 2``. The program uses the estimated convergence order to extrapolate to zero step-size, with reduced estimated error-bars. If ``in.order = 0``, no extrapolation is used, which is the most conservative input. The default order is usually acceptable, especially when combined with the default midpoint algorithm, see next section. While any non-negative order can be input, the theoretical orders of the four preset methods used *without* stochastic noise terms are: ``1`` for :func:`xEuler`; ``2`` for :func:`xRK2`; ``2`` for :func:`xMP`; ``4`` for :func:`xRK4`. Allowed values are:

    ::

        in.order >= 0

.. attribute:: in.seed

    *Default:* ``0``

    Random noise generation seed, for obtaining reproducible noise sequences. Only needed if ``in.noises > 0``

    ::

        in.seed >= 0

.. attribute:: in.graphs

    *Default:* number of observables

    This gives the number of observables or graphs computed. The default is the length of the cell array of observable functions. Normally, this is not initialized, as the default is typically used. Can be used to suppress data averaging.

    ::

        in.graphs >= 0

.. attribute:: in.print

    *Default:* ``1``

    Print flag for output information while running xSPDE. If ``print = 0``, most output is suppressed, while ``print = 1`` displays a progress report, and ``print = 2`` also generates a readable summary of the ``r`` lattice structure as a record.

    ::

        in.print >= 0

.. attribute:: in.raw

    *Default:* ``0``

    Flag for storing raw trajectory data. If this flag is turned on, raw trajectories are stored in memory and written to a file on completion. To make use of these, a file-name should be included!

    ::

        in.raw >= 0

.. attribute:: in.origin

    *Default:* ``[0, -in.ranges/2]``

    This displaces the graph origin for each simulation to a user-defined value. If omitted, all initial times in a sequence are zero, and the space origin is set to ``-in.ranges/2`` to give results that are symmetric about the origin:

    ::

        in.origin = [origin(1), ..., origin(4)]

.. attribute:: in.ipsteps

    *Default:* ``1`` for :func:`xEuler` and :func:`xRK2`, ``2`` for :func:`xMP` and :func:`xRK4`

    This specifies the number of interaction picture steps needed in a full propagation time-step. Default values are chosen according to the setting of :attr:`in.step`. Can be changed for custom integration methods.

    ::

        in.ipsteps = 1, 2, 3

.. attribute:: in.file

    *Default:* ``''``

    Matlab or *HDF5* file name for output data. Includes all data and parameter values, including raw trajectories if ``in.raw = 1``. If not needed just omit this. A Matlab filename should end in ``.mat``, while an HDF5 file requires the filename to end in ``.h5``.

    ::

        in.file = [origin(1), ..., origin(4)]

Advanced input functions
------------------------

Advanced input functions are user-definable functions which don’t usually need to be changed from default values. They allow customization and extension of xSPDE. These are as follows:

.. attribute:: in.grid(r)

    *Default:* :func:`xgrid`

    Initializes the grid of coordinates in space.

.. attribute:: in.noisegen

    *Default:* :func:`xnoisegen`

    Generates arrays of noise terms ``xi`` for each point in time.

.. attribute:: in.randomgen

    *Default:* :func:`xrandomgen`

    Generates a set of random fields ``rf`` to initialize the fields simulated.

.. attribute:: in.step

    *Default:* :func:`xMP`

    Specifies the stochastic integration routine for a step in time ``dt`` and noise ``xi``. It returns the new field ``a`` at space-time location ``r``, given the old field as input, and interaction-picture propagator :attr:`r.propagator` which is part of the lattice structure. This can be set to any of the predefined stochastic integration routines provided with xSPDE, described in the :ref:`chap-algorithms` chapter. User-written functions can also be used. The standard method, :func:`xMP`, is a midpoint integrator.

.. attribute:: in.prop

    *Default:* :func:`xprop`

    Returns the fields propagated in the interaction picture, depending on the propagator array :attr:`r.propagator`.

.. attribute:: in.propfactor

    *Default:* :func:`xpropfactor`

    Returns the transfer array :attr:`r.propagator`, used by the :attr:`in.prop` function. The time propagated is a fraction of the integration time-step, :attr:`r.dt`. It is equal to ``1 / in.ipsteps`` of the integration time-step.


Graphics inputs and functions
-----------------------------

The graphics parameters are also stored in the cell array ``input`` as a sequence of structures ``in``. This only need to be input when the graphs are generated, and can be changed at a later time to alter the graphics output. A sequence of simulations is graphed from ``input`` specifications.

If there is one simulation, just one structure can be input, without the sequence braces. The standard form of each parameter value, which should have the ``in.`` structure label added, is:

::

    in.label = parameter

If any inputs are omitted, there are default values which are set by the :func:`xgrpreferences` function, in all cases except for the comparison function :func:`in.compare`. The defaults can be changed by editing the :func:`xgrpreferences` function.

In the following descriptions, :attr:`in.graphs` is the total number of graphed variables of all types. The space coordinate, image, image-type and transverse data can be omitted if there is no spatial lattice, that is, if the dimension variable is set to one.

Graphics functions
~~~~~~~~~~~~~~~~~~

.. function:: in.compare(t,in)

    This is a cell array of functions. Each takes the time or frequency vector and returns comparison results for a graphed observable, as a function of real values versus time or frequency. Comparison results are graphed with a dashed line, for the two-dimensional graphs versus time. There is no default function handle.

Graphics parameters
~~~~~~~~~~~~~~~~~~~

For uniformity, the graphics parameters are cell arrays, indexed over the graph number using braces ``{}``. If a different type of input is used, like a scalar or matrix, xSPDE will attempt to convert the type. The axis labels are cell arrays, indexed over dimension.

Together with default values, they are:

.. attribute:: in.font

    *Default:* ``{18, ...}``

    This sets the default font size for the graph labels.

    ::

        in.font{n} > 0

.. attribute:: in.minbar

    *Default:* ``{0.01, ...}``

    This is the minimum relative error-bar that is plotted.

    ::

        in.minbar{n} >= 0

.. attribute:: in.images

    *Default:* ``{0, 0, 0, ...}``

    This is the number of 3D, transverse o-x-y images plotted as discrete time slices. Only valid if :attr:`in.dimension` is greater than 2. Note that, if present, the z-coordinate is set to its central value of ``z = 0``, when plotting the transverse images. This input should be from ``in.images(n) = 0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`in.graphs`:

    ::

        in.images{n} = 0 ... in.points(1)

.. attribute:: in.imagetype

    *Default:* ``{1, 1, ...}``

    This is the *type* of transverse image plotted. If an element is ``1``, a perspective surface plot is output, for ``2``, a gray plot with colours is output, or for ``3`` a contour plot with 10 equally spaced contours is generated. This has a vector length equal to :attr:`in.graphs`.

    ::

        in.imagetype{n} = 1, 2, 3

.. attribute:: in.transverse

    *Default:* ``{0, 0, ...}``

    This is the number of 2D, transverse o-x images plotted as discrete time slices. Only valid if :attr:`in.dimension` is greater than 2. Note that, if present, the y,z-coordinates are set to their central values, when plotting the transverse images. Each element should be from ``0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`in.graphs`:

    ::

        in.transverse{n}=0 ... in.points(1)

.. attribute:: in.headers

    *Default:* ``{1, 1, ...}``

    This is a Boolean variable with value ``true`` or ``1`` if graphs require headers giving the simulation name, and ``false`` or ``0`` with no headers. It is useful to include headings on graphs in preliminary stages, while they may not be needed in a published final result.

    ::

        in.headers{n} = 0, 1

.. attribute:: in.pdimension

    *Default:* ``{4, 4, ...}``

    This is the maximum plot dimension for each graphed quantity. The purpose is eliminate unwanted graphs. For example, it is useful to reduce the maximum dimension when averaging in space. Higher dimensional graphs are not needed, as the data is duplicated. Averaging can be useful for checking conservation laws, or for averaging over homogeneous data to reduce sampling errors.

    ::

        in.pdimension{n} = 1 ... 4

.. attribute:: in.xlabels

    *Default:* ``{'t', 'x', 'y', 'z'}``

    Labels for the graph axis independent variable labels, vector length of :attr:`in.dimension`. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.xlabels = {in.xlabels(1), ..., in.xlabels(in.dimension)}

.. attribute:: in.klabels

    *Default:* ``{'\\omega', 'k\_x', 'k\_y', 'k\_z'}``

    Labels for the graph axis Fourier transform labels, vector length of :attr:`in.dimension`. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.klabels = {in.klabels(1), ..., in.klabels(in.dimension)}

Graphics projections
~~~~~~~~~~~~~~~~~~~~

If there is a spatial lattice, the graphics program automatically generates several graphs for each observable, depending on space dimension. The maximum dimension that is plotted as set by :attr:`in.pdimension`. In the plots, the lattice is projected down to successively lower dimensions.

For each observable, the projection sequence is as follows:

-  If :attr:`in.dimension` is ``4``, a central :math:`z` point ``nz = 1 + floor(in.points(4)/2)`` is picked. For example, with 35 points, this gives the central point, ``nz = 18``.

-  This gives a three dimensional space-time lattice, which is treated the same as if :attr:`in.dimension` is ``3``.

-  If :attr:`in.images` are specified, two-dimensional :math:`x-y` plots are generated at equally spaced time intervals. If there is only one image, it is at the last time-point. Different plot-types are used depending on the setting of :attr:`in.imagetype`.

-  A central :math:`y` point ``ny = 1 + floor(in.points(3)/2)`` is picked. This gives a two dimensional space-time lattice, which is treated the same as if :attr:`in.dimension` is ``2``. If :attr:`in.transverse` is specified, one-dimensional :math:`x` plots are generated at equally spaced time intervals, as before.

-  A central :math:`x` point ``nx = 1 + floor(in.points(2)/2)`` is picked. This gives a one dimensional time lattice, which is treated the same as if :attr:`in.dimension` is ``1``.

-  Plots of observable vs time are obtained, including sampling errors and error bars. If comparison graphs are specified using :func:`in.compare` functions, they are plotted also, using a dotted line. A difference graph is also plotted when there is a comparison.


Averages and integrals
======================

Averages
--------

Lattice averages can allow one to extract stochastic results with reduced sampling errors. An average over the lattice is carried out using the :func:`xave` function, which is defined as follows:

.. function:: xave(o, [dx, r])

    This function takes a scalar observable ``o = [1, lattice]``, defined on the xSPDE lattice, and returns a space average with dimension ``[1, lattice]``. The input is an observable ``o``, and an optional lattices structure and vector switch ``dx``. If ``dx(j) > 0``, an average is taken over dimension ``j``. Dimensions are labelled from ``j = 1 ... 4`` as elsewhere. Time averages are ignored at present. Averages are returned at all lattice locations. To average over samples and all space dimensions, just use ``xave(o)``.

Higher dimensional graphs of lattice averages are generally not useful, as they are simply flat. The xSPDE program allows the user to remove unwanted higher dimensional graphs of average variables. This is achieved by setting the corresponding element of :attr:`in.pdimension` to the highest dimension required, which of course depends on which dimensions are averaged.

For example, to average over the entire space lattice and indicate that only time-dependent graphs are required, set ``dx = in.dx`` and:

::

    in.pdimension = 1

Note that :func:`xave` does not perform any average over ensembles, although this is done elsewhere for results calculated in any of the :attr:`in.observe` functions.

Integrals
---------

Integrals over the spatial lattice allow calculation of conserved or other global quantities. The code to take an integral over the lattice is carried out using the xSPDE :func:`xint` function:

.. attribute:: xint(o, dx, r)

    This function takes a scalar ``o``, and returns a space integral over selected dimensions with vector measure ``dx``. If ``dx(j) > 0`` an integral is taken over dimension ``j``. Dimensions are labelled from ``j = 1, ..., 4`` as in all xspde standards. Time integrals are ignored at present. Integrals are returned at all lattice locations. To integrate over an entire lattice, set ``dx = r.dx``, otherwise set ``dx(j) = r.dx(j)`` for selected dimensions ``j``.

As with averages, the xSPDE program allows the user to remove unwanted higher dimensional graphs when the integrated variable is used as an observable. For example, in a four dimensional simulation with integrals taken over the :math:`y` and :math:`z` coordinates, only :math:`t`- and :math:`x`-dependent graphs are required. Hence, set ``dx`` to ``[0, 0, r.dx(3), r.dx(4)]``, and:

::

    in.pdimension = 2

If momentum-space integrals are needed, use the transform switch to make sure that the field is Fourier transformed, and input :attr:`r.dk` instead of :attr:`r.dx`. Note that :func:`xint` returns a lattice observable, as required when used in the :attr:`in.observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.


.. _sec-lattice-structure:

Lattice structure
=================

Internally, xSPDE data is stored in a cell array, ``latt``, of structures ``r``, which is passed to functions. This includes all the data given above inside the ``in`` structure. In adition, it includes the table of computed parameters given below.

User application constants and parameters should not be reserved names; :attr:`in.c` and all names starting with ``in.c`` will always be available in all versions of xSPDE.

A lattice structure contains information about space-time grid and is passed to various functions, for instance :attr:`in.da` or :attr:`in.step`. The corresponding parameter is commonly marked as `r`.

.. attribute:: r.t

    Current value of time, :math:`t`.

.. attribute:: r.x

.. attribute:: r.y

.. attribute:: r.z

    Coordinate grids of :math:`x`, :math:`y`, :math:`z`.

.. attribute:: r.kx

.. attribute:: r.ky

.. attribute:: r.kz

    Grids in momentum space: :math:`k_x`, :math:`k_y`, :math:`k_z`.

.. attribute:: r.dt

    Output time-step.

.. attribute:: r.dx

    Steps in coordinate space: :math:`[t,x,y,z]`.

.. attribute:: r.dk

    Steps in momentum space: :math:`[\omega,k_{x},k_{y},k_{z}]`.

.. attribute:: r.propagator

    Contains the propagator array for the interaction picture.

.. attribute:: r.V

    Spatial lattice volume.

.. attribute:: r.K

    Momentum lattice volume.

.. attribute:: r.dV

    Spatial cell volume.

.. attribute:: r.dK

    Momentum cell volume.

.. attribute:: r.xc

    Space-time coordinate axes (vector cells).

.. attribute:: r.kc

    Computational axes in :math:`[\omega,k_{x},k_{y},k_{z}]` (vector cells).

.. attribute:: r.gk

    Graphics axes in :math:`[\omega,k_{x},k_{y},k_{z}]` (vector cells).

.. attribute:: r.kr

    Range in :math:`[\omega,k_{x},k_{y},k_{z}]` (vector).

.. attribute:: wtph

    Frequency phase-factors (vector).

.. attribute:: r.s.dx

    Initial stochastic normalization.

.. attribute:: r.s.dxt

    Propagating stochastic normalization.

.. attribute:: r.s.dk

    Initial :math:`k` stochastic normalization.

.. attribute:: r.s.dkt

    Propagating :math:`k` stochastic normalization.

.. attribute:: r.n.space

    Number of spatial lattice points.

.. attribute:: r.n.lattice

    Total lattice: ``in.ensembles(1) * r.n.space``.

.. attribute:: r.n.ensemble

    ``in.ensembles(2) * in.ensembles(3)``.

.. attribute:: r.n.random

    Number of initial random fields.

.. attribute:: r.n.noise

    Number of noise fields.

.. attribute:: r.d.int

    Dimensions for lattice integration (vector).

.. attribute:: r.d.a

    Dimensions for :math:`a` field (flattened, vector).

.. attribute:: r.d.r

    Dimensions for coordinates (flattened, vector).

.. attribute:: r.d.ft

    Dimensions for field transforms (vector).

.. attribute:: r.d.k

    Dimensions for noise transforms (vector).

.. attribute:: r.d.obs

    Dimensions for observations (vector).

.. attribute:: r.d.data

    Dimensions for average data (flattened, vector).

.. attribute:: r.d.raw

    Dimensions for raw data (flattened, vector).


Default functions
=================

These functions are used as defaults for simulations and can be overridden by the user.

.. function:: xinitial(~, r)

    Returns a field array filled with zeros.

.. function:: xtransfer(~, ~, a, ~)

    Returns the field ``a`` unchanged.

.. function:: xda(~, ~, r)

    Returns a derivative array filled with zeros.

.. function:: xlinear(~, r)

    Returns a linear response array filled with zeros.

.. function:: xobserve(a, ~)

    Returns the real part of ``a``.

.. function:: xrfilter(r)

    Returns an array of ones.

.. function:: xnfilter(r)

    Returns an array of ones.

.. function:: xgrid(r)

    Sets grid points in lattice from coordinate vectors. Returns the ``r`` structure with added grid points.

.. function:: xnoisegen(r)

    Generates random noise matrix :math:`\xi`.

.. function:: xrandomgen(r)

    Generates random field matrix :math:`w`.

.. function:: xpropfactor(nc, r)

    Returns the interaction picture propagation factor. ``nc`` is a check index, ``r`` is a lattice structure.


Frequently asked questions
==========================

Answers to some frequent questions, and reminders of points in this chapter are:

-  Can you average other stochastic quantities apart from the field?

   -  Yes: just specify this using the user function :attr:`in.observe`.

-  Can you have functions of the current time and space coordinate?

   -  Yes: xSPDE functions support this using the structure ``r``, as :attr:`r.t`, :attr:`r.x`, :attr:`r.y`, :attr:`r.z`.

-  Can you have several variables?

   -  Yes, input this using ``in.fields > 1``.

-  Are higher dimensional differential equations possible?

   -  Yes, this requires setting ``in.dimension > 1``.

-  Can you have spatial partial derivatives?

   -  Yes, provided they are linear in the fields; these are obtainable using the function :attr:`in.linear`.

-  Can you delete the graph heading?

   -  Yes, this is turned off if you set :attr:`in.headers` to ``0``.

-  Why are there two lines in the graphs sometimes?

   -  These are the one standard deviation sampling error limits, generated when ``in.ensembles(2,3) > 1``.

-  Why is there just one line in some graphs, with no sampling errors indicated?

   -  You need ``in.ensembles(2)`` or ``(3)`` for this; see previous question.

-  What are the error bars for?

   -  These are the estimated maximum errors due to finite step-sizes in time.
