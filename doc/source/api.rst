.. _chap-api:

**********
Public API
**********.. _chap-api:

**********
Public API
**********


High-level xSPDE functions and objects
======================================

The high-level xSPDE functions process the input parameters, producing simulation data and graphs. Function parameters are in the order of ``input``, which specifies the simulation, then ``data``, which is needed for graphics output.

.. function:: xspde(input)

    This is the combined xSPDE function. It accepts a simulation sequence, ``input``. This can be a single structure, ``in``, or else a cell array of structures, ``{in1,in2,..}``, for  sequences. Output graphs are displayed. It returns the output ``[error, input, data,raw]``, where ``error`` is the sum of simulation errors in :func:`xsim`, and difference errors found in the :func:`xgraph` comparisons. If a filename is specified in  the input, it writes an output data file including input and all output data. Raw data is stored on request. It calls the functions :func:`xsim` and :func:`xgraph`.


.. function:: xsim(input)

    This is the xSPDE simulation function. Like :func:`xspde`, it accepts input parameters in ``input``. It returns ``[maxerror, input, data, raw]``, where: ``maxerror`` is the sum of maximum step-size and maximum sampling errors, ``input`` is the full input structure or cell array for sequences, including default values, and ``data`` is a cell array of average observables. If the ``in.raw`` option is used, data for the actual trajectories is output in ``raw``. This can be run as a stand-alone function if no graphs are required.

.. function:: xgraph(data [,input])

    This is the xSPDE graphics function. It takes computed simulation  ``data`` and ``input``. It plots graphs, and returns the maximum difference ``diff`` from comparisons with user-specified comparison functions. The ``data`` should have as many cells as ``input`` cells, for sequences. 
    If ``data = 'filename.h5'`` or ``data= 'filename.mat'``, the specified file is read both for ``input`` and ``data``. Here ``.h5`` indicates an HDF5 file, and ``.mat`` indicates a Matlab file.
    When the``data`` input is given as a filename, input parameters in the file are replaced by any of the the new ``input`` parameters that are specified.  Any stored ``input`` can be overwritten, allowing graphs to be modified retrospectively.
    If the ``input`` is not present, a file-name must be specified in the ``data`` field. 



Open object-oriented architecture
----------------------------------

As well as extensibility through sequences, which was described in :ref:`chap-projects`, in the section :ref:`sec-sequential-integration`, the architecture of xSPDE allows functional extensions.

The input metadata includes both data and methods acting on the data, in the tradition of object-oriented programs. Yet there is no strict class typing. Users are encouraged to adapt the xSPDE program by adding input parameters and methods to the input structures.

*Such unorthodox object orientation is deliberate*.

This open extensibility permits arbitrary functions to be specified in the ``in`` structures. All functions and parameters have default values in xSPDE. It is also possible to include user defined functions, provided they satisfy the API definitions given below. This is achieved simply by including the relevant function handles and parameters in the input metadata.

Internal parameters and function handles from the ``in`` structures, together with any computed and default parameters and functions,  are stored in the :ref:`lattice structure <sec-lattice-structure>` ``r``. These are available to all user-defined functions. Use of pre-existing reserved names is not advisable, and the structure ``r.c`` is always reserved for user constants, parameters and functions if required.

The xSPDE software architecture is intended to be easily extended, and users are strongly encouraged to develop their own libraries and contribute to the xSPDE function pool. Because this generally requires new functions and parameters, the internal data architecture is as open as possible.

For example, to define your own integration function, include in the xSPDE input the line:

::

    in.step = @Mystep;

Next, include anywhere on your Matlab path, the function definition, for example:

::

    function a = Mystep(a,z,dt,r)
        % a = Mystep(a,z,dt,r) propagates a step my way.
        ...
        a = ...;
    end





Averages, integrals and derivatives
===================================

There are functions available in xSPDE for averages, integrals and derivatives. These can be used to calculate observables for plotting, but are also available for calculating stochastic derivatives as part of the stochastic equation. They operate in parallel over the local ensemble and lattice dimensions. They take a scalar quantity, for example a single field component, and return an average, a space integral, and a spatial derivative respectively. In each case the first argument is the field, the second argument is a vector defining the type of operation, and the last argument is the parameter structure, ``r``. If there are only two arguments, the operation vector is replaced by its default value.

Averages
--------

This function allows one to extract local ensemble averages which may be needed as part of an overall calculation involving a ratio or product of averages. In addition, spatial grid averages can be used to obtain stochastic results with reduced sampling errors if the overall grid is homogeneous.

An average is carried out using the :func:`xave` function, which is defined as follows:

.. function:: xave(o, [av, ] r)

    This function takes a scalar field or observable ``o = [1, n.lattice]``, defined on the xSPDE local lattice, and returns an average with dimension ``[1, n.lattice]``. The input is a field or observable ``o``, and an optional averaging switch ``av``. If ``av(j) > 0``, an average is taken over dimension ``j``. Dimensions are labelled from ``j = 1 ... 4`` as elsewhere. The first index indicates a local ensemble average, while subsequent indices indicate averages over the spatial grid. If the ``av`` vector is omitted, the average is only taken over the local ensemble. Averages are returned at all lattice locations. To average over both the local ensemble and all space dimensions, just use ``xave(o)``.

Higher dimensional graphs of grid averages are generally not useful, as they are simply flat. The xSPDE program allows the user to remove unwanted higher dimensional graphs of average variables. This is achieved by setting the corresponding element of :attr:`in.pdimension` to the highest dimension required, which depends on which dimensions are averaged.

For example, to average over the entire ensemble plus space lattice and indicate that only time-dependent graphs are required, set ``av = in.dx`` and:

::

    in.pdimension = 1

Note that :func:`xave` on its own gives identical results to those calculated in the :attr:`in.observe` functions. Its utility comes when more complex combinations or functions of ensemble averages are required.

Integrals
---------

Integrals over the spatial grid allow calculation of conserved or other global quantities. To take an integral over the spatial grid,  use the xSPDE :func:`xint` function:

.. attribute:: xint(o, [dx, ] r)

    This function takes a scalar ``o``, and returns a space integral over selected dimensions with vector measure ``dx``. If ``dx(j) > 0`` an integral is taken over dimension ``j``. Dimensions are labelled from ``j = 1, ...`` as in all xspde standards. Time integrals are ignored at present. Integrals are returned at all lattice locations. To integrate over an entire lattice, set ``dx = r.dx``, otherwise set ``dx(j) = r.dx(j)`` for selected dimensions ``j``.

As with averages, the xSPDE program allows the user to remove unwanted higher dimensional graphs when the integrated variable is used as an observable. For example, in a four dimensional simulation with integrals taken over the :math:`y` and :math:`z` coordinates, only :math:`t`- and :math:`x`-dependent graphs are required. Hence, set ``dx`` to ``[0, 0, r.dx(3), r.dx(4)]``, and:

::

    in.pdimension = 2

If momentum-space integrals are needed, use the transform switch to make sure that the field is Fourier transformed, and input :attr:`r.dk` instead of :attr:`r.dx`. Note that :func:`xint` returns a lattice observable, as required when used in the :attr:`in.observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.


Derivatives
-----------

The code to take a spatial derivative is carried out using the xSPDE :func:`xd` function:

.. attribute:: xd(o, [D, ] r)

This function takes a scalar ``o``, and returns a derivative over selected dimensions with a derivative ``D``.  Set ``D = r.Dx`` for a first order x-derivative, ``D = r.Dy`` for a first order y-derivative, and similarly ``D = r.Dz.*r.Dy`` for a cross-derivative in ``z`` and ``y``. Higher derivatives require powers of these. For higher dimensions use numerical labels, where ``D = r.Dx`` becomes ``D = r.D{1}``, and so on. Time derivatives are ignored at present. Derivatives are returned at all lattice locations.

If the derivative ``D`` is omitted, a first order x-derivative is returned.
Note that :func:`xd` returns a lattice observable, as required when used in the :attr:`in.observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.

.. _sec-lattice-structure:




Grid coordinates and time
=============================

Time and space
--------------

The default spatial grid
 for plotted output data is rectangular, with periodic boundary conditions in space, and

::

    dx(i) = in.ranges(i) / (in.points(i) - 1)

The time index is ``1``, and the space index ``i`` ranges from ``2`` to :attr:`in.dimension`. The maximum space-time dimension is ``in.dimension = 4``, while ``in.ranges(i)`` is the time and space duration of the simulation, and ``in.points(i)`` is the total number of plotted points in the ``i``-th direction.

Time is advanced in basic integration steps that are equal to or smaller than ``dx(1)``, for purposes of controlling and reducing errors:

::

    dt = dx(1) / (in.steps * nc)

Here, :attr:`in.steps` is the minimum number of steps used per plotted point, and ``nc = 1, 2`` is the check number. If ``nc = 1``, the run uses coarse time-divisions. If ``nc = 2`` the steps are halved in size for error-checking. Error-checking can be turned off if not required.

Low-level functions
-------------------

The xSPDE program is function oriented: low-level functions are used to define initial conditions, equations and observables, amongst many other things described below.

Functions of a single lattice have arguments in the following order:

-  the field ``a`` or initial random variable ``v``;
-  the stochastic noise ``z`` or other fields;
-  non-field arguments;
-  the grid structure ``r``.

The first argument, ``a``, is a real or complex vector field. This is a matrix whose first dimension is the field index. The second dimension is the lattice index.

The second argument, ``z``, if needed, is a real random noise, corresponding to :math:`\zeta` in the mathematical notation. This is a matrix whose first dimension is the noise index. The second dimension is the lattice index.

The last function argument is the  :ref:`lattice structure <sec-lattice-structure>`, ``r``. This contains data about the integration lattice. The most important constants are :attr:`r.t`, the current time, and the space coordinates, :attr:`r.x`, :attr:`r.y`, :attr:`r.z`. Other data stored in the lattice structure is explained in later chapters.

Functions of multiple lattice sequences take current arguments first, and the oldest arguments last.

Arrays
------

In all function calls, the variables used are matrices. The most important first dimension used is the field length :attr:`in.fields`. The second dimension in all arrays is the lattice index, with a length ``n.lattice = ensembles(1) * points(2) * ... * points(dimension)``. Here ``ensembles(1)`` is the number of stochastic samples integrated as an array.

For reference, the field dimensions are:

- ``a, da, L = [r.fields, r.nlattice]``;
- ``v = [r.randoms(1)+r.randoms(2), r.nlattice]``;
- ``z = [r.noises(1)+r.noises(2), r.nlattice]``;
- ``D.x, r.x, r.kx = [1, r.nlattice]``;
- ``o = [1, r.nlattice]``.

Each observable is defined by a function in a cell array with length :attr:`in.graphs`.


Simulation parameters
---------------------

For each simulation in the ``input`` sequence, the input parameters and functions are specified as a data structure, ``in``. These can be entered either interactively or as part of a simulation function file. The function file approach allows recycling and editing, so it is better for a large project.

There are extensive default preferences to simplify the inputs. If any inputs are omitted, there are default values which are set by inpreferences in all cases. These defaults are changed by editing the inpreferences function. The :func:`xgrpreferences` function is used to supply graphics default values.

**For vector or cell inputs, an input shorter than required is padded to the right using default values.**


.. _sec-input:

Input parameters and user functions
===================================

A sequence of simulations is obtained from inputs in a cell array, as ``input = {in1, in2, ...}``. The input parameters of each simulation in the sequence are specified in a Matlab structure. If there is one simulation, just one structure can be input, without the braces. This data is also passed to the :func:`xgraph` function. The inputs are numbers, vectors, strings, functions and cell arrays. All xSPDE metadata has preferred values, so only changes from the preferences need to be input. The resulting data is stored internally as a sequence of structures in a cell array, to describe the simulation sequence.

The standard way to input each parameter value is:

::

    in.label = parameter

The standard way to input each function is:

::

    in.label = @function-name

The inputs are scalar or vector parameters or function handles. Quantities relating to graphed averages are cell arrays, indexed by the graph number. The available inputs, with their default values in brackets, are given below.

Simulation metadata, including all preferred default values that were used in a particular simulation, is also stored for reference in any xSPDE output files. This is done in both the ``.mat`` and the ``.h5`` output files, so the entire simulation can be easily reconstructed or changed.

Note that inputs can be numbers, vectors, strings or cells arrays. To simplify the inputs, some conventions are used, as follows:

- All input data has default values
- Vector inputs of numbers are enclosed in square brackets, ``[...]``.
- Where multiple inputs of strings, functions or vectors are needed they should be enclosed in curly brackets, ``{...}``, to create a cell array.
- Vector or cell array inputs with only one member don’t require brackets.
- Incomplete or partial vector or cell array inputs are filled in with the last applicable default value.
- New function definitions can be just handles pointing elsewhere, or else defined inline.


Input parameters
----------------

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



.. attribute:: in.noises

    *Default:* :attr:`in.fields`

    This gives the number of stochastic noises generated per lattice point, in coordinate and momentum space respectively. Set to zero (``in.noises = 0``) for no noises. This is the number of *rows* in the noise-vector. Noises can be correlated either in ordinary or momentum spaces. The second input is the dimension of noises in k-space. It can be left out if zero.

    ::

        in.noises = [in.noises(1), in.noises(2)] >= 0.


.. attribute:: in.randoms

    *Default:* :attr:`in.noises`

    This gives the number of random fields generated per lattice point for the initial noise, in coordinate and momentum space. Set to zero (``in.randoms = 0``) for no random fields. Random fields can be correlated either in ordinary or momentum spaces. The second input is the dimension of random fields in momentum space. It can be left out if zero.

    ::

        in.randoms = [in.randoms(1), in.randoms(2)] >= 0

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

    There is one transform vector per observable. The ``j``-th index, ``t(j)``, indicates a Fourier transform on the ``j``-th axis. The normalization of the Fourier transform is such that the :math:`k=0` value in momentum space corresponds to the integral over space, with an additional factor of :math:`1/\sqrt{2\pi}`. This gives a Fourier integral which is symmetrically normalized in ordinary and momentum space. The Fourier transform that is graphed is such that
    :math:`k=0` is the *central* value.

.. attribute:: in.olabels

    *Default:* ``{'a_1', ...}``

    **Cell array** of labels for the graph axis observable functions. These are text labels that are used on the graph axes. The default value is ``'a_1'`` if the default observable is used, otherwise it is blank. This is overwritten by any subsequent label input when the graphics program is run:

    ::

        in.olabels{n} = 'string'

.. attribute:: in.c

    This starting letter is always reserved to store user-specified constants and parameters.  It is passed to user functions, and can be any data. All inputs --- including ``c`` data --- are copied into the data files and also the lattice structure ``r``.

    ::

        in.c = anything


Invariant inputs
----------------

The following can’t be changed during a sequence in the current xSPDE version --- the specified values for the first simulation will be used:

#. The extrapolation order

#. The number of ensembles (2)

#. The number of ensembles (3)

#. The output file-name


Input functions
---------------

A stochastic equation solver requires the definition of an initial distribution and a time derivative. In xSPDE, the time derivatives is divided up into a linear term including space derivatives, used to define an interaction picture, and the remaining derivatives. In addition, one must define quantities to be averaged over during the simulation, called graphs in xSPDE. These are all defined as functions, specified below.

.. attribute:: in.initial(v,r)

    *Default:* :func:`xinitial`

    Initializes the fields :math:`a` for the first simulation in a sequence. The initial Gaussian random field variable, ``v``, has unit variance if :attr:`in.dimension` is ``1`` or else is delta-correlated in space, with variance ``1/r.dV`` (:math:`\equiv 1/(dx_2...dx_d)`) for :math:`d` space-time dimensions. If :attr:`in.randoms` is specified in the input, ``v`` has a first dimension of ``in.randoms(1) + in.randoms(2)``. If not specified, the default for ``in.randoms`` is  ``in.noises``. If not specified, the default of :func:`in.initial` is ``a = 0``.

.. attribute:: in.transfer(v,r,a0,r0)

    *Default:* :func:`xtransfer`

    Initializes the fields :math:`a` for subsequent calculations in a sequence. Otherwise, this function behaves in a similar way to :attr:`in.initial`. The function includes the previous field ``a0`` and lattice ``r0``. The default set by :func:`xtransfer` is ``a = a0``.

.. attribute:: in.da(a,z,r)

    *Default:* :func:`xda`

    Calculates derivatives :math:`da` of the equation. The noise vector, ``z``, has variance :math:`1/(dx_{1}..dx_{d})`, for dimension :math:`d \le 4`, and a first dimension  whose default value is :attr:`in.fields` if :attr:`in.noises` are not given. Otherwise, it has a first dimension of ``in.noises(1) + in.noises(2)``. The second type of input noise allows for spatially correlated and filtered noise specified in momentum space.

.. attribute:: in.linear(D,r)

    *Default:* :func:`xlinear`

    A user-definable function which returns the linear coefficients :math:`L` in Fourier space. This is a function of the differential operator ``D``. The default is zero. Here ``D`` is a structure with components ``D.x``, ``D.y``, ``D.z``, which correspond to :math:`\partial / \partial x`, :math:`\partial / \partial y`, :math:`\partial / \partial z` respectively. Each component has an array dimension the same as the coordinate lattice.

.. attribute:: in.observe(a,r)

    *Default:* cell array of :func:`xobserve`

    **Cell array** of function handles that take the current field and returns a real observable ``o`` with dimension of ``[1, n.lattice]``. The default observable is the first real field amplitude. Note the use of braces for cell arrays! One can also input these individually as ``in.observe{1} = @(a,r) f(a,r)``, using an inline anonymous function. The total number of observe functions is stored internally as :attr:`in.graphs`. The fields ``a`` passed in the input are transformed according to the :attr:`in.transforms` metadata.

.. attribute:: in.rfilter(r)

    *Default:* :func:`xrfilter`

    Returns the momentum-space filters for the input random terms. Each component has an array dimension the same as the input random fields in momentum space, that is, the return dimension is ``[r.randoms(2), r.nlattice]``.

.. attribute:: in.nfilter(r)

    *Default:* :func:`xnfilter`

    Returns the momentum-space filters for the propagation noise terms. Each component has an array dimension the same as the random noises in momentum space, that is, the return dimension is ``[r.noises(2), r.nlattice]``.


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

    Flag for storing raw trajectory data. If this flag is turned on, raw trajectories are stored in memory. The raw data is returned in function calls and also written to a file on completion, if a file-name is included.

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

        in.ipsteps = 1, 2, 3, ..

.. attribute:: in.file

    *Default:* ``''``

    Matlab or *HDF5* file name for output data. Includes all data and parameter values, including raw trajectories if ``in.raw = 1``. If not needed just omit this. A Matlab filename should end in ``.mat``, while an HDF5 file requires the filename to end in ``.h5``. For a sequence of inputs, the filename should be given in the first structure of the sequence, and the entire sequence is stored.

    ::

        in.file = 'file-name'

Advanced input functions
------------------------

Advanced input functions are user-definable functions which don’t usually need to be changed from default values. They allow customization and extension of xSPDE. These are as follows:

.. attribute:: in.grid(r)

    *Default:* :func:`xgrid`

    Initializes the grid of coordinates in space.

.. attribute:: in.noisegen(r)

    *Default:* :func:`xnoisegen`

    Generates arrays of noise terms ``xi`` for each point in time.

.. attribute:: in.randomgen(r)

    *Default:* :func:`xrandomgen`

    Generates a set of initial random fields ``v`` to initialize the fields simulated.

.. attribute:: in.step(a,z,dt,r)

    *Default:* :func:`xMP`

    Specifies the stochastic integration routine for the field ``a``, given a step in time ``dt`` and noise ``z``, together with the interaction-picture propagator :attr:`r.propagator` which is part of the lattice structure. It returns the new field ``a``. This function can be set to any of the predefined stochastic integration routines provided with xSPDE, described in the :ref:`chap-algorithms` chapter. User-written functions can also be used. The standard method, :func:`xMP`, is a midpoint integrator.

.. attribute:: in.prop(a,r)

    *Default:* :func:`xprop`

    Returns the fields propagated for one step in the interaction picture, depending on the initial field ``a``, and the propagator array :attr:`r.propagator`. Note that the time-step used in :attr:`r.propagator` depends on the input time-step, the error-checking and the algorithm.

.. attribute:: in.propfactor(nc,r)

    *Default:* :func:`xpropfactor`

    Returns the transfer array :attr:`r.propagator`, used by the :attr:`in.prop` function. The time propagated is a fraction of the current integration time-step, :attr:`r.dt`. It is equal to ``1 / in.ipsteps`` of the integration time-step.


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

.. function:: in.function(data,in)

    This is a cell array of graphics function handles. Use when a graph is needed that is a function of the observed averages. The default value generates all the averages that are in the simulated data. The input is the data array of averages, and the output is the  data array that is plotted.
    
.. function:: in.xfunctions(x,in)

    This is a nested cell array of graphics axis transformations. Use when a graph is needed with an axis that is a function of the original axes. The default value generates all the averages that are in the simulated data. The input of the function is the original axis coordinates, and the output is the new coordinate set. Called as in.xfunctions{n}{nd} for the n-th graph and axis nd.
    

.. function:: in.compare(t,in)

    This is a cell array of comparative functions. Each takes the time or frequency vector - or whichever is the first dimension plotted - and returns comparison results for a graphed observable, as a function versus time or frequency (etc). Comparison results are graphed with a dashed line, for the two-dimensional graphs versus the first plotted dimension. There is no default function handle.

Graphics parameters
~~~~~~~~~~~~~~~~~~~

For uniformity, the graphics parameters that reference an individual data object are cell arrays, indexed over the graph number using braces ``{}``. If a different type of input is used, like a scalar or matrix, xSPDE will attempt to convert the type. The axis labels are cell arrays, indexed over dimension. The graph number used to index these cell arrays refers to the data object, and there can be multiple plots obtained, depending on the graphics input.

Together with default values, they are:

.. attribute:: in.gversion

    *Default:* ``'xGRAPH1.1'``

    This sets the version number of the graphics output. There is typically no need to input this.

    ::

        in.gversion = 'current version name'
        
    
.. attribute:: in.graphs

    *Default:* ``-1``

    If non-negative, this sets the maximum number of graphed datasets. Can be used to suppress unwanted graphs from an xSPDE graphics script. If omitted or set to the default value, all the graphs output from the in.grfunc graphics processing function are plotted.
    
    ::

        in.graphs = -1,..

.. attribute:: in.axes

    *Default:* ``{{0,0,0,..}}``

   Gives the axis and points plotted ``p`` for each plotted function. As special cases,  ``p = 0``, is the default value that gives the entire axis, while  ``p = -1`` generates one point on the axis, namely the last point for the time axis and the midpoint for the space axes. Other values are vector range indicators, for example ``p = 5`` plots the fifth point, while ``p = 1:4:41`` plots every fourth point. For each graph type ``n`` the axes can be individually specified. If more than three axes are specified, only the first three are used. The others are set to default values.

    ::

        in.axes{n} = {p1,p2,p3,..pd}

.. attribute:: in.font

    *Default:* ``{18, ...}``

    This sets the default font size for the graph labels. This can be changed per graph.

    ::

        in.font{n} > 0

.. attribute:: in.minbar

    *Default:* ``{0.01, ...}``

    This is the minimum relative error-bar that is plotted. Set to a large value to suppress unwanted error-bars, although its best not to ignore the error-bar information! 

    ::

        in.minbar{n} >= 0
        
        .. attribute:: in.esample

    *Default:* ``{1, ...}``

    This is the flag for plotting sampling error. Set to zero to suppress unwanted sampling error lines and just plot means, although its best not to ignore this information! 

    ::

        in.esample{n} >= 0

.. attribute:: in.images

    *Default:* ``{0, 0, 0, ...}``

    This is the number of 3D, transverse o-x-y movie images plotted as discrete time slices. Only valid if :attr:`in.dimension` is greater than 2. Note that, if present, the coordinates not plotted are set to their central value, for example ``z = 0``, when plotting the transverse images. This input should have a value from ``in.images(n) = 0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`in.graphs`:

    ::

        in.images{n} = 0 ... in.points(1)

.. attribute:: in.imagetype

    *Default:* ``{1, 1, ...}``

    This is the *type* of transverse o-x-y movie images plotted. If an element is ``1``, a perspective surface plot is output, for ``2``, a gray plot with colours is output, or for ``3`` a contour plot with 10 equally spaced contours is generated. This has a vector length equal to :attr:`in.graphs`.

    ::

        in.imagetype{n} = 1, 2, 3

.. attribute:: in.transverse

    *Default:* ``{0, 0, ...}``

    This is the number of 2D, transverse o-x images plotted as discrete time slices. Only valid if :attr:`in.dimension` is greater than 2. Note that, if present, the y,z-coordinates are set to their central values, when plotting the transverse images. Each element should be from ``0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`in.graphs`:

    ::

        in.transverse{n}=0 ... in.points(1)

.. attribute:: in.headers

    *Default:* ``{'head1', 'head2', ...}``

    This is a string variable giving the graph headers for each type of function plotted. The default value is an empty string ``''``, which gives the simulation heading. Use a space ``' '`` for no headers. It is useful to include simulation headers - which is the default - to identify graphs in preliminary stages, while they may not be needed in a published final result. 

    ::

        in.headers{n} = 'my_graph_header'

.. attribute:: in.pdimension

    *Default:* ``{3, 3, ...}``

    This is the maximum space-time grid dimension for each plotted quantity. The purpose is eliminate unwanted graphs. For example, it may be useful to reduce the maximum dimension when averaging in space. Higher dimensional graphs are not needed, as the data is duplicated. Averaging can be useful for checking conservation laws, or for averaging over homogeneous data to reduce sampling errors. All graphs are suppressed if it is set to zero. Any three dimensions can be chosen using the axes command.

    ::

        in.pdimension{n} \ge 0 
        
        

.. attribute:: in.xlabels

    *Default:* ``{'t', 'x', 'y', 'z'}`` or ``{'t', 'x_1', 'x_2', 'x_3'}``

    Labels for the graph axis independent variable labels, vector length of :attr:`in.dimension`. The numerical labeling default is used when the ``in.numberaxis`` option is set. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.xlabels = {in.xlabels(1), ..., in.xlabels(in.dimension)}

.. attribute:: in.klabels

    *Default:* ``{'\\omega', 'k\_x', 'k\_y', 'k\_z'}`` or ``{'\\omega', 'k\_1', 'k\_2', 'k\_3'}``

    Labels for the graph axis Fourier transform labels, vector length of :attr:`in.dimension`. The numerical labeling default is used when the ``in.numberaxis`` option is set. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.klabels = {in.klabels(1), ..., in.klabels(in.dimension)}

Graphics projections
~~~~~~~~~~~~~~~~~~~~

If there is a spatial grid, the graphics program automatically generates several graphs for each observable, depending on space dimension. The maximum dimension that is plotted as set by :attr:`in.pdimension`. In the plots, the lattice is projected down to successively lower dimensions.

For each observable, the projection sequence is as follows:

-  If :attr:`in.dimension` is ``4`` or greater, a central :math:`z` point ``nz = 1 + floor(in.points(4)/2)`` is picked. For example, with 35 points, the central point is ``nz = 18``.

-  This gives a three dimensional space-time lattice, which is treated the same as if :attr:`in.dimension` is ``3``.

-  If :attr:`in.images` are specified, two-dimensional :math:`x-y` plots are generated at equally spaced time intervals. If there is only one image, it is at the last time-point. Different plot-types are used depending on the setting of :attr:`in.imagetype`.

-  A central :math:`y` point ``ny = 1 + floor(in.points(3)/2)`` is picked. This gives a two dimensional space-time lattice, which is treated the same as if :attr:`in.dimension` is ``2``. If :attr:`in.transverse` is specified, one-dimensional :math:`x` plots are generated at equally spaced time intervals, as before.

-  A central :math:`x` point ``nx = 1 + floor(in.points(2)/2)`` is picked. This gives a one dimensional time lattice, which is treated the same as if :attr:`in.dimension` is ``1``.

-  Plots of observable vs time are obtained, including sampling errors and error bars. If comparison graphs are specified using :func:`in.compare` functions, they are plotted also, using a dotted line. A difference graph is also plotted when there is a comparison.

Parameter structure
===================

Internally, xSPDE parameters are stored in a cell array, ``latt``, of structures ``r``, which is passed to functions. This includes all the data given above inside the ``in`` structure. In addition, it includes the table of computed parameters given below.

User application constants and parameters should not be reserved names; :attr:`in.c` and all names starting with ``in.c`` will always be available in all versions of xSPDE.

A parameter structure contains information about the space-time grid and is passed to various functions, for instance :attr:`in.da` or :attr:`in.step`. The corresponding parameter is commonly marked as `r`.

.. attribute:: r.t

    Current value of time, :math:`t`.

.. attribute:: r.x

.. attribute:: r.y

.. attribute:: r.z

    Coordinate grids of :math:`x`, :math:`y`, :math:`z`.

.. attribute:: r.x{n}

    Higher dimensions are labeled numerically as :math:`x_1`,..  :math:`x_6`, and so on. This numerical axis convention can be set even for lower dimensions if ``in.numberaxis`` is set to 1.

.. attribute:: r.kx

.. attribute:: r.ky

.. attribute:: r.kz

    Grids in momentum space: :math:`k_x`, :math:`k_y`, :math:`k_z`.

.. attribute:: r.k5

    Higher dimensions are labeled numerically as :math:`k_5`,  :math:`k_6`, and so on.

.. attribute:: r.dt

    Output time-step between stored points for data averages.
    
.. attribute:: r.dtr

    Current reduced time-step used for integration.

.. attribute:: r.dx

    Steps in coordinate space: :math:`[t,x,y,z,x_5,..]`.

.. attribute:: r.dk

    Steps in momentum space: :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]`.

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

    Computational Fourier transform axes in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},.. ]` (vector cells).

.. attribute:: r.kg

    Graphics  Fourier transform axes in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]` (vector cells).

.. attribute:: r.kranges

    Range in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]` (vector).

.. attribute:: r.s.dx

    Initial stochastic normalization.

.. attribute:: r.s.dxt

    Propagating stochastic normalization.

.. attribute:: r.s.dk

    Initial :math:`k` stochastic normalization.

.. attribute:: r.s.dkt

    Propagating :math:`k` stochastic normalization.

.. attribute:: r.nspace

    Number of spatial lattice points: ``in.points(2) * .. * in.points(in.dimension)``.

.. attribute:: r.nlattice

    Total lattice: ``in.ensembles(1) * r.nspace``.

.. attribute:: r.ncopies

    Total copies of stochastic integrations: ``in.ensembles(2) * in.ensembles(3)``.

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

    Returns the real part of ``a(1,:)``.

.. function:: xrfilter(r)

    Returns an array of ones.

.. function:: xnfilter(r)

    Returns an array of ones.

.. function:: xgrid(r)

    Sets grid points in lattice from coordinate vectors. Returns the ``r`` structure with added grid points.

.. function:: xnoisegen(r)

    Generates random noise matrix :math:`z`.

.. function:: xrandomgen(r)

    Generates initial random field matrix :math:`v`.

.. function:: xpropfactor(nc, r)

    Returns the interaction picture propagation factor. ``nc`` is a check index, ``r`` is a lattice structure.


Frequently asked questions
==========================

Answers to some frequent questions, and reminders of points in this chapter are:

-  Can you average other stochastic quantities apart from the field?

   -  Yes: just specify the functions that need to be averaged using the user function :attr:`in.observe`.

-  Can you have functions of the current time and space coordinate?

   -  Yes: xSPDE functions support this using the structure ``r``, as :attr:`r.t`, :attr:`r.x`, :attr:`r.y`, :attr:`r.z`, or  :attr:`r.t`, ``r.x{1}``, and so on, for more than four space-time dimensions.

-  Can you have several independent stochastic variables?

   -  Yes, input this using ``in.fields > 1``.

-  Are higher dimensional differential equations possible?

   -  Yes, this requires setting ``in.dimension > 1``. This is essentially unlimited in xSPDE except for memory requirements.

-  Can you have spatial partial derivatives?

   -  Yes, provided they are linear in the fields; these are obtainable using the function :attr:`in.linear`.

-  Can you delete the graph heading?

   -  Yes, this is turned off if you set :attr:`in.headers` to ``0``.

-  Why are there two lines in the graphs sometimes?

   -  These are one standard deviation sampling error limits, generated when ``in.ensembles(2,3) > 1``.

-  Why is there just one line in some graphs, with no sampling errors indicated?

   -  You need ``in.ensembles(2)`` or ``(3)`` for this; see previous question.

-  What are the error bars for?

   -  These are the estimated maximum errors due to finite step-sizes in time.



High-level xSPDE functions and objects
======================================

The high-level xSPDE functions process the input parameters, producing simulation data and graphs. Function parameters are in the order of ``input``, which specifies the simulation, then ``data``, which is needed for graphics output.

.. function:: xspde(input)

    This is the combined xSPDE function. It accepts a simulation sequence, ``input``. This can be a single structure, ``in``, or else a cell array of structures, ``{in1,in2,..}``, for  sequences. Output graphs are displayed, and it returns the output ``[maxerror, input, data]``, where ``maxerror`` is the maximum error or difference found. If a filename is specified, it generates an output data file. It calls the functions :func:`xsim` and :func:`xgraph`.


.. function:: xsim(input)

    This is the xSPDE simulation function. Like :func:`xspde`, it accepts input parameters in ``input``. It returns ``[error, input, data, raw]``, where: ``error = [error(1),error(2)]`` is a vector of maximum step-size and sampling errors, ``input`` is the full input structure or cell array for sequences, including default values, and ``data`` is a cell array of average observables. If the ``in.raw`` option is used, data for the actual trajectories is output in ``raw``. This can be run as a stand-alone function if no graphs are required.

.. function:: xgraph(data [,input])

    This is the xSPDE graphics function. It takes computed simulation  ``data`` and ``input``. It plots graphs, and returns the maximum difference ``diff`` from comparisons with user-specified comparison functions. The ``data`` should have as many cells as ``input`` cells, for sequences. 
    If ``data = 'filename.h5'`` or ``data= 'filename.mat'``, the specified file is read both for ``input`` and ``data``. Here ``.h5`` indicates an HDF5 file, and ``.mat`` indicates a Matlab file.
    When the``data`` input is given as a filename, input parameters in the file are replaced by any of the the new ``input`` parameters that are specified.  Any stored ``input`` can be overwritten, allowing graphs to be modified with new labels.
    If the ``input`` field is not present, then a file-name must be specified in the ``data`` field. 



Open object-oriented architecture
----------------------------------

As well as extensibility through sequences, which was described in :ref:`chap-projects`, in the section :ref:`sec-sequential-integration`, the architecture of xSPDE allows functional extensions.

The input metadata includes both data and methods acting on the data, in the tradition of object-oriented programs. Yet there is no strict class typing. Users are encouraged to adapt the xSPDE program by adding input parameters and methods to the input structures.

*Such unorthodox object orientation is deliberate*.

This open extensibility permits arbitrary functions to be specified in the ``in`` structures. All functions and parameters have default values in xSPDE. It is also possible to include user defined functions, provided they satisfy the API definitions given below. This is achieved simply by including the relevant function handles and parameters in the input metadata.

Internal parameters and function handles from the ``in`` structures, together with any computed and default parameters and functions,  are stored in the :ref:`lattice structure <sec-lattice-structure>` ``r``. These are available to all user-defined functions. Use of pre-existing reserved names is not advisable, and the structure ``r.c`` is always reserved for user constants, parameters and functions if required.

The xSPDE software architecture is intended to be easily extended, and users are strongly encouraged to develop their own libraries and contribute to the xSPDE function pool. Because this generally requires new functions and parameters, the internal data architecture is as open as possible.

For example, to define your own integration function, include in the xSPDE input the line:

::

    in.step = @Mystep;

Next, include anywhere on your Matlab path, the function definition, for example:

::

    function a = Mystep(a,z,dt,r)
        % a = Mystep(a,z,dt,r) propagates a step my way.
        ...
        a = ...;
    end





Averages, integrals and derivatives
===================================

There are functions available in xSPDE for averages, integrals and derivatives. These can be used to calculate observables for plotting, but are also available for calculating stochastic derivatives as part of the stochastic equation. They operate in parallel over the local ensemble and lattice dimensions. They take a scalar quantity, for example a single field component, and return an average, a space integral, and a spatial derivative respectively. In each case the first argument is the field, the second argument is a vector defining the type of operation, and the last argument is the parameter structure, ``r``. If there are only two arguments, the operation vector is replaced by its default value.

Averages
--------

This function allows one to extract local ensemble averages which may be needed as part of an overall calculation involving a ratio or product of averages. In addition, spatial grid averages can be used to obtain stochastic results with reduced sampling errors if the overall grid is homogeneous.

An average is carried out using the :func:`xave` function, which is defined as follows:

.. function:: xave(o, [av, ] r)

    This function takes a scalar field or observable ``o = [1, n.lattice]``, defined on the xSPDE local lattice, and returns an average with dimension ``[1, n.lattice]``. The input is a field or observable ``o``, and an optional averaging switch ``av``. If ``av(j) > 0``, an average is taken over dimension ``j``. Dimensions are labelled from ``j = 1 ... 4`` as elsewhere. The first index indicates a local ensemble average, while subsequent indices indicate averages over the spatial grid. If the ``av`` vector is omitted, the average is only taken over the local ensemble. Averages are returned at all lattice locations. To average over both the local ensemble and all space dimensions, just use ``xave(o)``.

Higher dimensional graphs of grid averages are generally not useful, as they are simply flat. The xSPDE program allows the user to remove unwanted higher dimensional graphs of average variables. This is achieved by setting the corresponding element of :attr:`in.pdimension` to the highest dimension required, which depends on which dimensions are averaged.

For example, to average over the entire ensemble plus space lattice and indicate that only time-dependent graphs are required, set ``av = in.dx`` and:

::

    in.pdimension = 1

Note that :func:`xave` on its own gives identical results to those calculated in the :attr:`in.observe` functions. Its utility comes when more complex combinations or functions of ensemble averages are required.

Integrals
---------

Integrals over the spatial grid allow calculation of conserved or other global quantities. To take an integral over the spatial grid,  use the xSPDE :func:`xint` function:

.. attribute:: xint(o, [dx, ] r)

    This function takes a scalar ``o``, and returns a space integral over selected dimensions with vector measure ``dx``. If ``dx(j) > 0`` an integral is taken over dimension ``j``. Dimensions are labelled from ``j = 1, ...`` as in all xspde standards. Time integrals are ignored at present. Integrals are returned at all lattice locations. To integrate over an entire lattice, set ``dx = r.dx``, otherwise set ``dx(j) = r.dx(j)`` for selected dimensions ``j``.

As with averages, the xSPDE program allows the user to remove unwanted higher dimensional graphs when the integrated variable is used as an observable. For example, in a four dimensional simulation with integrals taken over the :math:`y` and :math:`z` coordinates, only :math:`t`- and :math:`x`-dependent graphs are required. Hence, set ``dx`` to ``[0, 0, r.dx(3), r.dx(4)]``, and:

::

    in.pdimension = 2

If momentum-space integrals are needed, use the transform switch to make sure that the field is Fourier transformed, and input :attr:`r.dk` instead of :attr:`r.dx`. Note that :func:`xint` returns a lattice observable, as required when used in the :attr:`in.observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.


Derivatives
-----------

The code to take a spatial derivative is carried out using the xSPDE :func:`xd` function:

.. attribute:: xd(o, [D, ] r)

This function takes a scalar ``o``, and returns a derivative over selected dimensions with a derivative ``D``.  Set ``D = r.D.x`` for a first order x-derivative, ``D = r.D.y`` for a first order y-derivative, and similarly ``D = r.D.x.*r.D.y`` for a cross-derivative in ``x`` and ``y``. Higher derivatives require powers of these. Time derivatives are ignored at present. Derivatives are returned at all lattice locations.

If the derivative ``D`` is omitted, a first order x-derivative is returned.
Note that :func:`xd` returns a lattice observable, as required when used in the :attr:`in.observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.

.. _sec-lattice-structure:




Grid coordinates and time
=============================

Time and space
--------------

The default spatial grid
 for plotted output data is rectangular, with periodic boundary conditions in space, and

::

    dx(i) = in.ranges(i) / (in.points(i) - 1)

The time index is ``1``, and the space index ``i`` ranges from ``2`` to :attr:`in.dimension`. The maximum space-time dimension is ``in.dimension = 4``, while ``in.ranges(i)`` is the time and space duration of the simulation, and ``in.points(i)`` is the total number of plotted points in the ``i``-th direction.

Time is advanced in basic integration steps that are equal to or smaller than ``dx(1)``, for purposes of controlling and reducing errors:

::

    dt = dx(1) / (in.steps * nc)

Here, :attr:`in.steps` is the minimum number of steps used per plotted point, and ``nc = 1, 2`` is the check number. If ``nc = 1``, the run uses coarse time-divisions. If ``nc = 2`` the steps are halved in size for error-checking. Error-checking can be turned off if not required.

Low-level functions
-------------------

The xSPDE program is function oriented: low-level functions are used to define initial conditions, equations and observables, amongst many other things described below.

Functions of a single lattice have arguments in the following order:

-  the field ``a`` or initial random variable ``v``;
-  the stochastic noise ``z`` or other fields;
-  non-field arguments;
-  the grid structure ``r``.

The first argument, ``a``, is a real or complex vector field. This is a matrix whose first dimension is the field index. The second dimension is the lattice index.

The second argument, ``z``, if needed, is a real random noise, corresponding to :math:`\zeta` in the mathematical notation. This is a matrix whose first dimension is the noise index. The second dimension is the lattice index.

The last function argument is the  :ref:`lattice structure <sec-lattice-structure>`, ``r``. This contains data about the integration lattice. The most important constants are :attr:`r.t`, the current time, and the space coordinates, :attr:`r.x`, :attr:`r.y`, :attr:`r.z`. Other data stored in the lattice structure is explained in later chapters.

Functions of multiple lattice sequences take current arguments first, and the oldest arguments last.

Arrays
------

In all function calls, the variables used are matrices. The most important first dimension used is the field length :attr:`in.fields`. The second dimension in all arrays is the lattice index, with a length ``n.lattice = ensembles(1) * points(2) * ... * points(dimension)``. Here ``ensembles(1)`` is the number of stochastic samples integrated as an array.

For reference, the field dimensions are:

- ``a, da, L = [r.fields, r.nlattice]``;
- ``v = [r.randoms(1)+r.randoms(2), r.nlattice]``;
- ``z = [r.noises(1)+r.noises(2), r.nlattice]``;
- ``D.x, r.x, r.kx = [1, r.nlattice]``;
- ``o = [1, r.nlattice]``.

Each observable is defined by a function in a cell array with length :attr:`in.graphs`.


Simulation parameters
---------------------

For each simulation in the ``input`` sequence, the input parameters and functions are specified as a data structure, ``in``. These can be entered either interactively or as part of a simulation function file. The function file approach allows recycling and editing, so it is better for a large project.

There are extensive default preferences to simplify the inputs. If any inputs are omitted, there are default values which are set by inpreferences in all cases. These defaults are changed by editing the inpreferences function. The :func:`xgrpreferences` function is used to supply graphics default values.

**For vector or cell inputs, an input shorter than required is padded to the right using default values.**


.. _sec-input:

Input parameters and user functions
===================================

A sequence of simulations is obtained from inputs in a cell array, as ``input = {in1, in2, ...}``. The input parameters of each simulation in the sequence are specified in a Matlab structure. If there is one simulation, just one structure can be input, without the braces. This data is also passed to the :func:`xgraph` function. The inputs are numbers, vectors, strings, functions and cell arrays. All xSPDE metadata has preferred values, so only changes from the preferences need to be input. The resulting data is stored internally as a sequence of structures in a cell array, to describe the simulation sequence.

The standard way to input each parameter value is:

::

    in.label = parameter

The standard way to input each function is:

::

    in.label = @function-name

The inputs are scalar or vector parameters or function handles. Quantities relating to graphed averages are cell arrays, indexed by the graph number. The available inputs, with their default values in brackets, are given below.

Simulation metadata, including all preferred default values that were used in a particular simulation, is also stored for reference in any xSPDE output files. This is done in both the ``.mat`` and the ``.h5`` output files, so the entire simulation can be easily reconstructed or changed.

Note that inputs can be numbers, vectors, strings or cells arrays. To simplify the inputs, some conventions are used, as follows:

- All input data has default values
- Vector inputs of numbers are enclosed in square brackets, ``[...]``.
- Where multiple inputs of strings, functions or vectors are needed they should be enclosed in curly brackets, ``{...}``, to create a cell array.
- Vector or cell array inputs with only one member don’t require brackets.
- Incomplete or partial vector or cell array inputs are filled in with the last applicable default value.
- New function definitions can be just handles pointing elsewhere, or else defined inline.


Input parameters
----------------

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



.. attribute:: in.noises

    *Default:* :attr:`in.fields`

    This gives the number of stochastic noises generated per lattice point, in coordinate and momentum space respectively. Set to zero (``in.noises = 0``) for no noises. This is the number of *rows* in the noise-vector. Noises can be correlated either in ordinary or momentum spaces. The second input is the dimension of noises in k-space. It can be left out if zero.

    ::

        in.noises = [in.noises(1), in.noises(2)] >= 0.


.. attribute:: in.randoms

    *Default:* :attr:`in.noises`

    This gives the number of random fields generated per lattice point for the initial noise, in coordinate and momentum space. Set to zero (``in.randoms = 0``) for no random fields. Random fields can be correlated either in ordinary or momentum spaces. The second input is the dimension of random fields in momentum space. It can be left out if zero.

    ::

        in.randoms = [in.randoms(1), in.randoms(2)] >= 0

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

    There is one transform vector per observable. The ``j``-th index, ``t(j)``, indicates a Fourier transform on the ``j``-th axis. The normalization of the Fourier transform is such that the :math:`k=0` value in momentum space corresponds to the integral over space, with an additional factor of :math:`1/\sqrt{2\pi}`. This gives a Fourier integral which is symmetrically normalized in ordinary and momentum space. The Fourier transform that is graphed is such that
    :math:`k=0` is the *central* value.


.. attribute:: in.c

    This starting letter is always reserved to store user-specified constants and parameters.  It is passed to user functions, and can be any data. All inputs --- including ``c`` data --- are copied into the data files and also the lattice structure ``r``.

    ::

        in.c = anything


Invariant inputs
----------------

The following can’t be changed during a sequence in the current xSPDE version --- the specified values for the first simulation will be used:

#. The extrapolation order

#. The number of ensembles (2)

#. The number of ensembles (3)

#. The output file-name


Input functions
---------------

A stochastic equation solver requires the definition of an initial distribution and a time derivative. In xSPDE, the time derivatives is divided up into a linear term including space derivatives, used to define an interaction picture, and the remaining derivatives. In addition, one must define quantities to be averaged over during the simulation, called graphs in xSPDE. These are all defined as functions, specified below.

.. attribute:: in.initial(v,r)

    *Default:* :func:`xinitial`

    Initializes the fields :math:`a` for the first simulation in a sequence. The initial Gaussian random field variable, ``v``, has unit variance if :attr:`in.dimension` is ``1`` or else is delta-correlated in space, with variance ``1/r.dV`` (:math:`\equiv 1/(dx_2...dx_d)`) for :math:`d` space-time dimensions. If :attr:`in.randoms` is specified in the input, ``v`` has a first dimension of ``in.randoms(1) + in.randoms(2)``. If not specified, the default for ``in.randoms`` is  ``in.noises``. If not specified, the default of :func:`in.initial` is ``a = 0``.

.. attribute:: in.transfer(v,r,a0,r0)

    *Default:* :func:`xtransfer`

    Initializes the fields :math:`a` for subsequent calculations in a sequence. Otherwise, this function behaves in a similar way to :attr:`in.initial`. The function includes the previous field ``a0`` and lattice ``r0``. The default set by :func:`xtransfer` is ``a = a0``.

.. attribute:: in.da(a,z,r)

    *Default:* :func:`xda`

    Calculates derivatives :math:`da` of the equation. The noise vector, ``z``, has variance :math:`1/(dx_{1}..dx_{d})`, for dimension :math:`d \le 4`, and a first dimension  whose default value is :attr:`in.fields` if :attr:`in.noises` are not given. Otherwise, it has a first dimension of ``in.noises(1) + in.noises(2)``. The second type of input noise allows for spatially correlated and filtered noise specified in momentum space.

.. attribute:: in.linear(r)

    *Default:* :func:`xlinear`

    A user-definable function which returns the linear coefficients :math:`L` in Fourier space. This is a function of the differential operator ``Dx``, ``Dy``, ``Dz``, which correspond to :math:`\partial / \partial x`, :math:`\partial / \partial y`, :math:`\partial / \partial z` respectively. Each component has an array dimension the same as the coordinate lattice. If axes are numbered, use  ``D{1}``, ``D{2}``, ``D{3}`` etc.

.. attribute:: in.observe(a,r)

    *Default:* cell array of :func:`xobserve`

    **Cell array** of function handles that take the current field and returns a real observable ``o`` with dimension of ``[1, n.lattice]``. The default observable is the first real field amplitude. Note the use of braces for cell arrays! One can also input these individually as ``in.observe{1} = @(a,r) f(a,r)``, using an inline anonymous function. The total number of observe functions is stored internally as :attr:`in.graphs`. The fields ``a`` passed in the input are transformed according to the :attr:`in.transforms` metadata.

.. attribute:: in.rfilter(r)

    *Default:* :func:`xrfilter`

    Returns the momentum-space filters for the input random terms. Each component has an array dimension the same as the input random fields in momentum space, that is, the return dimension is ``[r.randoms(2), r.nlattice]``.

.. attribute:: in.nfilter(r)

    *Default:* :func:`xnfilter`

    Returns the momentum-space filters for the propagation noise terms. Each component has an array dimension the same as the random noises in momentum space, that is, the return dimension is ``[r.noises(2), r.nlattice]``.


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

    Flag for storing raw trajectory data. If this flag is turned on, raw trajectories are stored in memory. The raw data is returned in function calls and also written to a file on completion, if a file-name is included.

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

        in.ipsteps = 1, 2, 3, ..

.. attribute:: in.file

    *Default:* ``''``

    Matlab or *HDF5* file name for output data. Includes all data and parameter values, including raw trajectories if ``in.raw = 1``. If not needed just omit this. A Matlab filename should end in ``.mat``, while an HDF5 file requires the filename to end in ``.h5``. For a sequence of inputs, the filename should be given in the first structure of the sequence, and the entire sequence is stored.

    ::

        in.file = 'file-name'

Advanced input functions
------------------------

Advanced input functions are user-definable functions which don’t usually need to be changed from default values. They allow customization and extension of xSPDE. These are as follows:

.. function:: in.function(data,in)

    This is a cell array of graphics function handles. Use when a graph is needed that is a function of the observed local averages over ``ensemble(1)``. The default value generates all the averages that are in the simulated data. The input is the data array of averages, and the output is another data array. This function can generate  error-bars and sampling errors in the local averages.

.. attribute:: in.grid(r)

    *Default:* :func:`xgrid`

    Initializes the grid of coordinates in space.

.. attribute:: in.noisegen(r)

    *Default:* :func:`xnoisegen`

    Generates arrays of noise terms ``xi`` for each point in time.

.. attribute:: in.randomgen(r)

    *Default:* :func:`xrandomgen`

    Generates a set of initial random fields ``v`` to initialize the fields simulated.

.. attribute:: in.step(a,z,dt,r)

    *Default:* :func:`xMP`

    Specifies the stochastic integration routine for the field ``a``, given a step in time ``dt`` and noise ``z``, together with the interaction-picture propagator :attr:`r.propagator` which is part of the lattice structure. It returns the new field ``a``. This function can be set to any of the predefined stochastic integration routines provided with xSPDE, described in the :ref:`chap-algorithms` chapter. User-written functions can also be used. The standard method, :func:`xMP`, is a midpoint integrator.

.. attribute:: in.prop(a,r)

    *Default:* :func:`xprop`

    Returns the fields propagated for one step in the interaction picture, depending on the initial field ``a``, and the propagator array :attr:`r.propagator`. Note that the time-step used in :attr:`r.propagator` depends on the input time-step, the error-checking and the algorithm.

.. attribute:: in.propfactor(nc,r)

    *Default:* :func:`xpropfactor`

    Returns the transfer array :attr:`r.propagator`, used by the :attr:`in.prop` function. The time propagated is a fraction of the current integration time-step, :attr:`r.dt`. It is equal to ``1 / in.ipsteps`` of the integration time-step.


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

.. function:: in.gfunction(data,in)

    This is a cell array of graphics function handles. Use when a graph is needed that is a function of the data generated by xsim, as a post-processing step. The default value is the simulated data. The input is the data array of averages from xsim, including error-estimates, and the output is the  data array that is plotted. Any handling of error estimates must be user-provided at this stage.
    
.. function:: in.xfunctions(x_nd,in)

    This is a nested cell array of graphics axis transformations. Use when a graph is needed with an axis that is a function of the original axes. The default value generates all the averages that are in the simulated data. The input of the function is the original axis coordinates, and the output is the new coordinate set. Called as in.xfunctions{n}{nd}(x_nd,in) for the n-th graph and axis nd, where x_nd is a vector of axis coordinate points for that axis dimension.
    

.. function:: in.compare(t,in)

    This is a cell array of comparative functions. Each takes the time or frequency vector - or whichever is the first dimension plotted - and returns comparison results for a graphed observable, as a function versus time or frequency (etc). Comparison results are graphed with a dashed line, for the two-dimensional graphs versus the first plotted dimension. There is no default function handle.

Graphics parameters
~~~~~~~~~~~~~~~~~~~

For uniformity, the graphics parameters that reference an individual data object are cell arrays, indexed over the graph number using braces ``{}``. If a different type of input is used, like a scalar or matrix, xSPDE will attempt to convert the type. The axis labels are cell arrays, indexed over dimension. The graph number used to index these cell arrays refers to the data object, and there can be multiple plots obtained, depending on the graphics input.

Together with default values, they are:

.. attribute:: in.gversion

    *Default:* ``'xGRAPH1.1'``

    This sets the version number of the graphics output. There is typically no need to input this.

    ::

        in.gversion = 'current version name'
        
    
.. attribute:: in.graphs

    *Default:* ``1``

    If non-negative, this sets the maximum number of graphed datasets. Can be used to suppress unwanted graphs from an xSPDE graphics script. If omitted, all the graphs output from the graphics processing function are plotted.
    
    ::

        in.graphs = 1,..
        
.. attribute:: in.olabels

    *Default:* ``{'a', ...}``

    **Cell array** of labels for the graph axis observables and functions. These are text labels that are used on the graph axes. The default value is ``'a_1'`` if the default observable is used, otherwise it is blank. This is overwritten by any subsequent label input when the graphics program is run:

    ::

        in.olabels{n} = 'string'


.. attribute:: in.axes

    *Default:* ``{{0,0,0,..}}``

   Gives the axis and points plotted ``p`` for each plotted function. As special cases,  ``p = 0``, is the default value that gives the entire axis, while  ``p = -1`` generates one point on the axis, namely the last point for the time axis and the midpoint for the space axes. Other values are vector range indicators, for example ``p = 5`` plots the fifth point, while ``p = 1:4:41`` plots every fourth point. For each graph type ``n`` the axes can be individually specified. If more than three axes are specified, only the first three are used. The others are set to default values.

    ::

        in.axes{n} = {p1,p2,p3,..pd}

.. attribute:: in.font

    *Default:* ``{18}``

    This sets the default font size for the graph labels. This can be changed per graph.

    ::

        in.font{n} > 0

.. attribute:: in.minbar

    *Default:* ``{0.01}``

    This is the minimum relative error-bar that is plotted. Set to a large value to suppress unwanted error-bars, although its best not to ignore the error-bar information! 

    ::

        in.minbar{n} >= 0
        
        .. attribute:: in.esample

    *Default:* ``{1}``

    This is the flag for plotting sampling error. Set to zero to suppress unwanted sampling error lines and just plot means, although its best not to ignore this information! 

    ::

        in.esample{n} >= 0

.. attribute:: in.images

    *Default:* ``{0, 0, 0, ...}``

    This is the number of 3D, transverse o-x-y movie images plotted as discrete time slices. Only valid if :attr:`in.dimension` is greater than 2. Note that, if present, the coordinates not plotted are set to their central value, for example ``z = 0``, when plotting the transverse images. This input should have a value from ``in.images(n) = 0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`in.graphs`:

    ::

        in.images{n} = 0 ... in.points(1)

.. attribute:: in.imagetype

    *Default:* ``{1, 1, ...}``

    This is the *type* of transverse o-x-y movie images plotted. If an element is ``1``, a perspective surface plot is output, for ``2``, a gray plot with colours is output, or for ``3`` a contour plot with 10 equally spaced contours is generated. This has a vector length equal to :attr:`in.graphs`.

    ::

        in.imagetype{n} = 1, 2, 3

.. attribute:: in.transverse

    *Default:* ``{0, 0, ...}``

    This is the number of 2D, transverse o-x images plotted as discrete time slices. Only valid if :attr:`in.dimension` is greater than 2. Note that, if present, the y,z-coordinates are set to their central values, when plotting the transverse images. Each element should be from ``0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`in.graphs`:

    ::

        in.transverse{n}=0 ... in.points(1)

.. attribute:: in.headers

    *Default:* ``{'', '', ...}``

    This is a string variable giving the graph headers for each type of function plotted. The default value is an empty string ``''``, which gives the simulation heading. Use a space ``' '`` for no headers. It is useful to include simulation headers - which is the default - to identify graphs in preliminary stages, while they may not be needed in a published final result. 

    ::

        in.headers{n} = 'my_graph_header'

.. attribute:: in.pdimension

    *Default:* ``{3, 3, ...}``

    This is the maximum space-time grid dimension for each plotted quantity. The purpose is eliminate unwanted graphs. For example, it may be useful to reduce the maximum dimension when averaging in space. Higher dimensional graphs are not needed, as the data is duplicated. Averaging can be useful for checking conservation laws, or for averaging over homogeneous data to reduce sampling errors. All graphs are suppressed if it is set to zero. Any three dimensions can be chosen using the axes command.

    ::

        in.pdimension{n} \ge 0 
        
        

.. attribute:: in.xlabels

    *Default:* ``{'t', 'x', 'y', 'z'}`` or ``{'t', 'x_1', 'x_2', 'x_3'}``

    Global labels for the graph axis independent variable labels, vector length of :attr:`in.dimension`. The numerical labeling default is used when the ``in.numberaxis`` option is set. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.xlabels = {in.xlabels(1), ..., in.xlabels(in.dimension)}

.. attribute:: in.klabels

    *Default:* ``{'\\omega', 'k\_x', 'k\_y', 'k\_z'}`` or ``{'\\omega', 'k\_1', 'k\_2', 'k\_3'}``

    Global labels for the graph axis Fourier transform labels, vector length of :attr:`in.dimension`. The numerical labeling default is used when the ``in.numberaxis`` option is set. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.klabels = {in.klabels(1), ..., in.klabels(in.dimension)}
        
        .. attribute:: in.glabels

    *Default:* ``{{'t', 'x', 'y', 'z'}}`` or ``{{'\omega', 'k_x', 'k_y', 'k_z'}}``

    Graph-dependent labels for the independent variable labels, nested cell array with first dimension  :attr:`in.graph`, second dimension :attr:`in.dimension`. 

    ::

        in.glabels{n} = {in.xlabels(1), ..., in.xlabels(in.dimension)}
        
        
         .. attribute:: in.lines

    *Default:* `` {{'-k','--k',':k','-.k','-ok','--ok',':ok','-.ok','-+k','--+k'}}``

    Line types for each line in every two-dimensional graph plotted.

    ::

        in.lines{n} = {linetype{1}, ..., linetype{nl}}
      

Graphics projections
~~~~~~~~~~~~~~~~~~~~

If there is a spatial grid, the graphics program automatically generates several graphs for each observable, depending on space dimension. The maximum dimension that is plotted as set by :attr:`in.pdimension`. In the plots, the lattice is projected down to successively lower dimensions.

For each observable, the projection sequence is as follows:

-  If :attr:`in.dimension` is ``4`` or greater, a central :math:`z` point ``nz = 1 + floor(in.points(4)/2)`` is picked. For example, with 35 points, the central point is ``nz = 18``.

-  This gives a three dimensional space-time lattice, which is treated the same as if :attr:`in.dimension` is ``3``.

-  If :attr:`in.images` are specified, two-dimensional :math:`x-y` plots are generated at equally spaced time intervals. If there is only one image, it is at the last time-point. Different plot-types are used depending on the setting of :attr:`in.imagetype`.

-  A central :math:`y` point ``ny = 1 + floor(in.points(3)/2)`` is picked. This gives a two dimensional space-time lattice, which is treated the same as if :attr:`in.dimension` is ``2``. If :attr:`in.transverse` is specified, one-dimensional :math:`x` plots are generated at equally spaced time intervals, as before.

-  A central :math:`x` point ``nx = 1 + floor(in.points(2)/2)`` is picked. This gives a one dimensional time lattice, which is treated the same as if :attr:`in.dimension` is ``1``.

-  Plots of observable vs time are obtained, including sampling errors and error bars. If comparison graphs are specified using :func:`in.compare` functions, they are plotted also, using a dotted line. A difference graph is also plotted when there is a comparison.

Parameter structure
===================

Internally, xSPDE parameters are stored in a cell array, ``latt``, of structures ``r``, which is passed to functions. This includes all the data given above inside the ``in`` structure. In addition, it includes the table of computed parameters given below.

User application constants and parameters should not be reserved names; :attr:`in.c` and all names starting with ``in.c`` will always be available in all versions of xSPDE.

A parameter structure contains information about the space-time grid and is passed to various functions, for instance :attr:`in.da` or :attr:`in.step`. The corresponding parameter is commonly marked as `r`.

.. attribute:: r.t

    Current value of time, :math:`t`.

.. attribute:: r.x

.. attribute:: r.y

.. attribute:: r.z

    Coordinate grids of :math:`x`, :math:`y`, :math:`z`.

.. attribute:: r.x{n}

    Higher dimensions are labeled numerically as :math:`x_1`,..  :math:`x_6`, and so on. This numerical axis convention can be set even for lower dimensions if ``in.numberaxis`` is set to 1.

.. attribute:: r.kx

.. attribute:: r.ky

.. attribute:: r.kz

    Grids in momentum space: :math:`k_x`, :math:`k_y`, :math:`k_z`.

.. attribute:: r.k5

    Higher dimensions are labeled numerically as :math:`k_5`,  :math:`k_6`, and so on.

.. attribute:: r.dt

    Output time-step between stored points for data averages.
    
.. attribute:: r.dtr

    Current reduced time-step used for integration.

.. attribute:: r.dx

    Steps in coordinate space: :math:`[t,x,y,z,x_5,..]`.

.. attribute:: r.dk

    Steps in momentum space: :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]`.

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

    Computational Fourier transform axes in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},.. ]` (vector cells).

.. attribute:: r.kg

    Graphics  Fourier transform axes in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]` (vector cells).

.. attribute:: r.kranges

    Range in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]` (vector).

.. attribute:: r.s.dx

    Initial stochastic normalization.

.. attribute:: r.s.dxt

    Propagating stochastic normalization.

.. attribute:: r.s.dk

    Initial :math:`k` stochastic normalization.

.. attribute:: r.s.dkt

    Propagating :math:`k` stochastic normalization.

.. attribute:: r.nspace

    Number of spatial lattice points: ``in.points(2) * .. * in.points(in.dimension)``.

.. attribute:: r.nlattice

    Total lattice: ``in.ensembles(1) * r.nspace``.

.. attribute:: r.ncopies

    Total copies of stochastic integrations: ``in.ensembles(2) * in.ensembles(3)``.

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

    Returns the real part of ``a(1,:)``.

.. function:: xrfilter(r)

    Returns an array of ones.

.. function:: xnfilter(r)

    Returns an array of ones.

.. function:: xgrid(r)

    Sets grid points in lattice from coordinate vectors. Returns the ``r`` structure with added grid points.

.. function:: xnoisegen(r)

    Generates random noise matrix :math:`z`.

.. function:: xrandomgen(r)

    Generates initial random field matrix :math:`v`.

.. function:: xpropfactor(nc, r)

    Returns the interaction picture propagation factor. ``nc`` is a check index, ``r`` is a lattice structure.


Frequently asked questions
==========================

Answers to some frequent questions, and reminders of points in this chapter are:

-  Can you average other stochastic quantities apart from the field?

   -  Yes: just specify the functions that need to be averaged using the user function :attr:`in.observe`.

-  Can you have functions of the current time and space coordinate?

   -  Yes: xSPDE functions support this using the structure ``r``, as :attr:`r.t`, :attr:`r.x`, :attr:`r.y`, :attr:`r.z`, or  :attr:`r.t`, ``r.x{1}``, and so on, for more than four space-time dimensions.

-  Can you have several independent stochastic variables?

   -  Yes, input this using ``in.fields > 1``.

-  Are higher dimensional differential equations possible?

   -  Yes, this requires setting ``in.dimension > 1``. This is essentially unlimited in xSPDE except for memory requirements.

-  Can you have spatial partial derivatives?

   -  Yes, provided they are linear in the fields; these are obtainable using the function :attr:`in.linear`.

-  Can you delete the graph heading?

   -  Yes, this is turned off if you set :attr:`in.headers` to ``0``.

-  Why are there two lines in the graphs sometimes?

   -  These are one standard deviation sampling error limits, generated when ``in.ensembles(2,3) > 1``.

-  Why is there just one line in some graphs, with no sampling errors indicated?

   -  You need ``in.ensembles(2)`` or ``(3)`` for this; see previous question.

-  What are the error bars for?

   -  These are the estimated maximum errors due to finite step-sizes in time.
