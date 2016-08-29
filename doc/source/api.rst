.. _chap-api:

**********
Public API
**********


High-level xSPDE functions and objects
======================================

The high-level xSPDE functions process the input parameters, producing simulation data and graphs. Function parameters are in the order of ``input``, which specifies the simulation, then ``data``, which is needed for graphics output.

.. function:: xspde (input)

    This is the combined xSPDE function. It accepts a simulation sequence, ``input``. This can be a single structure, ``in``, or else a cell array of structures, ``{in1,in2,..}``, for  sequences. Output graphs are displayed. It returns the output ``[error, input, data,rawdata]``, where ``error`` is the sum of simulation errors in :func:`xsim`, and difference errors found in the :func:`xgraph` comparisons. If a filename is specified in  the input, it writes an output data file including input and all output data. Raw data is stored on request. It calls the functions :func:`xsim` and :func:`xgraph`.


.. function:: xsim (input)

    This is the xSPDE simulation function. Like :func:`xspde`, it accepts input parameters in ``input``. It returns ``[maxerror, input, data, rawdata]``, where: ``maxerror`` is the sum of maximum step-size and maximum sampling errors, ``input`` is the full input structure or cell array for sequences, including default values, and ``data`` is a cell array of average observables. If the ``in.raw`` option is used, data for the actual trajectories is output in ``rawdata``. This can be run as a stand-alone function if no graphs are required.

.. function:: xgraph (data [,input])

    This is the xSPDE graphics function. It takes computed simulation  ``data`` and ``input``. It plots graphs, and returns the maximum difference ``diff`` from comparisons with user-specified comparison functions. The ``data`` should have as many cells as ``input`` cells, for sequences. 
    If ``data = 'filename.h5'`` or ``data= 'filename.mat'``, the specified file is read both for ``input`` and ``data``. Here ``.h5`` indicates an HDF5 file, and ``.mat`` indicates a Matlab file.
    When the ``data`` input is given as a filename, input parameters in the file are replaced by any of the the new ``input`` parameters that are specified.  Any stored ``input`` can be overwritten when the graphs are generated. This allows graphs of data to be modified retrospectively.



Open object-oriented architecture
----------------------------------

As well as extensibility through sequences, which was described in :ref:`chap-projects`, in the section :ref:`sec-sequential-integration`, the architecture of xSPDE allows functional extensions.

The input metadata includes both data and methods acting on the data, in the tradition of object-oriented programs. Yet there is no strict class typing. Users are encouraged to adapt the xSPDE program by adding input parameters and methods to the input structures.

This open object orientation is deliberate. An open extensibility permits arbitrary functions to be specified in the ``in`` structures. All functions and parameters have default values in xSPDE. It is also possible to include user defined functions, provided they satisfy the API definitions given below. This is achieved simply by including the relevant function handles and parameters in the input metadata.

Internal parameters and function handles from the ``in`` structures, together with any computed and default parameters and functions,  are stored in the :ref:`lattice structure <sec-parameter-structure>` ``r``. These are available to all user-defined functions. Use of pre-existing reserved names is not advisable, and the structure ``r.c`` is always reserved for user constants, parameters and functions if required.

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
    
These low-level functions are described in detail in the next sections of this chapter.


Low-level xSPDE functions
=========================

The xSPDE program is function oriented: low-level functions are used to define initial conditions, equations and observables, amongst many other things given in detail below. To make this process simpler, argument passing conventions are used. Common parameters are always passed using a lattices structure variable ``r`` as the last argument.

Functions of a single lattice have arguments in the following order:

-  the field ``a`` or initial random variable ``v``;
-  the stochastic noise ``z`` or other fields;
-  non-field arguments;
-  the lattice structure ``r``.

The first argument, ``a``, is a real or complex vector field. This is a matrix whose first dimension is the field index, with a range of 1 to :attr:`fieldsplus`, with :attr:`fieldsplus` = :attr:`fields` (1) + :attr:`fields` (2). The second dimension is the lattice index.

The second argument, ``z``, if needed, is a real random noise, corresponding to :math:`\zeta` in the mathematical notation. This is a matrix whose first dimension is the noise index. The second dimension is the lattice index.

The last function argument is the  :ref:`lattice structure <sec-parameter-structure>`, ``r``. This contains data about the integration details and lattice, stored as ``r.name``. Important constants in the structure are :attr:`t`, the current time, and  space coordinates, :attr:`x`, :attr:`y`, :attr:`z`. Other data stored in the structure is explained in later chapters.

Functions of multiple lattice sequences take current arguments first, and the oldest arguments last.

Integration arrays
------------------

In all integration function calls, the variables used are matrices. The first dimension used is the stochastic field length :attr:`fields` (1). The second dimension in all field arrays is the lattice index, with a length ``n.lattice = ensembles(1) * points(2) * ... * points(dimension)``. Here ``ensembles(1)`` is the number of stochastic samples integrated as an array.

The field dimensions for the flattened arrays passed to xSIM integration functions are:

- ``a = [r.fieldsplus, r.nlattice]``
- ``rv = [r.randoms(1)+r.randoms(2), r.nlattice]``
- ``w = [r.noises(1)+r.noises(2), r.nlattice]``
- ``r.Dx, r.x, r.kx = [1, r.nlattice]``

Data arrays
-----------

Each observable used to generate graph data is defined by a function in a cell array with length :attr:`graphs`. There are two stages of averaging. First, an average over a local ensemble at a single time-point is performed using the  :func:`observe` function. Next, if more sophisticated data is required, an optional  :func:`function` is used to transform data.

The first dimension ``lines`` is initially determined by the :func:`observe` function. This can be modifed if required  by the data transformation  :func:`function`. It is typically one for a single-line graph, but can be greater. The last dimensions in all data arrays is the vector of time-space dimensions: ``points = [points(1), ... ,points(dimension)]``. 

- ``d{n} = [lines,1, points]``.

If the optional :func:`function` method is used to transform data within xSIM, the entire average data cell array from every :func:`observe` function is passed after local averaging, to allow all transformations. On output from xSIM to xGRAPH, the data arrays are augmented by the addition of error estimates, addressed using the second index. 


Simulation parameters
---------------------

For each simulation in the ``input`` sequence, the input parameters and functions are specified as a data structure, ``in``. These can be entered either interactively or as part of a simulation function file. The function file approach allows recycling and editing, so it is better for a large project.

There are extensive default preferences to simplify the inputs. If any inputs are omitted, there are default values which are set by inpreferences in all cases. These defaults are changed by editing the inpreferences function. The :func:`xgpreferences` function is used to supply graphics default values.

**For vector or cell inputs, an input shorter than required is padded to the right using default values.**



.. _sec-parameters:

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


xSIM parameters
----------------



.. attribute:: version

    *Default:* ``'xSIM2.3'``

    This sets the current version number of the  simulation program. There is typically no need to input this.

    ::

        in.version = 'current version name'
        

.. attribute:: name

    *Default:* ``' '``

    Name used to label simulation, usually corresponding to the equation or problem solved. This can be added or removed from graphs using the :attr:`headers` Boolean variable, as explained in the section on graphics parameters.

    ::

        in.name = 'your project name'

.. attribute:: dimension

    *Default:* ``1``

    The total space-time dimension is labelled, unsurprisingly,

    ::

        in.dimension = 1...4

.. attribute:: fields

    *Default:* ``[1,0]``

    These are real or complex variables stored at each lattice point, and are the independent variables for integration. The fields are vectors that can have any dimension. The first number is the number of real or complex fields that are initialized by the :func:`initial` function and integrated using the :func:`da` derivative. The optional second number is the number of real or complex auxiliary fields specified with the :func:`define` function.

    ::

        in.fields(1,2) = 0, 1, 2, ...
        
.. attribute:: fieldsplus

    *Default:* ``1``

    This is the total of stochastic plus defined fields. This is calculated internally: :attr:`fieldsplus` = :attr:`fields` (1) + :attr:`fields` (2).

    ::

        in.fieldsplus = 0, 1, 2, ...




.. attribute:: noises

    *Default:* :attr:`fields` (1)

    This gives the number of stochastic noises generated per lattice point, in coordinate and momentum space respectively. Set to zero (``in.noises = 0``) for no noises. This is the number of *rows* in the noise-vector. Noises can be delta-correlated or correlated in space. The second input is the dimension of noises in k-space. It can be left out if zero. This allows use of finite correlation lengths when needed, by including a frequency filter function that is used to multiply the noise in Fourier-space. The Fourier-space noise variance is the square of the filter function. Note that the first noise index, noises(1), indicates how many independent noise fields are generated, while noises(2) indicates how many of these are are fourier-transformed, filtered and then inverse fourier transformed to give correlations. These appear as extra noises, so the total is noises(1)+noises(2). The filtered noises have a finite correlation length. They are also correlated with the first noises(2) noises they are generated from. 

    ::

        in.noises = [in.noises(1), in.noises(2)] >= 0.


.. attribute:: randoms

    *Default:* :attr:`noises`

    This gives the number of random fields generated per lattice point for the initial noise, in coordinate and momentum space. Set to zero (``in.randoms = 0``) for no random fields. Random fields can be delta-correlated or correlated in space. The second input is the dimension of random fields in momentum space. It can be left out if zero. The Fourier-space random variance is the square of the filter function. Note that the first noise index, in.randoms(1), indicates how many independent random fields are generated, while in.randoms(2) indicates how many of these are are fourier-transformed, filtered and then inverse fourier transformed. These appear as additional random fields, so the total is in.randoms(1)+in.randoms(2). The filtered noises have a finite correlation length. They are correlated with the first in.randoms(2) random fields they are generated from, just as with the noise terms. 

    ::

        in.randoms = [in.randoms(1), in.randoms(2)] >= 0

.. attribute:: ranges

    *Default:* ``[10, 10, ...]``

    Each lattice dimension has a coordinate range, given by:

    ::

        in.ranges = [in.ranges(1), ..., in.ranges(dimension)]

    In the temporal graphs, the first coordinate is plotted over ``0:in.ranges(1)``. All other coordinates are plotted over ``-in.ranges(n)/2:in.ranges(n)/2``. The default value is ``10`` in each dimension.

.. attribute:: points

    *Default:* ``[51, 35, ..., 35]``

    The rectangular lattice of points plotted for each dimension are defined by a vector giving the number of points in each dimension:

    ::

        in.points = [in.points(1), ..., in.points(in.dimension)]

    The default values are simply given as a rough guide for initial calculations. Large, high dimensional lattices take more time to integrate. Increasing :attr:`points` improves graphics resolution, and gives better accuracy in each relevant dimension as well, but requires more memory. Speed is improved when the lattice points are a product of small prime factors.

.. attribute:: steps

    *Default:* ``1``

    Number of time-steps per plotted point. The total number of integration steps in a simulation is therefore ``in.steps * (in.points(1)-1)``. Thus, :attr:`steps` can be increased to improve the accuracy, but gives no change in graphics resolution. **Increase** steps to give a **lower** time-discretization error:

    ::

        in.steps = 1, 2, ...

.. attribute:: ensembles

    *Default:* ``[1, 1, 1]``

    Number of independent stochastic trajectories simulated. This is specified in three levels to allow maximum parallelism. The first gives within-thread parallelism, allowing vector instructions. The second gives a number of independent trajectories calculated serially. The third gives multi-core parallelism, and requires the Matlab parallel toolbox. Either ``in.ensembles(2)`` or ``in.ensembles(3)`` are required if sampling error-bars are to be calculated.

    ::

        in.ensembles = [in.ensembles(1), in.ensembles(2), in.ensembles(3)] >= 1

    The *total* number of stochastic trajectories or samples is ``ensembles(1) * ensembles(2) * ensembles(3)``.

    
.. attribute:: boundaries

    *Default:* ``[0, 0, ...]``

    Type of spatial boundary conditions used, set for each dimension independently, and used in the partial differential equation solutions. The default option, or ``0``, is periodic. If ``1``,  Neumann boundaries are used, with normal derivatives set to zero.  If ``2``,  Dirichlet boundaries are used, with field values set to zero. Note that in the current xSPDE code, setting non-periodic boundaries requires the use of finite difference type derivatives, without the option of an interaction picture derivative. Using Fourier derivatives will automatically make the boundary conditions periodic.

    ::

        in.boundaries = [0, in.boundaries(2), in.boundaries(3)..] >= 0

    Dimensions for setting the boundary conditions are numbered starting from the time dimension, for consistency with numbering conventions elsewhere. However, only the space dimension boundaries are used here, for :math:`j > 1`.

.. attribute:: transforms

    *Default:* ``{0}``

    **Cell array** that defines the different transform spaces used to calculate field observables. This has the structure

    ::

        in.transforms{n} = [t(1), ..., t(4)] >= 0

    There is one transform vector per observable. The ``j``-th index, ``t(j)``, indicates a Fourier transform on the ``j``-th axis. The normalization of the Fourier transform is such that the :math:`k=0` value in momentum space corresponds to the integral over space, with an additional factor of :math:`1/\sqrt{2\pi}`. This gives a Fourier integral which is symmetrically normalized in ordinary and momentum space. The Fourier transform that is graphed is such that
    :math:`k=0` is the *central* value.

.. attribute:: olabels

    *Default:* ``{'a_1', ...}``

    **Cell array** of labels for the graph axis observable functions. These are text labels that are used on the graph axes. The default value is ``'a_1'`` if the default observable is used, otherwise it is blank. This is overwritten by any subsequent label input when the graphics program is run:

    ::

        in.olabels{n} = 'string'

.. attribute:: c

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


Advanced input parameters
-------------------------

More advanced input parameters, which don’t usually need to be changed from default values, are as follows:

.. attribute:: iterations

    *Default:* ``4``

    For iterative algorithms like the implicit midpoint method, the iteration count is set here, typically around 3-4. Will increase the integration accuracy if set higher, but it may be better to increase :attr:`steps` if this is needed. With non-iterated algorithms, this input is not used:

    ::

        in.iterations = 1, 2, ...

.. attribute:: checks

    *Default:* ``1``

    This defines how many times the integration is carried out for error-checking purposes. If :attr:`checks` is `0`, there is one integration, but no checking at smaller time-steps. For error checking, set ``in.checks = 1``, which repeats the calculation at a shorter time-step --- but with identical noise --- to obtain the error bars, taking three times longer overall:

    ::

        in.checks = 0, 1

.. attribute:: order

    *Default:* ``1``

    This is the extrapolation order, which is **only** used if ``in.checks = 1``. The program uses the estimated convergence order to extrapolate to zero step-size, with reduced estimated error-bars. If ``in.order = 0``, no extrapolation is used, which is the most conservative input. The default order is usually acceptable, especially when combined with the default midpoint algorithm, see next section. While any non-negative order can be input, the theoretical orders of the four preset methods used *without* stochastic noise terms are: ``1`` for :func:`xEuler`; ``2`` for :func:`xRK2`; ``2`` for :func:`xMP`; ``4`` for :func:`xRK4`. Allowed values are:

    ::

        in.order >= 0

.. attribute:: seed

    *Default:* ``0``

    Random noise generation seed, for obtaining reproducible noise sequences. Only needed if ``in.noises > 0``

    ::

        in.seed >= 0

.. attribute:: graphs

    *Default:* number of observables

    This gives the number of observables computed. The default is the length of the cell array of observe functions. Normally, this is not initialized, as the default is typically used. Can be used to suppress data averaging.

    ::

        in.graphs >= 0
        
.. attribute:: functions

    *Default:* number of functional transformations

    This gives the number of graphs computed, which are functions of the observables. The default is the length of the cell array of observe functions. Normally, this is not initialized, as the default is typically used. 

    ::

        in.functions >= 0


.. attribute:: print

    *Default:* ``1``

    Print flag for output information while running xSPDE. If ``print = 0``, most output is suppressed, while ``print = 1`` displays a progress report, and ``print = 2`` also generates a readable summary of the ``r`` lattice structure as a record.

    ::

        in.print >= 0

.. attribute:: raw

    *Default:* ``0``

    Flag for storing raw trajectory data. If this flag is turned on, raw trajectories are stored in memory. The raw data is returned in function calls and also written to a file on completion, if a file-name is included.

    ::

        in.raw >= 0

.. attribute:: origin

    *Default:* ``[0, -in.ranges/2]``

    This displaces the graph origin for each simulation to a user-defined value. If omitted, all initial times in a sequence are zero, and the space origin is set to ``-in.ranges/2`` to give results that are symmetric about the origin:

    ::

        in.origin = [origin(1), ..., origin(4)]

.. attribute:: ipsteps

    *Default:* ``1`` for :func:`xEuler` and :func:`xRK2`, ``2`` for :func:`xMP` and :func:`xRK4`

    This specifies the number of interaction picture steps needed in a full propagation time-step. Default values are chosen according to the setting of :func:`step`. Can be changed for custom integration methods.

    ::

        in.ipsteps = 1, 2, 3, ..

.. attribute:: file

    *Default:* ``''``

    Matlab or *HDF5* file name for output data. Includes all data and parameter values, including raw trajectories if ``in.raw = 1``. If not needed just omit this. A Matlab filename should end in ``.mat``, while an HDF5 file requires the filename to end in ``.h5``. For a sequence of inputs, the filename should be given in the first structure of the sequence, and the entire sequence is stored.

    ::

        in.file = 'file-name'


.. _sec-functions:

xSIM functions
===============

The structure of xsim makes use of many functions, some of which are internal, and some user supplied. This the the main mechanism for extensibility.

Input functions
---------------

A stochastic equation solver requires the definition of an initial distribution and a time derivative. In xSPDE, the time derivatives is divided up into a linear term including space derivatives, used to define an interaction picture, and the remaining derivatives. In addition, one must define quantities to be averaged over during the simulation, called graphs in xSPDE. These are all defined as functions, specified below.

.. function::  initial (rv,r)

    *Default:* :func:`xinitial`

    Initializes the fields :math:`a` for the first simulation in a sequence. The initial Gaussian random field variable, ``rv``, has unit variance if :attr:`dimension` is ``1`` or else is delta-correlated in space, with variance ``1/r.dv`` (:math:`\equiv 1/(dx_2...dx_d)`) for :math:`d` space-time dimensions. If :attr:`randoms` is specified in the input, ``rv`` has a first dimension of ``randoms(1) + randoms(2)``. If not specified, the default for ``randoms`` is  ``noises``, and the default of :func:`initial` is ``a = 0``.

.. function:: transfer(rv,r,a0,r0)

    *Default:* :func:`xtransfer`

    Initializes the fields :math:`a` for subsequent calculations in a sequence. Otherwise, this function behaves in a similar way to :func:`initial`. The function includes the previous field ``a0`` and lattice ``r0``. The default set by :func:`xtransfer` is ``a = a0``.

.. function::  da (a,w,r)

    *Default:* :func:`xda`

    Calculates derivatives :math:`da` of the equation. The noise vector, ``w``, has variance :math:`1/(dx_{1}..dx_{d})`, for dimension :math:`d \le 4`, and a first dimension  whose default value is :attr:`fields` if :attr:`noises` are not given. Otherwise, it has a first dimension of ``in.noises(1) + in.noises(2)``. The second type of input noise allows for spatially correlated and filtered noise specified in momentum space.
    
    
    .. function::  define (a,w,r)

    *Default:* :func:`xdefine`

    Calculates auxiliary field values during propagation.

.. function:: linear (r)

    *Default:* :func:`xlinear`

    A user-definable function which returns the linear coefficients :math:`L` in Fourier space. This is a function of the differential operator ``Dx``, ``Dy``, ``Dz``, which correspond to :math:`\partial / \partial x`, :math:`\partial / \partial y`, :math:`\partial / \partial z` respectively. Each component has an array dimension the same as the coordinate lattice. If axes are numbered, use  ``D{1}``, ``D{2}``, ``D{3}`` etc.

.. function:: observe (a,r)

    *Default:* cell array of :func:`xobserve`

    **Cell array** of function handles that take the current field and returns a real observable ``o`` with dimension of ``[1, n.lattice]``. The default observable is the first real field amplitude. Note the use of braces for cell arrays! One can also input these individually as ``in.observe{1} = @(a,r) f(a,r)``, using an inline anonymous function. The total number of observe functions is stored internally as :attr:`graphs`. The fields ``a`` passed in the input are transformed according to the :attr:`functions` metadata.

.. function::  rfilter (r)

    *Default:* :func:`xrfilter`

    Returns the momentum-space filters for the input random terms. Each component has an array dimension the same as the input random fields in momentum space, that is, the return dimension is ``[r.randoms(2), r.nlattice]``.

..function:: nfilter (r)

    *Default:* :func:`xnfilter`

    Returns the momentum-space filters for the propagation noise terms. Each component has an array dimension the same as the random noises in momentum space, that is, the return dimension is ``[r.noises(2), r.nlattice]``.


Advanced input functions
------------------------

Advanced input functions are user-definable functions which don’t usually need to be changed from default values. They allow customization and extension of xSPDE. These are as follows:

.. function:: xave (o, [av, ] r)

    This function takes a vector or scalar field or observable, for example ``o = [1, n.lattice]``, defined on the xSPDE local lattice, and returns an average over the spatial lattice with the same dimension. The input is a field or observable ``o``, and an optional averaging switch ``av``. If ``av(j) > 0``, an average is taken over dimension ``j``. Space dimensions are labelled from ``j = 2 ... `` as elsewhere.  If the ``av`` vector is omitted, the average is taken over all space directions.  To average over the local ensemble and all space dimensions, use ``xave(o)``. Averages are returned at all lattice locations.
    
.. function:: xint (o, [dx, ] r)

    This function takes a scalar or vector quantity ``o``, and returns a  space integral over selected dimensions with vector measure ``dx``. If ``dx(j) > 0`` an integral is taken over dimension ``j``. Space dimensions are labelled from ``j = 2, ...`` as elsewhere. Time integrals are ignored at present.  To integrate over an entire lattice, set ``dx = r.dx``, otherwise set ``dx(j) = r.dx(j)`` for selected dimensions ``j``.  Integrals are returned at all lattice locations.


.. function:: xd (o, [D, ] r)

    This function takes a scalar or vector quantity ``o``, and returns a spectral derivative over selected dimensions with a derivative ``D``, by Fourier transforming the data.  Set ``D = r.Dx`` for a first order x-derivative, ``D = r.Dy`` for a first order y-derivative, and similarly ``D = r.Dz.*r.Dy`` for a cross-derivative in ``z`` and ``y``. Higher derivatives require powers of these, for example `D = r.Dz.^4``. For higher dimensions use numerical labels, where ``D = r.Dx`` becomes ``D = r.D{2}``, and so on. If the derivative ``D`` is omitted, a first order x-derivative is returned.

.. function:: xd1 (o, [dir, ] r)

    This takes a scalar or vector ``o``, and returns a first derivative with an axis direction ``dir`` using finite differences.  Set ``dir = 2`` for an x-derivative, ``dir = 3`` for a y-derivative.  Time derivatives are ignored at present. Derivatives are returned at all lattice locations. The boundary condition is set by the in.boundaries input. It can be made periodic, which is the default, or Neumann with zero derivative, or Dirichlet with zero field.

.. function:: xd2 (o, [dir, ] r)

	This takes a scalar or vector ``o``, and returns the second  derivative in axis direction ``dir``.  Set ``dir = 2`` for an x-derivative, ``dir = 3`` for a y-derivative.  All other properties are exactly the same as :func:`xd1`.


.. function:: function (data,in)

    This is a cell array of data function handles. Use when simulation data is needed that is a function of the :func:`observe` local averages over ``ensemble(1)``. The default value simply generates all the averages that are in the simulated data. The input to the ``n``-th function is the cell array of averages, and the output is a data array for the ``n``-th graph. This function is used at simulation time, and  generates both  error-bars and sampling errors in the graphed results.

.. function:: grid (r)

    *Default:* :func:`xgrid`

    Initializes the grid of coordinates in space.

.. function:: noisegen (r)

    *Default:* :func:`xnoisegen`

    Generates arrays of noise terms ``xi`` for each point in time.

.. function:: randomgen (r)

    *Default:* :func:`xrandomgen`

    Generates a set of initial random fields ``v`` to initialize the fields simulated.

.. function:: step (a,w,dt,r)

    *Default:* :func:`xRK4`

    Specifies the stochastic integration routine for the field ``a``, given a step in time ``dt`` and noise ``w``, together with the interaction-picture propagator :attr:`propagator` which is part of the lattice structure. It returns the new field ``a``. This function can be set to any of the predefined stochastic integration routines provided with xSPDE, described in the :ref:`chap-algorithms` chapter. User-written functions can also be used. The standard method, :func:`xRK4`, is a fourth-order Runge-Kutta. Another very useful alternative, :func:`xMP`, is a midpoint integrator.

.. function:: prop (a,r)

    *Default:* :func:`xprop`

    Returns the fields propagated for one step in the interaction picture, depending on the initial field ``a``, and the propagator array :attr:`propagator`. Note that the time-step used in :attr:`propagator` depends on the input time-step, the error-checking and the algorithm.

.. function:: propfactor (nc,r)

    *Default:* :func:`xpropfactor`

    Returns the transfer array :attr:`propagator`, used by the :attr:`prop` function. The time propagated is a fraction of the current integration time-step, :attr:`dt`. It is equal to ``1 / ipsteps`` of the integration time-step.



.. _sec-gparameters:

xGRAPH parameters
=================

The graphics parameters are also stored in the cell array ``input`` as a sequence of structures ``in``. This only need to be input when the graphs are generated, and can be changed at a later time to alter the graphics output. A sequence of simulations is graphed from ``input`` specifications.

If there is one simulation, just one structure can be input, without the sequence braces. The standard form of each parameter value, which should have the ``in.`` structure label added, is:

::

    in.label = parameter

If any inputs are omitted, there are default values which are set by the :func:`xgpreferences` function, in all cases except for the comparison function :func:`compare`. The defaults can be changed by editing the :func:`xgpreferences` function.

In the following descriptions, :attr:`graphs` is the total number of graphed variables of all types. The space coordinate, image, image-type and transverse data can be omitted if there is no spatial lattice, that is, if the dimension variable is set to one.

For uniformity, the graphics parameters that reference an individual data object are cell arrays, indexed over the graph number using braces ``{}``. If a different type of input is used, like a scalar or matrix, xSPDE will attempt to convert the type. The axis labels are cell arrays, indexed over dimension. The graph number used to index these cell arrays refers to the data object, and there can be multiple plots obtained, depending on the graphics input.

Together with default values, they are:

.. attribute:: gversion

    *Default:* ``'xGRAPH2.3'``

    This sets the current version number of the graphics program. There is typically no need to input this.

    ::

        in.gversion = 'current version name'
        
    
.. attribute:: graphs

    *Default:* ``in.functions``

    If specified, this sets the maximum number of graphed datasets. Can be used to suppress unwanted graphs from an xSPDE graphics script. If omitted, all the data output from the in.functions data processing functions are plotted.
    
    ::

        in.graphs = 1,..
        
.. attribute:: olabels

    *Default:* ``{'a', ...}``

    **Cell array** of labels for the graph axis observables and functions. These are text labels that are used on the graph axes. The default value is ``'a_1'`` if the default observable is used, otherwise it is blank. This is overwritten by any subsequent label input when the graphics program is run:

    ::

        in.olabels{n} = 'string'

.. attribute:: axes

    *Default:* ``{{0,0,0,..}}``

   Gives the axis and points plotted ``p`` for each plotted function. As special cases,  ``p = 0``, is the default value that gives the entire axis, while  ``p = -1`` generates one point on the axis, namely the last point for the time axis and the midpoint for the space axes. Other values are vector range indicators, for example ``p = 5`` plots the fifth point, while ``p = 1:4:41`` plots every fourth point. For each graph type ``n`` the axes can be individually specified. If more than three axes are specified, only the first three are used. The others are set to default values.

    ::

        in.axes{n} = {p1,p2,p3,..pd}

.. attribute:: font

    *Default:* ``{18, ...}``

    This sets the default font sizes for the graph labels, indexed by graph. This can be changed per graph.

    ::

        in.font{n} > 0

.. attribute:: minbar

    *Default:* ``{0.01, ...}``

    This is the minimum relative error-bar that is plotted. Set to a large value to suppress unwanted error-bars, although its best not to ignore the error-bar information! This can be changed per graph.

    ::

        in.minbar{n} >= 0
        
        .. attribute:: esample

    *Default:* ``{1, ...}``

    This is the flag for plotting sampling error. Set to zero to suppress unwanted sampling error lines and just plot means, although its best not to ignore this information! This can be changed per graph.

    ::

        in.esample{n} >= 0

.. attribute:: images

    *Default:* ``{0, 0, 0, ...}``

    This is the number of 3D, transverse o-x-y movie images plotted as discrete time slices. Only valid if :attr:`dimension` is greater than 2. Note that, if present, the coordinates not plotted are set to their central value, for example ``z = 0``, when plotting the transverse images. This input should have a value from ``in.images(n) = 0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`graphs`:

    ::

        in.images{n} = 0 ... in.points(1)

.. attribute:: imagetype

    *Default:* ``{1, 1, ...}``

    This is the *type* of transverse o-x-y movie images plotted. If an element is ``1``, a perspective surface plot is output, for ``2``, a gray plot with colours is output, or for ``3`` a contour plot with 10 equally spaced contours is generated. This has a vector length equal to :attr:`graphs`.

    ::

        in.imagetype{n} = 1, 2, 3

.. attribute:: transverse

    *Default:* ``{0, 0, ...}``

    This is the number of 2D, transverse o-x images plotted as discrete time slices. Only valid if :attr:`dimension` is greater than 2. Note that, if present, the y,z-coordinates are set to their central values, when plotting the transverse images. Each element should be from ``0`` up to a maximum value of the number of plotted time-points. It has a vector length equal to :attr:`graphs`:

    ::

        in.transverse{n}=0 ... in.points(1)

.. attribute:: headers

    *Default:* ``{'head1', 'head2', ...}``

    This is a string variable giving the graph headers for each type of function plotted. The default value is an empty string ``''``, which gives the overall simulation heading. Use a space ``' '`` to suppress graphics headers entirely. It is useful to include simulation headers - which is the default - to identify graphs in preliminary stages, while they may not be needed in a final result. 

    ::

        in.headers{n} = 'my_graph_header'

.. attribute:: pdimension

    *Default:* ``{3, 3, ...}``

    This is the maximum space-time grid dimension for each plotted quantity. The purpose is eliminate unwanted graphs. For example, it may be useful to reduce the maximum dimension when averaging in space. Higher dimensional graphs are not needed, as the data is duplicated. Averaging can be useful for checking conservation laws, or for averaging over homogeneous data to reduce sampling errors. All graphs are suppressed if it is set to zero. Any three dimensions can be chosen using the axes command.

    ::

        in.pdimension{n} \ge 0 
        
        

.. attribute:: xlabels

    *Default:* ``{'t', 'x', 'y', 'z'}`` or ``{'x_1', 'x_2', 'x_3', 'x_4'}``

    Labels for the graph axis independent variable labels, vector length of :attr:`dimension`. The numerical labeling default is used when the ``in.numberaxis`` option is set. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.xlabels = {in.xlabels(1), ..., in.xlabels(in.dimension)}

.. attribute:: klabels

    *Default:* ``{'\\omega', 'k\_x', 'k\_y', 'k\_z'}`` or ``{'k\_1', 'k\_2', 'k\_3', 'k\_4'}``

    Labels for the graph axis Fourier transform labels, vector length of :attr:`dimension`. The numerical labeling default is used when the ``in.numberaxis`` option is set. *Note, these are typeset in Latex mathematics mode!*

    ::

        in.klabels = {in.klabels(1), ..., in.klabels(in.dimension)}

.. attribute:: glabels

    *Default:* ``{{'t', 'x', 'y', 'z'}}`` or ``{{'\omega', 'k_x', 'k_y', 'k_z'}}``

    Graph-dependent labels for the independent variable labels, nested cell array with first dimension  :attr:`graphs`, second dimension :attr:`dimension`. 

    ::

        in.glabels{n} = {in.xlabels(1), ..., in.xlabels(in.dimension)}
        
        
         .. attribute:: lines

    *Default:* `` {{'-k','--k',':k','-.k','-ok','--ok',':ok','-.ok','-+k','--+k'}}``

    Line types for each line in every two-dimensional graph plotted.

    ::

        in.lines{n} = {linetype{1}, ..., linetype{nl}}
      


.. _sec-gfunctions:

xGRAPH functions
================

.. function:: gfunction (data,in)

    This is a cell array of graphics function handles. Use when a graph is needed that is a functional transformation of the observed averages. The default value generates all the averages that are in the simulated data. The input is the data cell array of averages, and the output is the  data array that is plotted. Note that in general the cell index is used to describe a given graph, while the first vector index in the graphed data indexes a line in the graph. For multidimensional data, the graphics program automatically generates several different projections of a given graph to allow a complete picture.
    
.. function:: xfunctions (x_nd,in)

    This is a nested cell array of graphics axis transformations. Use when a graph is needed with an axis that is a function of the original axes.  The input of the function is the original axis coordinates, and the output is the new coordinate set. The default value generates the input axes. Called as in.xfunctions{n}{nd}(x_nd,in) for the n-th graph and axis nd, where x_nd is a vector of axis coordinate points for that axis dimension.
    

.. function:: compare (t,in)

    This is a cell array of comparative functions. Each takes the time or frequency vector - or whichever is the first dimension plotted - and returns comparison results for a graphed observable, as a function versus time or frequency (etc). Comparison results are graphed with a dashed line, for the two-dimensional graphs versus the first plotted dimension. There is no default function handle.



.. _sec-parameter-structure:

Parameter structure
===================

Internally, xSPDE parameters and function handles are stored in a cell array, ``latt``, of structures ``r``, which is passed to functions. This includes all the data given above inside the ``in`` structure. In addition, it includes the table of computed parameters given below.

User application constants and parameters should not be reserved names. No reserved name uses capitals (except ``D``), special symbols, or starts with :attr:`c`. Therefore, all names starting with ``in.c``, special symbols or capitals - (except ``D`` for derivatives) -  will be available to the user name-space in future versions of xSPDE.

A parameter structure contains information about the space-time grid and is passed to various functions, for instance :func:`da` or :func:`step`. The corresponding parameter is accessed in the structure `r`, for example, `r.x`.

.. attribute:: t

    Current value of time, :math:`t`.

.. attribute:: x

.. attribute:: y

.. attribute:: z

    Coordinate grids of :math:`x`, :math:`y`, :math:`z`.

.. attribute:: x{n}

    Higher dimensions are labeled numerically as :math:`x_1`,..  :math:`x_6`, and so on. This numerical axis convention can be set even for lower dimensions if ``in.numberaxis`` is set to 1.

.. attribute:: kx

.. attribute:: ky

.. attribute:: kz

    Grids in momentum space: :math:`k_x`, :math:`k_y`, :math:`k_z`.

.. attribute:: k{n}

    Higher dimensions are labeled numerically as :math:`k_5`,  :math:`k_6`, and so on.

.. attribute:: dt

    Output time-step between stored points for data averages.
    
.. attribute:: dtr

    Current reduced time-step used for integration.

.. attribute:: dx

    Steps in coordinate space: :math:`[t,x,y,z,x_5,..]`.

.. attribute:: dk

    Steps in momentum space: :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]`.

.. attribute:: propagator

    Contains the propagator array for the interaction picture.

.. attribute:: v

    Spatial lattice volume.

.. attribute:: kv

    Momentum lattice volume.

.. attribute:: dv

    Spatial cell volume.

.. attribute:: dkv

    Momentum cell volume.

.. attribute:: xc

    Space-time coordinate axes (vector cells).

.. attribute:: kc

    Computational Fourier transform axes in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},.. ]` (vector cells).

.. attribute:: kg

    Graphics  Fourier transform axes in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]` (vector cells).

.. attribute:: kranges

    Range in :math:`[\omega,k_{x},k_{y},k_{z},k_{5},..]` (vector).

.. attribute:: s.dx

    Initial stochastic normalization.

.. attribute:: s.dxt

    Propagating stochastic normalization.

.. attribute:: s.dk

    Initial :math:`k` stochastic normalization.

.. attribute:: s.dkt

    Propagating :math:`k` stochastic normalization.

.. attribute:: nspace

    Number of spatial lattice points: ``in.points(2) * .. * in.points(in.dimension)``.

.. attribute:: nlattice

    Total lattice: ``in.ensembles(1) * nspace``.

.. attribute:: ncopies

    Total copies of stochastic integrations: ``in.ensembles(2) * in.ensembles(3)``.

.. attribute:: d.int

    Dimensions for lattice integration (vector).

.. attribute:: d.a

    Dimensions for :math:`a` field (flattened, vector).

.. attribute:: d.r

    Dimensions for coordinates (flattened, vector).

.. attribute:: d.ft

    Dimensions for field transforms (vector).

.. attribute:: d.k

    Dimensions for noise transforms (vector).

.. attribute:: d.obs

    Dimensions for observations (vector).

.. attribute:: d.data

    Dimensions for average data (flattened, vector).

.. attribute:: d.raw

    Dimensions for raw data (flattened, vector).


Default functions
=================

These functions are used as defaults for simulations and can be overridden by the user.

.. function:: xinitial (~, r)

    Returns a field array filled with zeros.

.. function:: xtransfer (~, ~, a, ~)

    Returns the field ``a`` unchanged.

.. function:: xda (~, ~, r)

    Returns a derivative array filled with zeros.
    
.. function:: xdefine (~, ~, r)

    Returns a define array filled with zeros.

.. function:: xlinear (~, r)

    Returns a linear response array filled with zeros.

.. function:: xobserve (a, ~)

    Returns the real part of ``a(1,:)``.

.. function:: xrfilter (r)

    Returns an array of ones.

.. function:: xnfilter (r)

    Returns an array of ones.

.. function:: xgrid (r)

    Sets grid points in lattice from coordinate vectors. Returns the ``r`` structure with added grid points.

.. function:: xnoisegen (r)

    Generates random noise matrix :math:`z`.

.. function:: xrandomgen (r)

    Generates initial random field matrix :math:`v`.

.. function:: xpropfactor (nc, r)

    Returns the interaction picture propagation factor. ``nc`` is a check index, ``r`` is a lattice structure.


Frequently asked questions
==========================

Answers to some frequent questions, and reminders of points in this chapter are:

-  Can you average other stochastic quantities apart from the field?

   -  Yes: just specify the functions that need to be averaged using the user function :func:`observe`.

-  Can you have functions of the current time and space coordinate?

   -  Yes: xSPDE functions support this using the structure ``r``, as :attr:`t`, :attr:`x`, :attr:`y`, :attr:`z`, or  :attr:``x{1}``, and so on, for more than four space-time dimensions.

-  Can you have several independent stochastic variables?

   -  Yes, input this using ``in.fields(1) > 1``.
   
- I need to have auxiliary fields are defined as functions of other fields.

   -  No problem, input this specification using ``in.fields(2) > 0``, and define the extra fields with the :func:`define` function.

-  Are higher dimensional differential equations possible?

   -  Yes, this requires setting ``in.dimension > 1``. This is essentially unlimited in xSPDE except for memory requirements.

-  Can you have spatial partial derivatives?

   -  Yes, provided they are linear in the fields; these are obtainable using the function :func:`linear`. If they are more general,  use the derivative functions :func:`xd` or or you want special, nonperiodic boundary conditions, then use the finite difference methods,  :func:`xd1`  and  :func:`xd2` .

-  Can you delete the graph heading?

   -  Yes, this is turned off if you set :attr:`headers` to ``0``.

-  Why are there two nearly parallel lines in the graphs sometimes?

   -  These are one standard deviation sampling error limits, generated when ``in.ensembles(2,3) > 1``.

-  Why is there just one line in some graphs, with no sampling errors indicated?

   -  You need ``in.ensembles(2)`` or ``(3)`` for this; see previous question.

-  What are the error bars for?

   -  These are the estimated maximum errors due to finite step-sizes in time.


