****
xSIM
****

The simulation program logic is straightforward. The main code is a very compact function called :func:`xspde`. This calls :func:`xsim`, for the simulation, then :func:`xgraph` for the graphics. Most of the work is done by other specialized functions. Input parameters come from an input array, output is saved either in a ``data`` array, or else in a specified file. When completed, timing and error results are printed. In this chapter, we go into the workings of the simulation program, xsim.


How xSIM works
==============

To summarize the previous chapters, xSIM will solve stochastic partial different equations for a vector field :math:`\boldsymbol{a}(t,\boldsymbol{x})` and vector noise :math:`\boldsymbol{\zeta}(t,\boldsymbol{x})`, of form:

.. math::

    \frac{\partial}{\partial t}\boldsymbol{a}(t,\boldsymbol{x})=\mathbf{A}\left[\boldsymbol{a}\right]+\underline{\mathbf{B}}\left[\boldsymbol{a}\right]\cdot\boldsymbol{\zeta}(t,\boldsymbol{x})+\underline{\mathbf{L}}\left[\boldsymbol{\nabla},\boldsymbol{a}\right]

It can also solve ordinary stochastic equations, or partial differential equations without noise. Extensive error checking outputs are available. Both initial stochastic conditions and noise can have nonlocal spatial filters applied. All inputs are entered as part of an object-oriented structure. This includes the functions used to specify the equations and the quantities to average. The outputs can be either stored, or graphed interactively.

xSPDE includes a built-in multidimensional graphics tool, xGRAPH, treated in the next chapter.

Sequences
---------

In many types of application, a sequence of stochastic equations requires simulation. In these cases the final field value after integration of one equation becomes the initial value of the next equation in the sequence. At the end of each simulation loop, global averages and error-bars are calculated and stored for output.

Sequences therefore are the basic concept used in both the input of parameters to xSPDE, and the storage of data generated.

Input and data arrays
---------------------

To explain xSIM in full detail,

-  Simulation inputs are stored in the ``input`` cell array.

-  This describes a *sequence* of simulations, so that ``input = {in1, in2, ...}``.

-  Each structure ``in`` describes a simulation, whose output fields are the input of the next.

-  The main function is called using ``[error, input, data, rawdata] = xsim([rawdata,] input)``.

-  Averages are recorded sequentially in the ``data`` cell array.

-  Raw trajectory data is stored in the ``raw`` cell array if required.

The sequence ``input`` has a number of individual simulation objects ``in``. Each includes parameters that specify the simulation, with functions that give the equations and observables. If there is only one simulation, just one individual specification ``in`` is needed. All outputs can be saved to disk storage if required.

The optional [data,] input is only used when there is previous raw data that needs analysis. If this is present, no new simulation takes place. Any ``observe`` functions in the new ``input`` will be employed to take further averages over the existing raw data. This allows re-analysis of large simulation data-sets without more simulations.

The returned ``error`` is a vector: the first component is the maximum error found, the second component is the elapsed time. The returned ``input`` structure is available to the user to give the data file-name, in case xSIM needs to store data with a new file-name. For data security, it will not overwrite existing data.

If xSIM is called within xSPDE, it will generate graphs with its own graphics program xGRAPH. Otherwise, data can be stored then graphed later using xGRAPH.

Customization options
---------------------

There are a wide range of customization options available for those who wish to have the very own xSIM version.

Customization options include functions the allow user definition of:

- inputs    
- interfaces
- stochastic equations
- boundary values   
- mean observables
- linear propagators
- coordinate grids
- noise correlations
- integration methods

**There are four internal options for stochastic integration methods, but arbitrary user specification is also possible.**

User-defined functions have to return specified array sizes, compatible with the internal arrays in xSIM. These sizes are checked by xSIM prior to a simulation. The xSIM program will print out a record of its progress.



Averaged data
================

Observables and functions
--------------------------

To allow options for taking averages, these are carried out in two stages. The first type of average is a local average, taken over any function of the locally stored ensemble of trajectories. These use the :func:observe functions, specified by the user. The default is the real values of each of the fields, stored as a vector. Multiple observe functions can be used, and they are defined as a cell array of functions.

Next, any function can be taken of these local averages, using the :func:function transformations, again specified by the user. The default is the original set of local averages. This is useful if different combinations, such as normalised ratios are needed, or to combine the averages at different times. These second level function outputs are then averaged again over a second level of ensemble averaging, if specified. This is used to obtain estimates of sampling and step-size errors in the final data outputs.

This is explained below in more detail.

Observe functions
-----------------

During the calculation, observables are calculated and averaged over the ``ensembles(1)`` parallel trajectories in the :func:`xpath` function. These are determined by the functions in the :func:`observe` cell array.

The number of :func:`observe` functions may be smaller or larger than the number of vector fields. The observable may be a scalar or vector. These include the averages over the ensembles, and can be visualized as a single graph with one or more lines. The :func:`observe` functions use for input and output the flat or
matrix type internal arrays.

Next, arbitrary functional transforms can be taken, using the :attr:`function` cell array. These functions can use as their input the full set of :func:`observe` output data cell arrays, including a time index. They default to the original :func:`observe` data if they are not user-defined. Functional transforms are most useful if one wishes to use functions which require knowledge of normalization or ensemble averages of lower-level data. There can be more :attr:`function` definitions than :func:`observe` functions if needed.

Each :func:`observe` function or transformation in :func:`xsim` defines a single logical  ``graph`` for the simulation output. However, the graphics function :func:`xgraph` can generate  several projections or views of the same dataset, as explained below.

Combined observables: ``data``
-------------------------------

These results are added to the earlier results in the cell array ``data``, to create a combined set of graphs for the simulation. Initially, both the first and second moment is stored, in order to allow calculation of the sampling error in each quantity.  These are averaged over the higher level ensembles, to allow estimates of sampling errors. Each resulting graph or average data is each stored  in an array of size

.. data:: data  -  all graphics datasets from one sequence member collected in a cell array

    **Cell Array**, has dimension: ``data{graphs}``, made up of a collection of arrays:

#.  graph: observable or function making up a single graph with

    **Array**, has dimension: ``(lines, points(1),..,points(d), errorchecks)``.

In the simplest case, there is just one vector component per average. More generally, the number of components is larger than this if there is a requirement to compare different lines in one graph. Note that, unlike the propagating field, the time dimension is fully expanded.  This is necessary in order to generate outputs at each of the ``points(1)`` time slices. Thus, all the space-time dimensions are stored.

When step-size checking is turned on using the :attr:`checks` flag set to ``1``, a low resolution field is stored for comparison with a high-resolution field of half the step-size, to obtain the time-step error. The observables which are stored have three check indices which are all included in the array. These are the high resolution means, together with error-bars due to time-steps, and estimates of high-resolution standard deviations due to sampling statistics.

Because of the error-checking, the last data dimension, errorchecks, is the total number of components in the data array due to error-checking.  After ensemble averaging, this index is typically ``c = 1, 2, 3``, which is used to index over the:

#. mean value,

#. time-step error-bars and

#. sampling errors

respectively for each space-time point and each graphed function. As a result, the output ``data`` ready for graphing with xGRAPH includes step-size error bars and plotted lines for the two estimated upper and lower standard deviations, obtained from the statistical moments.


Stochastic flowchart
====================

The main program logic is nearly self-explanatory. It has four functions
and two main arrays that store results.

.. _fig-flowchart:
.. figure:: Figures/Flowchart.*

   xSPDE flowchart, showing the data, lattice and
   field processing.

There are also two important computational routines behind the scenes, which need to be kept in mind. These are :func:`da`, which is short for difference in :math:`a`. This is completely user specified, and gives a local step in time. The next workhorse routine is :func:`xprop`. This is not a beefy Rugby forward, but calculates spatial propagation.

The logical order is as follows:

:func:`xsim` decides the overall workflow, and parallel operation at a high level. Here, ``in.ensembles(3)`` is used to specify parallel integration, with a ``parfor`` loop. The random seeds include data from the loop index to make sure the noise is independent for each ensemble member, including parallel ensembles.

.. function:: xinpreferences

    is called by :func:`xlattice` to set the defaults that are not already entered.

.. function:: xlattice

    creates a space-time lattice from the input data, which is a data-structure. This also initializes the actual ``data`` array for averaging purposes. Next, a loop is initiated over an ensemble of fields for checking and ensemble averaging. The calculations inside the loop can all be carried our in parallel, if necessary. These internal steps are actually relatively simple.

.. function:: xensemble

    repeats each stochastic path for the check/ensemble loop. It is important to notice that the random seed is reset at the start of each ensemble loop. The seed has a unique value that is different for each ensemble member. Note that for successive simulations that are **not** stored in the same data array, the seed should ideally be manually chosen differently for inputs to successive integration blocks, in order to guarantee independent noise sequences. The check variable can be set to ``in.checks = 0,1``. The integration is executed only once with ``in.checks = 0``. With ``in.checks = 1``, there is another error-checking integration, using half the step-size the second time. This takes three times as long overall. The matrices used to define the interaction picture transformations are stored **for each check loop,** as they vary with step-size.

.. function:: xpath

    propagates the field ``a`` over a path in time. There are :attr:`steps` time-steps for each point stored in time, to allow for greater accuracy without excessive data storage, where needed. This integrates the equations for a predetermined time duration. Note that the random seed has the same value for **both** the check loops. This is because the same number of random variates must be generated in the same order to allow accurate extrapolation. The two loops must use the same random numbers, or else the check is not accurate. For random numbers generated during the integration, the coarse step will add two fine step random noises together, to achieve the goal of identical noise behavior. Results of any required averages, variances and checks are accumulated in the ``data`` array.

.. function:: xprop

    uses either Fourier space or finite differences to calculate a step in the interaction picture, using linear transformations that are pre-calculated. There are both linear transformations and momentum dependent terms available. These are pre-calculated by the :func:`xlattice` function, and stored in the ``prop`` arrays.

Simulation user functions
-------------------------

:func:`initial`

    is used to initialize each integration in time. This is a user-defined function, which can involve random numbers if there is an initial probability distribution. This creates a stochastic field on the lattice, called ``a``. Initialization functions can use coordinates, ``r.t``,``r.x``, ``r.y``, ``r.z``, or for larger dimensions, using numerical lattice labels ``r.x{1}``, ``r.x{2}``, ``r.x{3}``, ``r.x{4}``. Numerical labels can be used for any number of dimension if the switch ``numberaxis=1``. The default is :func:`xinitial`, which sets fields to zero.

:func:`step`

    is the algorithm or method computes each space-time point in the lattice. This also generates the random numbers fields at each time-step. It can be user-modified by setting the handle in.step. The default is ``in.step = xRK4``.

:func:`observe`

    is a cell array of observation functions whose output is averaged over the ensembles, called from :func:`xpath`. In general, this returns an array whose first coordinate is the line-number of the n-th graph. The default, :func:`xobserve`, returns the real amplitudes. The return value is averaged over the local ensemble and stored as data, ``d{n}``. Note that the input of :func:`observe` is the complete field array.

:func:`function`

    is a cell array of functions used when graphs are needed that are functions of the observed averages. The default value is simply ``d{n}``. This is further averaged over higher ensembles to obtain sampling error estimates. Note that the input of :func:`function` is the complete data cell array, ``d``, which includes all the space-time averages for all the observe functions available.


:func:`linear`

    is the linear response, including transverse derivatives in space. The default, :func:`xlinear`, sets this to zero. Derivatives are specified using arrays ``r.Dx``, ``r.Dy``, ``r.Dz``, or for larger dimensions, using numerical lattice labels ``r.D{2}``, ``r.D{3}``, ``r.D{4}``, ``r.D{5}``.

:func:`da`

    is called by :func:`step` to calculate derivatives at every step in the process, including the stochastic terms. Returns a vector with ``in.fields(1)`` first components.

:func:`define`

    is called by :func:`step` to calculate auxiliary fields at every step in the process. Returns a vector with ``in.fields(2)`` first components.

Details of the different parts of the program are given below. Note that the functions ``tic()`` and ``toc()`` are called to time each simulation.

The xSPDE data and arrays that are user accessible are parameters ``r``, fields ``a``,  average observables ``data``, and raw trajectories ``rawdata``. Apart from the parameters, which are Matlab structures, all fields and data are arrays.

Data arrays and ensembles
=========================

There is a unified index model in all xSPDE arrays. However, in the internal calculations of derivatives and observables, these indices are flattened to give a matrix, as explained below. In all cases, the underlying  xSPDE array index ordering is kept exactly the same:

#. field index :math:`i`

#. time, t index :math:`j_1`

#. x index :math:`j_2`

#. y index :math:`j_3`

#. z index :math:`j_4` ..

#. ensemble or error-checking index :math:`e` or :math:`c`

The number of space dimensions is arbitrary. To conserve storage, one time - the current one - is stored for propagating fields. The ensemble index can be adjusted to increase or decrease local memory usage. If needed, all data generated can be saved in ``rawdata`` arrays.

The fields ``a`` are complex arrays stored discretely on space or momentum grids. Internally, the fields are matrices stored on the flattened xSPDE internal lattice, with just two indices only. Transformations to Fourier space are used both for interaction picture propagation [Caradoc-Davies2000]_ and for averages over Fourier space.

Two different types of Fourier representations are used. In xsim, Fourier transformations are for propagation, which requires the fastest possible methods, and uses :math:`k=0` as the first or lowest index. In xgraph, Fourier transformations are for graphical representations. Hence, the  indices are re-ordered to a conventional index ordering, with negative momentum values in the first index position.

The :ref:`parameters <sec-parameters>` are stored in a structure called, simply, ``r``. It is available to all user-definable routines. The label ``r`` is chosen here for no special reason, and can be changed by the user. These structures reside in a static internal cell array that combines both input and lattice parameters, including the interaction picture transformations. The data is generally different for each simulation in a sequence.

Averaged results are called observables in xSPDE. For each sequence, these are stored in either space or Fourier domains, in the array ``data``, as determined by the :attr:`transforms` vector for each observable. This is a vector of switches for each of the space-time coordinates. The ``data`` arrays obtained in the program as calculations progress are stored in cell arrays, ``cdata``, indexed by a sequence index.

If required, ``rawdata`` ensemble data consisting of all the trajectories ``a`` developing in time can be stored and output. This is memory intensive, and is only done if the :attr:`raw` option is set to ``1``.

All calculated data, including fields, observables and graphics results, is stored in arrays of implicit or explicit rank (2+d), where d is the space-time dimension given in the input. The first index is a field index :math:`(i)`,  while the next indices :math:`j\equiv j_{1},\ldots j_{d}\equiv \mathbf{j}` are for time and space, and the last is a statistics/errors index :math:`(e)`. The space-time dimension d is unlimited.

xSPDE flattened arrays
----------------------

When the fields, noises or coordinates are integrated by the xSPDE integration functions, they are flattened to a matrix. The first index is the field index, and the combined second index covers all the rest. It is more convenient when calculating derivatives and observables in xSIM, to use these flattened arrays or matrices. They are obtained by combining indices :math:`(\mathbf{j},e)` into a flattened second index :math:`J`. This is faster and more compact notationally. Hence, when used in xSPDE functions, the fields are indexed as :math:`a(i,J)`.

xSPDE array types
-----------------

There are several different types of arrays used. Note that for the field, noise and coordinate arrays, only one time index is stored, so :math:`j_1=1`. The stored ensemble index is for the lowest level statistical ensemble, :math:`e_1`. These arrays are as follows:

• Field arrays,   :math:`a(i,\mathbf{j},e_1)` - these have an ensemble index of up to :math:`e_1=ensembles(1)`, but just a single point in time for efficiency.  The fields are flattened to give :math:`a(i,J)`.

• Random and noise arrays,  :math:`w(n,\mathbf{j},e_1)` - these are like field arrays, except that they contain random numbers for the stochastic equations. Random and noise fields are flattened to give :math:`w(n,J)`, where `n` ranges over the available number of noise variables.

• Coordinate arrays :math:`r.x\{l\}(1,\mathbf{j},e_1)` - these store the values of coordinates at grid-points, depending on the axis :math:`l=2,\ldots d` , and are part of the main internal data structure, `r`. These only have a single first index. Coordinates are flattened to give :math:`r.x\{l\}(1,J)`. For less than four total dimensions, this notation is replaced by :math:`r.t`,:math:`r.x(1,J)`,:math:`r.y(1,J)`,:math:`r.z(1,J)`. There is a similar array in momentum space, :math:`k\{l\}(1,J)`.

• Raw arrays,  :math:`r\{s,c,e_2\}(i,\mathbf{j},e_1)` - like fields, but with all points stored. Use with care, as they take up large amounts of memory! Here, we use the notation that :math:`\mathbf{j} = j_1,j_2,\ldots j_d` for :math:`d` space-time dimensions. Note that when output or saved, these have additional cell indices: :math:`s=1,\ldots S` is the sequence number, :math:`c=1,2` for  error-checking the time-step, and  :math:`h=1,\ldots e_2*e_3` for the combined serial and parallel ensemble index. To keep track of all data, an error-check and  ensemble index are needed here.

• Data arrays,  :math:`d\{g\}(\ell,\mathbf{j},c)` - these store the averages, or arbitrary functions of them, with an error-checking index :math:`c=1,2,3`, to store checking data at all time points. Here :math:`g` is the graph index, :math:`\ell' is the line  index. No ensemble index is needed, as these are already ensemble averaged at the first level, so the last index is used to store the checking data at this stage in the code. Here :math:`\mathbf{j} = j_1,j_2,\ldots j_d` space-time points. If the field is transformed, the :math:`\mathbf{j} ` index gives the index in either normal or Fourier space-time, as indicated by the :attr:`transforms` flag.

• Graphics data arrays,  :math:`gd\{s\}\{g\}(\ell,\mathbf{j},c)`  - these store the data that is actually plotted, and can include further functional transformations if required.

The first index :math:`\ell` in a graphics or data array describes different lines on a graph. There can be different first dimensions between fields, noises and output data, as they are specified using different parameters. For only a single output graph, the cell index is not needed.

All outputs have an extra high-level cell index :math:`\{g\}` called the graph or function index. This corresponds to the index :math:`\{g\}` of the observe function used to generate averages. One can have several data arrays in a larger cell arrays to make a number of distinct output graphs labelled :math:`g`, each with multiple averages. Sequences generate separate graphics arrays in sequence, labelled by the first graphics cell index.

More details of ensembles, grids and the internal lattice are given below. Note that the term ``lattice`` is used to refer to the total internal field storage. This combines the local ensemble and the spatial grid together.




Ensembles
---------

Ensembles are used for averaging over stochastic trajectories. They come in three layers: local, serial and parallel, in order to optimize simulations for memory and for parallel operations. The ``in.ensembles(1)`` local  trajectories are used for array-based parallel ensemble averaging, indexed by :math:`e_1`. These trajectories are stored in one array, to allow fast on-chip parallel processing.

Distinct stochastic trajectories are also organized at a higher level into a set of ``in.ensembles(2)`` serial ensembles for statistical purposes, which allows a more precise estimate of sampling error bars. For greater speed, these can  be integrated using ``in.ensembles(3)`` parallel threads. In raw data, these are combined and indexed by the :math:`e_2` cell index.

This hierarchical organization allows allows flexibility in allocating memory and optimizing parallel processing. It is usually faster to have larger values of ``in.ensembles(1)``, but more memory intensive. Using larger values of ``in.ensembles(2)`` is slower, but requires less memory.  Using larger values of ``in.ensembles(3)`` is fast, but requires the Matlab parallel toolbox, and uses both threads and memory resources. It is generally not effective to increase ``in.ensembles(3)`` above the maximum number of available computational cores.

In summary, the stochastic ensembles are defined as follows:

#. Local ensemble: The first or local ensemble contains ``ensembles(1)`` trajectories stored on the xSPDE internal lattice and processed using matrix operations. These are averaged using vector instructions, and indexed locally with the :math:`e_1` index.

#. Serial ensemble: The second or serial ensemble contains ``ensembles(2)`` of the local ensembles, processed in a sequence to conserve memory.

#. Parallel ensemble: The third or parallel ensemble contains ``ensembles(3)`` of the serial ensembles processed in parallel using different threads to allow multi-core and multi-CPU parallel operations. The serial and parallel ensembles are logically equivalent, and give identical results. They are indexed by the combined :math:`e_2` cell index in raw data.


Coordinates, integrals and derivatives
================================================


Time and space
--------------

The default space-time grid for plotted output data is rectangular, with

::

    dx(i) = in.ranges(i) / (in.points(i) -1)

The time index is ``1``, and the space index ``i`` ranges from ``2`` to :attr:`dimension`. The maximum space-time dimension is  unlimited, while ``in.ranges(i)`` is the time and space duration of the simulation, and ``in.points(i)`` is the total number of sampled points available in the ``i``-th direction. The input ``in.boundaries=-1,0,1`` changes the space boundary condition, and is given independently for each field, dimension and boundary. The inputs are ``-1`` for Neumann or specified derivative  boundaries (also used for time), ``0`` for periodic boundaries (the default value) and ``1`` for Dirichlet or vanishing field  boundaries.



Time is advanced in basic integration steps that are equal to or smaller than ``dx(1)``, for purposes of controlling and reducing errors:

::

    dt = dx(1) / (in.steps * nc)

Here, :attr:`steps` is the minimum number of steps used per plotted point, and ``nc = 1, 2`` is the check number. If ``nc = 1``, the run uses coarse time-divisions. If ``nc = 2`` the steps are halved in size for error-checking. Error-checking can be turned off if not required.

The xSPDE space and momentum grid can have any dimension, provided there is enough memory. However, default label values are limited to ten, since more than ten total dimensions will require very large time and storage requirements, and is seldom practical unless the grid is extremely coarse.






Space grid
-------------

We define the grid cell size :math:`dx_{j}` in the :math:`j`-th dimension in terms of maximum range :math:`r_{j}`, the number of points :math:`n_{j}:`, and the boundary value :math:`r_{j}`, as:

.. math::

    dx_{j}=\frac{r_{j}}{n_{j}+b_{j}}.

Each grid starts at a value defined by the vector :attr:`origin`. Using the default values, the time grid starts at :math:`t=0` and ends at :math:`t=T=r_{1}`, for :math:`n=1,\ldots N_{j}`:

.. math::

    t\left(n\right)=(n-1)dt.

Unless there is an offset origin , the :math:`j`-th coordinate grid starts at :math:`-r_{j}/2` and ends at :math:`r_{j}/2` , so that, for :math:`n=1,\ldots n_{j}`:

.. math::

    x_{j}\left(n\right)=-r_{j}/2+(n-1)dx_{j}.

Momentum grid
--------------

All fields can be transformed into Fourier space for taking averages in the :func:`observe` function. This is achieved with the user-defined :attr:`transforms` cell array. This is a cell array of vector switches. For any graph and dimension where :attr:`transforms` is set to unity, the corresponding Fourier transform is taken.

The momentum space graphs and spectral methods all use a Fourier transform definition so that, for :math:`d` dimensions:

.. math::

    \tilde{\boldsymbol{a}}\left(\boldsymbol{k},\omega\right)=\frac{1}{\left(2\pi\right)^{d/2}}\int d\boldsymbol{x}e^{i(\omega t-\boldsymbol{k}\cdot\boldsymbol{x})}\boldsymbol{a}\left(\boldsymbol{x},t\right)

In order to match this to the standard definition of a discrete FFT, the :math:`j`-th momentum lattice cell size :math:`dk_{j}` in the :math:`j`-th dimension is defined in terms of the number of points :math:`N_{j}:`

.. math::

    dk_{j}=\frac{2\pi}{dx_{j}N_{j}}.

The momentum range is therefore

.. math::

    K_{j}=\left(N_{j}-1\right)dk_{j},

while the momentum lattice starts at :math:`-K_{j}/2` and ends at :math:`K_{j}/2` , so that when graphing the data:

.. math::

    k_{j}\left(n\right)=-K_{j}/2+(j-1)dk_{j}.

However, due to the standard definitions of discrete Fourier transforms, the order used during computation and stored in the data arrays is different, namely:

.. math::

    k_{j}\left(n\right)=0..(N_{j}-1)/2)dk_{j},-(N_{j}-1)/2)dk_{j},.-dk_{j}






Averages
--------

There are functions available in xSPDE for grid averages, spatial integrals and derivatives to handle the spatial grid. These can be used to calculate observables for plotting, but are also available for calculating stochastic derivatives as part of the stochastic equation. They operate in parallel over the local ensemble and lattice dimensions. They take a vector or scalar quantity, for example a single field component, and return an average, a space integral, and a spatial derivative respectively. In each case the first argument is the field, the second argument is a vector defining the type of operation, and the last argument is the parameter structure, ``r``. If there are only two arguments, the operation vector is replaced by its default value.

Spatial grid averages can be used to obtain stochastic results with reduced sampling errors if the overall grid is homogeneous. An average is carried out using the builtin xSPDE function :func:`xave` with arguments ``(o, [av, ] r)``.

This takes a vector or scalar field or observable, for example ``o = [1, n.lattice]``, defined on the xSPDE local lattice, and returns an average over the spatial lattice with the same dimension. The input is a field or observable ``o``, and an optional averaging switch ``av``. If ``av(j) > 0``, an average is taken over dimension ``j``. Space dimensions are labelled from ``j = 2 ... 4`` as elsewhere.  If the ``av`` vector is omitted, the average is taken over all space directions.  To average over the local ensemble and all space dimensions, use ``xave(o)``. Averages are returned at all lattice locations.

Higher dimensional graphs of grid averages are generally not useful, as they are simply flat. The xSPDE program allows the user to remove unwanted higher dimensional graphs of average variables. This is achieved by setting the corresponding element of :attr:`pdimension` to the highest dimension required, which depends on which dimensions are averaged.

For example, to average over the entire ensemble plus space lattice and indicate that only time-dependent graphs are required, set ``av = in.dx`` and:

::

    in.pdimension = 1

Note that :func:`xave` on its own gives identical results to those calculated in the :func:`observe` functions. Its utility comes when more complex combinations or functions of ensemble averages are required. If the :attr:`transforms` switch is set, then momentum space averages are returned.

Integrals
---------

Integrals over the spatial grid allow calculation of conserved or other global quantities. To take an integral over the spatial grid,  use the xSPDE function :func:`xint` with arguments ``(o, [dx, ] r)``.

    This function takes a scalar or vector quantity ``o``, and returns a trapezoidal space integral over selected dimensions with vector measure ``dx``. If ``dx(j) > 0`` an integral is taken over dimension ``j``. Dimensions are labelled from ``j = 1, ...`` as in all xSPDE standards. Time integrals are ignored at present. Integrals are returned at all lattice locations. To integrate over an entire lattice, set ``dx = r.dx``, otherwise set ``dx(j) = r.dx(j)`` for selected dimensions ``j``.

As with averages, the xSPDE program allows the user to remove unwanted higher dimensional graphs when the integrated variable is used as an observable. For example, in a four dimensional simulation with integrals taken over the :math:`y` and :math:`z` coordinates, only :math:`t`- and :math:`x`-dependent graphs are required. Hence, set ``dx`` to ``[0, 0, r.dx(3), r.dx(4)]``, and:

::

    in.pdimension = 2

If momentum-space integrals are needed, use the :attr:`transforms` switch to make sure that the field is Fourier transformed, and input :attr:`dk` instead of :attr:`dx`. Note that :func:`xint` returns a lattice observable, as required when used in the :func:`observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.




Derivatives in equations
------------------------

xSPDE can support either spectral or finite difference methods for derivatives. The default spectral method used is a discrete Fourier transform, but other methods can be added, as the code is inherently extensible. These derivatives are obtained through function calls.

The code to take a spectral derivative, using spatial Fourier transforms, is carried out using the xSPDE :func:`xd` function with arguments ``(o, [D, ] r)``. This can be used both in calculating derivatives for equations, and for averages or observables if they are needed.

This function takes a scalar or vector quantity ``o``, and returns a spectral derivative over selected dimensions with a derivative ``D``, by Fourier transforming the data.  Set ``D = r.Dx`` for a first order x-derivative, ``D = r.Dy`` for a first order y-derivative, and similarly ``D = r.Dz.*r.Dy`` for a cross-derivative in ``z`` and ``y``. Higher derivatives require powers of these, for example `D = r.Dz.^4``. For higher dimensions use numerical labels, where ``D = r.Dx`` becomes ``D = r.D{2}``, and so on. Time derivatives are ignored at present. Derivatives are returned at all lattice locations.

If the derivative ``D`` is omitted, a first order x-derivative is returned.
Note that :func:`xd` returns a lattice observable, as required when used in the :func:`observe` function. If the integral is used in another function, it returns a matrix of dimension ``[1, lattice]``.

Finite difference first derivatives
-----------------------------------

The code to take a first order spatial derivative with finite difference methods is carried out using the xSPDE function :func:`xd1` with arguments ``(o, [dir, ] r)``.

This takes a scalar or vector ``o``, and returns a first derivative with an axis direction ``dir``.  Set ``dir = 2`` for an x-derivative, ``dir = 3`` for a y-derivative.  Time derivatives are ignored at present. Derivatives are returned at all lattice locations.

If the direction ``dir`` is omitted, an x-derivative is returned. These derivatives can be used both in calculating propagation, and in calculating observables. The boundary condition is set by the in.boundaries input. It can be made periodic, which is the default, or Neumann with zero derivative, or Dirichlet with zero field.

Finite difference second derivatives
------------------------------------

The code to take a second order spatial derivative with finite difference methods is carried out using the xSPDE :func:`xd2` with arguments ``(o, [dir, ] r)`` function.

This takes a scalar or vector ``o``, and returns the second  derivative in axis direction ``dir``.  Set ``dir = 2`` for an x-derivative, ``dir = 3`` for a y-derivative.  All other properties are exactly the same as :func:`xd1`.




Interaction picture and Fourier transforms
==========================================

The xSPDE algorithms all allow the use of a sequence of interaction pictures. Each successive interaction picture is referenced to :math:`t=t_{n}`, for the n-th step starting at :math:`t=t_{n}`, so :math:`\boldsymbol{a}_{I}(t_{n})=\boldsymbol{a}(t_{n})\equiv\boldsymbol{a}_{n}`. It is possible to solve stochastic partial differential equations in xSPDE using explicit derivatives, but this is often less efficient.

A conventional discrete Fourier transform (DFT) using a fast Fourier transform method is employed for the interaction picture (IP) transformations used in computations, as this is fast and simple. In one dimension, this is given by a sum over indices starting with zero, rather than the Matlab convention of one. Hence, if  :math:`\tilde{m}=m-1`:

.. math::
 A_{\tilde{n}}=\mathcal{F}\left(a\right)=\sum_{\tilde{m}=0}^{N-1}a_{\tilde{m}}\exp\left[-2\pi i\tilde{m}\tilde{n}/N\right]

Suppose the spatial grid spacing is :math:`dx`, and the number of grid points is :math:`N`, then the maximum range from the first to last point is:

.. math::

    R=(N-1)dx

We note that the momentum grid spacing is

.. math::

    dk=\frac{2\pi}{Ndx}

The IP Fourier transform can be written in terms of an FFT as

.. math::

    \boldsymbol{A}\left(\boldsymbol{k}_{\boldsymbol{n}}\right)=\prod_{j}\left[\sum_{\tilde{m}_{j}}\exp\left[-i\left(dk_{j}dx_{j}\right)\tilde{m}_{j}\tilde{n}_{j}\right]\right]

The inverse FFT Fourier transforms automatically divide by the correct factors of :math:`\prod_{j}N_{j}` to ensure invertibility. Note also that due to the periodicity of the exponential function, negative momenta are obtained if we consider an ordered lattice such that:

.. math::

    \begin{aligned}
    k_{j} & = (j-1)dk\,\,\,(j\le N/2)\\
    k_{j} & = (j-1-N)dk\,\,(j>N/2)
    \end{aligned}

This Fourier transform is multiplied by an appropriate factor to propagate in the interaction picture, than an inverse Fourier transform is applied. While it is for interaction picture transforms, an additional scaling factor is applied to obtain transformed fields in averages.

In other words, in the averages

.. math::

 \tilde{a}_{n} = \frac{dt}{\sqrt{2\pi}} A_{\tilde{n}'}

where the indexing change indicates that graphed momenta are stored from negative to positive values. Note also that for frequency spectra a positive sign is used in the frequency exponent, to agree with physics conventions.


Interaction picture derivatives
-------------------------------

For calculating derivatives in the interaction picture, the notation :math:`D` indicates a derivative. To explain, one integrates by parts:

.. math::

    D^{p}\tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right)=\left[ik_{x}\right]^{p}\tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right)=\frac{1}{\left(2\pi\right)^{d/2}}\int d\boldsymbol{x}e^{-i\boldsymbol{k}\cdot\boldsymbol{x}}\left[\frac{\partial}{\partial x}\right]^{p}\boldsymbol{a}\left(\boldsymbol{x}\right)\label{eq:Fourier derivative}

This means, for example, that to calculate a one dimensional space derivative in the Linear interaction picture routine, one uses:

- :math:`\nabla_{x}\rightarrow` ``r.Dx``

Here ``r.Dx`` returns an array of momenta in cyclic order in dimension :math:`d` as defined above, suitable for an FFT calculation. The imaginary :math:`i` is not needed to give the correct sign, from the equation above. Instead, it is included in the D array. In two dimensions, the code to return a full two-dimensional Laplacian is:

- :math:`\boldsymbol{\nabla}^{2}=\nabla_{x}^{2}+\nabla_{y}^{2}\rightarrow` ``r.Dx.^2+r.Dy.^2``

Note that the dot in the notation of ``.^`` is needed to take the square of each element in the array.

Spectra in the time-domain
--------------------------

For calculating a spectrum in the time-domain, the method of inputting a :attr:`transforms` switch is used, with ``transforms{n}(1) = 1`` to turn on Fourier transforms in the time domain for the n-th observable. This requires much more dedicated internal memory.

To conserve memory, one can use more internal :attr:`steps` combined with less :attr:`points`. In order to ensure that spectral results are independent of memory conservation strategies, xSPDE uses a technique of trapezoidal averaging when calculating frequency spectra.

With this method, all fields are averaged internally using trapezoidal integration in time over any internal steps, to give the average midpoint value.  After this, the resulting step-averaged fields are then Fourier transformed.

For example, in the simplest case of just one internal step, with no error-checking, this means that the field used to calculate a spectrum is:

.. math::

    \bar{a}_{j}=\left({a}_{j-1}+{a}_{j}\right)/2,

which corresponds to the time in the spectral Fourier transform, of:

.. math::

    \bar{t}_{j}=\left({t}_{j-1}+{t}_{j}\right)/2.

For an error-checking calculation with two internal :attr:`steps`, there are four successive valuations: :math:`a_{j1}`, :math:`a_{j2}`, :math:`a_{j3}`, :math:`a_{j}`, with the last value the one plotted at :math:`t_{j}`. In this case, for spectral calculations one averages according to:

.. math::

 \bar{a}_{j}=\left({a}_{j-1}+2({a}_{j1}+{a}_{j2}+{a}_{j3})+{a}_{j}\right)/8.

When there are even larger numbers of internal steps, either from error-checking or from using the internal :attr:`steps` parameter, one proceeds similarly by carrying out a trapezoidal average over all internal steps.

In addition, one must define the first field :math:`\bar{a}_{1}`. Due to the cyclic nature of discrete Fourier transforms, this is also logically the last field value.  Hence, this is set equal to the corresponding cyclic average of the first and last field value, in order to reduce aliasing errors at high frequencies in the resulting spectrum:

.. math::

    \bar{a}_{1}=\frac{1}{2} \left({a}_{N}+{a}_{1}\right),

which corresponds to a time in the spectral Fourier transform of:

.. math::

    \bar{t}_{1} = {t}_{1}-dt/2 \equiv {t}_{N}+dt/2.

This aliasing of virtual times, one higher and one lower than any integration time, is a consequence of the discrete Fourier transform method. It also means that the effective total integration time in the Fourier transform definition is :math:`T_{eff} = T+dt = 2\pi/d\omega`, where :math:`T` is the total integration time, and :math:`dt` is the time interval between integration points.



Fields
======

In the xSIM code, the complex vector field ``a`` is generally stored as a compressed or flattened matrix with dimensions ``[fields, lattice]``. Here ``nlattice`` is the total number of lattice points including an ensemble dimension, to increase computational efficiency:

::

    nlattice =  nspace * ensembles (1)

The total number of space points ``r.nspace`` is given by:

::

    nspace = points (2) * ... * points (in.dimension)

The use of a matrix for the fields is convenient in that fast matrix operations are possible in a high-level language.

In different subroutines it may be necessary to expand out this array to more easily reference the array structure. However, the internal field structure ``a`` at a single time-point is as follows

.. data:: a

** - Array** of dimension: (:attr:`fieldsplus`,  ``nlattice``)

Note: Here, :attr:`fieldsplus` = :attr:`fields` (1) + :attr:`fields` (2) is the total number of field components. Here :attr:`fields` (1) are the dynamical fields, while :attr:`fields` (2) are defined or auxiliary fields that are sometimes necessary.  Also, :attr:``in.ensembles`` (1) is the number of statistical samples processed as a parallel vector. This can be set to one to save data space, or increased to improve parallel efficiency. The time dimension :attr:`points` (1) is always compressed to one during calculations. During spectral calculations, and for raw output, the full length of the time lattice, :attr:`points` (1), is stored, which increases memory requirements.


Raw output
---------------

If required, by using the switch :attr:`raw` set to one.  xSPDE  will then store every trajectory generated. This is raw, unprocessed data, so there is no graph index. The raw data output is stored in an output cell array :data:`raw`. The array is written to disk using the Matlab file-name, on completion, provided a file name is input, and is also available as an xSIM function output.

The cell indices are: sequence index, error-checking index, ensemble index.

.. data:: raw

    ** - Cell Array**, has dimension: ``raw{sequence, check, in.ensemble(2)*in.ensemble(3)}``

If thread-level parallel processing is used, these are also stored in the cell array, which is indexed over both the parallel and serial ensemble. Inside each raw cell is at least one complete space-time :data:`a` stored as a complex array, with indices for the field index, the time-space lattice, and the samples.

Each location in the cell array stores one time-space-sample trajectory in xSPDE, which is a real or complex array with (:attr:`dimension` + 2) indices, noting that :attr:`points` is a vector with :attr:`dimension` indices. Here the dynamical fields are  expanded to more easily reference the array structure. The full, expanded field structure ``a`` at a single time-point is as follows

.. data:: a

  ** - Array** of dimension: (:attr:`fieldsplus`,  :attr:`points`, :attr:`ensembles` (1))

Here, :attr:`fieldsplus` = :attr:`fields` (1) + :attr:`fields` (2) is the total number of field components, where :attr:`fields` (1) are the dynamical fields, while :attr:`fields` (2) are defined or auxiliary fields.  Also, :attr:``in.ensembles`` (1) is the number of statistical samples processed as a parallel vector. This can be set to one to save data space, or increased to improve parallel efficiency. Provided no frequency information is needed, the time dimension :attr:`points` (1) is compressed to one during calculations. During spectral calculations, and for raw output, the full length of the time lattice, :attr:`points` (1), is stored, which increases memory requirements.


The main utility of raw data is for storing data-sets from large simulations for later re-analysis. It is also a platform for further development of analytic tools for third party developers, to treat statistical features not included in the functional tools provided. For example, one might need to plot histograms of distributions from this.
