*******************
Logic and sequences
*******************

The simulation program logic is straightforward. It is a very compact function called :func:`xspde`. This calls :func:`xsim`, for the simulation, then :func:`xgraph` for the graphics. Most of the work is done by other specialized functions. Input parameters come from an input array, output is saved either in a ``data`` array, or else in a specified file. When completed, timing and error results are printed.


How it works
============

To summarize the previous chapters, xSPDE will solve stochastic partial different equations for a vector field :math:`\boldsymbol{a}(t,\boldsymbol{x})` and vector noise :math:`\boldsymbol{\zeta}(t,\boldsymbol{x})`, of form:

.. math::

    \frac{\partial}{\partial t}\boldsymbol{a}(t,\boldsymbol{x})=\mathbf{A}\left[\boldsymbol{a}\right]+\underline{\mathbf{B}}\left[\boldsymbol{a}\right]\cdot\boldsymbol{\zeta}(t,\boldsymbol{x})+\underline{\mathbf{L}}\left[\boldsymbol{\nabla}\right]\cdot\boldsymbol{a}

It can also solve ordinary stochastic equations, or partial differential equations without noise. Extensive error checking outputs are available. Both initial stochastic conditions and noise can have nonlocal spatial filters applied. All inputs are entered as part of an object-oriented structure. This includes the functions used to specify the equations and the quantities to average. The outputs can be either stored, or graphed interactively. xSPDE includes built-in multidimensional graphics tools.

Sequences
---------

In many types of application, a sequence of stochastic equations requires simulation. In these cases the final field value after integration of one equation becomes the initial value of the next equation in the sequence. 

Sequences therefore are the basic concept used in both the input of parameters to xSPDE, and the storage of data generated.

Input and data arrays
---------------------

To explain xSPDE in full detail,

-  Simulation inputs are stored in the ``input`` cell array.

-  This describes a *sequence* of simulations, so that ``input = {in1, in2, ...}``.

-  Each structure ``in`` describes a simulation (see :ref:`sec-input` for details), whose output is the input of the next.

-  The main function is called using ``data = xspde(input)``.

-  Averages are recorded sequentially in the ``data`` cell array.

-  Raw trajectory data is stored in the ``raw`` cell array if required.

The sequence ``input`` has a number of individual simulation objects ``in``. Each includes parameters that specify the simulation, with functions that give the equations and observables. If there is only one simulation, just one individual specification ``in`` is needed. In addition, xSPDE generates graphs with its own graphics program.

Customization options
---------------------

There are a wide range of customization options available for those who wish to have the very own xSPDE version.

Customization options include functions the allow user definition of:

- inputs    
- interfaces
- stochastic equations   
- mean observables
- linear propagators
- coordinate grids
- noise correlations
- integration methods

**There are four internal options for stochastic integration methods, but arbitrary user specification is also possible.**

The program will print out a record of its progress, then generate the specified graphs.


Stochastic flowchart
====================

The main program logic is nearly self-explanatory. It has four functions
and two main arrays that store results.

.. _fig-flowchart:
.. figure:: Figures/Flowchart.*

   xSPDE flowchart, showing the data, lattice and
   field processing.

There are also two important computational routines behind the scenes, which need to be kept in mind. These are :attr:`in.da`, which is short for difference in :math:`a`. This is completely user specified, and gives a local step in time. The next workhorse routine is :func:`xprop`. This is not a beefy Rugby forward, but calculates spatial propagation.

The logical order is as follows:

:func:`xsim` decides the overall workflow, and parallel operation at a high level. Here, ``in.ensembles(3)`` is used to specify parallel integration, with a ``parfor`` loop. The random seeds include data from the loop index to make sure the noise is independent for each ensemble member, including parallel ensembles.

.. function:: xinpreferences

    is called by :func:`xlattice` to set the defaults that are not already entered.

.. function:: xlattice

    creates a space-time lattice from the input data, which is a data-structure. This also initializes the actual ``data`` array for averaging purposes. Next, a loop is initiated over an ensemble of fields for checking and ensemble averaging. The calculations inside the loop can all be carried our in parallel, if necessary. These internal steps are actually relatively simple.

.. function:: xensemble

    repeats each stochastic path for the check/ensemble loop. It is important to notice that the random seed is reset at the start of each ensemble loop. The seed has a unique value that is different for each ensemble member. Note that for successive simulations that are **not** stored in the same data array, the seed should ideally be manually chosen differently for inputs to successive integration blocks, in order to guarantee independent noise sequences. The check variable can be set to ``in.errorchecks = 1, 2``. This is the total number of integrations carried out. The integration is executed once with ``in.errorchecks = 1``. With ``in.errorchecks = 2``, there are two integrations, using half the step-size the second time. This takes three times as long overall. The matrices used to define the interaction picture transformations are stored **for each check loop,** as they vary with step-size.

.. function:: xpath

    propagates the field ``a`` over a path in time. There are :attr:`in.steps` time-steps for each point stored in time, to allow for greater accuracy without excessive data storage, where needed. This integrates the equations for a predetermined time duration. Note that the random seed has the same value for **both** the check loops. This is because the same number of random variates must be generated in the same order to allow accurate extrapolation. The two loops must use the same random numbers, or else the check is not accurate. For random numbers generated during the integration, the coarse step will add two fine step random noises together, to achieve the goal of identical noise behavior. Results of any required averages, variances and checks are accumulated in the ``data`` array.

.. function:: xprop

    uses Fourier space to calculate a step in the interaction picture, using linear transformations that are pre-calculated. There are both linear transformations and momentum dependent terms available. These are pre-calculated by the :func:`xlattice` function, and stored in the ``prop`` arrays.

Simulation user functions
-------------------------

:attr:`in.initial`

    is used to initialize each integration in time. This is a user-defined function, which can involve random numbers if there is an initial probability distribution. This creates a stochastic field on the lattice, called ``a``. Initialization functions can use coordinates, ``r.x``, ``r.y``, ``r.z``, or for larger dimensions, using numerical lattice labels ``r.x{1}``, ``r.x{2}``, ``r.x{3}``, ``r.x{4}``. Numerical labels can be used for any number of dimension if the switch ``numberaxis=1``. The default is :func:`xinitial`, which sets fields to zero.

:attr:`in.step`

    is the algorithm or method computes each space-time point in the lattice. This also generates the random numbers fields at each time-step. It can be user-modified by setting the handle in.step.

:attr:`in.observe{n}`

    is the n-th observation function whose output is averaged over the ensembles, called from :func:`xpath`. In general, this returns an array whose first coordinate is the line-number of the n-th graph. The default, :func:`xobserve`, returns the real amplitudes. The return value is averaged over the local ensemble and stored as data, ``d{n}``.
    
:attr:`in.function{n}`

    is used when a graph is needed that is a function of the observed averages. The default value is simply ``d{n}``. This is further averaged over higher ensembles to obtain sampling error estimates.


:attr:`in.linear`

    is the linear response, including transverse derivatives in space. The default, :func:`xlinear`, sets this to zero. Derivatives are specified using arrays ``r.Dx``, ``r.Dy``, ``r.Dz``, or for larger dimensions, using numerical lattice labels ``r.D{1}``, ``r.D{2}``, ``r.D{3}``, ``r.D{4}``.

:attr:`in.da`

    is called by :attr:`in.step` to calculate derivatives at every step in the process, including the stochastic terms.

Details of the different parts of the program are given below. Note that the functions ``tic()`` and ``toc()`` are called to time each simulation.


Graphics function, xgraph
=========================

At the end of the loop, global averages and error-bars are calculated. The main functions involved are:

:func:`xgraph` is called by xSPDE when the ensemble loops finished. The results are graphed and output if required.

.. function:: xgpreferences

    is called by :func:`xgraph` to set the graphics defaults that are not already entered.

Comparison results are calculated if available from the user-specified :attr:`in.compare`, an error summary is printed, and the results plotted using the :func:`xgraph` routine, which is a function that graphs the observables. It is prewritten to cover a range of useful graphs, but can be modified to suit the user. The code is intended to cascade down from higher to lower dimension, generating different types of user-defined graphs. Each type of graph is generated once for each specified graphics function.

The code is intended to cascade down from higher to lower dimension, generating different types of user-defined graphs. Each type of graph is generated once for each specified graphics function. The graphics axes that are used for plotting, and the points plotted, are defined using the optional axes input parameters, where :attr:`in.axes{n}` indicates the n-th specified graphics function or set of generated graphical data.

If there are no axes inputs, or the inputs are zero - for example,
``in.axes{1} = {0,0,0}``, then only the lowest dimensions are plotted, up to 3. If the axes inputs project out a single point in a given dimension, - for example, ``axes{1}={0,31,-1,0}``, these axes are suppressed in the plots. This reduces the effective dimension of the data - in this case to two dimensions. 

Examples:

• ``axes{1}={0}``
  - For function 1, plot all the time points; higher dimensions get defaults.

• ``axes{2}={-1,0}``
  - For function 2, plot the maximum time (the default), and all x-points.

• ``axes{3}={1:4:51,32,64}``
  - For function 3, plot every 4-th time point at x point 32, y point 64

• ``axes{4}={0,1:4:51,0}``
  - For function 4, plot all time points, every 4-th x point, and all y-points.

Note that -1 indicates a default point, which is the last point on the time axis, and the midpoint on the other axes. 

The pdimension input can also be used to reduce dimensionality, as this sets the maximum effective plotted dimension. For example, ``pdimension{1}=1`` means that only plots vs time are output for the first function plotted. Note that in the following, t,x,y,z may be replaced by corresponding higher dimensions if there are axes that are suppressed. Slices can be taken at any desired point, not just the midpoint. Using the standard notation of, for example, ``axes{1}={6:3:81}``, can be used to modify the starting, interval, and finishing points for complete control on the plot points.

Results depend on the value of :attr:`in.dimension`, or else the effective graphics dimension if axes are suppressed:

- ``4``: for the highest space dimension, only a slice through :math:`z=0` is available. This is then graphed as if it was in three dimensions.

- ``3``: for two dimensions, distinct graphic images of observable *vs x,y* are plotted at :attr:`in.images` time slices. Otherwise, only a slice through :math:`y=0` is available. This is then treated as if it was in two dimensions.

- ``2``: for two dimensions, one three-dimensional image of observable *vs x,t* is plotted. Otherwise, only a slice through :math:`x=0` is available. This is otherwise treated as in one dimension.

- ``1``: for one dimensions, one image of observable *vs* :math:`t` is plotted, with data at each lattice point in time. Exact results, error bars and sampling error bounds are included if available.

In addition to time-dependent graphs, the :func:`xgraph` function can generate :attr:`in.images` (3D) and :attr:`in.transverse` (2D) plots at specified points in time, up to a maximum given by the number of time points specified. The number of these can be individually specified for each graphics output. The images available are specified in :attr:`in.imagetype`: 3D perspective plots, grey-scale colour plots and contour plots.

Graphics user functions
-----------------------

:attr:`in.gfunction`

    This is used when a graph is needed that is a function of the data coming from the simulation package, since this data can be analysed at a later time. Error estimates are less accurate when this function is used, due to error-propagation effects that may occur after averaging, unless corrected for explicitly in the graphics function. 

:attr:`in.xfunctions`

    This is used when a graph is needed whose axes are a function of the original axes. 

:attr:`in.compare`

    This is used when a two-dimensional graph is needed with a comparison line.

Error control
=============

The final 2D output graphs will have error-bars if :attr:`in.errorchecks` is set to ``2``, which is also the default parameter setting. This is to make sure the final results are accurate. Error-bars below a minimum relative size compared to the vertical range of the plot, specified by the graphics variable :attr:`in.minbar`, are not plotted. There is a clear strategy if the errors are too large.

Either increase the :attr:`in.points`, which gives more plotted points and lower errors, or increase the :attr:`in.steps`, which reduces the step size without changing the graphical resolution. The default algorithm and extrapolation order can be changed, read the xSPDE manual when doing this. Error bars on the graphs can be removed by setting ``in.errorchecks = 1`` or increasing :attr:`in.minbar` in final graphs.

If ``in.ensembles(2) > 1`` or ``in.ensembles(3) > 1``, which allows xSPDE to calculate sampling errors, it will plot upper and lower limits of one standard deviation. If the sampling errors are too large, try increasing ``in.ensembles(1)``, which increases the trajectories in a single thread. An alternative is to increase ``in.ensembles(2)``. This is slower, but is only limited by the compute time, or else to increase ``in.ensembles(3)``, which gives higher level parallelization. Each is limited in different ways; the first by memory, and the second by time, the third by the number of available cores. Sampling error control helps ensures accuracy.

Note that error bars and sampling errors are only graphed for 2D graphs of results vs time. The error-bars are not plotted when they are below a user-specified size, to improve graphics quality. Higher dimensional graphs do not include this, for visibility reasons, but they are still recorded in the data files. Errors caused by the spatial lattice are not checked automatically in the xSPDE code. They must be checked by manually, by comparing results with different transverse lattice ranges and step-size.
