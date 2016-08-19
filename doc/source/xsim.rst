****
xSIM
****

The simulation program logic is straightforward. The main code is a very compact function called :func:`xspde`. This calls :func:`xsim`, for the simulation, then :func:`xgraph` for the graphics. Most of the work is done by other specialized functions. Input parameters come from an input array, output is saved either in a ``data`` array, or else in a specified file. When completed, timing and error results are printed. In this chapter, we go into the workings of the simulation program, xsim.


How xSIM works
==============

To summarize the previous chapters, xSIM will solve stochastic partial different equations for a vector field :math:`\boldsymbol{a}(t,\boldsymbol{x})` and vector noise :math:`\boldsymbol{\zeta}(t,\boldsymbol{x})`, of form:

.. math::

    \frac{\partial}{\partial t}\boldsymbol{a}(t,\boldsymbol{x})=\mathbf{A}\left[\boldsymbol{a}\right]+\underline{\mathbf{B}}\left[\boldsymbol{a}\right]\cdot\boldsymbol{\zeta}(t,\boldsymbol{x})+\underline{\mathbf{L}}\left[\boldsymbol{\nabla}\right]\cdot\boldsymbol{a}

It can also solve ordinary stochastic equations, or partial differential equations without noise. Extensive error checking outputs are available. Both initial stochastic conditions and noise can have nonlocal spatial filters applied. All inputs are entered as part of an object-oriented structure. This includes the functions used to specify the equations and the quantities to average. The outputs can be either stored, or graphed interactively. 

xSPDE includes a built-in multidimensional graphics tool, xGRAPH, treated in the next chapter.

Sequences
---------

In many types of application, a sequence of stochastic equations requires simulation. In these cases the final field value after integration of one equation becomes the initial value of the next equation in the sequence. 

Sequences therefore are the basic concept used in both the input of parameters to xSPDE, and the storage of data generated.

Input and data arrays
---------------------

To explain xSPDE in full detail,

-  Simulation inputs are stored in the ``input`` cell array.

-  This describes a *sequence* of simulations, so that ``input = {in1, in2, ...}``.

-  Each structure ``in`` describes a simulation, whose output fields are the input of the next.

-  The main function is called using ``data = xspde(input)``.

-  Averages are recorded sequentially in the ``data`` cell array.

-  Raw trajectory data is stored in the ``raw`` cell array if required.

The sequence ``input`` has a number of individual simulation objects ``in``. Each includes parameters that specify the simulation, with functions that give the equations and observables. If there is only one simulation, just one individual specification ``in`` is needed. In addition, xSPDE generates graphs with its own graphics program.

Customization options
---------------------

There are a wide range of customization options available for those who wish to have the very own xSIM version.

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

There are also two important computational routines behind the scenes, which need to be kept in mind. These are :func:`da`, which is short for difference in :math:`a`. This is completely user specified, and gives a local step in time. The next workhorse routine is :func:`xprop`. This is not a beefy Rugby forward, but calculates spatial propagation.

The logical order is as follows:

:func:`xsim` decides the overall workflow, and parallel operation at a high level. Here, ``in.ensembles(3)`` is used to specify parallel integration, with a ``parfor`` loop. The random seeds include data from the loop index to make sure the noise is independent for each ensemble member, including parallel ensembles.

.. function:: xinpreferences

    is called by :func:`xlattice` to set the defaults that are not already entered.

.. function:: xlattice

    creates a space-time lattice from the input data, which is a data-structure. This also initializes the actual ``data`` array for averaging purposes. Next, a loop is initiated over an ensemble of fields for checking and ensemble averaging. The calculations inside the loop can all be carried our in parallel, if necessary. These internal steps are actually relatively simple.

.. function:: xensemble

    repeats each stochastic path for the check/ensemble loop. It is important to notice that the random seed is reset at the start of each ensemble loop. The seed has a unique value that is different for each ensemble member. Note that for successive simulations that are **not** stored in the same data array, the seed should ideally be manually chosen differently for inputs to successive integration blocks, in order to guarantee independent noise sequences. The check variable can be set to ``in.checks = 0,1``. The integration is executed once with ``in.checks = 0``. With ``in.checks = 1``, there is another error-checking integration, using half the step-size the second time. This takes three times as long overall. The matrices used to define the interaction picture transformations are stored **for each check loop,** as they vary with step-size.

.. function:: xpath

    propagates the field ``a`` over a path in time. There are :attr:`steps` time-steps for each point stored in time, to allow for greater accuracy without excessive data storage, where needed. This integrates the equations for a predetermined time duration. Note that the random seed has the same value for **both** the check loops. This is because the same number of random variates must be generated in the same order to allow accurate extrapolation. The two loops must use the same random numbers, or else the check is not accurate. For random numbers generated during the integration, the coarse step will add two fine step random noises together, to achieve the goal of identical noise behavior. Results of any required averages, variances and checks are accumulated in the ``data`` array.

.. function:: xprop

    uses either Fourier space or finite differences to calculate a step in the interaction picture, using linear transformations that are pre-calculated. There are both linear transformations and momentum dependent terms available. These are pre-calculated by the :func:`xlattice` function, and stored in the ``prop`` arrays.

Simulation user functions
-------------------------

:func:`initial`

    is used to initialize each integration in time. This is a user-defined function, which can involve random numbers if there is an initial probability distribution. This creates a stochastic field on the lattice, called ``a``. Initialization functions can use coordinates, ``r.x``, ``r.y``, ``r.z``, or for larger dimensions, using numerical lattice labels ``r.x{1}``, ``r.x{2}``, ``r.x{3}``, ``r.x{4}``. Numerical labels can be used for any number of dimension if the switch ``numberaxis=1``. The default is :func:`xinitial`, which sets fields to zero.

:func:`step`

    is the algorithm or method computes each space-time point in the lattice. This also generates the random numbers fields at each time-step. It can be user-modified by setting the handle in.step. The default is ``in.step = xRK4``.

:func:`observe`

    is a cell array of observation functions whose output is averaged over the ensembles, called from :func:`xpath`. In general, this returns an array whose first coordinate is the line-number of the n-th graph. The default, :func:`xobserve`, returns the real amplitudes. The return value is averaged over the local ensemble and stored as data, ``d{n}``. Note that the input of :func:`observe` is the complete field array.
    
:func:`function`

    is a cell array of functions used when graphs are needed that are functions of the observed averages. The default value is simply ``d{n}``. This is further averaged over higher ensembles to obtain sampling error estimates. Note that the input of :func:`function` is the complete data cell array, ``d``, which includes all the averages for all the observe functions available.


:func:`linear`

    is the linear response, including transverse derivatives in space. The default, :func:`xlinear`, sets this to zero. Derivatives are specified using arrays ``r.Dx``, ``r.Dy``, ``r.Dz``, or for larger dimensions, using numerical lattice labels ``r.D{2}``, ``r.D{3}``, ``r.D{4}``, ``r.D{5}``.

:func:`da`

    is called by :func:`step` to calculate derivatives at every step in the process, including the stochastic terms. Returns a vector with ``in.fields(1)`` first components.
    
:func:`define`

    is called by :func:`step` to calculate auxiliary fields at every step in the process. Returns a vector with ``in.fields(2)`` first components.

Details of the different parts of the program are given below. Note that the functions ``tic()`` and ``toc()`` are called to time each simulation.

The xSPDE data and arrays that are user accessible are parameters ``r``, fields ``a``,  average observables ``data``, and raw trajectories ``rawdata``. Apart from the parameters, which are Matlab structures, all fields and data are arrays. 

Indices in arrays
=================

There is a unified index model in all xSPDE arrays. However, in the internal calculations of derivatives and observables, these indices are flattened to give a matrix, as explained below. In all cases, the underlying  xSPDE array index ordering is kept exactly the same:

#. field index :math:`i`

#. ensemble or error-checking index :math:`e` or :math:`c`

#. time, t index :math:`j_1`

#. x index :math:`j_2`

#. y index :math:`j_3`

#. z index :math:`j_4` ..

The number of space dimensions is arbitrary. To conserve storage, one time - the current one - is stored for propagating fields. The ensemble index can be adjusted to increase or decrease local memory usage. If needed, all data generated can be saved in ``rawdata`` arrays.

The fields ``a`` are complex arrays stored discretely on space or momentum grids. Internally, the fields are matrices stored on the flattened xSPDE internal lattice, with just two indices only. Transformations to Fourier space are used both for interaction picture propagation [Caradoc-Davies2000]_ and for averages over Fourier space. 

Two different types of Fourier representations are used. In xsim, Fourier transformations are for propagation, which requires the fastest possible methods, and uses :math:`k=0` as the first or lowest index. In xgraph, Fourier transformations are for graphical representations. Hence, the  indices are re-ordered to a conventional index ordering, with negative momentum values in the first index position.


The :ref:`parameters <sec-parameters>` are stored in a structure called, simply, ``r``. It is available to all user-definable routines. The label ``r`` is chosen because the parameters include the grid coordinates in space and time. These structures reside in a static internal cell array that combines both input and lattice parameters, including the interaction picture transformations, called :data:`latt`. The data in :data:`latt` is different for each simulation in a sequence.

Averaged results are called observables in xSPDE. For each sequence, these are stored in either space or Fourier domains, in the array ``data``, as determined by the :attr:`transforms` vector for each observable. This is a vector of switches for each of the space-time coordinates. The ``data`` arrays obtained in the program as calculations progress are stored in cell arrays, ``cdata``, indexed by a sequence index.

If required, ``rawdata`` ensemble data consisting of all the trajectories ``a`` developing in time can be stored and output. This is memory intensive, and is only done if the :attr:`raw` option is set to ``1``.

All calculated data, including fields, observables and graphics results, is stored in arrays of implicit or explicit rank (2+d), where d is the space-time dimension given in the input. The first index is a field index :math:`(i)`, the second a statistics/errors index :math:`(e)`, while the remaining indices :math:`j\equiv j_{1},\ldots j_{d}\equiv j_{1},\mathbf{j}` are for time and space. The space-time dimension d is unlimited. 

xSPDE array types
-----------------

There are several different types of arrays used. These are as follows:

• Field arrays,   :math:`a(i,e_1,1,\mathbf{j})` - these have an ensemble index of up to :math:`e_1=ensembles(1)`, but just a single, i.e., present time-point for efficiency.

• Random and noise arrays,  :math:`w(i,e_1,1,\mathbf{j})` - these are like field arrays, except that they contain random numbers for the stochastic equations.

• Coordinate arrays :math:`x\{l\}(1,e_1,1,\mathbf{j})` - these store the values of coordinates at grid-points, depending on the axis :math:`l=2,\ldots d` .

• Raw arrays,  :math:`r\{s,c,e_2\}(i,e_1,j)` - like fields, but with all points stored. Use with care, as they take up large amounts of memory. Note that when output or saved, these have additional cell indices: :math:`s=1,\ldots S` is the sequence number, :math:`c=1,2` for the error-checking of the time-step :math:`e_2=1,2` for the combined serial and parallel ensemble index. To keep track of all the data, an error-check and two ensemble indices are needed here.

• Data arrays,  :math:`d\{n\}(i,c,j)` - these store the averages, or arbitrary functions of them, with an error-checking index :math:`c=1,2,3`, to store checking data at all time points. No ensemble index is needed, as these are ensemble averages, so the second index is used to store the checking data at this stage in the code. If the data is transformed, the :math:`j` index gives the index in Fourier space.

• Graphics arrays,  :math:`g\{n\}(i,c,j)`  - these store the data that is plotted, and can include further functional transformations if required.

The field index :math:`i` in a graphics or data array describes different lines on a graph. There can be quite different first dimensions between fields, noises and output data, as they are specified using different parameters. If only a single output graph is wanted, the observe cell index is not needed.

All outputs have an extra high-level cell index :math:`\{n\}` called the graph or function index. This corresponds to the index :math:`\{n\}` of the observe function used to generate averages. One can have several data arrays in a larger cell arrays to make a number of distinct output graphs labelled :math:`n`, each with multiple averages. Sequences generate separate graphics arrays.

More details of ensembles, grids and the internal lattice are given below. Note that the term ``lattice`` is used to refer to the total internal field storage. This combines the local ensemble and the spatial grid together. 


xSPDE flattened arrays
----------------------

When the fields, noises or coordinates are integrated by the xSPDE integration functions, they are flattened to a matrix. The first index is the field index, and the combined second index covers all the rest. It is simply more convenient when calculating derivatives and observables in xSIM, to use these flattened arrays or matrices. They are obtained by combining indices :math:`(e,j)` into a flattened second index :math:`J`. This is faster and more compact notationally. Hence, when used in xSPDE functions, the fields are indexed as :math:`a(i,J)`. 

Ensembles
================

Ensembles are used for averaging over stochastic trajectories. They come in three layers: local, serial and parallel, in order to optimize simulations for memory and for parallel operations. The ``in.ensembles(1)`` local  trajectories are used for array-based parallel ensemble averaging, indexed by :math:`e_1`. These trajectories are stored in one array, to allow fast on-chip parallel processing. 

Distinct stochastic trajectories are also organized at a higher level into a set of ``in.ensembles(2)`` serial ensembles for statistical purposes, which allows a more precise estimate of sampling error bars. For greater speed, these can  be integrated using ``in.ensembles(3)`` parallel threads. In raw data, these are combined and indexed by the :math:`e_2` cell index. 

This hierarchical organization allows allows flexibility in allocating memory and optimizing parallel processing. It is usually faster to have larger values of ``in.ensembles(1)``, but more memory intensive. Using larger values of ``in.ensembles(2)`` is slower, but requires less memory.  Using larger values of ``in.ensembles(3)`` is fast, but requires the Matlab parallel toolbox, and uses both threads and memory resources. It is generally not effective to increase ``in.ensembles(3)`` above the maximum number of available computational cores.

In summary, the ensembles are defined as follows:

Local ensemble
--------------

The first or local ensemble contains ``ensembles(1)`` trajectories stored on the xSPDE internal lattice and processed using matrix operations. These are averaged using vector instructions, and indexed locally with the :math:`e_1` index.

Serial and parallel ensembles
-----------------------------

The second or serial ensemble contains ``ensembles(2)`` of the local ensembles, processed in a sequence to conserve memory. 
 
The third or parallel ensemble contains ``ensembles(3)`` of the serial ensembles processed in parallel using different threads to allow multi-core and multi-CPU parallel operations.

The serial and parallel ensembles are logically equivalent, and give identical results. They are indexed by the combined :math:`e_2` cell index in raw data.


Coordinates, integrals and derivatives
================================================


Time and space
--------------



Time is advanced in basic integration steps that are equal to or smaller than ``dx(1)``, for purposes of controlling and reducing errors:

::

    dt = dx(1) / (in.steps * nc)

Here, :attr:`steps` is the minimum number of steps used per plotted point, and ``nc = 1, 2`` is the check number. If ``nc = 1``, the run uses coarse time-divisions. If ``nc = 2`` the steps are halved in size for error-checking. Error-checking can be turned off if not required.

The xSPDE space and momentum grid can have any dimension, provided there is enough memory. Using more than six to ten total dimensions causes large time requirements and is not very practical.

The default spatial grid
 for plotted output data is rectangular, with 

::

    dx(i) = in.ranges(i) / (in.points(i) - 1)

The time index is ``1``, and the space index ``i`` ranges from ``2`` to :attr:`dimension`. The maximum space-time dimension is ``in.dimension = 4``, while ``in.ranges(i)`` is the time and space duration of the simulation, and ``in.points(i)`` is the total number of plotted points in the ``i``-th direction.




Space grid
-------------

We define the grid cell size :math:`dx_{j}` in the :math:`j`-th dimension in terms of maximum range :math:`R_{j}` and the number of points :math:`N_{j}:`

.. math::

    dx_{j}=\frac{R_{j}}{N_{j}-1}.

Each grid starts at a value defined by the vector :attr:`origin`. Using the default values, the time grid starts at :math:`t=0` and ends at :math:`t=T=r_{1}`, for :math:`n=1,\ldots N_{j}`:

.. math::

    t\left(n\right)=(n-1)dt.

The :math:`j`-th coordinate grid starts at :math:`-r_{j}/2` and ends at :math:`r_{j}/2` , so that, for :math:`n=1,\ldots N_{j}`:

.. math::

    x_{j}\left(n\right)=-R_{j}/2+(n-1)dx_{j}.

Momentum grid
--------------

The momentum space graphs and spectral methods all use a Fourier transform definition so that, for :math:`d` dimensions:

.. math::

    \tilde{\boldsymbol{a}}\left(\boldsymbol{k},\omega\right)=\frac{1}{\left(2\pi\right)^{d/2}}\int d\boldsymbol{x}e^{i(\omega t-\boldsymbol{k}\cdot\boldsymbol{x})}\boldsymbol{a}\left(\boldsymbol{x},t\right)

In order to match this to the standard definition of a discrete FFT, the :math:`j`-th momentum lattice cell size :math:`dk_{j}` in the :math:`j`-th dimension is defined in terms of the number of points :math:`N_{j}:`

.. math::

    dk_{j}=\frac{2\pi}{dx_{j}N_{j}}.

The momentum range is therefore

.. math::

    K_{j}=\left(N_{j}-1\right)dk_{j},

while the momentum lattice starts at :math:`-k_{j}/2` and ends at :math:`k_{j}/2` , so that when graphing the data:

.. math::

    k_{j}\left(n\right)=-K_{j}/2+(N_{j}-1)dk_{j}.
    
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

Note that :func:`xave` on its own gives identical results to those calculated in the :func:`observe` functions. Its utility comes when more complex combinations or functions of ensemble averages are required.

Integrals
---------

Integrals over the spatial grid allow calculation of conserved or other global quantities. To take an integral over the spatial grid,  use the xSPDE function :func:`xint` with arguments ``(o, [dx, ] r)``.

    This function takes a scalar or vector quantity ``o``, and returns a trapezoidal space integral over selected dimensions with vector measure ``dx``. If ``dx(j) > 0`` an integral is taken over dimension ``j``. Dimensions are labelled from ``j = 1, ...`` as in all xSPDE standards. Time integrals are ignored at present. Integrals are returned at all lattice locations. To integrate over an entire lattice, set ``dx = r.dx``, otherwise set ``dx(j) = r.dx(j)`` for selected dimensions ``j``.

As with averages, the xSPDE program allows the user to remove unwanted higher dimensional graphs when the integrated variable is used as an observable. For example, in a four dimensional simulation with integrals taken over the :math:`y` and :math:`z` coordinates, only :math:`t`- and :math:`x`-dependent graphs are required. Hence, set ``dx`` to ``[0, 0, r.dx(3), r.dx(4)]``, and:

::

    in.pdimension = 2

If momentum-space integrals are needed, use the transform switch to make sure that the field is Fourier transformed, and input :attr:`r.dk` instead of :attr:`r.dx`. Note that :func:`xint` returns a lattice observable, as required when used in the :func:`observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.


Spectral derivatives
--------------------

xSPDE can support either spectral or finite difference methods for derivatives. The default spectral method used is a discrete Fourier transform, but other methods can be added, as the code is inherently extensible.

The code to take a spectral derivative, using spatial Fourier transforms, is carried out using the xSPDE :func:`xd` function with arguments ``(o, [D, ] r)``.

This function takes a scalar or vector quantity ``o``, and returns a spectral derivative over selected dimensions with a derivative ``D``, by Fourier transforming the data.  Set ``D = r.Dx`` for a first order x-derivative, ``D = r.Dy`` for a first order y-derivative, and similarly ``D = r.Dz.*r.Dy`` for a cross-derivative in ``z`` and ``y``. Higher derivatives require powers of these, for example `D = r.Dz.^4``. For higher dimensions use numerical labels, where ``D = r.Dx`` becomes ``D = r.D{2}``, and so on. Time derivatives are ignored at present. Derivatives are returned at all lattice locations.

If the derivative ``D`` is omitted, a first order x-derivative is returned.
Note that :func:`xd` returns a lattice observable, as required when used in the :func:`observe` function. If the integral is used in another function, note that it returns a matrix of dimension ``[1, lattice]``.

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


A conventional fast Fourier transform (FFT) is employed for the interaction picture (IP) transformations used in computations, as this is fast and simple. In one dimension, this is given by a sum over indices starting with zero, rather than the Matlab convention of one. Hence, if :math:`\tilde{m}=m-1`:

.. math::

    \tilde{a}_{\tilde{n}}=\mathcal{F}\left(a\right)=\sum_{\tilde{m}=0}^{N-1}a_{\tilde{m}}\exp\left[-2\pi i\tilde{m}\tilde{n}/N\right]

Suppose the spatial grid spacing is :math:`dx`, and the number of grid points is :math:`N`, then the maximum range from the first to last point is:

.. math::

    R=(N-1)dx

We note that the momentum grid spacing is

.. math::

    dk=\frac{2\pi}{Ndx}

The IP Fourier transform can be written in terms of an FFT as

.. math::

    \tilde{\boldsymbol{a}}\left(\boldsymbol{k}_{\boldsymbol{n}}\right)=\prod_{j}\left[\sum_{\tilde{m}_{j}}\exp\left[-i\left(dk_{j}dx_{j}\right)\tilde{m}_{j}\tilde{n}_{j}\right]\right]

The inverse FFT Fourier transforms automatically divide by the correct factors of :math:`\prod_{j}N_{j}` to ensure invertibility. Note also that due to the periodicity of the exponential function, negative momenta are obtained if we consider an ordered lattice such that:

.. math::

    \begin{aligned}
    k_{j} & = (j-1)dk\,\,\,(j\le N/2)\\
    k_{j} & = (j-1-N)dk\,\,(j>N/2)
    \end{aligned}
    

    
Derivatives
-----------

For calculating derivatives in the interaction picture, the notation :math:`D` indicates a derivative. To explain, one integrates by parts:

.. math::

    D^{p}\tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right)=\left[ik_{x}\right]^{p}\tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right)=\frac{1}{\left(2\pi\right)^{d/2}}\int d\boldsymbol{x}e^{-i\boldsymbol{k}\cdot\boldsymbol{x}}\left[\frac{\partial}{\partial x}\right]^{p}\boldsymbol{a}\left(\boldsymbol{x}\right)\label{eq:Fourier derivative}

This means, for example, that to calculate a one dimensional space derivative in the Linear interaction picture routine, one uses:

- :math:`\nabla_{x}\rightarrow` ``r.Dx``

Here ``r.Dx`` returns an array of momenta in cyclic order in dimension :math:`d` as defined above, suitable for an FFT calculation. The imaginary :math:`i` is not needed to give the correct sign, from the equation above. Instead, it is included in the D array. In two dimensions, the code to return a full two-dimensional Laplacian is:

- :math:`\boldsymbol{\nabla}^{2}=\nabla_{x}^{2}+\nabla_{y}^{2}\rightarrow` ``r.Dx.^2+r.Dy.^2``

Note that the dot in the notation of ``.^`` is needed to take the square of each element in the array.

Fields
======

In the xSIM code, the complex vector field ``a`` is generally stored as a compressed or flattened matrix with dimensions ``[fields, lattice]``. Here ``lattice`` is the total number of lattice points including an ensemble dimension, to increase computational efficiency:

::

    lattice = in.ensembles(1) * r.nspace

The total number of space points ``r.nspace`` is given by:

::

    r.nspace = in.points(2) * ... * in.points(in.dimension)

The use of a matrix for the fields is convenient in that fast matrix operations are possible in a high-level language.



In different subroutines it may be necessary to expand out this array to more easily reference the array structure. The full, expanded field structure ``a`` at a single time-point is as follows

:: data:: a

    [in.fieldsplus, in.ensembles(1), 1, in.points(2) ,... , in.points(dimension)] 

Note: Here, :attr:`fieldsplus` = :attr:`fields` (1) + :attr:`fields` (2) is the total number of field components and ``in.ensembles(1)`` is the number of statistical samples processed as a parallel vector. This can be set to one to save data space, or increased to improve parallel efficiency. Provided no frequency information is needed, the time dimension ``in.points(1)`` is compressed to one during calculations. During spectral calculations, the full length of the time lattice, ``in.points(1)``, is stored, which increases memory requirements.

.. data:: latt

    This includes a propagation array :attr:`r.propagator`, used in the interaction picture calculations. There are two momentum space propagators, for coarse and fine steps respectively, which are computed when they are needed.

Raw data
--------

If required, by using the switch :attr:`raw` set to one,  xSPDE can store every trajectory generated. This is raw, unprocessed data, so there is no graph index. This raw data is stored in a cell array :data:`rawdata`. The array is written to disk using the Matlab file-name, on completion, provided a file name is input.

The cell indices are: sequence index, error-checking index, ensemble index.

.. data:: rawdata

    **Cell Array**, has dimension: ``rawdata{sequence, check, in.ensemble(2)*in.ensemble(3)}``

If thread-level parallel processing is used, these are also stored in the cell array, which is indexed over both the parallel and serial ensemble. Inside each raw cell is at least one complete space-time :data:`field` stored as a complex array, with indices for the field index, the samples, and the time-space lattice. 

Each location in the cell array stores one sample-time-space trajectory in xSPDE, which is a real or complex array with (:attr:`dimension` + 2) indices, noting that :attr:`points` is a vector with :attr:`dimension` indices :

.. data:: field

    **Array**, has dimension: ``(:attr:`fieldsplus`, :attr:`ensemble(1)`, :attr:`points`)``

The main utility of raw data is for storing data-sets from large simulations for later re-analysis. It is also a platform for further development of analytic tools for third party developers, to treat statistical features not included in the functional tools provided. For example, one might need to plot histograms of distributions from this.


