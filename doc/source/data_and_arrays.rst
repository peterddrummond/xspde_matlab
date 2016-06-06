***************
Data and arrays
***************

The xSPDE data and arrays that are user accessible are parameters ``r``, fields ``a``,  average observables ``data``, and raw trajectories ``raw``. Apart from the parameters, which are Matlab structures, all fields and data are arrays. There is a unified index ordering of:


#. field index 

#. ensemble or error index

#. time index

#. space index 1

#. space index 2

#. space index 3 ..

The number of space coordinates is arbitrary. However, to conserve storage, time dimensions are reduced to a single, current, index for propagating fields. The ensemble index can be adjusted to increase or decrease local memory usage. If needed, all data can be stored in `raw' arrays.

The fields ``a`` are complex arrays stored discretely on spatial or momentum grids. Internally, the fields are matrices stored on the flattened xSPDE internal lattice. This combines a local ensemble and a spatial grid. Temporary transformations to the Fourier domain are used for calculating interaction picture propagation [Caradoc-Davies2000]_ and averages over Fourier space. 

Two different types of Fourier representations are used. In xsim, Fourier transformations are for propagation, which requires the fastest possible methods, and uses :math:`k=0` as the first or lowest index. In xgraph, Fourier transformations are for graphical representations. Hence, the  indices are re-ordered to a conventional index ordering, with negative momentum values in the first index position.


The :ref:`parameters <sec-parameters>` are stored in a structure called, simply, ``r``. It is available to all user-definable routines. The label ``r`` is chosen because the parameters include the grid coordinates in space and time. These structures reside in a static internal cell array that combines both input and lattice parameters, including the interaction picture transformations, called :data:`latt`. The data in :data:`latt` is different for each simulation in a sequence.

Averaged results are called observables in xSPDE. For each sequence, these are stored in either space or Fourier domains, in the array ``data``, as determined by the :attr:`in.transforms` vector for each observable. This is a vector of switches for each of the space-time coordinates. The ``data`` arrays obtained in the program as calculations progress are stored in cell arrays, ``cdata``, indexed by a sequence index.

If required, ``raw`` ensemble data consisting of all the trajectories ``a`` developing in time can be stored and output. This is memory intensive, and is only done if the :attr:`in.raw` option is set to ``1``.

More details of ensembles, grids and the internal lattice are given below. Note that the term ``lattice`` is used to refer to the total internal field storage. This combines the local ensemble and the spatial grid together. 

Ensembles
================

Ensembles are used for averaging over stochastic trajectories. They come in three layers: local, serial and parallel, in order to optimize simulations for memory and for parallel operations. The ``in.ensembles(1)`` local  trajectories are used for array-based parallel ensemble averaging. These trajectories are stored in one array, to allow fast on-chip parallel processing. Distinct stochastic trajectories are also organized at a higher level into a set of ``in.ensembles(2)`` serial ensembles for statistical purposes, which allows a more precise estimate of sampling error bars. For greater speed, these can  be integrated using ``in.ensembles(3)`` parallel threads.

This hierarchical organization allows allows flexibility in allocating memory and optimizing parallel processing. It is usually faster to have larger values of ``in.ensembles(1)``, but more memory intensive. Using larger values of ``in.ensembles(2)`` is slower, but requires less memory.  Using larger values of ``in.ensembles(3)`` is fast, but requires the Matlab parallel toolbox, and uses both threads and memory resources. It is generally not effective to increase ``in.ensembles(3)`` above the maximum number of available computational cores.

In summary, the ensembles are defined as follows:

Local ensemble
--------------

The first or local ensemble contains ``ensembles(1)`` trajectories stored on the xSPDE internal lattice and processed using vector or matrix operations. 

Serial ensemble
--------------

The second or serial ensemble contains ``ensembles(2)`` of the local ensembles, processed in a sequence to conserve memory. 

Parallel ensemble
--------------
 
The third or parallel ensemble contains ``ensembles(3)`` of the serial ensembles processed in parallel using different threads to allow multi-core and multi-CPU parallel operations.


Grids in x and k
================

The xSPDE space and momentum grid can have any dimension, provided there is enough memory. Using more than six to ten total dimensions causes large time requirements and is not very practical.

The xSPDE algorithms all use a sequence of interaction pictures. Each successive interaction picture is referenced to :math:`t=t_{n}`, for the n-th step starting at :math:`t=t_{n}`, so :math:`\boldsymbol{a}_{I}(t_{n})=\boldsymbol{a}(t_{n})\equiv\boldsymbol{a}_{n}`. It is possible to solve stochastic partial differential equations in xSPDE using explicit derivatives, but this is generally less efficient. To understand spatial discretization and the interaction picture, we first must understand the xSPDE spatial grids.



Space grid
-------------

We define the grid cell size :math:`dx_{j}` in the :math:`j`-th dimension in terms of maximum range :math:`R_{j}` and the number of points :math:`N_{j}:`

.. math::

    dx_{j}=\frac{R_{j}}{N_{j}-1}.

Each grid starts at a value defined by the vector :attr:`in.origin`. Using the default values, the time grid starts at :math:`t=0` and ends at :math:`t=T=r_{1}`, for :math:`n=1,\ldots N_{j}`:

.. math::

    t\left(n\right)=(n-1)dt.

The :math:`j`-th coordinate grid starts at :math:`-r_{j}/2` and ends at :math:`r_{j}/2` , so that, for :math:`n=1,\ldots N_{j}`:

.. math::

    x_{j}\left(n\right)=-R_{j}/2+(n-1)dx_{j}.

Momentum grid
----------------

The momentum space graphs use a Fourier transform definition so that, for :math:`d` dimensions:

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
    


Computational Fourier transforms
================================

A conventional fast Fourier transform (FFT) is used for the interaction picture (IP) transformations used in computations, as this is fast and simple. In one dimension, this is given by a sum over indices starting with zero, rather than the Matlab convention of one. Hence, if :math:`\tilde{m}=m-1`:

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

For calculating derivatives and propagating in the interaction picture, the notation :math:`D` indicates a derivative. To explain, one integrates by parts:

.. math::

    D^{p}\tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right)=\left[ik_{x}\right]^{p}\tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right)=\frac{1}{\left(2\pi\right)^{d/2}}\int d\boldsymbol{x}e^{-i\boldsymbol{k}\cdot\boldsymbol{x}}\left[\frac{\partial}{\partial x}\right]^{p}\boldsymbol{a}\left(\boldsymbol{x}\right)\label{eq:Fourier derivative}

This means, for example, that to calculate a one dimensional space derivative in the Linear routine, one uses:

- :math:`\nabla_{x}\rightarrow` ``r.Dx``

Here ``r.Dx`` returns an array of momenta in cyclic order in dimension :math:`d` as defined above, suitable for an FFT calculation. The imaginary :math:`i` is not needed to give the correct sign, from the equation above. Instead, it is included in the D array. In two dimensions, the code to return a full two-dimensional Laplacian is:

- :math:`\boldsymbol{\nabla}^{2}=\nabla_{x}^{2}+\nabla_{y}^{2}\rightarrow` ``r.Dx.^2+r.Dy.^2``

Note that the dot in the notation of ``.^`` is needed to take the square of each element in the array.


Graphics transforms
===================

All transforms defined in the observables are obtained from a cell array of vectors called :attr:`in.transforms`, which determines if a given coordinate axis is transformed prior to a given observable being measured. This can be turned on and off independently for each observable and axis. The coordinate axes are specified in the order of ``t,x,y,z,..``.

The index ordering and normalization used in the standard discrete FFT approach is efficient for interaction picture propagation, but not useful for graphing, since graphics routines prefer the momenta to be monotonic, i.e. in the order:

.. math::

    k_{j}\left(n\right)=-K_{j}/2+(n-1)dk_{j}.

Accordingly, all momentum indices for observable data and axes are re-ordered when graphing, although they are initially stored in the computational order.


Fields
======

In the xSPDE code, the complex vector field ``a`` is genrally stored as a compressed or flattened matrix with dimensions ``[fields, lattice]``. Here ``lattice`` is the total number of lattice points including an ensemble dimension, to increase computational efficiency:

::

    lattice = in.ensembles(1) * r.nspace

The total number of space points ``r.nspace`` is given by:

::

    r.nspace = in.points(2) * ... * in.points(in.dimension)

The use of a matrix for the fields is convenient in that fast matrix operations are possible in a high-level language.



In different subroutines it may be necessary to expand out this array to more easily reference the array structure. The expanded field structure ``a`` is as follows

::

    [in.fields, in.ensembles(1), 1, in.points(2) ,... , in.points(dimension)] 

Note: Here, :attr:`in.fields` is the number of field components and ``in.ensembles(1)`` is the number of statistical samples processed as a parallel vector. This can be set to one to save data space, or increased to improve parallel efficiency. Provided no frequency information is needed, the time dimension ``in.points(1)`` is compressed to one during calculations. During spectral calculations, the full length of the time lattice, ``in.points(1)``, is stored, which increases memory requirements.

.. data:: latt

    This includes a propagation array :attr:`r.propagator`, used in the interaction picture calculations. There are two momentum space propagators, for coarse and fine steps respectively, which are computed when they are needed.


Data
====

Observables: ``data``
---------------------

During the calculation, observables are calculated and averaged over the ``ensembles(1)`` parallel trajectories in the :func:`xpath` function. These are determined by the functions in the :attr:`in.observe` cell array.

The number of :attr:`in.observe` functions may be smaller or larger than the number of vector fields. The observable may be a scalar or vector. These include the averages over the ensembles, and can be visualized as a single graph with one or more lines.

Next, arbitrary functional transforms can be taken, using the :attr:`in.function` cell array. These functions can use as their input any of the :attr:`in.observe` output data arrays. They default to the original :attr:`in.observe` data if they are not user-defined. The results are added to the earlier results in the array ``data``, to create graphs for each function. At this stage, both the first and second moment is stored, in order to allow calculation of the sampling error in each quantity.  These are averaged over the higher level ensembles, to allow estimates of sampling errors.

Functional transforms are most useful if one wishes to use functions which require knowledge of normalization or ensemble averages of lower-level data.

Each resulting graph or average data is each stored  in an array of size

::

    [components, errors, in.points(1), in.points(2), ... , in.points(dimension)] 

In the simplest case, there is just one vector component per average. More generally, the number of components is larger than this if there is a requirement to compare different results in one graph. Note that, unlike the propagating field, the time dimension is fully expanded.  This is necessary in order to generate outputs at each of the ``in.points(1)`` time slices. 

When step-size checking is turned on using the :attr:`in.errorchecks` flag set to ``2``, a low resolution field is stored for comparison with a high-resolution field of half the step-size, to obtain the time-step error.

The observables which are stored have three error indices which are all included in the array. These are the high resolution means, together with error-bars due to time-steps, and estimates of high-resolution standard deviations due to sampling statistics.

The observable ``data`` which is plotted automatically includes step-size error bars and plotted lines for the two estimated upper and lower standard deviations, obtained from the statistical moments.

In summary, after ensemble averaging, the second index is ``errors = 1, 2, 3``, which is used to index over the

#. mean value,

#. time-step error-bars and

#. sampling errors

respectively for each space-time point and each graphed function.

Data from each simulation in a sequence is packed into successive cells of an overall cell array :data:`cdata`. This is used to store the total data in a sequence of simulations.

All these fields are resident in memory. They can be re-accessed and replotted, using the :func:`xgraph` function, if required. In summary:

.. data:: cdata

    **Cell Array**, has dimension: ``cdata{sequence}{graph}``.

.. data:: observable or function

    **Array**, has dimension: ``(components, errors, in.points(1), ... in.points(in.dimension))``.

The cell index enumerates first the sequence number and then the graph number. The second array index (``1``, ``2``, ``3``) give the error-checking status of the data. If there is no error-bar checking, the second data array is zero. If there is no sampling error checking, the third data array is zero.

Graphics Data
=============

Observables: ``data``
---------------------

During the calculation, observables are calculated and averaged over the ``ensembles(1)`` parallel trajectories in the :func:`xpath` function. The results are added to the earlier results in the array ``data``, to create graphs for each observable. At this stage, both the first and second moment is stored, in order to allow calculation of the sampling error in each quantity.

There are :attr:`in.graphs` real observables, which are typically determined by the number of functions defined in the :attr:`in.observe` cell array, unless there are further definitions of functional transformations. The number of :attr:`in.graphs` may be smaller or larger than the number of vector fields. The observable field includes all the necessary averages over the ensembles.

Raw data
========

If required, xSPDE can store every trajectory generated.

This raw data is stored in a cell array :data:`raw`. The array is written to disk using the Matlab file-name, on completion, provided a file name is input.

The cell indices are: sequence index, error-checking index, ensemble index.

.. data:: raw

    **Cell Array**, has dimension: ``raw{seq, errcheck, in.ensemble(2)*in.ensemble(3)}``

If thread-level parallel processing is used, these are also stored in the cell array, which is indexed over both the parallel and serial ensemble. Inside each raw cell is at least one complete space-time :data:`field` stored as a complex array, with indices for the field index, the samples, the time index and the  space lattice. 

The sample-time-space trajectory in xSPDE  is a real or complex array with (in.dimension+2) indices:

.. data:: field

    **Array**, has dimension: ``(in.fields, in.ensemble(1), in.points)``

The main utility of raw data is for storing data-sets from large simulations for later re-analysis. It is also a platform for further development of analytic tools for third party developers, to treat statistical features not included in the functional tools provided. For example, one might need to plot histograms of distributions from this.
