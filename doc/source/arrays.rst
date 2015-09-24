******
Arrays
******

The two important internal types of data that are user accessible are: ``a`` and ``r``.

The fields ``a`` are complex arrays defined on spatial or momentum grids. Internally, the fields ``a`` are just matrices stored on a flattened space lattice, except for temporary transformations to the Fourier domain for calculating interaction picture propagation [Caradoc-Davies2000]_. Two different types of Fourier representations are used, depending on whether the transformations are for propagation, which requires the fastest possible methods, or whether the transformation is for graphing, which uses a conventional index ordering.

If required, ``raw`` ensemble data consisting of all the trajectories ``a`` developing in time can be stored and output. This is memory intensive, and is only done if the :attr:`in.raw` option is set to ``1``.

The :ref:`lattice data <sec-lattice-structure>` is a structure called, simply, ``r``. It is available to all user-definable routines. The label ``r`` is chosen because the lattice parameters have the important function of storing the grid coordinates in space and time. These structures reside in a static internal cell array of inputs and lattice parameters, including the interaction picture transformations, called :data:`latt`. The data in :data:`latt` is different for each simulation in a sequence.

Averaged results for each sequence are stored in either space or Fourier domains, in the array ``data``, as determined by the :attr:`in.transforms` vector for the observable. The ``data`` arrays obtained in the program as calculations progress are stored in cell arrays, ``cdata``, indexed by a sequence index.

The internal spatial grid definitions are as follows:

Grids in x and k
================

The algorithms all use a sequence of interaction pictures. Each successive interaction picture is referenced to :math:`t=t_{n}`, for the n-th step starting at :math:`t=t_{n}`, so :math:`\boldsymbol{a}_{I}(t_{n})=\boldsymbol{a}(t_{n})\equiv\boldsymbol{a}_{n}`. To understand the interaction picture, we first must understand the xSPDE lattice.

Space lattice
-------------

We define the lattice cell size :math:`dx_{j}` in the :math:`j`-th dimension in terms of maximum range :math:`R_{j}` and the number of points :math:`N_{j}:`

.. math::

    dx_{j}=\frac{R_{j}}{N_{j}-1}.

Each lattice starts at a value defined by the vector :attr:`in.origin`. Using the default values, the time lattice starts at :math:`t=0` and ends at :math:`t=T=r_{1}`, for :math:`n=1,\ldots N_{j}`:

.. math::

    t\left(n\right)=(n-1)dt.

The :math:`j`-th coordinate lattice starts at :math:`-r_{j}/2` and ends at :math:`r_{j}/2` , so that, for :math:`n=1,\ldots N_{j}`:

.. math::

    x_{j}\left(n\right)=-R_{j}/2+(n-1)dx_{j}.

Momentum lattice
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

while the momentum lattice starts at :math:`-k_{j}/2` and ends at :math:`k_{j}/2` , so that:

.. math::

    k_{j}\left(n\right)=-K_{j}/2+(n-1)dk_{j}.


Computational Fourier transforms
================================

A conventional fast Fourier transform (FFT) is used for the interaction picture (IP) transformations used in computations, as this is fast and simple. In one dimension, this is given by a sum over indices starting with zero, rather than the Matlab convention of one. Hence, if :math:`\tilde{m}=m-1`:

.. math::

    \tilde{a}_{\tilde{n}}=\mathcal{F}\left(a\right)=\sum_{\tilde{m}=0}^{N-1}a_{\tilde{m}}\exp\left[-2\pi i\tilde{m}\tilde{n}/N\right]

Suppose the lattice spacing is :math:`dx`, and the number of lattice points is :math:`N`, then the maximum range from the first to last point is:

.. math::

    R=(N-1)dx

We note that the momentum lattice spacing is

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

- :math:`\nabla_{x}\rightarrow` ``D.x``

Here ``D.x`` returns an array of momenta in cyclic order in dimension :math:`d` as defined above, suitable for an FFT calculation. The imaginary :math:`i` is not needed to give the correct sign, from the equation above. Instead, it is included in the D array. In two dimensions, the code to return a full two-dimensional Laplacian is:

- :math:`\boldsymbol{\nabla}^{2}=\nabla_{x}^{2}+\nabla_{y}^{2}\rightarrow` ``D.x.^2+D.y.^2``

Note that the dot in the notation of ``.^`` is needed to take the square of each element in the array.


Graphics transforms
===================

The index ordering and normalization used in the standard discrete FFT approach is efficient for interaction picture propagation, but not useful for graphing, since graphics routines prefer the momenta to be monotonic, i.e. in the order:

.. math::

    k_{j}\left(n\right)=-K_{j}/2+(n-1)dk_{j}.

All transforms defined in the observables are obtained from a vector called :attr:`in.transforms`, which determines if a given coordinate axis is transformed prior to a given observable being measured. This can be turned on and off independently for each observable. The space and time transforms are defined in the next sub-sections.

Time transforms
---------------

We define a Fourier transform in time as:

.. math::

   \tilde{\boldsymbol{a}}\left(\omega\right) = \int_{0}^{T}\frac{dt}{\left(2\pi\right)^{1/2}}\boldsymbol{a}\left(t\right)\exp\left[i\omega t\right]

To achieve this in one dimension, note that:

.. math::

    dtd\omega=\frac{2\pi}{N}.

Defining :math:`\kappa=\omega^{\max}/2`:

.. math::

   \begin{aligned}
   \tilde{\boldsymbol{a}}\left(\omega_{\tilde{n}}\right) & = \frac{dt}{\sqrt{2\pi}}\sum_{\tilde{m}=0}^{N-1}\exp\left[i\left(\tilde{n}d\omega-\kappa\right)\left(\tilde{m}dt\right)\right]a_{\tilde{m}} \\
    & = \frac{dt}{\sqrt{2\pi}}\sum_{\tilde{m}=0}^{N-1}\exp\left[i\left(\tilde{n}\tilde{m}dtd\omega-\kappa\tilde{m}dt\right)\right]a_{\tilde{m}} \\
    & = \frac{dt}{\sqrt{2\pi}}\sum_{\tilde{m}=0}^{N-1}\exp\left[2\pi i\tilde{m}\tilde{n}/N\right]e^{-i\kappa\tilde{m}dt}a_{\tilde{m}} \\
    & = \frac{\sqrt{2\pi}}{d\omega}\times\frac{1}{N}\sum_{\tilde{m}=0}^{N-1}\exp\left[2\pi i\tilde{m}\tilde{n}/N\right]e^{-i\kappa t}a_{\tilde{m}}\end{aligned}

Hence to get an ordered Fourier transform for graphing data, with the usual physics and mathematics definitions, we must **premultiply** by :math:`e^{-i\kappa t}Ndt/\sqrt{2\pi}`, then take a discrete IFFT. This is taken care of internally in the xSPDE transform routines.

Space transforms
----------------

We define a Fourier transform in space as:

.. math::

   \tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right) = \int\frac{dV}{\left(2\pi\right)^{d/2}}\boldsymbol{a}\left(\boldsymbol{x}\right)\exp\left[-i\boldsymbol{k}\cdot\boldsymbol{x}\right]

To achieve this with an FFT in one dimension, let :math:`\kappa=K/2`,
and :math:`\rho=R/2`, and define:

.. math::

   \begin{aligned}
   \tilde{\boldsymbol{a}}\left(k_{\tilde{n}}\right) & = \frac{dx}{\sqrt{2\pi}}\sum_{\tilde{m}=0}^{N-1}\exp\left[-i\left(\tilde{n}dk-K/2\right)\left(\tilde{m}dx-R/2\right)\right]a_{\tilde{m}} \\
    & = \frac{dx}{\sqrt{2\pi}}\sum_{\tilde{m}=0}^{N-1}\exp\left[-i\left(\tilde{n}\tilde{m}dxdk-K\tilde{m}dx/2-R\tilde{n}dk/2+KR/4\right)\right]a_{\tilde{m}} \\
    & = e^{-iR\left(K/2-\tilde{n}dk\right)/2}\frac{dx}{\sqrt{2\pi}}\sum_{\tilde{m}=0}^{N-1}\exp\left[-2\pi i\tilde{m}\tilde{n}/N\right]e^{iK\tilde{m}dx/2}a_{\tilde{m}} \\
    & = e^{-iRk/2}\frac{dx}{\sqrt{2\pi}}\sum_{\tilde{m}=0}^{N-1}\exp\left[-2\pi i\tilde{m}\tilde{n}/N\right]e^{iKx/2+iKR/4}a_{\tilde{m}}\end{aligned}

Hence to get an ordered Fourier transform for graphing data, with the usual mathematical definitions, we must **premultiply** by a phase factor, take a discrete FFT, then **post-multiply** by *another* phase factor. This is taken care of internally in the xSPDE transform routines.


Fields
======

In the xSPDE code, the complex vector field ``a`` is stored as a complex matrix with dimensions ``[fields, lattice]``. Here ``lattice`` is the total number of lattice points including an ensemble dimension, to increase computational efficiency:

::

    lattice = in.ensembles(1) * r.n.space

The total number of space points ``r.n.space`` is given by:

::

    r.n.space = in.points(2) * ... * in.points(in.dimension)

The use of a matrix for the fields is convenient in that fast matrix operations are possible in a high-level language.

The ``in.ensembles(1)`` trajectories are used for array-based parallel ensemble averaging. These trajectories are stored in parallel in one array, to allow fast on-chip parallel processing. Distinct stochastic trajectories are also organized at a higher level into a set of ``in.ensembles(2)`` ensembles for statistical purposes, which allows a more precise estimate of sampling error bars. These can also be integrated in parallel using ``in.ensembles(3)`` parallel threads.

This hierarchical organization allows allows flexibility in allocating memory and optimizing parallel processing. It is usually faster to have larger values of ``in.ensembles(1)``, but more memory intensive. Using larger values of ``in.ensembles(2)`` is slower, but requires less memory.

In different subroutines it maybe necessary to expand out this array to more easily reference the array structure. The expanded structure is as follows

**Array** ``a`` has dimension: ``(in.fields, in.ensembles(1), in.points(2), ..., in.points(in.dimension))``.

Note: Here, :attr:`in.fields` is the number of field components and ``in.ensembles(1)`` is the number of statistical samples processed as a parallel vector. This can be set to one to save data space, or increased to improve parallel efficiency. Provided no frequency information is needed, the time dimension ``in.points(1)`` is compressed to one during calculations. During spectral calculations, the full length of the time lattice, ``in.points(1)``, is stored, which increases memory requirements.

.. data:: latt

    This includes a propagation array :attr:`r.propagator`, used in the interaction picture calculations. There are two momentum space propagators, for coarse and fine steps respectively, which are computed when they are needed.


Data
====

Observables: ``data``
---------------------

During the calculation, observables are calculated and averaged over the ``ensembles(1)`` parallel trajectories in the :func:`xpath` function. The results are added to the earlier results in the array ``data``, to create graphs for each observable. At this stage, both the first and second moment is stored, in order to allow calculation of the sampling error in each quantity.

There are :attr:`in.graphs` real observables, which are determined by the number of functions defined in the :attr:`in.observe` cell array. The number of :attr:`in.graphs` may be smaller or larger than the number of vector fields. The observable field includes all the necessary averages over the ensembles.

When step-size checking is turned on using the :attr:`in.errorchecks` flag set to ``2``, a low resolution field is stored for comparison with a high-resolution field of half the step-size, to obtain the time-step error.

The observable ``data`` which is stored therefore involves three arrays which are all included in the data array. These are the high resolution means, together with error-bars due to time-steps, and estimates of high-resolution standard deviations due to sampling statistics.

The observable ``data`` which is plotted therefore includes step-size error bars and plotted lines for the two estimated upper and lower standard deviations, obtained from the statistical moments.

In summary, data from each simulation is stored internally in an array of size

::

    errors * in.points(1) * n.space * in.graphs.

This is a flattened version of the full data dimension, which is logically

::

    errors * in.points(1) * ... * in.points(dimension) * in.graphs

This is necessary in order to generate outputs at each of the ``in.points(1)`` time slices. Here ``errors = 1, 2, 3`` is used to index over the

#. mean value,

#. time-step error-bars and

#. sampling errors

respectively for each space-time point and each graphed function.

Data from each simulation in a sequence is packed into successive cells of an overall cell array :data:`cdata`. This is used to store the total data in a sequence of simulations.

All these fields are resident in memory. They can be re-accessed and replotted, using the :func:`xgraph` function, if required. In summary:

.. data:: cdata

    **Cell Array**, has dimension: ``cdata{sequence}``.

.. data:: data

    **Array**, has dimension: ``(errors, in.points(1), ... in.points(in.dimension), in.graphs)``.

The cell index enumerates the sequence number. The first array index (``1``, ``2``, ``3``) give the error-checking status of the data. If there is no error-bar checking, the second data array is zero. If there is no sampling error checking, the third data array is zero.

Raw data
========

Although the quantity of data generated can be overwhelming, xSPDE can store every trajectory generated if asked to do so.

This raw data is stored in a cell array :data:`raw`. The array is written to disk using the Matlab file-name, on completion, provided a file name is input.

The cell indices are: the ensemble index, the error-checking index and the sequence index.

.. data:: raw

    **Cell Array**, has dimension: ``raw{ensemble, err, seq}``

Inside each cell is at least one complete space-time :data:`field` stored as a complex array, with indices for the field index, the sample-space lattice, and the time index. The sample-space lattice structure internal to xSPDE means that a subensemble of individual stochastic fields is integrated in parallel. These are defined as a real or complex array:

.. data:: field

    **Array**, has dimension: ``(in.fields, in.lattice, in.points(1))``

While this is a lengthy description, and an even larger array, it is also necessary if all the raw data needs to be extracted.

The main utility of the raw data is to provide a platform for further development of analytic tools for third party developers, to treat statistical features not included in the functional tools provided. For example, the basic xSPDE package does not provide histograms of distributions.
