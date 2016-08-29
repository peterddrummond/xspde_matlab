*******
xGRAPH
*******


The graphics program :func:`xgraph`  inputs data from :func:`xsim` simulations, then  graphs them in a variety of multidimensional graphics formats. This is the xSPDE graphics function, which uses inputs  ``data`` and, optionally, ``input``. It plots graphs, and returns the maximum difference ``diff`` from comparisons with user-specified comparison functions. 
    If required, the first argument can be a data file-name. The specified file is then read both for ``input`` and ``data``. The stored input parameters in the file can be replaced by any of the the new ``input`` parameters that are specified.


Input and data arrays
---------------------

To explain xGRAPH in full detail,

-  Data inputs are stored in the ``data`` cell array.

-  This describes a *sequence* of simulations, so that ``data = {dat1, dat2, ...}``.

-  Each structure ``data`` describes a simulation, with separate graphs for each one.

-  The main function is called using ``diff = xgraph(data[,input])``.

-  The optional second input argument is to allow modification of the graphs.

-  The diff output is available when there are comparisons made on the graphed data.

The input data sequence ``cdata`` is a cell array with a number of individual simulation objects ``dat1,...``. Each includes all the parameters that specify the simulation, together with the generated data. If there is only one simulation, just one individual ``dat`` is needed. 

Customization options
---------------------

There are a wide range of customization options available for those who wish to have the very own xGRAPH version.

Customization options include functions and parameters to allow user definition of:

- graphed functions of the input data
- plot ranges and axes for each dimension   
- plot dimensions, ie, two or thee dimensional
- plot types and line types   
- comparisons with user-specified functions
- error bars

The program will print out a record of its progress, then generate the specified graphs.

Graphics functions
--------------------------

To allow options for taking averages, these are carried out in two stages. The first type of average is a local average, taken over any function of the locally stored ensemble of trajectories. These use the :func:observe functions, specified by the user. The default is the real values of each of the fields, stored as a vector. Multiple observe functions can be used, and they are defined as a cell array of functions.

Next, any function can be taken of these local averages, using the :func:function transformations, again specified by the user. The default is the original set of local averages. This is useful if different combinations are needed of the local averages. These second level function outputs are then averaged again over a second level of ensemble averaging, if specified. This is used to obtain estimates of sampling and step-size errors in the final data outputs.



Graphics transforms
===================

All transforms defined in the observables are obtained from a cell array of vectors called :attr:`transforms`, which determines if a given coordinate axis is to be transformed prior to a given observable being measured. This can be turned on and off independently for each observable and axis. The coordinate axes are specified in the order of ``t,x,y,z,..``. This must match the :attr:`transforms` attributes used to generate the data, as no additional transform takes place.

The index ordering and normalization used in the standard discrete FFT approach is efficient for interaction picture propagation, but not useful for graphing, since graphics routines prefer the momenta to be monotonic, i.e. in the order:

.. math::

    k_{j}\left(n\right)=-K_{j}/2+(n-1)dk_{j}.

Accordingly, as explained above, all momentum indices for observable data and axes are re-ordered when graphing, although they are initially stored in the computational order.



Sequenced observables: ``data``
--------------------------------

For graphics input, cell data from each simulation in a sequence is packed into successive cells of an overall cell array :data: `data` . This is used to store the total graphics data in a sequence of simulations. All these fields are resident in memory, and can be stored for re-use. They can be re-accessed and replotted, using the :func:`xgraph` function, if required, with array dimensions:

.. :data:: data - combined graphics data from entire sequence

    **Cell Array**, has dimension: ``dat{sequence}{graph}``.

#.  graph: observable or function making up a single graph

    **Array**, has dimension: ``(components, checks, in.points(1), ... in.points(in.dimension))``.

The cell index enumerates first the sequence number and then the graph number. The second array index (``1``, ``2``, ``3``) give the error-checking status of the data. If there is no error-bar checking, the second data array is zero. If there is no sampling error checking, the third data array is zero.


In summary, observables are calculated and averaged over the ``ensembles(1)`` parallel trajectories in the :func:`xpath` function. The results are added to the earlier results in the array ``data``, to create graphs for each observable. 
There are :attr:`graphs` real observables, which are determined by the number of functions defined in the :func:`observe` cell array, unless there are additional functional transformations. The number of :attr:`graphs` may be smaller or larger than the number of vector fields. The stored cdata includes all the necessary averages over the ensembles in a complete sequence.


Graphics function
=================



:func:`xgraph` is called by xSPDE when the ensemble loops finished. The results are graphed and output if required.

.. function:: xgpreferences

    is called by :func:`xgraph` to set the graphics defaults that are not already entered.

Comparison results are calculated if available from the user-specified :attr:`compare`, an error summary is printed, and the results plotted using the :func:`xgraph` routine, which is a function that graphs the observables. It is prewritten to cover a range of useful graphs, but can be modified to suit the user. The code is intended to cascade down from higher to lower dimension, generating different types of user-defined graphs. Each type of graph is generated once for each specified graphics function.

The code is intended to cascade down from higher to lower dimension, generating different types of user-defined graphs. Each type of graph is generated once for each specified graphics function. The graphics axes that are used for plotting, and the points plotted, are defined using the optional axes input parameters, where :attr:`axes` indicates the axes preferences for n-th graph or set of generated graphical data.

If there are no :attr:`axes` inputs, or the inputs are zero - for example,
``in.axes{1} = {0,0,0}``, then only the lowest dimensions are plotted, up to 3. If the axes inputs project out a single point in a given dimension, - for example, ``axes{1}={0,31,-1,0}``, these axes are suppressed in the plots. This reduces the effective dimension of the data - in this case to two dimensions. 

Examples:

• ``axes{1}={0}``
  - For function 1, plot all the time points; higher dimensions get defaults.

• ``axes{2}={-1,0}``
  - For function 2, plot the maximum time (the default), and all x-points. The first or time axis is suppressed. 

• ``axes{3}={1:4:51,32,64}``
  - For function 3, plot every 4-th time point at x point 32, y point 64

• ``axes{4}={0,1:4:51,0}``
  - For function 4, plot all time points, every 4-th x point, and all y-points.

Note that -1 indicates a default point, which is the last point on the time axis, and the midpoint on the other axes. This has the effect of suppressing this dimension in any plots.

The pdimension input can also be used to reduce dimensionality, as this sets the maximum effective plotted dimension. For example, ``pdimension{1}=1`` means that only plots vs time are output for the first function plotted. This is equivalent to setting ``axes{1}={0,-1,-1,-1,-1}``. Note that in the following, t,x,y,z are replaced by corresponding higher dimensions if there are axes that are suppressed. Slices can be taken at any desired point, not just the midpoint. Using the standard notation of, for example, ``axes{1}={6:3:81}``, can be used to modify the starting, interval, and finishing points for complete control on the plot points.

Results depend on the value of :attr:`dimension`, or else the effective graphics dimension if axes are suppressed:

- ``4``: for the highest space dimension, only a slice through :math:`z=0` is available. This is then graphed as if it was in three dimensions.

- ``3``: for two dimensions, distinct graphic images of observable *vs x,y* are plotted at :attr:`images` time slices. Otherwise, only a slice through :math:`y=0` is available. This is then treated as if it was in two dimensions.

- ``2``: for two dimensions, one three-dimensional image of observable *vs x,t* is plotted. Otherwise, only a slice through :math:`x=0` is available. This is otherwise treated as in one dimension.

- ``1``: for one dimensions, one image of observable *vs* :math:`t` is plotted, with data at each lattice point in time. Exact results, error bars and sampling error bounds are included if available.

In addition to time-dependent graphs, the :func:`xgraph` function can generate :attr:`images` (3D) and :attr:`transverse` (2D) plots at specified points in time, up to a maximum given by the number of time points specified. The number of these can be individually specified for each graphics output. The images available are specified in :attr:`imagetype`: 3D perspective plots, grey-scale colour plots and contour plots.

Graphics user functions
=======================

:attr:`gfunction`

    This is used when a graph is needed that is a function of the data coming from the simulation package, since this data can be analysed at a later time. Error estimates are less accurate when this function is used, due to error-propagation effects that may occur after averaging, unless corrected for explicitly in the graphics function. 

:attr:`xfunctions`

    This is used when a graph is needed whose axes are a function of the original axes. 

:attr:`compare`

    This is used when a two-dimensional graph is needed with a comparison line.

Error checks
============

The final 2D output graphs will have error-bars if :attr:`checks` is set to ``1``, which is also the default parameter setting. This is to make sure the final results are accurate. Error-bars below a minimum relative size compared to the vertical range of the plot, specified by the graphics variable :attr:`minbar`, are not plotted. There is a clear strategy if the errors are too large.

Either increase the :attr:`points`, which gives more plotted points and lower errors, or increase the :attr:`steps`, which reduces the step size without changing the graphical resolution. The default algorithm and extrapolation order can be changed, read the xSPDE manual when doing this. Error bars on the graphs can be removed by setting ``in.checks = 0`` or increasing :attr:`minbar` in final graphs.

If ``in.ensembles(2) > 1`` or ``in.ensembles(3) > 1``, which allows xSPDE to calculate sampling errors, it will plot upper and lower limits of one standard deviation. If the sampling errors are too large, try increasing ``in.ensembles(1)``, which increases the trajectories in a single thread. An alternative is to increase ``in.ensembles(2)``. This is slower, but is only limited by the compute time, or else to increase ``in.ensembles(3)``, which gives higher level parallelization. Each is limited in different ways; the first by memory, and the second by time, the third by the number of available cores. Sampling error control helps ensures accuracy.

Note that error bars and sampling errors are only graphed for 2D graphs of results vs time. The error-bars are not plotted when they are below a user-specified size, to improve graphics quality. Higher dimensional graphs do not include this, for visibility reasons, but they are still recorded in the data files. Errors caused by the spatial lattice are not checked automatically in the xSPDE code. They must be checked by manually, by comparing results with different transverse lattice ranges and step-size.


Graphics projections
====================

If there is a spatial grid, the graphics program automatically generates several graphs for each observable, depending on space dimension. The maximum dimension that is plotted as set by :attr:`pdimension`. In the plots, the lattice is projected down to successively lower dimensions.

For each observable, the projection sequence is as follows:

-  If :attr:`dimension` is ``4`` or greater, a central :math:`z` point ``nz = 1 + floor(in.points(4)/2)`` is picked. For example, with 35 points, the central point is ``nz = 18``.

-  This gives a three dimensional space-time lattice, which is treated the same as if :attr:`dimension` is ``3``.

-  If :attr:`images` are specified, two-dimensional :math:`x-y` plots are generated at equally spaced time intervals. If there is only one image, it is at the last time-point. Different plot-types are used depending on the setting of :attr:`imagetype`.

-  A central :math:`y` point ``ny = 1 + floor(in.points(3)/2)`` is picked. This gives a two dimensional space-time lattice, which is treated the same as if :attr:`dimension` is ``2``. If :attr:`transverse` is specified, one-dimensional :math:`x` plots are generated at equally spaced time intervals, as before.

-  A central :math:`x` point ``nx = 1 + floor(in.points(2)/2)`` is picked. This gives a one dimensional time lattice, which is treated the same as if :attr:`dimension` is ``1``.

-  Plots of observable vs time are obtained, including sampling errors and error bars. If comparison graphs are specified using :func:`compare` functions, they are plotted also, using a dotted line. A difference graph is also plotted when there is a comparison.

