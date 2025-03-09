.. _chap-interactive:

*****************
Interactive xSPDE
*****************

All xSPDE simulations require parameters stored in an input structure used by the xSPDE toolbox. Inputs have default values, which can be changed by creating fields in the input structure. A complete list of parameters and their uses is given in :ref:`chap-api`. 

The simplest way to use xSPDE is via the interactive Matlab command window, illustrated in this chapter.  In this mode of operation, any required parameters are entered into the input structure  ``p``, and then the command ``xspde(p)`` runs and graphs the simulation.
Note that one must enter ``clear`` first to erase previous data in the Matlab workspace, unless the previous data is being recycled. 

Stochastic equations
====================

An ordinary stochastic equation  [Kloeden1995]_ for a real or complex vector :math:`\boldsymbol{a}` is:

.. math::

    \frac{\partial\boldsymbol{a}}{\partial t}=\dot{\boldsymbol{a}}=\boldsymbol{A}\left(\boldsymbol{a}\right)+\underline{\mathbf{B}}\left(\boldsymbol{a}\right)\boldsymbol{\zeta}(t).

Here :math:`\boldsymbol{A}` is a vector function, :math:`\underline{\mathbf{B}}` a matrix, and :math:`\boldsymbol{\zeta}` is a real Gaussian distributed noise vector such that :math:`\left\langle \boldsymbol{\zeta}\right\rangle = 0`, and :

.. math::

    \left\langle \zeta_{i}\left(t\right)\zeta_{j}\left(t\right)\right\rangle = \delta\left(t-t^\prime\right)\delta_{ij}.

xSPDE simulations
=================


To simulate a stochastic equation like this interactively with xSPDE, first make sure Matlab path is pointing to the xSPDE folder and type ``clear`` to clear old data.

Next, enter the xspde parameters (see :ref:`chap-api`) into the command window, as follows:

::

    p.label1 = <parameter1>
    p.label2 = ...
    p.deriv = @(a,w,p) <expression for da/dt>
    xspde(p)

- The notation ``p.label = parameter`` creates a field in the structure ``p`` (which is created dynamically, if it has not been defined before).
- The notation ``@(..)`` is the Matlab shorthand for an anonymous function.
- The parameters passed to the :func:`deriv` function are: ``a``, the stochastic variable; ``w``,  the random noise; and ``p``, the input structure with additional coordinates and parameters.
- :func:`xspde` is called with the input structure as the only argument.
- parameters or functions that are omitted are replaced with default values.
- a sequence of simulations requires an input list: `p1,p2..`.

Once the simulation is completed, xSPDE will generate graphs of averages of any required observables or moments. If needed, simulated data is stored in specified files for later use, which requires a file-name to be entered.


Interactive examples
====================

The random walk
---------------

The first example is the simplest possible stochastic equation:

.. math::

    \dot{a}=w(t),

with a complete xSPDE script in Matlab below, and output in :numref:`fig-simplest-case-wiener`.

::

    p.deriv = @(a,w,p) w; xspde(p);

- Here :func:`deriv` defines the derivative function. The notation @(a,w,p) defines an inline function. In this case the derivative equals the noise ``w``. Other parameters have their default values.
- The last argument of xspde user functions is  ``p``, containing the parameters required for the simulation.

.. _fig-simplest-case-wiener:
.. figure:: Figures/Wiener_1.*

   The simplest case: a random walk.

Laser quantum noise
-------------------

Next we treat a model for the quantum noise of a single mode laser:

.. math::

    \dot{a}=\left(1-\left|a\right|^{2}\right)a+b\zeta(t),

where :math:`\zeta=\left(w_{1}+iw_{2}\right)`, so that:

.. math::

    \left\langle \zeta(t)\zeta^{*}(t^\prime)\right\rangle =2\delta\left(t-t^\prime\right).

Here the coefficient :math:`b` describes the quantum noise of the laser, and is inversely proportional to the equilibrium photon number.  An xSPDE script in Matlab is given below, for the case of :math:`b=0.01`, with an output graph in :numref:`fig-laser`. Note the use of the clear command here to clean up the Matlab workspace before the start.

::

    clear
    p.noises = 2;
    p.observe = @(a,p) abs(a)^2;
    p.olabels = '|a|^2';
    p.deriv = @(a,w,p) (1-abs(a)^2).*a+0.01*(w(1)+i*w(2));
    xspde(p)

.. _fig-laser:
.. figure:: Figures/Laser.*

   Simulation of the stochastic equation describing a laser turning on.

Note that:

- :attr:`noises` is the number of noises,
- :func:`observe` is the graphed function,
- :attr:`olabels` gives the axis label.


Ito and Stratonovich equations
==============================

The xSPDE toolbox is primarily designed to treat Stratonovich equations [Gardiner2004]_, which are the broad-band limit of a finite band-width random noise equation, with derivatives evaluated at the midpoint in time of a time-step.


An equivalent type of stochastic equation is the Ito form. This is written in a similar way to a Stratonovich equation, except that this corresponds to a limit where derivatives are evaluated at the start of each step. To avoid confusion, we can write an Ito equation as a difference equation:

.. math::

    d\boldsymbol{a}=\boldsymbol{A}^{I}\left[\boldsymbol{a}\right]+\underline{\mathbf{B}}\left[\boldsymbol{a}\right]\cdot d\boldsymbol{w}(t).

Here:

.. math:: 

 \left\langle dw_{i}\left(\boldsymbol{x}\right) dw_{j}\left(\boldsymbol{x}^\prime\right)\right\rangle =\delta_{ij}dt. 

When :math:`\mathbf{\mathsf{B}}` is not a constant, the Ito drift term is different to the Stratonovich one. This difference occurs because the noise term is non-differentiable. The relationship is that

.. math::

    A_{i} = A_{i}^{I}-\frac{1}{2}\sum_{j,m}\frac{\partial B_{ij}}{\partial a_{m}}B_{mj}.
    
Provided the noise coefficient :math:`B` is constant - which is called additive noise - there is no real difference between the two types of equation. Otherwise, it is essential to know which type of stochastic equation it is, in order to get unambiguous results!

Financial calculus
------------------

The Black-Scholes equation is a well-known Ito stochastic equation, used to price financial options. It describes the fluctuations in a stock value:

.. math::

    da=\mu a\,dt+\sigma a\,dw,

where :math:`\left\langle dw^{2}\right\rangle =dt`. Since the noise is multiplicative, the equation is different in Ito and Stratonovich forms of stochastic calculus. The corresponding Stratonovich equation, as used in xSPDE is:

.. math::

    \dot{a}=\left(\mu-\sigma^{2}/2\right)a+\sigma a\,w(t).

An interactive xSPDE script in Matlab is given below with an output graph in :numref:`fig-black-scholes`, for the case of a volatile stock with :math:`\mu=0.1`, :math:`\sigma=1`. Note the spiky behaviour, typical of multiplicative noise, and also of the risky stocks in the small capitalization portions of the stock market.

::

    clear
    p.initial = @(rv,p) 1;
    p.deriv = @(a,w,p) -0.4*a+a*w;
    xspde(p)

.. _fig-black-scholes:
.. figure:: Figures/Black-Scholes.*

   Simulation of the Black-Scholes equation describing stock prices.

-  Here :func:`initial` describes the initialization function.
-  The first argument of ``@(v,p)`` is ``v``, an initial random variable.
-  The error-bars are estimates of step-size error.
-  Errors can be reduced by using more time-steps: see :ref:`chap-projects`.

This graph is of a single stochastic realisation. Generation of averages is also straightforward. This is described in :ref:`chap-projects`.


Stochastic partial differential equations
=========================================

More generally, xSPDE solves [Werner1997]_ a stochastic partial differential equation for a complex vector field defined with arbitrary transverse dimension. The total dimensions is then input as :math:`d=2,3,...`. Equations of this type occur in many disciplines, including biology, chemistry, engineering and physics. They are in differential form as

.. math::

    \frac{\partial\boldsymbol{a}}{\partial t}=\boldsymbol{A}\left[\boldsymbol{a}\right]+\underline{\mathbf{B}}\left[\boldsymbol{a}\right]\cdot\boldsymbol{\zeta}(t)+\underline{\mathbf{L}}\left[\boldsymbol{\nabla},
    \boldsymbol{a}\right]

Here :math:`\boldsymbol{a}` is a real or complex vector or vector field. The initial conditions are arbitrary functions. :math:`\boldsymbol{A}\left[\boldsymbol{a}\right]` and :math:`\underline{\mathbf{B}}\left[\boldsymbol{a}\right]` are vector and matrix functions of :math:`\boldsymbol{a}`, :math:`\underline{\mathbf{L}}\left[\boldsymbol{\nabla},\boldsymbol{a}\right]` is a function of the fields and their derivatives, and :math:`\mathbf{\boldsymbol{\zeta}}=\left[\boldsymbol{\zeta}^{x},\boldsymbol{\zeta}^{k}\right]` are real delta-correlated noise fields such that:

.. math::

    \begin{split}
    \left\langle \zeta_{i}^{x}\left(t,\boldsymbol{x}\right)\zeta_{j}^{x}\left(t,\boldsymbol{x}^\prime\right)\right\rangle  & = \delta\left(\boldsymbol{x}-\boldsymbol{x}^\prime\right)\delta\left(t-t^\prime\right)\delta_{ij}\nonumber \\
    \left\langle \zeta_{i}^{k}\left(t,\boldsymbol{k}\right)\zeta_{j}^{k}\left(t,\boldsymbol{k}^\prime\right)\right\rangle  & = f(\boldsymbol{k})\delta\left(\boldsymbol{k}-\boldsymbol{k}^\prime\right)\delta\left(t-t^\prime\right)\delta_{ij}.\end{split}

Note that the x and k noise term for each value of the index are generated from the same underlying stochastic process. This is necessary because there are some equations that require both a filtered and unfiltered noise generated from the same underlying random number distribution. If these correlations are not wanted, and the noises are required to be independent, then different noise indices must be used. 

Transverse boundary conditions are assumed periodic as the default option, which allows the use of efficient spectral Fourier transform propagation codes. Other types of boundary conditions available are Neumann boundaries with zero normal derivatives, and Dirichlet boundaries with zero fields at the boundary. These require the use of finite difference methods. The boundary type can be individually specified in each axis direction. The term :math:`\underline{\mathbf{L}}\left[\boldsymbol{\nabla}\right]` may be omitted, as space derivatives can also be treated directly in the derivative function, and this is necessary with Neumann or Dirichlet boundaries. The momentum filter :math:`f(\boldsymbol{k})` is an arbitrary user-specified function, allowing for spatially correlated noise.

To treat stochastic partial differential equations or SPDEs, in the most efficient  way, the equations are divided into the first two terms, which are essentially an ordinary stochastic equation, and the last term which is assumed linear, and therefore gives a linear partial differential equation:

.. math::

    \frac{\partial\boldsymbol{a}}{\partial t}=\underline{\mathbf{L}}\left[\boldsymbol{\nabla}\right]\cdot\boldsymbol{a}

The *interaction picture* is a moving reference frame used to solve the linear part of the equation exactly, defined by an exponential transformation. This is carried out internally by matrix multiplications and Fourier transforms.

In more detail, in Fourier space, if :math:`\tilde{\boldsymbol{a}}\left(\boldsymbol{k}\right)=\mathcal{F}\left[\boldsymbol{a}\left(\mathbf{x}\right)\right]` is the Fourier transform of :math:`\boldsymbol{a}`, we simply define:

.. math::

    \tilde{\boldsymbol{a}}(\boldsymbol{k},dt)=\mathcal{P}\left(\boldsymbol{k},dt\right)\mathbf{\tilde{a}}_{I}\left(\boldsymbol{k},dt\right)

where the propagation function can be written intuitively as :math:`\mathcal{P}=\exp\left[\underline{\mathbf{L}}(\mathbf{D})dt\right]`, where :math:`\mathbf{D}=i\boldsymbol{k}\sim\nabla`. The function :math:`\underline{\mathbf{L}}(\mathbf{D})` is input using the xSPDE linear response function :func:`linear`. With this definition, at each step the equation that is solved can be re-written in a more readily soluble form as:

.. math::

    \frac{\partial\boldsymbol{a}_{I}}{\partial t}=\mathcal{D}\left[\mathcal{F}^{-1}\mathcal{P}\left(\mathcal{F}\boldsymbol{a}_{I}\right)\right]

The total derivative in the interaction picture is the xSPDE derivative function :func:`deriv`:

.. math:: \dot{\boldsymbol{a}}_{I}=\boldsymbol{A}+\underline{\mathbf{B}}\,\boldsymbol{\zeta}

where usually :math:`\boldsymbol{A}`, :math:`\underline{\mathbf{B}}` are evaluated at the midpoint which is the origin in the interaction picture.  For convenience, the final output is calculated in the original picture, with at least two interaction picture (IP) transformations per time-step.


Note that there are many other types of partial differential equation that can be treated with xSPDE, even if the interaction picture method doesn't apply. This occurs when there are nonlinear functions with arbitrary derivatives, or derivatives that are non-diagonal in the vector indices. For these cases, the space derivatives are evaluated inside the derivative term :func:`deriv`:. If there are higher order time derivatives as well, these can be re-expressed as a set of first-order time derivatives, provided the problem is an initial-value problem.


Symmetry breaking
-----------------

An example of a SPDE with space-time dimensions of :math:`d=3`  is the stochastic Ginzburg-Landau equation. This describes symmetry breaking, in which the system develops a spontaneous phase which can vary spatially. The model is widely used in fields ranging from lasers to magnetism, superconductivity, superfluidity and even particle physics:

.. math::

    \dot{a}=\left(1-\left|a\right|^{2}\right)a+b\zeta(t)+ic\nabla^{2}a

where

.. math::

    \left\langle \zeta(x)\zeta^{*}(x^\prime)\right\rangle =2\delta\left(t-t^\prime\right)\delta\left(x-x^\prime\right).

An xSPDE script is given below, for parameter values of :math:`b=0.001` and :math:`c=0.01`, with the output graphed in :numref:`fig-symmetry-breaking`. Note that in the graph, the range ``-5<x<5`` is the default xSPDE coordinate range, while
the ``.*`` notation is used in functions here, as fields require element-wise multiplication.

::

    clear
    p.noises = 2;
    p.dimensions = 3;
    p.steps = 10;
    p.linear = @(p) i*0.01*(p.Dx.^2+p.Dy.^2);
    p.observe = @(a,~) abs(a).^2;
    p.olabels = '|a|^2';
    p.deriv = @(a,w,~) (1-abs(a).^2).*a+0.001*(w(1,:,:)+i*w(2,:,:));
    xspde(p)

Here:

- :attr:`dimensions` is the space-time dimensions, with an :math:`x-t` plot given here.
- :attr:`steps` gives the integration steps per plot-point, for improved accuracy.
- :attr:`linear` is the linear operator --- an imaginary laplacian
- ``p.Dx`` indicates a derivative operation, :math:`\partial/\partial x`. See the reference entry for :attr:`linear` for more information.

.. _fig-symmetry-breaking:
.. figure:: Figures/GinzLand.*

   Simulation of the stochastic equation describing symmetry breaking in two dimensions. Spatial fluctuations are caused by the different phase-domains that interfere. The graph obtained here is projected onto the :math:`y=0` plane.

