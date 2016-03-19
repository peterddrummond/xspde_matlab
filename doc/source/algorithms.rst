.. _chap-algorithms:

**********
Algorithms
**********

Stochastic, partial and ordinary differential equations are central to numerical mathematics. Certainly, ordinary differential equations have been known in some form ever since calculus was invented. There are a truly extraordinary number of algorithms used to solve these equations. One program cannot possibly provide all of them.

xSPDE currently provides four built-in choices of algorithm. All built-in methods are defined in an interaction picture. All can be used with any space dimension, including ``in.dimension = 1``, which gives an ordinary stochastic equation. All can be used either with stochastic or with non-stochastic equations. When applied to stochastic equations, the Euler method requires an Ito form of stochastic equation, while the others should be used with the Stratonovich form of these equations. Each uses the interaction picture to take care of exactly soluble linear terms.

If you have a favorite integration method that isn’t here, don’t panic. User-defined algorithms can be added freely. You can easily add your own. The existing methods are listed below, and the corresponding ``.m``-files can be used as a model. Call the routine, for example ``"myalgorithm.m"``, set ``in.step = @myalgorithm``, then adjust the value of :attr:`in.ipsteps` if the interaction-picture transform length must be changed to a new value.

Similarly, the interaction-picture transformation, :attr:`in.prop`, can also be changed if the built-in choice is not adequate for your needs.


Notation
========

The general equation treated here is given in differential form as

.. math::

    \frac{\partial\boldsymbol{a}}{\partial t} =\boldsymbol{A}\left[\boldsymbol{a}, t \right]+\underline{\mathbf{B}}\left[\boldsymbol{a}, t \right] \cdot\boldsymbol{\zeta}(t)+ \underline{\mathbf{L}}\left[\boldsymbol{\nabla}\right]\cdot\boldsymbol{a}.


It is convenient for the purposes of describing interaction picture methods,  to introduce an abbreviated notation as:
        
.. math::
    
    \begin{aligned} \mathcal{D}\left[\mathbf{a}, t \right]=\boldsymbol{A}\left[\boldsymbol{a},t \right]+\underline{\mathbf{B}}\left[\boldsymbol{a},t \right]\cdot\boldsymbol{\zeta}(t).
    \end{aligned}


Hence, we can rewrite the differential equation in the form:
    
    
.. math::

    \frac{\partial\boldsymbol{a}}{\partial t}=\mathcal{D}\left[\mathbf{a}, t \right]+\underline{\mathbf{L}}\left[\boldsymbol{\nabla}\right]\cdot\boldsymbol{a}.


Next, we define a linear propagator. This is given formally by:
    
.. math::

  \mathcal{P}\left(\Delta t \right) = \exp \left( \Delta t \underline{\mathbf{L}}\left[\boldsymbol{\nabla}\right] \right)


Typically, but not necessarily, this is evaluated in Fourier space, where it should be just a diagonal term in the momentum vector conjugate to the transverse space coordinate. It will then involve a Fourier transform, multiplication by an appropriate function of the momentum, and then an inverse Fourier transform afterwards.


xSPDE algorithms
================

The four built-in algorithms provided are:

.. function:: xEuler

    The first-order Euler method, a simple first-order explicit approach.

.. function:: xRK2

    A second order Runge-Kutta method.

.. function:: xMP

    The midpoint method: a semi-implicit, second order algorithm.

.. function:: xRK4

    A fourth order Runge-Kutta method, which is a popular ODE solver.

For simplicity, the stochastic noise is assumed constant throughout the interval :math:`dt`. The reader is referred to the literature ([Drummond1991]_, [Kloeden1995]_, [Werner1997]_, [Higham2001]_) for more details.

However, a word of caution is in order. For stochastic equations, which are non-differentiable, the classifications of convergence order should be taken *cum grano salis.* In other words, don’t believe it. Stochastic convergence is a complex issue, and the usual rules of calculus don’t apply. This is because stochastic noise is non-differentiable. It has relative fluctuations proportional to :math:`1/\sqrt{dtdV}`, for noise defined on a lattice with temporal cell-size :math:`dt` and spatial cell-size :math:`dV`. Hence the usual differentiability and smoothness properties required to give high-order convergence for standard Runge-Kutta methods are simply not present.

Higher order, more complex algorithms for stochastic integration do exist, but they are not included in the current xSPDE distribution. The reason for this is simply that stochastic integration errors are often dominated by the sampling error, which makes the practical advantage of using high-order algorithms less significant in most calculations.  

All is not completely lost however, since xSPDE will attempt to estimate both the step-size and the sampling error, so you can check convergence yourself.


Euler
=====

This is an explicit Ito-Euler method using an interaction picture. While very traditional, it is not generally recommended except for testing purposes. If it is used, very small step-sizes will generally be necessary to reduce errors to a usable level.

This is because it is is only convergent to first order, and therefore tends to have large errors. It is designed for use with an Ito form of stochastic equation. It requires one IP transform per step (``in.ipsteps = 1``). Starting from time :math:`t=t_{n}`, to get the next time point at :math:`t=t_{n+1}=t_{n}+\Delta t`,  one calculates:

.. math::

    \begin{aligned}
    \Delta\mathbf{a}_{n} & = \Delta t\mathcal{D}\left[\mathbf{a}_{n}, t_{n}\right] \\
    \mathbf{a}_{n+1} & = \mathcal{P}\left(\Delta t\right)\cdot\left[\mathbf{a}_{n}+\Delta\mathbf{a}_{n}\right]\end{aligned}


Second order Runge-Kutta
========================

This is a second order Runge-Kutta method using an interaction picture [Caradoc-Davies2000]_. It is convergent to second order in time for non-stochastic equations, but for stochastic equations it can be more slowly convergent than the midpoint method. It requires two IP transforms per step, but each is a full time-step long (``in.ipsteps = 1``).

To get the next time point, one calculates:

.. math::

    \begin{aligned}
    \bar{\mathbf{a}} & = \mathcal{P}\left(\Delta t\right)\cdot\left[\mathbf{a}_{n}\right] \\
    \mathbf{d}^{(1)} & = \Delta t\mathcal{P}\left(\Delta t\right)\cdot\mathcal{D}\left[\mathbf{a}_{n},  t_{n} \right] \\
    \mathbf{d}^{(2)} & = \Delta t\mathcal{D}\left[\bar{\mathbf{a}}+\mathbf{d}^{(1)}, t_{n+1} \right] \\
    \mathbf{a}_{n+1} & = \bar{\mathbf{a}}+\left(\mathbf{d}^{(1)}+\mathbf{d}^{(2)}\right)/2\end{aligned}


Midpoint
========

This is an implicit midpoint method using an interaction picture. It gives good results for stochastic [Drummond1991]_ and stochastic partial differential equations [Werner1997]_. While it is only convergent to second order in time for non-stochastic equations, it is strongly convergent and robust. It requires two half-length IP transforms per step (``in.ipsteps = 2``).

To get the next time point, one calculates a midpoint derivative iteratively at time to get the next time point at :math:`t=t_{n+1/2}=t_{n}+\Delta t/2`,  to give an estimated midpoint field :math:`\bar{\mathbf{a}}^{(i)}`, usually with three iterations:

.. math::

    \begin{aligned}
    \bar{\mathbf{a}}^{(0)} & = \mathcal{P}\left(\frac{\Delta t}{2}\right)\cdot\left[\mathbf{a}_{n}\right] \\
    \bar{\mathbf{a}}^{(i)} & = \bar{\mathbf{a}}^{(0)}+\frac{\Delta t}{2}\mathcal{D}\left[\bar{\mathbf{a}}^{(i-1)}, t_{n+1/2} \right] \\
    \mathbf{a}_{n+1} & = \mathcal{P}\left(\frac{\Delta t}{2}\right)\cdot\left[2\bar{\mathbf{a}}^{(i)}-\bar{\mathbf{a}}^{(0)}\right]
    \end{aligned}


Fourth order Runge-Kutta
========================

This is a fourth order Runge-Kutta method using an interaction picture [Caradoc-Davies2000]_. It is convergent to fourth order in time for non-stochastic equations, but for stochastic equations it can be more slowly convergent than the midpoint method. It requires four half-length IP transforms per step (``in.ipsteps = 2``). To get the next time point, one calculates four derivatives sequentially:

.. math::

    \begin{aligned}
    \bar{\mathbf{a}} & = \mathcal{P}\left(\frac{\Delta t}{2}\right)\cdot\left[\mathbf{a}_{n}\right] \\
    \mathbf{d}^{(1)} & = \frac{\Delta t}{2}\mathcal{P}\left(\frac{\Delta t}{2}\right)\cdot\mathcal{D}\left[\mathbf{a}_{n}, t_{n}\right] \\
    \mathbf{d}^{(2)} & = \frac{\Delta t}{2}\mathcal{D}\left[\bar{\mathbf{a}}+\mathbf{d}^{(1)}, t_{n+1/2} \right] \\
    \mathbf{d}^{(3)} & = \frac{\Delta t}{2}\mathcal{D}\left[\bar{\mathbf{a}}+\mathbf{d}^{(2)}, t_{n+1/2} \right] \\
    \mathbf{d}^{(4)} & = \frac{\Delta t}{2}\mathcal{D}\left[\mathcal{P}\left(\frac{\Delta t}{2}\right)\left[\bar{\mathbf{a}}+2\mathbf{d}^{(3)}, t_{n+1} \right]\right] \\
    \mathbf{a}_{n+1} & = \mathcal{P}\left(\frac{\Delta t}{2}\right)\cdot\left[\bar{\mathbf{a}}+\left(\mathbf{d}^{(1)}+2\left(\mathbf{d}^{(2)}+\mathbf{d}^{(3)}\right)\right)/3\right]+\mathbf{d}^{(4)}/3
    \end{aligned}

This might seem like the obvious choice, having the highest order. However, it can actually converge at a range of apparent rates, depending on the relative importance of stochastic and non-stochastic terms. Due to its reliance on differentiability, it may converge more slowly than the midpoint method with stochastic terms present.

The actual error is best judged by measuring it, as explained next.


Convergence checks
==================

To check convergence, xSPDE repeats the calculations at least twice for checking step-sizes, and many times more in stochastic cases. *If you think this is too boring and slow, turn it off.* However, you won’t know your errors!

Whatever the application, you will find the error-estimates useful. If the errors are too large, and this is relative to the application, you should decrease the time-steps or increase the number of samples. Which to do entirely depends on the type of error. In xSPDE, the step-size error due to finite time-step sizes is called the "step" error. The sampling error due to finite samples of trajectories is called the "sample" error. The maximum value of each of these, calculated over the set of all computed observables, is printed out at the end of the run.

Where there is 2D graphical output, the error bars give the step-size error, if you have ``in.check = 2``. To distinguish the error types, two lines are graphed for an upper and lower standard deviation departure from the mean, indicating the sampling error. This is only plotted if the total number of ensembles is greater than one, preferably at least 10--20 to give reliable estimates.

Note that the sample error is usually reasonably accurate. It occasionally may underestimate errors for pathological distributions. The step error is generally the more cautious of the two, and tends to overestimate errors. Neither should be relied as more than a rough guide.

As a check, the code allows users to graph a defined 2D exact result, if known, for comparison and testing purposes. These are graphed using dashed lines. This facility can be turned on or off for each observable using Boolean variables. This can be useful even if no exact result is known, but there is a known conservation law.

In summary, there are three types of convergence checks, all of which appear in the output as printed maximum values and projected two-dimensional graphs:

-  Error bars indicate the error due to finite step-size
-  Upper and lower solid lines indicate the :math:`\pm\sigma` sampling error bounds
-  Dashed lines indicate comparison values, which are useful when there are exact results for testing


Extrapolation order and error bars
==================================

For checking step-size errors, xSPDE allows the user to specify ``in.errorchecks = 2``, which is the default option. This gives one integration at the specified step-size, and one at half the specified step-size. The data is plotted at the fine step-size. The standard error-bar, with no extrapolation, has a half-size equal to the difference of fine and coarse step graphed results.

To allow for extrapolation, xSPDE allows user input of an assumed extrapolation order called :attr:`in.order`. If this is done, and errorchecks are set to 2 to allow successive integration with two different step-sizes, the output of all data graphed will be extrapolated to the specified order. In this case, the error bar half-size is set to the difference of the fine estimate and the *extrapolated* estimate.

Extrapolation is a well-known technique for improving the accuracy of a differential equation solver. Suppose an algorithm has a result with a known convergence order :math:`n`. This means that for small enough step-size, integration results :math:`R\left(dt\right)` with step-size :math:`dt` have an error of size :math:`dt^{n}`, that is:

.. math::

    R\left(dt\right)=R_{0}+E\left(dt\right)=R_{0}+k.dt^{n}.

Hence, from two results at different values of :math:`dt,` differing by a factor of :math:`2`, one would obtain

.. math::

    \begin{aligned}
    R_{1} & = R\left(dt\right)=R_{0}+k.dt^{n} \\
    R_{2} & = R\left(2dt\right)=R_{0}+2^{n}k.dt^{n}.
    \end{aligned}

The true result, extrapolated to the small-step size limit, is therefore given by giving more weight to the fine step-size result, while *subtracting* from this a correction due to the coarse step-size calculation:

.. math::

    R_{0}=\frac{\left[R_{1}-R_{2}2^{-n}\right]}{\left[1-2^{-n}\right]}.

Thus, for example, if we define a factor :math:`\epsilon` as

.. math::

    \epsilon\left(n\right)=\frac{1}{\left[2^{n}-1\right]}=\left(1,\frac{1}{3},\frac{1}{7}\ldots\right),

then the true results are obtained from extrapolation to zero step-size as:

.. math::

    R_{0}=\left(1+\epsilon\right)R_{1}-\epsilon R_{2}.

The built-in algorithms have convergence order as ordinary differential equation integrators of 1, 2, 2, 4 respectively, and should converge to this order at small step-sizes.

However, the situation is not as straightforward for stochastic equations. First order convergence is always obtainable stochastically. In addition, second order convergence is generally obtainable with the midpoint algorithm, although this is not guaranteed: it depends on the precise noise term. However, the Runge-Kutta algorithms used do **not** converge to the standard ODE order for stochastic equations. Hence extrapolation should be used with extreme caution in stochastic calculations.

While extrapolated results are usually inside those given by the default error-bars, **extrapolation with too high an order can under-estimate the resulting error bars.** Therefore, xSPDE assumes a cautious default order of ``in.order = 1``. Note that one can set ``in.order = 0`` to obtain fine resolution values and error bars without extrapolation, but this is generally less accurate.


Sampling errors
===============

Sampling error estimation in xSPDE uses sub-ensemble averaging. Ensembles are specified in three levels. The first, ``in.ensemble(1)``, is called the number of samples for brevity. All computed quantities returned by the :attr:`in.observe` functions are first averaged over the samples, which are calculated efficiently using a parallel vector of trajectories. By the central limit theorem, these sample averages are distributed as a normal distribution at large sample number.

Next, the sample averages are averaged **again** over the two higher level ensembles, if specified. This time, the variance is accumulated. The variance of these distributions is used to estimate a standard deviation in the mean, since each computed quantity is now a normally distributed result. This method is applied to all the :attr:`in.graphs` observables. The two lines generated represent :math:`\bar{o}\pm\sigma`, where :math:`o` is the observe function output, and :math:`\sigma` is the standard deviation in the mean.

The highest level ensemble, ``in.ensemble(3)``, is used for parallel simulations. This requires the Matlab parallel toolbox. Either type of high-level ensemble, or both together, can be used to calculate sampling errors.

Note that one standard deviation is not a strong bound; errors are expected to exceed this value in 32% of observed measurements. Another point to remember is that stochastic errors are often correlated, so that a group of points may all have similar errors due to statistical sampling.
