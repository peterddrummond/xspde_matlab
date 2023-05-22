README

xSPDE

The xSPDE code is an extensible Stochastic Partial Differential Equation solver.

It is a stochastic toolbox for constructing simulations, which is applicable to many stochastic problems. It has a modular design which can be changed to suit different applications, and includes strategies for calculating errors. At a basic level just one or two lines of input are enough to specify the equation. For advanced users, the architecture is open and extensible in numerous ways.

Stochastic equations are equations with random noises. They occur in many fields of science, engineering, economics and other disciplines. xSPDE can solve SDEs, ODEs and PDEs as well as SPDEs: partial differential stochastic equations. These include partial spatial derivatives, like the Maxwell or Schrödinger equations.

There are many equations of this type, and xSPDE can treat a wide range. It has a configurable design. Different simulations can be carried out sequentially. Averages, moments, scatter-plots, probabilities and graphics code is included, together with chi-square tests and comparisons with analytic theory or data.

The code supports parallelism at both the vector instruction level and at the thread level, using the parallel toolbox if available. It can integrate any number of complex or real vector fields in any number of space dimensions. It uses sub-ensemble averaging and extrapolation to obtain error estimates, with up to second order (SDE) or fifth order (ODE) global convergence accuracy in time.

All algorithms used are readily extensible, if different methods are required.


xSPDE

The xSPDE folder should contain the xSIM, xGRAPH, xMANUAL and xAMPLES folder. It requires Matlab or Octave.

Set your path to include this folder and subfolders, or use an xSPDE toolbox. 

Run xAMPLES/Batchtest.m to test operation.

XSPDE runs fastest combined with the Matlab parallel toolbox.

If you have no parallel toolbox, set ensembles(3) = 1 in the examples.

Happy xspde-ing.


