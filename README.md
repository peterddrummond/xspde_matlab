# xspde_matlab

The xSPDE code is an extensible Stochastic Partial Differential Equation solver. It is a stochastic toolbox for constructing simulations, which is applicable to many stochastic problems. It has a modular design which can be changed to suit different applications, and includes strategies for calculating errors. At a basic level just one or two lines of input are enough to specify the equation. For advanced users, the entire architecture is open and extensible in numerous ways.

Stochastic equations are equations with random noise terms. They occur in in many fields of science, engineering, economics and other disciplines. xSPDE can solve both ordinary and partial differential stochastic equations. These include partial spatial derivatives, like the Maxwell or Schrödinger equations. There are many equations of this type, and xSPDE can treat a wide range. It has a configurable functional design. The general structure permits drop-in replacements of the functions provided. Different simulations can be carried out sequentially, to simulate the various stages in an experiment or other process.

The code supports parallelism at both the vector instruction level and at the thread level, using Matlab matrix instructions and the parallel toolbox. It calculates averages of arbitrary functions of any number of complex or real fields. It uses sub-ensemble averaging and extrapolation to obtain accurate error estimates.


xSPDE

The xSPDE folder contains the xSIM, xGRAPH and xAMPLES folder. It requires Matlab or Octave.

An extensive user's guide is in provided in the doc subfolder in pdf and html form.

Set your path to include this folder and subfolders, or use an xSPDE toolbox. 

Run xAMPLES/Batchtest.m to test operation.

XSPDE runs fastest combined with the Matlab parallel toolbox.

If you have no parallel toolbox, set ensembles(3) = 1 in the examples.

Happy xspde-ing.


