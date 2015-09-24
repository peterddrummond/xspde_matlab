###############################################################
xSPDE: **extensible** Stochastic Partial Differential Equations
###############################################################

*Peter D. Drummond, Simon Kiesewetter*

xSPDE was developed at the Centre for Quantum and Optical Science, Swinburne University of Technology, Melbourne, Victoria, Australia.

Thanks for much valuable feedback and many useful suggestions from: Bogdan Opanchuk, Rodney Polkinghorne, Simon Kiesewetter, Laura Rosales-Zarate, King Ng and Run Yan Teh.

.. rubric:: xSPDE: a stochastic toolbox

The xSPDE code is an extensible Stochastic Partial Differential Equation solver. It is a **stochastic toolbox** for constructing simulations, which is applicable to many stochastic problems [Gardiner2004]_. It has a modular design which can be changed to suit different applications, and includes strategies for calculating errors. At a basic level just one or two lines of input are enough to specify the equation. For advanced users, the entire architecture is open and extensible in numerous ways.

Versions with an `.m` ending are written in Matlab, an interpreted scientific language of The Mathworks Inc. This version is best regarded as a prototyping platform. A code for new applications can be quickly developed and tested. This will not be quite as fast as a dedicated code but can be written easily and *understandably*.

The xSPDE logo is a three-pointed star that symbolizes that the code is suitable for all three domains: ordinary, partial or stochastic equations. Coincidentally, the logo is also similar to the three-pointed star from a Mercedes 280SE, renowned for its advanced automotive engineering.

Readers of this document may also wish to try XMDS [Colecutt2001]_, and its successor, `XMDS2 <http://sourceforge.net/projects/xmds/>`_ [Dennis2013]_, which are similar programs using XML input files.


.. rubric:: Applications

Stochastic equations are equations with random noise terms [Gardiner2004]_. They occur in in many fields of science, engineering, economics and other disciplines. xSPDE can solve both ordinary and partial differential stochastic equations. These include partial spatial derivatives, like the Maxwell or Schrödinger equations.

There are many equations of this type, and xSPDE can treat a wide range. It has a configurable functional design. The general structure permits drop-in replacements of the functions provided. Different simulations can be carried out sequentially, to simulate the various stages in an experiment or other process.

The code supports parallelism at both the vector instruction level and at the thread level, using Matlab matrix instructions and the parallel toolbox. It calculates averages of arbitrary functions of any number of complex or real fields. It uses sub-ensemble averaging and extrapolation to obtain accurate error estimates.

*Note: xSPDE is distributed without any guarantee, under the MIT open-source license. Please test it yourself before use.*


.. rubric:: Source code

The source code for xSPDE is available at `github.com/peterddrummond/xspde <https://github.com/peterddrummond/xspde>`_.


.. toctree::
   :maxdepth: 2

   interactive
   projects
   tutorial
   logic_and_data
   arrays
   algorithms
   api
   references


******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
