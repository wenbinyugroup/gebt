.. gebt documentation master file, created by
   sphinx-quickstart on Wed Nov 20 10:02:56 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GEBT
==================

GEBT (Geometrically Exact Beam Theory) is a code implementing the mixed variational formulation of the geometrically exact intrinsic beam theory developed by Prof. Hodges of Georgia Institute of Technology [hodges1990mixed]_, which captures all the geometrical nonlinearities obtainable by a beam model.
It is a companion code of VABS, a general-purpose cross-sectional analysis tool, to enable efficient yet high-fidelity analysis of slender structures, whether they are made of composite materials or not.

For effective design space explorations, we need to simplify the original nonlinear three-dimensional (3D) analysis of slender structures into a two-dimensional (2D) cross-sectional analysis and a one-dimensional (1D) nonlinear beam analysis.
This has been achieved using VABS along with GEBT.
VABS takes a finite element mesh of the cross-section including all the details of geometry and material as inputs to calculate the sectional properties including structural and inertial properties.
These properties are needed for GEBT to predict the global behavior of the slender structure.
The 3D pointwise displacement/strain/stress distribution within the structure can also be recovered by providing outputs from GEBT into VABS.




Features
--------

GEBT adopts the mixed variational formulation of the geometrically exact beam theory of Prof. Hodges.
It uses the lowest possible shape functions and the element matrices are calculated exactly without numerical integration.
Since it is a mixed variational formulation, the complete set of variables can be directly calculated.
For example, for static analysis, we can obtain three displacements, three rotations, three forces, and three moments.
GEBT can treat an arbitrary assembly of beams made of arbitrary material and oriented arbitrarily in the 3D space.
GEBT uses dynamic link libraries (DLLs) to encapsulate the analysis capability so that it has true plug-n-play capability which is convenient for integration into other environments.
Now GEBT can be used both as a standalone application and a callable library.
A GEBT manual for developers can be requested separately.




Functionalities
---------------

The most recent version of GEBT, GEBT 4.0, has the following functionalities:

1. Static, both linear and nonlinear, analysis of beam assemblies with straight and/or initially curved/twisted members under prescribed concentrated or distributed forces/moments or displacements/rotations.
2. Static, both linear and nonlinear, analysis of beam assemblies with straight members and/or initially curved/twisted with linearly varying sectional properties.
3. Static, nonlinear, analysis of beam assemblies with straight and/or initially curved/twisted members under prescribed distributed or concentrated follower forces/moments.
4. Steady state of the dynamic response of beam assemblies with straight and/or initially curved/twisted members under prescribed distributed or concentrated dead or follower forces/moments.
5. Transient dynamic response of beam assemblies with straight and/or initially curved/twisted members under prescribed distributed or concentrated dead or follower forces/moments.
6. Sensitivities for the global beam responses for all the aforementioned analysis capabilities.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   convention
   input
   output
   bib
