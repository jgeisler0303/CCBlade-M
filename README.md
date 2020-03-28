# CCBlade-M
This a basic reimplementation for MATLAB of the python code from here: https://github.com/WISDEM/CCBlade

It is a solver for the blade element momentum (BEM) method for analyzing wind turbine aerodynamic forces along the length of a blade. It uses a special residual function and zero finding to achieve fast and robust (guaranteed) convergence. More information can be found in this paper: https://onlinelibrary.wiley.com/doi/abs/10.1002/we.1636

The code is basic insofar as it does not implement the following features
* Tilt and cone angle, blade pre-curve and pre-sweep
* Yaw misalignment
* Blade azimuth angle
* Wind shear
* Spline smoothing of cl and cd (linear interpolation is used)
* Glauert Buhl correction (Glauert correction is used instead)
* Interpolation of Reynolds number
* Span-wise interpolation of chord to thickness ratio
* Inverse analysis
* Calculation of derivatives

Deviating from the original code, the loads (forces and moments) are calculated assuming a linear variation of the line loads at the nodes.

The code expects the aerodynamic parameters in a structure that is closely resembling the structure of the FAST input files, but could be easily adapted to other sources.

The code is meant as a basis for further development in MATLAB and especially in C++ as mex functions or other standalone applications.

The code uses an implementation of Brent's Method from the German Wikipedia that seems to be slightly faster than the built-in MATLAB ``fzero`` function.

# Super fast MATLAB mex (C++) implemenation 
In addition to the standard ``CCBlade.m`` function there is an identical C++ implementaion of the code (``CCBlade_mex.cpp``) that can be compiled into a MATLAB Mex function (``CCBlade_mex``) with identical arguments and return value as ``CCBlade.m`` but several orders of magnitude faster.

The C++ code is self-sufficient and does not rely on any additional libraries. Especially, it contains an implementation of the interpolation function and Brent's Method.

In order to compile the mex function, a C++ compiler needs to be installed in MATLAB (tested with gcc). In Octave, install the ``liboctave-dev`` package. The function ``makeCCBlade_mex`` can be used to automatically compile the soure code with all necessary compiler flags.

# Test function
There is a test function ``testCCBlade`` that will:
* Load the aerodynamic parameters from the supplied FAST 5MW turbine configuration into a suitable structure.
* Pass the data to ``CCBlade`` to calculate the values of a single operating point.
* Compile the mex function.
* Call the mex implementation of CCBlade in a double loop to calculate an entire aerodynamic field. 
* Plot a surface plot of the cp field and the operating point from ``CCBlade``

The test function assumes that the FAST MATLAB utilities from this repository: https://github.com/OpenFAST/matlab-toolbox are located in a directory parallel to CCBlade.
