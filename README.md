# SGLDBench
A Suite of Benchmarks for Stress-guided Lightweight Design

This suite is the associated software of the paper: 
	SGLDBench: A Benchmark Suite for Stress-Guided Lightweight 3D Designs
	by Junpeng Wang, Simon Niedermayr, Christoph Neuhauser, Dennis R. Bukenberger, Jun Wu, and Rüdiger Westermann
	arXiv: 
	
# Overview
SGLDBench is developed based on MatLab, but has extentions to call external modules implemented using Python or binary executables (on Windows).

At the bottom of SGLDBench is an efficient Geometric Multigrid Solver for high-resolution Static Finite Element Analysis Simulation on Cartesian grid
This solver is implemented in MatLab, the MEX function is used for higher efficiency.

Based on this solver, six different lightweight design strategies are included:

Method 1 --- Topology Optimization
Method 2 --- Porous Infill Optimizaiton
Method 3 --- Stress-aware Graded Voronoi Diagram
Method 4 --- Principal Stress Lines-guided Infill
Method 5 --- Stress-aligned Conforming Lattice
Method 6 --- Stress-aligned Volumetric Michell's Trusses

An easy-to-use WebGL-based render is provided to compare and analyze different designs.

# Dependence
To run SGLDBench, a MatLab release is needed at least, and its "Image Processing Toolbox" needs to be installed as well. 
Recommending versions R2022b and newer.

All MEX functions are already compiled on Windows and included in the repository. In case one needs to re-compile them for reason, a C/C++ compiler will
be needed as well. To compile them, one can directly run the MatLab script './SGLDBench/src/MEXfuncs/Run2CompileMEXfiles.m', this works for both Window and Linux OS.

To run Method 3 (Stress-aware Graded Voronoi Diagram), Python needs to be installed!

Method 1, 2, 4 can be seamlessly used on both Windows and Linux.

Method 3, 5, 6 only supports Windows for now since we use the Windows-compiled executable of TetGen (https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) 
to create the gateway tet-mesh for all these three methods.

# Usages
SGLDBench follows the pipeline of
Modeling -> Applying Boundary Conditions -> Numerical Simulation -> Visualization
