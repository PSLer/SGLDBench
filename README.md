# SGLDBench
A Suite of Benchmarks for Stress-guided Lightweight Design

This suite is the associated software of the paper: 
	SGLDBench: A Benchmark Suite for Stress-Guided Lightweight 3D Designs
	by Junpeng Wang, Simon Niedermayr, Christoph Neuhauser, Dennis R. Bukenberger, Jun Wu, and Rüdiger Westermann
	arXiv: 
	
# 1. Overview
SGLDBench is developed based on MatLab, but has extentions to call external modules implemented using Python or binary executables (on Windows).

At the bottom of SGLDBench is an efficient geometric multigrid solver for high-resolution static Finite Element Analysis simulation on Cartesian grid
This solver is implemented in MatLab, the MEX functionality is used for higher efficiency.

Built upon this solver, six different lightweight design strategies are provided:

Method 1 --- Topology Optimization
Method 2 --- Porous Infill Optimizaiton
Method 3 --- Stress-aware Graded Voronoi Diagram
Method 4 --- Principal Stress Lines-guided Infill
Method 5 --- Stress-aligned Conforming Lattice
Method 6 --- Stress-aligned Volumetric Michell's Trusses

An easy-to-use WebGL-based render is provided to compare and analyze the generated different designs.

# 2. Dependence
To run SGLDBench, a MatLab release is needed at least, and its "Image Processing Toolbox" and "Parallel Computing Toolbox" need to be installed as well. 
Recommending versions R2022b and newer.

All MEX functions are already compiled on Windows and included in the repository. In case one needs to re-compile them for some reasons, a C/C++ compiler will
be needed. To compile them, one can directly run the MatLab script "./SGLDBench/src/MEXfuncs/Run2CompileMEXfiles.m", this works for both Window and Linux.

To run Method 3 (Stress-aware Graded Voronoi Diagram), Python needs to be installed! The README in directory "./SGLDBench/externalModules/GradedVoronoiDiagram"
instructs how to automatically install dependencies.

Method 1, 2, 4 can be seamlessly used on both Windows and Linux.

Method 3, 5, 6 only supports Windows for now since we use the Windows-compiled executable of TetGen (https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) 
to create the gateway tet-mesh for all these three methods.

For Method 6, the repository is not directly included in SGLDBench, but one can download it from https://github.com/rarora7777/VolumetricTruss and its 
dependency "gptoolbox" from https://github.com/alecjacobson/gptoolbox/. Then, putting them into the directory './SGLDBench/externalModules/', SGLDBench takes
the default folder names "VolumetricTruss-master" and "gptoolbox-master" and the suggested directory to set the search path.

# 3. Usages
Upon launching SGLDBench by running "./SGLDBench/SGLDBench_Main.m", the entire process follows the pipeline of
"Modeling -> Applying Boundary Conditions -> Numerical Simulation -> Visualization"
# 3.1 Modeling

# 3.2 Applying Boundary Conditions

# 3.3 Numerical Simulation

# 3.4 Visualization

# 4. Citation
If you use the code and data of SGLDBench, please cite it as 

to-dos
1. Initialize U_ in Stiffness Evaluation instead of Stress Analysis>>>>>>>>>>>
2. include density values into .TopVoxel
3. Supports to evaluate external density layout or graph
3. nodeCoords_, eleCentroids_