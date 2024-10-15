# SGLDBench

A Suite of Benchmarks for Stress-guided Lightweight Design

This repository is the associated software of the paper: 
> "SGLDBench: A Benchmark Suite for Stress-Guided Lightweight 3D Designs"
> by Junpeng Wang, Simon Niedermayr, Christoph Neuhauser, Dennis R. Bukenberger, Jun Wu, and Rüdiger Westermann.
> arXiv: 


## 1. Overview

SGLDBench is developed based on MatLab, but has extensions to call external modules implemented using Python or binary executables (on Windows).

At the bottom of SGLDBench is an efficient geometric multigrid solver for high-resolution static Finite Element Analysis (FEA) simulations on Cartesian grids.
This solver is implemented in MatLab, the MEX functionality is used for higher efficiency.

Built upon this solver, six different lightweight design strategies are provided:
- Method 1 --- Topology Optimization
- Method 2 --- Porous Infill Optimizaiton
- Method 3 --- Stress-aware Graded Voronoi Diagram
- Method 4 --- Principal Stress Lines-guided Infill
- Method 5 --- Stress-aligned Conforming Lattice
- Method 6 --- Stress-aligned Volumetric Michell's Trusses

An easy-to-use WebGL-based render is provided to compare and analyze the generated different designs.


## 2. Dependencies

To run SGLDBench, a MatLab release is needed at least, and its "Image Processing Toolbox" and "Parallel Computing Toolbox" need to be installed as well. 
Recommending versions R2022b and newer.

All MEX functions are already compiled on Windows and included in the repository. In case one needs to re-compile them for some reasons, a C/C++ compiler will
be needed. To compile them, one can directly run the MatLab script "./SGLDBench/src/MEXfuncs/Run2CompileMEXfiles.m", which works for both Window and Linux.

To run Method 3 (Stress-aware Graded Voronoi Diagram), **Python needs to be installed**. The README in directory `./SGLDBench/externalModules/GradedVoronoiDiagram`
instructs how to automatically install dependencies.

Methods 1, 2, 4 can be seamlessly used on both Windows and Linux.

Methods 3, 5, 6 only support Windows for now since we use the Windows-compiled executable of TetGen (https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) 
to create the gateway tet-mesh for all these three methods.

For Method 6, the repository is not directly included in SGLDBench, but one can be either fetched by calling
```shell
git submodule update --init --recursive
```
... or manually downloaded from https://github.com/rarora7777/VolumetricTruss together with its dependency "gptoolbox" from https://github.com/alecjacobson/gptoolbox.
In case a manual download is preferred, these repositories need to be put into the directory `./SGLDBench/externalModules/`.
SGLDBench takes the default folder names `VolumetricTruss` and `gptoolbox` and the suggested work directory to set the search path.


## 3. Usages

The entire process follows the pipeline of "Modeling -> Specifying Boundary Conditions -> Stress Analysis -> Structural Generation -> and Visual Analysis".

### 3.1 via GUI

One can launch the GUI of SGLDBench by running `./SGLDBench/SGLDBench_Main.m`. The demo video shows some guidelines on how to use it. 

### 3.2 via Script

In the directory `./SGLDBench/QuickAccess`, we also provide the scripts to quickly access the targeted methods.


## 4. Data

### 4.1 Input

SGLDBench primarily takes triangle surface meshes in ".obj" or ".ply" format as input to describe the design domain;

Besides, SGLDBench can also take the tailored FEM voxel model file (".TopVoxel") as input for Methods 1, 2, 4. This file includes information of voxel volume, 
boundary conditions, passive elements (optional), and the density value of each voxel (optional). 
Please refer to the demo dataset `./SGLDBench/data/Part_R256.TopVoxel` for details.

Three common shapes (cuboid, L-shape, and cylinder) in structural design and optimization are integrated in SGLDBench for testing. For the cuboid domain, several
built-in boundary conditions are also provided (through the GUI).

### 4.2 Output

All the (intermediate) output data is placed in the directory `./SGLDBench/out/`.
SGLDBench stores the design as volumetric data in the NIFTI format (".nii"). All the targeted results of these 6 methods are named `DesignVolume.nii`.
The corresponding stress-to-stress alignment metric is named `alignmentMetricVolume_byStress.nii`.

One can also save the created voxel model (incl. voxel volume, boundary conditions, ...) as an ASCII file (".TopVoxel") for repeated use.

## 5. Cite

If you use the code and data of SGLDBench, please cite it as

```bibtex
TODO - add .bib entry
```
