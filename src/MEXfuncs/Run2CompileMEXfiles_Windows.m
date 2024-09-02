%% In case needed (only tested on Windows 10, 11)
clear all; clc;
mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" AssembleCmptStencilFromFinestLevel.cpp
mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" AssembleCmptStencilFromNonFinestLevel.cpp
mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_Interpolation_MatrixFree_mex.cpp
mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_KbyU_MatrixFree_mex.cpp
mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_KbyU_MatrixFree8x8_mex.cpp
mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_Restriction_MatrixFree_mex_advanced.cpp