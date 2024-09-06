%% In case needed (only tested on Windows 10, 11)
clear all; clc;
if ispc
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_AssembleCmptStencilFromFinestLevel.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_AssembleCmptStencilFromNonFinestLevel.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_Interpolation_MatrixFree_mex.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_KbyU_MatrixFree_mex.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_KbyU_MatrixFree8x8_mex.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_Restriction_MatrixFree_mex_advanced.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" TopOpti_CmptUnitCompliance_mex.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" TopOpti_SetupDensityFilter_mex.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" TopOpti_PerformDensityFiltering_mex.cpp
	mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" MMA_mex.cpp MMAseq.cpp
elseif isunix
	mex -setup C++
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_AssembleCmptStencilFromFinestLevel.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_AssembleCmptStencilFromNonFinestLevel.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_Interpolation_MatrixFree_mex.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_KbyU_MatrixFree_mex.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_KbyU_MatrixFree8x8_mex.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_Restriction_MatrixFree_mex_advanced.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TopOpti_CmptUnitCompliance_mex.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TopOpti_SetupDensityFilter_mex.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" TopOpti_PerformDensityFiltering_mex.cpp
	mex -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" MMA_mex.cpp MMAseq.cpp
else
	disp('Platform not supported');
end