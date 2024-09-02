
clear all; clc;
mex -setup C++
mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" AssembleCmptStencilFromFinestLevel.cpp
mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" AssembleCmptStencilFromNonFinestLevel.cpp
mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_Interpolation_MatrixFree_mex.cpp
mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_KbyU_MatrixFree_mex.cpp
mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_KbyU_MatrixFree8x8_mex.cpp
mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_Restriction_MatrixFree_mex_advanced.cpp
