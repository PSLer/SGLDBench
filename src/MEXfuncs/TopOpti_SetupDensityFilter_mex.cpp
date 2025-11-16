#include "mex.h"
#include "matrix.h"
#include <algorithm>
#include <omp.h>
#include <cmath>

// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" TopOpti_SetupDensityFilter_mex.cpp"
/* Main MEX function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double rMin = (double) mxGetScalar(prhs[0]);
	int numElements = (int) mxGetScalar(prhs[1]);
	int32_T *eleMapForward = (int32_T *) mxGetData(prhs[2]);
	
	int resX = (int) mxGetScalar(prhs[3]);
	int resY = (int) mxGetScalar(prhs[4]);
	int resZ = (int) mxGetScalar(prhs[5]);
	
    plhs[0] = mxCreateDoubleMatrix(numElements, 1, mxREAL);
    double *Y = mxGetPr(plhs[0]);
	
	const int rMinCeilMinus1 = (int)std::ceil(rMin)-1;
	const int resXY = resX*resY;

	#pragma omp parallel for collapse(3) schedule(static)
	for (int kk=1; kk<=resZ; ++kk) {
		for (int ii=1; ii<=resX; ++ii) {
			for (int jj=1; jj<=resY; ++jj) {
				int e1MapBack = (kk-1)*resXY+(ii-1)*resY+jj;
				int32_T e1 = eleMapForward[e1MapBack-1];
				if (e1) {
					double weightsSum = 0;
					for (int kk2=std::max(kk-rMinCeilMinus1,1); kk2<=std::min(kk+rMinCeilMinus1, resZ); ++kk2) {
						for (int ii2=std::max(ii-rMinCeilMinus1,1); ii2<=std::min(ii+rMinCeilMinus1,resX); ++ii2) {
							for (int jj2=std::max(jj-rMinCeilMinus1,1); jj2<=std::min(jj+rMinCeilMinus1,resY); ++jj2) {
								int e2MapBack = (kk2-1)*resXY+(ii2-1)*resY+jj2;
								int32_T e2 = eleMapForward[e2MapBack-1];
								if (e2) {
									double distance = rMin - std::sqrt((ii-ii2)*(ii-ii2) + (jj-jj2)*(jj-jj2) + (kk-kk2)*(kk-kk2));
									weightsSum += std::max(0.0, distance);
								}	
							}
						}
					}
					Y[e1-1] = weightsSum;	
				}
			}
		}
	}
}

