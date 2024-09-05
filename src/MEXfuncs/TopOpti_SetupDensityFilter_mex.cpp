#include "mex.h"
#include "matrix.h"
#include <algorithm>
#include <omp.h>
#include <cmath>
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))
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
	
	int kk;
    // Initialize Y to zeros
    #pragma omp parallel for schedule(static)	
    for (kk = 0; kk < numElements; ++kk) {
        Y[kk] = 0.0;
    }
	
	#pragma omp parallel
	{
		double weightsSum, e2e1Weight;
		#pragma omp for
		for (kk=1; kk<=resZ; ++kk) {
			for (int ii=1; ii<=resX; ++ii) {
				for (int jj=1; jj<=resY; ++jj) {
					int e1MapBack = (kk-1)*resX*resY+(ii-1)*resY+jj;
					int32_T e1 = eleMapForward[e1MapBack-1];
					if (e1) {
						weightsSum = 0;
						for (int kk2=MAX(kk-(ceil(rMin)-1),1); kk2<=MIN(kk+(ceil(rMin)-1), resZ); ++kk2) {
							for (int ii2=MAX(ii-(ceil(rMin)-1),1); ii2<=MIN(ii+(ceil(rMin)-1),resX); ++ii2) {
								for (int jj2=MAX(jj-(ceil(rMin)-1),1); jj2<=MIN(jj+(ceil(rMin)-1),resY); ++jj2) {
									int e2MapBack = (kk2-1)*resX*resY+(ii2-1)*resY+jj2;
									int32_T e2 = eleMapForward[e2MapBack-1];
									if (e2) {
										double distance = rMin - std::sqrt((ii-ii2)*(ii-ii2) + (jj-jj2)*(jj-jj2) + (kk-kk2)*(kk-kk2));
										e2e1Weight = MAX(0.0, distance);
										weightsSum += e2e1Weight;
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
}

