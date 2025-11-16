#include "mex.h"
#include <omp.h>
#include <stdlib.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp:experimental /std:c++20" TopOpti_CmptUnitCompliance_mex.cpp"
inline double compute_ce(const double* __restrict Ke, const double* __restrict Ue);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "3 inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "One output required.");
    }

    // Retrieve the inputs
    const mxArray *U_mx        = prhs[0];
    const mxArray *eNodMat_mx  = prhs[1];
    const mxArray *Ke_mx       = prhs[2];
	
    double *U = mxGetPr(U_mx);
    int32_T *eNodMat = (int32_T *)mxGetData(eNodMat_mx);
    double *Ke = mxGetPr(Ke_mx);
    mwSize numElements = mxGetM(eNodMat_mx);


    // Validate input dimensions
    if (!mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input U must be a column vector.");
    }
    if (!mxIsInt32(prhs[1]) || mxGetN(prhs[1]) != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input eNodMat must be an Mx8 int32 matrix.");
    }
    if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 24 || mxGetN(prhs[2]) != 24) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNot24x24", "Input Ke must be a 24x24 double matrix.");
    }

    plhs[0] = mxCreateDoubleMatrix(numElements, 1, mxREAL);
    double *Y = mxGetPr(plhs[0]);
	
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < numElements; ++i) {
		int32_T ieNod[8];
        double  Ue[24];
		for (int j = 0; j < 8; ++j) {
			ieNod[j] = eNodMat[i + j * numElements];
		}
			
		// Create Ue matrix with dimensions [24,1]
		for (int j = 0; j < 8; ++j) {
            mwIndex node = (mwIndex)ieNod[j] - 1;
            mwIndex baseIndex = 3 * node;
			Ue[3 * j + 0] = U[baseIndex + 0];
			Ue[3 * j + 1] = U[baseIndex + 1];
			Ue[3 * j + 2] = U[baseIndex + 2];			
		}
	
		Y[i] = compute_ce(Ke, Ue);
    }
}

inline double compute_ce(const double* __restrict Ke,
                         const double* __restrict Ue)
{
    double y[24];
#if 1
    // y = Ke * Ue
    #pragma omp simd
    for (int i = 0; i < 24; ++i) {
        y[i] = 0.0;
    }

    // Accumulate column by column (contiguous in memory)
    for (int j = 0; j < 24; ++j) {
        double uj = Ue[j];
        const double* col = Ke + j*24;  // j-th column, 24 contiguous entries

        #pragma omp simd
        for (int i = 0; i < 24; ++i) {
            y[i] += col[i] * uj;
        }
    }

    // q = Ue^T * y
    double q_local = 0.0;
    #pragma omp simd reduction(+:q_local)
    for (int i = 0; i < 24; ++i) {
        q_local += Ue[i] * y[i];
    }
#else
	for (int i = 0; i < 24; ++i) {
		double sum = 0.0;
		#pragma omp simd reduction(+:sum)
		for (int j = 0; j < 24; ++j) {
			sum += Ke[i + j*24] * Ue[j];
		}
		y[i] = sum;
	}
	
	// q = Ue^T * y
	double q_local = 0.0;
	#pragma omp simd reduction(+:q_local)
	for (int i = 0; i < 24; ++i) {
		q_local += Ue[i] * y[i];
	}	
#endif
    return q_local;
}
